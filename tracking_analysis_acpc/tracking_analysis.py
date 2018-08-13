import matplotlib
matplotlib.use('Agg')
import logging
import matplotlib.pyplot as plt
#import wrfcube
#import ramscube
#from wrfmpdiag import calculate_wrf_mp_path,calculate_rams_mp_path,lump_processes
from cloudtrack import maketrack,segmentation_3D, segmentation_2D
from cloudtrack import plot_tracks_mask_field_loop
from cloudtrack import plot_mask_cell_track_follow,plot_mask_cell_track_static,plot_mask_cell_track_static_timeseries
from .cell_analysis import extract_cell_cubes_subset
import iris
from iris.analysis import MAX
import os
import pandas as pd
import numpy as np

from matplotlib import colors
#from analysis_pathways.extract_cubes_pathways import track_cell_orientation_subset,load_AUX

from html_animate import make_animation

from multiprocessing.pool import ThreadPool
import dask

from wrfmpdiag import processes_colors

processes_colors, processes_names = processes_colors(microphysics_scheme='RAMS', colors_processes='lumped')

def add_geopotential_height(cube):
    from .make_geopotential_height_coord import geopotential_height_coord,geopotential_height_coord_stag
    coord_names=[coord.name() for coord in cube.coords()]
    if ('model_level_number' in coord_names) and ('geopotential_height' not in coord_names):
        if cube.coord('model_level_number').shape[0]==95:
            cube.add_aux_coord(geopotential_height_coord_stag,data_dims=cube.coord_dims('model_level_number'))
        if cube.coord('model_level_number').shape[0]==94:
            cube.add_aux_coord(geopotential_height_coord,data_dims=cube.coord_dims('model_level_number'))
    return cube



def tracking_analysis(filenames,
                      mode,
                      top_savedir_data,
                      top_savedir_tracking,
                      top_plotdir,
                      load_function,
                      cell_selection=None,
                      parallel=False,
                      loading_period=None,
                      tracking_period=None,
                      parameters_tracking=None,
                      parameters_segmentation_TWC=None,
                      parameters_segmentation_w=None,
                      plotting_parameters=None,
                      plotting_period=None,
                      n_core=1):
    
    print('n_core: '+str(n_core))
    logging.debug('ncore: '+str(n_core) )
    
    if n_core==1:
        dask.config.set(scheduler='single-threaded')
    elif n_core>1:
        dask.set_options(pool=ThreadPool(n_core-1))
        
    logging.debug('start tracking analysis')
    logging.debug('top_savedir_data: '+ top_savedir_data)
    logging.debug('top_savedir_tracking: '+ top_savedir_tracking)
    logging.debug('top_plotdir: '+ top_plotdir)

    logging.debug('number of processors used :'+ str(n_core))
   
    if loading_period:
        constraint_loading_period=iris.Constraint(time = lambda cell: loading_period[0] <= cell <  loading_period[1])
    elif loading_period is None:
        constraint_loading_period=iris.Constraint(None)

    if tracking_period:
        constraint_tracking_period=iris.Constraint(time = lambda cell: tracking_period[0] <= cell <  tracking_period[1])
    elif tracking_period is None:
        constraint_tracking_period=iris.Constraint(None)

    #Tracking setup:


    if 'load_data' in mode:
        
        logging.info('start loading data for tracking')

        savedir=os.path.join(top_savedir_data)
        os.makedirs(os.path.join(savedir),exist_ok=True)
        
        logging.debug(' start loading w')
        W=load_function(filenames,'w').extract(constraint_loading_period)                 
        iris.save([W],os.path.join(savedir,'Data_w.nc'))    
        logging.debug(' w loaded and saved')
        
        W_max=W.collapsed('model_level_number',MAX)
        logging.debug('w_max calculated')
        iris.save([W_max],os.path.join(savedir,'Data_w_max.nc'))    
        logging.debug('w max loaded and saved')

        # Create 2D field of maximum midlevel vertical velocity
        # Constraint for choosing midlevel heights for w mask
        constraint_midlevel=iris.Constraint(model_level_number= lambda cell: 30 < cell < 48)
        W_mid = W.extract(constraint_midlevel) # get midlevel data
        W_mid_max = W_mid.collapsed('model_level_number',MAX) # get Maximum column updraft in mid-level data
        iris.save([W_mid_max],os.path.join(savedir,'Data_w_mid_max.nc'))    
        logging.debug('w_mid_max loaded and saved')
        
        del(W,W_max,W_mid,W_mid_max)

        
        logging.debug(' start loading water content')
        IWC=load_function(filenames,'IWC').extract(constraint_loading_period)  
        LWC=load_function(filenames,'LWC').extract(constraint_loading_period)  
        TWC=IWC+LWC
        TWC.rename('total water content')
        iris.save([TWC,LWC,IWC],os.path.join(savedir,'Data_WC.nc'))#                
        logging.debug('water content loaded and saved')
        del(TWC,IWC,LWC)


        logging.debug(' start loading water path')
        Airmass=load_function(filenames,'airmass').extract(constraint_loading_period) 
        logging.debug('airmass loaded')

        Airmass_path=load_function(filenames,'airmass_path').extract(constraint_loading_period) 
        logging.debug('airmass path loaded')

        IWP=load_function(filenames,'IWP').extract(constraint_loading_period) 
        logging.debug('ice water path loaded')

        
        LWP=load_function(filenames,'LWP').extract(constraint_loading_period)      
        logging.debug('liquid water path loaded')

        TWP=IWP+LWP
        TWP.rename('TWP')

        iris.save([TWP,LWP,IWP,Airmass,Airmass_path],os.path.join(savedir,'Data_WP.nc'))
        del(TWP,IWP,LWP,Airmass_path,Airmass)

        logging.debug('water path loaded and saved')
        
        logging.debug(' start process rates')

        processes_lumped=load_function(filenames,'processes_lumped').extract(constraint_tracking_period)                 
        iris.save(processes_lumped,os.path.join(savedir,'Data_Processes_lumped.nc'))
        logging.debug('loaded and saved processes')
        del(processes_lumped)

        logging.debug('data loaded and saved')

        
                
    if 'tracking' in mode:
        
        logging.info('start tracking')

        logging.debug('start loading data for tracking')

        loaddir=os.path.join(top_savedir_data)
        savedir=os.path.join(top_savedir_tracking)
        os.makedirs(savedir,exist_ok=True)

        W_max=iris.load_cube(os.path.join(loaddir,'Data_w_max.nc')).extract(constraint_tracking_period)  
        W_mid_max=iris.load_cube(os.path.join(loaddir,'Data_w_mid_max.nc')).extract(constraint_tracking_period)  

        
        logging.debug('data loaded for tracking')

        logging.info('start tracking based on midlevel maximum vertical velocity:')
        Track,Features,Features_unfiltered,Track_unfilled=maketrack(W_mid_max,return_intermediate=True,**parameters_tracking)
        
        Track.to_hdf(os.path.join(savedir,'Track.hdf'),'table')
        Features.to_hdf(os.path.join(savedir,'Features.hdf'),'table')
        Features_unfiltered.to_hdf(os.path.join(savedir,'Features_unfiltered.hdf'),'table')
        Track_unfilled.to_hdf(os.path.join(savedir,'Track_unfilled.hdf'),'table')

        logging.info('tracking based on midlevel maximum vertical velocity finished and saved')
        
        logging.info('start tracking based on column maximum vertical velocity:')

        Track_w_max,Features_w_max,Features_unfiltered_w_max,Track_unfilled_w_max=maketrack(W_max,return_intermediate=True,**parameters_tracking)
    
        Track_w_max.to_hdf(os.path.join(savedir,'Track_w_max.hdf'),'table')
        Features_w_max.to_hdf(os.path.join(savedir,'Features_w_max.hdf'),'table')
        Features_unfiltered_w_max.to_hdf(os.path.join(savedir,'Features_unfiltered_w_max.hdf'),'table')
        Track_unfilled_w_max.to_hdf(os.path.join(savedir,'Track_unfilled_w_max.hdf'),'table')
        
        logging.info('tracking based on column maximum vertical velocity finished and saved')

        logging.info('tracks calculated and saved')

        logging.info('tracks calculated and saved')

            
    if 'segmentation' in mode:
                
        logging.info('start watershedding')
        # Set up loading and saving directories:
        loaddir=os.path.join(top_savedir_data)
        savedir=os.path.join(top_savedir_tracking)
        os.makedirs(savedir,exist_ok=True)
        
         # Read Trajectories:
        Track=pd.read_hdf(os.path.join(savedir,'Track.hdf'),'table')

        # perform 3D segmentation based on total condenate
        logging.info('start watershedding TWC')
        TWC=iris.load_cube(os.path.join(loaddir,'Data_WC.nc'),'total water content').extract(constraint_tracking_period)        
        Mask_TWC,Track_TWC=segmentation_3D(Track,TWC,**parameters_segmentation_TWC)
        iris.save([Mask_TWC],os.path.join(savedir,'Watershed_TWC.nc'))                
        Track_TWC.to_hdf(os.path.join(savedir,'Track_TWC.hdf'),'table')
        logging.debug('segmentation TWC performed and saved')
        del(TWC,Mask_TWC,Track_TWC)
        
        # perform 3D segmentation based on updraft velocity
        logging.info('start watershedding w')
        W=iris.load_cube(os.path.join(loaddir,'Data_w.nc')).extract(constraint_tracking_period)             
        Mask_w,Track_w=segmentation_3D(Track,W,**parameters_segmentation_w)
        iris.save([Mask_w],os.path.join(savedir,'Watershed_w.nc'))
        Track_w.to_hdf(os.path.join(savedir,'Track_w.hdf'),'table')
        logging.debug('segmentation w performed and saved')
        del(W,Mask_w,Track_w)

        # perform 2D segmentation based maximum midlevel updraft velocity (used for tracking)
        logging.debug('start watershedding w_mid_max')
        W_mid_max=iris.load_cube(os.path.join(loaddir,'Data_w_mid_max.nc')).extract(constraint_tracking_period)         
        Mask_w_mid_max,Track_w_mid_max=segmentation_2D(Track,W_mid_max,**parameters_segmentation_w)
        iris.save([Mask_w_mid_max],os.path.join(savedir,'Watershed_w_mid_max.nc'))                
        Track_w_mid_max.to_hdf(os.path.join(savedir,'Track_w_mid_max.hdf'),'table')
        logging.debug('segmentation w_mid_max performed and saved')
        del(W_mid_max,Mask_w_mid_max,Track_w_mid_max)

    

    
        #%%        
    if 'plot' in mode:
        savedir_data=os.path.join(top_savedir_data)
        savedir_tracking=os.path.join(top_savedir_tracking)
        plotdir=os.path.join(top_plotdir)
    
        Track=pd.read_hdf(os.path.join(savedir_tracking,'Track.hdf'),'table')
        Features=pd.read_hdf(os.path.join(savedir_tracking,'Features.hdf'),'table')

        W=iris.load_cube(os.path.join(savedir_data,'Data_w.nc')).extract(constraint_tracking_period)     
        W_mid_max=iris.load_cube(os.path.join(savedir_data,'Data_w_mid_max.nc')).extract(constraint_tracking_period)     
        W_max=iris.load_cube(os.path.join(savedir_data,'Data_w_max.nc')).extract(constraint_tracking_period)     

#        # TWP=iris.load_cube(os.path.join(savedir_data,'Data.nc'),'total water path')    
#        # TWC=iris.load_cube(os.path.join(savedir_data,'Data.nc'),'total water content')    
#    
        Mask_TWC=iris.load_cube(os.path.join(savedir_tracking,'Watershed_TWC.nc'),'watershedding_output_mask')
#        Mask_w=iris.load_cube(os.path.join(savedir_tracking,'Watershed_w.nc'),'watershedding_output_mask')
        Mask_w_mid_max=iris.load_cube(os.path.join(savedir_tracking,'Watershed_w_mid_max.nc'),'watershedding_output_mask')
#       
        logging.debug('data loaded from savefile')

        axis_extent=plotting_parameters['axis_extent']
    #%%

#        plot_tracks_mask_field_loop(Track,W_max,Mask_w,
#                                    plot_dir=os.path.join(plotdir,'w_max_Mask_w'),name='w_max_Mask_w',
#                                     axis_extent=axis_extent,#figsize=figsize,orientation_colorbar='horizontal',pad_colorbar=0.2,
#                                     vmin=0,vmax=30,
#                                     plot_outline=True,plot_marker=True,marker_track=',',plot_number=True)
#        logging.debug('W max Tracks plotted')
        
        logging.debug('start plotting W max Tracks')

        plot_tracks_mask_field_loop(track=Track,field=W_max,mask=Mask_TWC,features=Features,
                                    plot_dir=os.path.join(plotdir,'w_max_Mask_TWC'),name='w_max_Mask_TWC',
                                    axis_extent=axis_extent,#figsize=figsize,orientation_colorbar='horizontal',pad_colorbar=0.2,
                                    vmin=0,vmax=30,
                                    plot_outline=True,plot_marker=True,marker_track='x',plot_number=True,plot_features=True)
        logging.debug('W max Tracks plotted')
        
        make_animation(input_dir=os.path.join(plotdir,'w_max_Mask_TWC'),
                       output=os.path.join(plotdir,'Animation_w_max_Mask_TWC','w_max_Mask_TWC.html'))
        
        logging.debug('start plotting W mid max Tracks')

        plot_tracks_mask_field_loop(track=Track,field=W_max,mask=Mask_w_mid_max,features=Features,
                                    plot_dir=os.path.join(plotdir,'w_max_Mask_w_mid_max'),name='w_max_Mask_w_mid_max',
                                    axis_extent=axis_extent,#figsize=figsize,orientation_colorbar='horizontal',pad_colorbar=0.2,
                                    vmin=0,vmax=30,
                                    plot_outline=True,plot_marker=True,marker_track='x',plot_number=True,plot_features=True)
        logging.debug('W mid max Tracks plotted')
        
        make_animation(input_dir=os.path.join(plotdir,'w_max_Mask_w_mid_max'),
                       output=os.path.join(plotdir,'Animation_w_max_Mask_w_mid_max','w_max_Mask_w_mid_max.html'))


    if 'plot_mask' in mode:
        
        
        logging.debug('Start plotting tracks and masks for individual cells')

        savedir_data=os.path.join(top_savedir_data)
        savedir_tracking=os.path.join(top_savedir_tracking)
        plotdir=os.path.join(top_plotdir)
        
        features=pd.read_hdf(os.path.join(savedir_tracking,'Features.hdf'),'table')
        track_TWC=pd.read_hdf(os.path.join(savedir_tracking,'Track_TWC.hdf'),'table')
        track_w=pd.read_hdf(os.path.join(savedir_tracking,'Track_TWC.hdf'),'table')

        W_max=iris.load_cube(os.path.join(savedir_data,'Data_w_max.nc')).extract(constraint_tracking_period)    
        W_mid_max=iris.load_cube(os.path.join(savedir_data,'Data_w_mid_max.nc')).extract(constraint_tracking_period)    
        TWP=iris.load_cube(os.path.join(savedir_data,'Data_WP.nc'),'TWP').extract(constraint_tracking_period)  

        mask_TWC=iris.load_cube(os.path.join(savedir_tracking,'Watershed_TWC.nc'),'watershedding_output_mask')
        mask_w=iris.load_cube(os.path.join(savedir_tracking,'Watershed_w.nc'),'watershedding_output_mask')

        logging.debug('data loaded from savefile')    
                                 
#        if cell_selection is None:
#            cell_selection=np.unique(track['particle'])
        if plotting_period:
            track_TWC_time=track_TWC.loc[(track_TWC['time']>plotting_period[0]) & (track_TWC['time']<plotting_period[1])]
            track_w_time=track_TWC.loc[(track_w['time']>plotting_period[0]) & (track_w['time']<plotting_period[1])]

#        cell_selection=track_TWC_time.groupby('particle').agg({'ncells':'max'}).nlargest(columns='ncells',n=50).index.values
        if cell_selection is None:
            cell_selection=track_TWC_time['particle'].unique()

        cmap_TWP = colors.LinearSegmentedColormap.from_list("", ["white","lightblue",'lightsteelblue','cornflowerblue','royalblue','rebeccapurple','indigo'])
#        bounds_TWP = [2, 3, 7, 9, 15]
#        norm_TWP = colors.BoundaryNorm(bounds_TWP, cmap_TWP.N)
        levels_TWP=[0,0.5,1,1.5,2,3,4,5,10,15,20,30,40,50]
        
        cmap_w = colors.LinearSegmentedColormap.from_list("", ["forestgreen",'olivedrab','yellowgreen','khaki','gold', ])
#        cmap_w = colors.LinearSegmentedColormap.from_list("", ['orange',"forestgreen",'red','yellowgreen','purple','gold', ])

        levels_w=[1,3,5,10,15]
#        bounds_w=np.concat(np.array(levels_w)-np.diff(levels_w)
        bounds_w=[0.5, 2, 4, 7, 12, 18]
        norm_w = colors.BoundaryNorm(bounds_w, cmap_w.N)

        for particle in cell_selection:  
            logging.debug('Start plotting mask for cell'+ str(particle))
            # Make plots centred around and moving with the tracked cell (10km each side)


#            plot_mask_cell_track_follow(particle=particle,
#                                           track=track_TWC, 
#                                           cog=None,
#                                           features=features,
#                                           mask_total=mask_TWC,
#                                           field_contour=W_mid_max, 
#                                           label_field_contour='midlevel max w (m/s)',
#                                           cmap_field_contour=cmap_w,
#                                           vmin_field_contour=0,vmax_field_contour=15,levels_field_contour=[0,1,2,5,10,15],
#                                           field_filled=TWP, 
#                                           label_field_filled='TWP (kg/m^2)',
#                                           cmap_field_filled=cmap_TWP,
#                                           vmin_field_filled=0,vmax_field_filled=50,levels_field_filled=[0,0.5,1,1.5,2,3,4,5,10,15,20,30,40,50],
#                                           width=10000,
#                                           name= 'Masking_TWC_follow_'+str(int(particle)), 
#                                           plotdir=os.path.join(plotdir,'Masking_TWC_follow'),
#                                           )
#            
#            make_animation(input_dir=os.path.join(plotdir,'Masking_TWC_follow','Masking_TWC_follow_'+str(int(particle))),
#                           output=os.path.join(plotdir,'Animations_Masking_TWC_follow','Masking_TWC_follow_'+str(int(particle)))
#                            ,delay=200)
#            
#             Make plots static and covering the whole path of the tracked cell:
#
#            plot_mask_cell_track_static(particle=particle,
#                                           track=track_TWC, 
#                                           cog=None,
#                                           features=features,
#                                           mask_total=mask_TWC,
#                                           field_contour=W_mid_max, 
#                                           label_field_contour='midlevel max w (m/s)',
#                                           cmap_field_contour=cmap_w,
#                                           vmin_field_contour=0,vmax_field_contour=15,levels_field_contour=[0,1,2,5,10,15],
#                                           field_filled=TWP, 
#                                           label_field_filled='TWP (kg/m^2)',
#                                           cmap_field_filled=cmap_TWP,
#                                           vmin_field_filled=0,vmax_field_filled=50,levels_field_filled=[0,0.5,1,1.5,2,3,4,5,10,15,20,30,40,50],
#                                           width=10000,n_extend=2,
#                                           name= 'Masking_TWC_static_'+str(int(particle)), 
#                                           plotdir=os.path.join(plotdir,'Masking_TWC_static'),
#                                           )
#            make_animation(input_dir=os.path.join(plotdir,'Masking_TWC_static','Masking_TWC_static_'+str(int(particle))),
#                           output=os.path.join(plotdir,'Animations_Masking_TWC_static','Masking_TWC_static_'+str(int(particle))))
#
            
            plot_mask_cell_track_static_timeseries(particle=particle,
                                           track=track_TWC, 
                                           cog=None,
                                           features=features,
                                           mask_total=mask_TWC,
                                           field_contour=W_mid_max, 
                                           label_field_contour='midlevel max w (m/s)',
                                           cmap_field_contour=cmap_w,norm_field_contour=norm_w,
                                           vmin_field_contour=0,vmax_field_contour=15,levels_field_contour=levels_w,
                                           field_filled=TWP, 
                                           label_field_filled='TWP (kg/m^2)',
                                           cmap_field_filled=cmap_TWP,
                                           vmin_field_filled=0,vmax_field_filled=50,levels_field_filled=levels_TWP,
                                           track_variable=track_TWC,variable='ncells',variable_label='number of cells (watershedding)',
                                           width=10000,n_extend=2,
                                           name= 'Masking_TWC_static_ncells_'+str(int(particle)), 
                                           plotdir=os.path.join(plotdir,'Masking_TWC_static_ncells'),
                                           )
            make_animation(input_dir=os.path.join(plotdir,'Masking_TWC_static_ncells','Masking_TWC_static_ncells_'+str(int(particle))),
                           output=os.path.join(plotdir,'Animations_Masking_TWC_static_ncells','Masking_TWC_static_ncells_'+str(int(particle))),
                           delay=200)

        if cell_selection is None:
            cell_selection_w=track_w_time['particle'].unique()

        for particle in cell_selection_w:  
            
#            logging.debug('Start plotting mask for cell'+ str(particle))
#
#             #Make plots static and covering the whole path of the tracked cell:
#            plot_mask_cell_track_static(particle=particle,
#                                           track=track_TWC, 
#                                           cog=None,
#                                           features=features,
#                                           mask_total=mask_w,
#                                           field_contour=W_mid_max, 
#                                           label_field_contour='midlevel max w (m/s)',
#                                           cmap_field_contour=cmap_w,
#                                           vmin_field_contour=0,vmax_field_contour=15,levels_field_contour=[0,1,2,5,10,15],
#                                           field_filled=W_max, 
#                                           label_field_filled='max w (kg/m^2)',
#                                           cmap_field_filled=cmap_TWP,
#                                           vmin_field_filled=0,vmax_field_filled=15,levels_field_filled=[0,1,2,5,10,15],
#                                           width=10000,
#                                           name= 'Masking_w_static_'+str(int(particle)), 
#                                           plotdir=os.path.join(plotdir,'Masking_w_static'),
#                                           )
#            make_animation(input_dir=os.path.join(plotdir,'Masking_w_static','Masking_w_static_'+str(int(particle))),
#                           output=os.path.join(plotdir,'Animations_Masking_w_static','Masking_w_static_'+str(int(particle)))
#                            ,delay=200)
#
            
            logging.debug('Start plotting mask for cell'+ str(particle))

            plot_mask_cell_track_static_timeseries(particle=particle,
                                           track=track_w, 
                                           cog=None,
                                           features=features,
                                           mask_total=mask_w,
                                           field_contour=W_mid_max, 
                                           label_field_contour='midlevel max w (m/s)',
                                           cmap_field_contour=cmap_w,norm_field_contour=norm_w,
                                           vmin_field_contour=0,vmax_field_contour=15,levels_field_contour=levels_w,
                                           field_filled=W_max, 
                                           label_field_filled='max w (kg/m^2)',
                                           cmap_field_filled=cmap_TWP,
                                           vmin_field_filled=0,vmax_field_filled=15,levels_field_filled=levels_w,
                                           track_variable=track_TWC,variable='ncells',variable_label='number of cells (watershedding)',
                                           width=10000,n_extend=2,
                                           name= 'Masking_w_static_ncells_'+str(int(particle)), 
                                           plotdir=os.path.join(plotdir,'Masking_w_static_ncells')
                                           )
            make_animation(input_dir=os.path.join(plotdir,'Masking_w_static_ncells','Masking_w_static_ncells_'+str(int(particle))),
                           output=os.path.join(plotdir,'Animations_Masking_w_static_ncells','Masking_w_static_ncells_'+str(int(particle))),
                           delay=200)

    if 'interpolation_mass' in mode:     
        
        # Calculated summed up profiles and integrated values (hydrometeor mass and process rates) following each cell:
        logging.debug('Start calculations for individual cells')
            
        logging.debug('loading tracks and Mask')

        tracks=pd.read_hdf(os.path.join(top_savedir_tracking,'Track.hdf'),'table')
        mask=iris.load_cube(os.path.join(top_savedir_tracking,'Watershed_TWC.nc'),'watershedding_output_mask')
        
        TWC=iris.load_cube(os.path.join(top_savedir_data,'Data_WC.nc'),'total water content').extract(constraint_tracking_period)        
        logging.debug('loaded total water content from file')
        for coord in TWC.coords():
            coord.var_name=coord.name()
        
#        WC=iris.load(os.path.join(top_savedir_data,'Data_WC.nc')).extract(constraint_tracking_period)        
#        logging.debug('loaded water content from file')
#        for cube in WC:
#            for coord in cube.coords():
#                coord.var_name=coord.name()
        
        # load airmass (used to calculate total values from mixing ratios)
        airmass=iris.load_cube(os.path.join(top_savedir_data,'Data_WP.nc'),'airmass').extract(constraint_tracking_period) 
        for coord in airmass.coords():
            coord.var_name=coord.name()
        logging.debug('airmass loaded from file')

        
        if cell_selection is None:
            cell_selection=tracks['particle'].unique()

        # create empty list to store DataFrames with integrated values for each individual cell so that these can joined together in one at the end
        list_track=[]
        
        
        # loop over individual cells for analysis
        for particle in cell_selection:                        
            
            logging.debug( 'Start calculating Data for cell '+ str(int(particle)))
            savedir_cell=os.path.join(top_savedir_tracking,'cells',str(int(particle)))
            os.makedirs(savedir_cell,exist_ok=True)

            height_levels=np.arange(0,20001,1000.)
            
            mass_sum=iris.cube.CubeList()
            
            for cube in [TWC]:
                cube_sum=cube*airmass
                cube_sum.rename('total_water_mass')
                cube_sum=add_geopotential_height(cube_sum)                
                mass_sum.append(cube_sum)
                
            mass_cell_sum, mass_cell_integrated,track_mass_integrated=extract_cell_cubes_subset(cubelist_in=mass_sum,
                                                                                                mask=mask,
                                                                                                track=tracks,
                                                                                                particle=particle,
                                                                                                z_coord='model_level_number',
                                                                                                height_levels=height_levels)
            
            iris.save(mass_cell_sum,os.path.join(savedir_cell,'Mass_profile_'+str(int(particle))+'.nc'))
            iris.save(mass_cell_integrated,os.path.join(savedir_cell,'Mass_integrated_'+str(int(particle))+'.nc'))
            track_mass_integrated.to_hdf(os.path.join(savedir_cell,'track_mass_integrated_'+str(int(particle))+'.h5'),'table')
            
            logging.debug( 'mass calculated and saved to '+ savedir_cell)

            logging.debug( 'Data calculated and saved to '+ savedir_cell)
            
            list_track.append(track_mass_integrated)
        
        # concatenate DataFrames from intividual cells into one DataFrame containing all integrated values for all cells and save to file
        track_mass=pd.concat(list_track)
        track_mass.to_hdf(os.path.join(top_savedir_tracking,'track_mass_integrated.h5'),'table')

        logging.debug( 'Data calculated and saved to file for all cells')


    if 'plot_mask_mass' in mode:
        
        
        logging.debug('Start plotting tracks and masks for individual cells')

        savedir_data=os.path.join(top_savedir_data)
        savedir_tracking=os.path.join(top_savedir_tracking)
        plotdir=os.path.join(top_plotdir)
        
        features=pd.read_hdf(os.path.join(savedir_tracking,'Features.hdf'),'table')
        track_TWC=pd.read_hdf(os.path.join(savedir_tracking,'Track_TWC.hdf'),'table')
        track_w=pd.read_hdf(os.path.join(savedir_tracking,'Track_TWC.hdf'),'table')

        W_max=iris.load_cube(os.path.join(savedir_data,'Data_w_max.nc')).extract(constraint_tracking_period)    
        W_mid_max=iris.load_cube(os.path.join(savedir_data,'Data_w_mid_max.nc')).extract(constraint_tracking_period)    
        TWP=iris.load_cube(os.path.join(savedir_data,'Data_WP.nc'),'TWP').extract(constraint_tracking_period)  

        mask_TWC=iris.load_cube(os.path.join(savedir_tracking,'Watershed_TWC.nc'),'watershedding_output_mask')
        mask_w=iris.load_cube(os.path.join(savedir_tracking,'Watershed_w.nc'),'watershedding_output_mask')

        logging.debug('data loaded from savefile')    
                                 
#        if cell_selection is None:
#            cell_selection=np.unique(track['particle'])
        if plotting_period:
            track_TWC_time=track_TWC.loc[(track_TWC['time']>plotting_period[0]) & (track_TWC['time']<plotting_period[1])]
            track_w_time=track_TWC.loc[(track_w['time']>plotting_period[0]) & (track_w['time']<plotting_period[1])]

#        cell_selection=track_TWC_time.groupby('particle').agg({'ncells':'max'}).nlargest(columns='ncells',n=50).index.values
        if cell_selection is None:
            cell_selection=track_TWC_time['particle'].unique()

        cmap_TWP = colors.LinearSegmentedColormap.from_list("", ["white","lightblue",'lightsteelblue','cornflowerblue','royalblue','rebeccapurple','indigo'])
        levels_TWP=[0,0.5,1,1.5,2,3,4,5,10,15,20,30,40,50]

        cmap_w = colors.LinearSegmentedColormap.from_list("", ["forestgreen",'olivedrab','yellowgreen','khaki','gold'])
        levels_w=[1,2,5,10,15]


        track_mass=pd.read_hdf(os.path.join(savedir_tracking,'track_mass_integrated.h5'),'table')

        processes=list(processes_names.keys())
        colors_processes=processes_colors
 
        for particle in cell_selection:  
            logging.debug('Start plotting mask for cell'+ str(particle))
            
            plot_mask_cell_track_static_timeseries(particle=particle,
                                           track=track_mass, 
                                           cog=None,
                                           features=features,
                                           mask_total=mask_TWC,
                                           field_contour=W_mid_max, 
                                           label_field_contour='midlevel max w (m/s)',
                                           cmap_field_contour=cmap_w,
                                           vmin_field_contour=0,vmax_field_contour=15,levels_field_contour=levels_w,
                                           field_filled=TWP, 
                                           label_field_filled='TWP (kg/m^2)',
                                           cmap_field_filled=cmap_TWP,
                                           vmin_field_filled=0,vmax_field_filled=50,levels_field_filled=levels_TWP,
                                           track_variable=track_mass,variable='total_water_mass',
                                           variable_ylabel='total condensate mass (kg)',
                                           width=10000,n_extend=2,
                                           name= 'Masking_TWC_static_mass_'+str(int(particle)), 
                                           plotdir=os.path.join(plotdir,'Masking_TWC_static_mass'),
                                           )
            make_animation(input_dir=os.path.join(plotdir,'Masking_TWC_static_mass','Masking_TWC_static_mass_'+str(int(particle))),
                           output=os.path.join(plotdir,'Animations_Masking_TWC_static_mass','Masking_TWC_static_mass_'+str(int(particle))),
                           delay=200)

    if 'interpolation_processes' in mode:     
        
        # Calculated summed up profiles and integrated values (hydrometeor mass and process rates) following each cell:
        logging.debug('Start calculations for individual cells')
            
        logging.debug('loading tracks and Mask')

        tracks=pd.read_hdf(os.path.join(top_savedir_tracking,'Track.hdf'),'table')
        mask=iris.load_cube(os.path.join(top_savedir_tracking,'Watershed_TWC.nc'),'watershedding_output_mask')

        logging.debug('start loading processes form preprocessed file:')
        processes_lumped=iris.load(os.path.join(top_savedir_data,'Data_Processes_lumped.nc')).extract(constraint_tracking_period)        
##        
#        logging.debug(' start process rates from data files')
#
#        processes_lumped=load_function(filenames,'processes_lumped').extract(constraint_tracking_period)                 
#        logging.debug('loaded process rate')
#        
        for cube in processes_lumped:
            for coord in cube.coords():
                coord.var_name=coord.name()

        
        logging.debug('loaded processes')

        logging.debug('start loading auxillary variables')
        
        
        # load airmass (used to calculate total values from mixing ratios)
        airmass=iris.load_cube(os.path.join(top_savedir_data,'Data_WP.nc'),'airmass').extract(constraint_tracking_period) 
        for coord in airmass.coords():
            coord.var_name=coord.name()
        logging.debug('airmass loaded from file')

        
        if cell_selection is None:
            cell_selection=tracks['particle'].unique()

        # create empty list to store DataFrames with integrated values for each individual cell so that these can joined together in one at the end
        list_track=[]
        
        # loop over individual cells for analysis
        for particle in cell_selection:                        
            
            logging.debug( 'Start calculating Data for cell '+ str(int(particle)))
            savedir_cell=os.path.join(top_savedir_tracking,'cells',str(int(particle)))
            os.makedirs(savedir_cell,exist_ok=True)

            height_levels=np.arange(0,20001,1000.)
                            
            # Sum up process rates for individual cells
            processes_lumped_sum=iris.cube.CubeList()
            
            
            for cube in processes_lumped:
                for coord in cube.coords():
                    coord.var_name=coord.name()
                cube_sum=cube*airmass
                cube_sum.rename(cube.name())
                cube_sum=add_geopotential_height(cube_sum)
                processes_lumped_sum.append(cube_sum)
#            
#
            processes_lumped_cell_sum, processes_lumped_cell_integrated,track_processes_lumped_integrated=extract_cell_cubes_subset(cubelist_in=processes_lumped_sum,
                                                                                                                                 mask=mask,
                                                                                                                                 track=tracks,
                                                                                                                                 particle=particle,
                                                                                                                                 z_coord='model_level_number',
                                                                                                                                 height_levels=height_levels)
            
            iris.save(processes_lumped_cell_sum,os.path.join(savedir_cell,'Processes_lumped_profile_'+str(int(particle))+'.nc'))
            iris.save(processes_lumped_cell_integrated,os.path.join(savedir_cell,'Processes_lumped_integrated_'+str(int(particle))+'.nc'))
            track_processes_lumped_integrated.to_hdf(os.path.join(savedir_cell,'track_processes_lumped_integrated_'+str(int(particle))+'.h5'),'table')
            
            logging.debug( 'processes calculated and saved to '+ savedir_cell)
#
            list_track.append(track_processes_lumped_integrated)
        
        # concatenate DataFrames from intividual cells into one DataFrame containing all integrated values for all cells and save to file
        track_processes=pd.concat(list_track)
        track_processes.to_hdf(os.path.join(top_savedir_tracking,'track_processes_integrated.h5'),'table')

        logging.debug( 'Data calculated and saved to file for all cells')

    if 'plot_mask_processes' in mode:
        
        
        logging.debug('Start plotting tracks and masks for individual cells')

        savedir_data=os.path.join(top_savedir_data)
        savedir_tracking=os.path.join(top_savedir_tracking)
        plotdir=os.path.join(top_plotdir)
        
        features=pd.read_hdf(os.path.join(savedir_tracking,'Features.hdf'),'table')
        track_TWC=pd.read_hdf(os.path.join(savedir_tracking,'Track_TWC.hdf'),'table')
        track_w=pd.read_hdf(os.path.join(savedir_tracking,'Track_TWC.hdf'),'table')

        W_max=iris.load_cube(os.path.join(savedir_data,'Data_w_max.nc')).extract(constraint_tracking_period)    
        W_mid_max=iris.load_cube(os.path.join(savedir_data,'Data_w_mid_max.nc')).extract(constraint_tracking_period)    
        TWP=iris.load_cube(os.path.join(savedir_data,'Data_WP.nc'),'TWP').extract(constraint_tracking_period)  

        mask_TWC=iris.load_cube(os.path.join(savedir_tracking,'Watershed_TWC.nc'),'watershedding_output_mask')
        mask_w=iris.load_cube(os.path.join(savedir_tracking,'Watershed_w.nc'),'watershedding_output_mask')

        logging.debug('data loaded from savefile')    
                                 
#        if cell_selection is None:
#            cell_selection=np.unique(track['particle'])
        if plotting_period:
            track_TWC_time=track_TWC.loc[(track_TWC['time']>plotting_period[0]) & (track_TWC['time']<plotting_period[1])]
            track_w_time=track_TWC.loc[(track_w['time']>plotting_period[0]) & (track_w['time']<plotting_period[1])]

#        cell_selection=track_TWC_time.groupby('particle').agg({'ncells':'max'}).nlargest(columns='ncells',n=50).index.values
        if cell_selection is None:
            cell_selection=track_TWC_time['particle'].unique()

        cmap_TWP = colors.LinearSegmentedColormap.from_list("", ["white","lightblue",'lightsteelblue','cornflowerblue','royalblue','rebeccapurple','indigo'])
        levels_TWP=[0,0.5,1,1.5,2,3,4,5,10,15,20,30,40,50]

        cmap_w = colors.LinearSegmentedColormap.from_list("", ["forestgreen",'olivedrab','yellowgreen','khaki','gold'])
        levels_w=[0,1,2,5,10,15]

        track_processes=pd.read_hdf(os.path.join(savedir_tracking,'track_processes_integrated.h5'),'table')

        processes=list(processes_names.keys())
        colors_processes=processes_colors
 
        for particle in cell_selection:  
            logging.debug('Start plotting mask for cell'+ str(particle))
            
            plot_mask_cell_track_static_timeseries(particle=particle,
                                           track=track_processes, 
                                           cog=None,
                                           features=features,
                                           mask_total=mask_TWC,
                                           field_contour=W_mid_max, 
                                           label_field_contour='midlevel max w (m/s)',
                                           cmap_field_contour=cmap_w,
                                           vmin_field_contour=0,vmax_field_contour=15,levels_field_contour=levels_w,
                                           field_filled=TWP, 
                                           label_field_filled='TWP (kg/m^2)',
                                           cmap_field_filled=cmap_TWP,
                                           vmin_field_filled=0,vmax_field_filled=50,levels_field_filled=levels_TWP,
                                           track_variable=track_processes,variable=processes,
                                           variable_label=processes,variable_color=colors_processes,
                                           variable_ylabel='integrated process rate (kg/s)',variable_legend=True,
                                           width=10000,n_extend=2,
                                           name= 'Masking_TWC_static_processes_'+str(int(particle)), 
                                           plotdir=os.path.join(plotdir,'Masking_TWC_static_processes'),
                                           )
            make_animation(input_dir=os.path.join(plotdir,'Masking_TWC_static_processes','Masking_TWC_static_processes_'+str(int(particle))),
                           output=os.path.join(plotdir,'Animations_Masking_TWC_static_processes','Masking_TWC_static_processes_'+str(int(particle))),
                           delay=200)



    if 'plot_lifetime' in mode:
        savedir_tracking=os.path.join(top_savedir_tracking)
        track_TWC=pd.read_hdf(os.path.join(savedir_tracking,'Track_TWC.hdf'),'table')
        plot_dir=os.path.join(top_plotdir,'Lifetime')
        os.makedirs(plot_dir,exist_ok=True)
        plot_dir_cells=os.path.join(plot_dir,'Cells')
        os.makedirs(plot_dir_cells,exist_ok=True)
        fig4,ax4=plt.subplots(nrows=1,ncols=1,figsize=(15/2.54,15/2.54))
        Track_particle=track_TWC.groupby('particle')
        for particle, Track in Track_particle:
            fig5,ax5=plt.subplots(nrows=1,ncols=1,figsize=(15/2.54,15/2.54))
            ax4.plot(Track['time_cell'].minutes,Track['ncells'],ls='-',lw=0.3,legend=False)
            ax5.plot(Track['time_cell'].minutes,Track['ncells'],ls='-',lw=1,legend=False,color='navy')
            ax5.set_xlabel('cell lifetime (min)')
            ax5.set_ylabel('ncells')
            ax5.set_xlim([0,2*1e9*3600])
            ax4.set_ylim([0,max(10,1.1*Track['ncells'].max())])
            ax5.set_xticks(1e9*3600*np.arange(0,2,0.25))
#                    ax5.xaxis.set_major_locator(hours)
#                    ax5.xaxis.set_major_formatter(formatter)
            ax5.set_title('cell: '+str(particle))
            filename=os.path.join(plot_dir_cells,'lifetime'+'_'+str(int(particle))+'.png')
            fig5.savefig(filename,dpi=300)
            plt.close(fig5)
        ax4.set_xlabel('cell lifetime (min)')
        ax4.set_ylabel('ncells')
        ax4.set_xlim([0,2*1e9*3600])
        ax4.set_xticks(1e9*3600*np.arange(0,2,0.25))
#                ax4.xaxis.set_major_locator(hours)
#                ax4.xaxis.set_major_formatter(formatter)

        ax4.set_ylim([0,5000])
        plt.close(fig4)
        filename=os.path.join(plot_dir,'lifetimes'+'_'+'all'+'.png')
        fig4.savefig(filename,dpi=300)
        logging.info('lifetimes ncell plotted')

        make_animation(input_dir=plot_dir_cells,
                       output=os.path.join(plot_dir,'Animation_Lifetime_cells'),delay=200)

