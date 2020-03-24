import logging
import iris
import os
import pandas as pd

from .plot import color_TWP,color_w
from .plot import plot_mask_cell_track_static_timeseries
from .core import add_geopotential_height
import numpy as np

def interpolation_plot_mask_mass_w_TWC(cell=None,savedir_tracking=None,savedir_data=None,plotdir=None,cell_selection=None,constraint_tracking_period=None):
    from html_animate import make_animation
    logging.debug('Start plotting tracks and masks for individual cells')
        
    features=pd.read_hdf(os.path.join(savedir_tracking,'Features.h5'),'table')
    track=pd.read_hdf(os.path.join(savedir_tracking,'Track.h5'),'table')

    W_mid_max=iris.load_cube(os.path.join(savedir_data,'Data_w_mid_max.nc')).extract(constraint_tracking_period)    
    TWP=iris.load_cube(os.path.join(savedir_data,'Data_WP.nc'),'TWP').extract(constraint_tracking_period)  

    mask=iris.load_cube(os.path.join(savedir_tracking,'Mask_Segmentation_w_TWC.nc'),'segmentation_mask')

    cmap_TWP,levels_TWP,bounds_TWP,norm_TWP=color_TWP()
    cmap_w,levels_w,bounds_w,norm_w=color_w()

    logging.debug('data loaded from savefile')    

    if cell_selection is None:
        cell_selection=track['cell'].dropna().unique()

    for cell in cell_selection:                        
        logging.debug('Start plotting mask for cell'+ str(cell))
        
        savedir_cell=os.path.join(savedir_tracking,'cells',str(int(cell)))
        track_mass=pd.read_hdf(os.path.join(savedir_cell,'track_mass_integrated_w_TWC_'+str(int(cell))+'.h5'),'table')

        plot_mask_cell_track_static_timeseries(cell=cell,
                                       track=track_mass, 
                                       cog=None,
                                       features=features,
                                       mask_total=mask,
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
                                       name= 'Masking_TWC_static_mass_'+str(int(cell)), 
                                       plotdir=os.path.join(plotdir,'Masking_TWC_static_mass'),
                                       )
        make_animation(input_dir=os.path.join(plotdir,'Masking_TWC_static_mass','Masking_TWC_static_mass_'+str(int(cell))),
                       output=os.path.join(plotdir,'Animations_Masking_TWC_static_mass','Masking_TWC_static_mass_'+str(int(cell))),
                       delay=200,make_gif=False)

def interpolation_plot_processes_w_TWC(savedir_tracking=None,savedir_data=None,plotdir=None,cell_selection=None,constraint_tracking_period=None,plotting_period=None):
    from mpdiag import processes_colors
    from html_animate import make_animation
    logging.debug('Start plotting tracks and masks for individual cells')
    features=pd.read_hdf(os.path.join(savedir_tracking,'Features.h5'),'table')
    track=pd.read_hdf(os.path.join(savedir_tracking,'Track.h5'),'table')
    
    W_mid_max=iris.load_cube(os.path.join(savedir_data,'Data_w_mid_max.nc')).extract(constraint_tracking_period)    
    TWP=iris.load_cube(os.path.join(savedir_data,'Data_WP.nc'),'TWP').extract(constraint_tracking_period)  
    mask=iris.load_cube(os.path.join(savedir_tracking,'Mask_Segmentation_w_TWC.nc'),'segmentation_mask')
    logging.debug('data loaded from savefile')    
    cmap_TWP,levels_TWP,bounds_TWP,norm_TWP=color_TWP()
    cmap_w,levels_w,bounds_w,norm_w=color_w()
    if plotting_period:
            track_time=track.loc[(track['time']>plotting_period[0]) & (track['time']<plotting_period[1])]

    logging.debug('loading tracks and Mask')

    tracks=pd.read_hdf(os.path.join(savedir_tracking,'Track.h5'),'table')
    mask=iris.load_cube(os.path.join(savedir_tracking,'Mask_Segmentation_w_TWC.nc'),'segmentation_mask')

    processes_file=os.path.join(savedir_data,'Data_Processes_lumped.nc')
    
    if os.path.exists(processes_file):       
        logging.debug('start loading processes form preprocessed file:')
        processes_lumped=iris.load(processes_file).extract(constraint_tracking_period)        
#        else:#        
#            logging.debug(' start process rates from data files')
#            processes_lumped=load_function(filenames,'processes_lumped').extract(constraint_tracking_period)                 

    for cube in processes_lumped:
        for coord in cube.coords():
            coord.var_name=coord.name()
            coord.coord_system=None

    logging.debug('loaded processes')

    logging.debug('start loading auxillary variables')
    
    
    # load airmass (used to calculate total values from mixing ratios)
    airmass=iris.load_cube(os.path.join(savedir_data,'Data_Airmass.nc'),'airmass').extract(constraint_tracking_period) 
    for coord in airmass.coords():
        coord.var_name=coord.name()
    logging.debug('airmass loaded from file')


    height_levels=np.arange(0,20001,1000.)
                    
    # Sum up process rates for individual cells
    processes_lumped_sum=iris.cube.CubeList()
    
    for cube in processes_lumped:
        cube_sum=cube*airmass
        cube_sum.rename(cube.name())
        cube_sum=add_geopotential_height(cube_sum)
        processes_lumped_sum.append(cube_sum)
#            
    if cell_selection is None:
        cell_selection=track_time['cell'].dropna().unique()
        
    logging.debug(f'cell_selection: {cell_selection}')

    # loop over individual cells for analysis
    for cell in cell_selection:                        
        logging.debug( 'Start calculating Data for cell '+ str(int(cell)))
        savedir_cell=os.path.join(savedir_tracking,'cells',str(int(cell)))
        os.makedirs(savedir_cell,exist_ok=True)

#
        processes_lumped_cell_sum, processes_lumped_cell_integrated,track_processes_lumped_integrated=extract_cell_cubes_subset(cubelist_in=processes_lumped_sum,
                                                                                                                             mask=mask,
                                                                                                                             track=tracks,
                                                                                                                             cell=cell,
                                                                                                                             z_coord='model_level_number',
                                                                                                                             height_levels=height_levels)
        
        iris.save(processes_lumped_cell_sum,os.path.join(savedir_cell,'Processes_lumped_profile_w_TWC.nc'),zlib=True,complevel=4)
        iris.save(processes_lumped_cell_integrated,os.path.join(savedir_cell,'Processes_lumped_integrated_w_TWC.nc'),zlib=True,complevel=4)
        track_processes_lumped_integrated.to_hdf(os.path.join(savedir_cell,'track_processes_lumped_integrated_w_TWC.h5'),'table')
        
        logging.debug( 'processes calculated and saved to '+ savedir_cell)

        savedir_cell=os.path.join(savedir_tracking,'cells',str(int(cell)))

        colors_processes, processes_names = processes_colors(microphysics_scheme='RAMS', colors_processes='lumped')
        processes=list(processes_names.keys())
        track_processes=pd.read_hdf(os.path.join(savedir_cell,'track_processes_lumped_integrated_w_TWC_'+str(int(cell))+'.h5'),'table')

        plot_mask_cell_track_static_timeseries(cell=cell,
                                       track=track_processes, 
                                       cog=None,
                                       features=features,
                                       mask_total=mask,
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
                                       name= 'Masking_w_TWC_static_processes_'+str(int(cell)), 
                                       plotdir=os.path.join(plotdir,'Masking_w_TWC_static_processes'),
                                       )
        make_animation(input_dir=os.path.join(plotdir,'Masking_w_TWC_static_processes','Masking_w_TWC_static_processes_'+str(int(cell))),
                       output=os.path.join(plotdir,'Animations_Masking_w_TWC_static_processes','Masking_w_TWC_static_processes_'+str(int(cell))),
                       delay=200,make_gif=False)

def interpolation_plot_processes_TWC(savedir_tracking=None,savedir_data=None,plotdir=None,cell_selection=None,constraint_tracking_period=None,plotting_period=None):
    from mpdiag import processes_colors
    from html_animate import make_animation
    logging.debug('Start plotting tracks and masks for individual cells')
    features=pd.read_hdf(os.path.join(savedir_tracking,'Features.h5'),'table')
    track=pd.read_hdf(os.path.join(savedir_tracking,'Track.h5'),'table')
    
    W_mid_max=iris.load_cube(os.path.join(savedir_data,'Data_w_mid_max.nc')).extract(constraint_tracking_period)    
    TWP=iris.load_cube(os.path.join(savedir_data,'Data_WP.nc'),'TWP').extract(constraint_tracking_period)  
    mask=iris.load_cube(os.path.join(savedir_tracking,'Mask_Segmentation_w_TWC.nc'),'segmentation_mask')
    logging.debug('data loaded from savefile')    
    cmap_TWP,levels_TWP,bounds_TWP,norm_TWP=color_TWP()
    cmap_w,levels_w,bounds_w,norm_w=color_w()
    if plotting_period:
            track_time=track.loc[(track['time']>plotting_period[0]) & (track['time']<plotting_period[1])]

    logging.debug('loading tracks and Mask')

    tracks=pd.read_hdf(os.path.join(savedir_tracking,'Track.h5'),'table')
    mask=iris.load_cube(os.path.join(savedir_tracking,'Mask_Segmentation_TWC.nc'),'segmentation_mask')

    processes_file=os.path.join(savedir_data,'Data_Processes_lumped.nc')
    
    if os.path.exists(processes_file):       
        logging.debug('start loading processes form preprocessed file:')
        processes_lumped=iris.load(processes_file).extract(constraint_tracking_period)        
#        else:#        
#            logging.debug(' start process rates from data files')
#            processes_lumped=load_function(filenames,'processes_lumped').extract(constraint_tracking_period)                 

    for cube in processes_lumped:
        for coord in cube.coords():
            coord.var_name=coord.name()
            coord.coord_system=None

    logging.debug('loaded processes')

    logging.debug('start loading auxillary variables')
    
    
    # load airmass (used to calculate total values from mixing ratios)
    airmass=iris.load_cube(os.path.join(savedir_data,'Data_Airmass.nc'),'airmass').extract(constraint_tracking_period) 
    for coord in airmass.coords():
        coord.var_name=coord.name()
    logging.debug('airmass loaded from file')


    height_levels=np.arange(0,20001,1000.)
                    
    # Sum up process rates for individual cells
    processes_lumped_sum=iris.cube.CubeList()
    
    for cube in processes_lumped:
        cube_sum=cube*airmass
        cube_sum.rename(cube.name())
        cube_sum=add_geopotential_height(cube_sum)
        processes_lumped_sum.append(cube_sum)
#            
    if cell_selection is None:
        cell_selection=track_time['cell'].dropna().unique()

    logging.debug(f'cell_selection: {cell_selection}')
        
    # loop over individual cells for analysis
    for cell in cell_selection:                        
        logging.debug( 'Start calculating Data for cell '+ str(int(cell)))
        savedir_cell=os.path.join(savedir_tracking,'cells',str(int(cell)))
        os.makedirs(savedir_cell,exist_ok=True)

#
        processes_lumped_cell_sum, processes_lumped_cell_integrated,track_processes_lumped_integrated=extract_cell_cubes_subset(cubelist_in=processes_lumped_sum,
                                                                                                                             mask=mask,
                                                                                                                             track=tracks,
                                                                                                                             cell=cell,
                                                                                                                             z_coord='model_level_number',
                                                                                                                             height_levels=height_levels)
        
        iris.save(processes_lumped_cell_sum,os.path.join(savedir_cell,'Processes_lumped_profile_TWC.nc'),zlib=True,complevel=4)
        iris.save(processes_lumped_cell_integrated,os.path.join(savedir_cell,'Processes_lumped_integrated.nc'),zlib=True,complevel=4)
        track_processes_lumped_integrated.to_hdf(os.path.join(savedir_cell,'track_processes_lumped_integrated_TWC.h5'),'table')
        
        logging.debug( 'processes calculated and saved to '+ savedir_cell)

        savedir_cell=os.path.join(savedir_tracking,'cells',str(int(cell)))

        colors_processes, processes_names = processes_colors(microphysics_scheme='RAMS', colors_processes='lumped')
        processes=list(processes_names.keys())
        track_processes=pd.read_hdf(os.path.join(savedir_cell,'track_processes_lumped_integrated_TWC_'+str(int(cell))+'.h5'),'table')

        plot_mask_cell_track_static_timeseries(cell=cell,
                                       track=track_processes, 
                                       cog=None,
                                       features=features,
                                       mask_total=mask,
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
                                       name= 'Masking_w_TWC_static_processes_'+str(int(cell)), 
                                       plotdir=os.path.join(plotdir,'Masking_w_TWC_static_processes'),
                                       )
        make_animation(input_dir=os.path.join(plotdir,'Masking_w_TWC_static_processes','Masking_w_TWC_static_processes_'+str(int(cell))),
                       output=os.path.join(plotdir,'Animations_Masking_w_TWC_static_processes','Masking_w_TWC_static_processes_'+str(int(cell))),
                       delay=200,make_gif=False)

def interpolation_plot_mask_mass_TWC(cell=None,savedir_tracking=None,savedir_data=None,plotdir=None,cell_selection=None,constraint_tracking_period=None):
    from html_animate import make_animation
    logging.debug('Start plotting tracks and masks for individual cells')
        
    features=pd.read_hdf(os.path.join(savedir_tracking,'Features.h5'),'table')
    track=pd.read_hdf(os.path.join(savedir_tracking,'Track.h5'),'table')

    W_mid_max=iris.load_cube(os.path.join(savedir_data,'Data_w_mid_max.nc')).extract(constraint_tracking_period)    
    TWP=iris.load_cube(os.path.join(savedir_data,'Data_WP.nc'),'TWP').extract(constraint_tracking_period)  

    mask=iris.load_cube(os.path.join(savedir_tracking,'Mask_Segmentation_TWC.nc'),'segmentation_mask')

    cmap_TWP,levels_TWP,bounds_TWP,norm_TWP=color_TWP()
    cmap_w,levels_w,bounds_w,norm_w=color_w()

    logging.debug('data loaded from savefile')    

    if cell_selection is None:
        cell_selection=track['cell'].dropna().unique()

    for cell in cell_selection:                        
        logging.debug('Start plotting mask for cell'+ str(cell))
        
        savedir_cell=os.path.join(savedir_tracking,'cells',str(int(cell)))
        track_mass=pd.read_hdf(os.path.join(savedir_cell,'track_mass_integrated_'+str(int(cell))+'.h5'),'table')

        plot_mask_cell_track_static_timeseries(cell=cell,
                                       track=track_mass, 
                                       cog=None,
                                       features=features,
                                       mask_total=mask,
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
                                       name= 'Masking_TWC_static_mass_'+str(int(cell)), 
                                       plotdir=os.path.join(plotdir,'Masking_TWC_static_mass'),
                                       )
        make_animation(input_dir=os.path.join(plotdir,'Masking_TWC_static_mass','Masking_TWC_static_mass_'+str(int(cell))),
                       output=os.path.join(plotdir,'Animations_Masking_TWC_static_mass','Masking_TWC_static_mass_'+str(int(cell))),
                       delay=200,make_gif=False)
