import iris
import logging
import os
import matplotlib.pyplot as plt
import pandas as pd

def plot_mask_processes_TWC(savedir_tracking=None,savedir_data=None,plotdir=None,cell_selection=None,constraint_tracking_period=None,plotting_period=None):
    from mpdiag import processes_colors
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
    if plotting_period:
            track_time=track.loc[(track['time']>plotting_period[0]) & (track['time']<plotting_period[1])]
    
    if cell_selection is None:
        cell_selection=np.unique(track_time['cell'])
        
    logging.debug(f'cell_selection: {cell_selection}')

    # loop over individual cells for analysis
    for cell in cell_selection:                        
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
                                           name= 'Masking_TWC_static_processes_'+str(int(cell)), 
                                           plotdir=os.path.join(plotdir,'Masking_TWC_static_processes'),
                                           )
        make_animation(input_dir=os.path.join(plotdir,'Masking_TWC_static_processes','Masking_TWC_static_processes_'+str(int(cell))),
                       output=os.path.join(plotdir,'Animations_Masking_TWC_static_processes','Masking_TWC_static_processes_'+str(int(cell))),
                       delay=200,make_gif=False)

def plot_mask_processes_w_TWC(savedir_tracking=None,savedir_data=None,plotdir=None,cell_selection=None,constraint_tracking_period=None,plotting_period=None):
    from mpdiag import processes_colors
    from html_animate import make_animation
    logging.debug('Start plotting tracks and masks for individual cells')
    features=pd.read_hdf(os.path.join(savedir_tracking,'Features.h5'),'table')
    track=pd.read_hdf(os.path.join(savedir_tracking,'Track_w_TWC.h5'),'table')
    W_mid_max=iris.load_cube(os.path.join(savedir_data,'Data_w_mid_max.nc')).extract(constraint_tracking_period)    
    TWP=iris.load_cube(os.path.join(savedir_data,'Data_WP.nc'),'TWP').extract(constraint_tracking_period)  
    mask=iris.load_cube(os.path.join(savedir_tracking,'Mask_Segmentation_w_TWC.nc'),'segmentation_mask')
    logging.debug('data loaded from savefile')    
    cmap_TWP,levels_TWP,bounds_TWP,norm_TWP=color_TWP()
    cmap_w,levels_w,bounds_w,norm_w=color_w()
    if plotting_period:
            track_time=track.loc[(track['time']>plotting_period[0]) & (track['time']<plotting_period[1])]
    
    if cell_selection is None:
        cell_selection=np.unique(track_time['cell'])
        
    logging.debug(f'cell_selection: {cell_selection}')

    # loop over individual cells for analysis
    for cell in cell_selection:                        
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

def plots_maps_loop(savedir_data=None,savedir_tracking=None,plotdir=None,plotting_parameters=None,constraint_tracking_period=None):
    from html_animate import make_animation

    track=pd.read_hdf(os.path.join(savedir_tracking,'Track.h5'),'table')
    features=pd.read_hdf(os.path.join(savedir_tracking,'Features.h5'),'table')

    W_max=iris.load_cube(os.path.join(savedir_data,'Data_w_max.nc')).extract(constraint_tracking_period)
    TWP=iris.load_cube(os.path.join(savedir_data,'Data_WP.nc'),'TWP')    

    Mask_TWC=iris.load_cube(os.path.join(savedir_tracking,'Mask_Segmentation_TWC.nc'),'segmentation_mask')
    Mask_w=iris.load_cube(os.path.join(savedir_tracking,'Mask_Segmentation_w.nc'),'segmentation_mask')
    Mask_w_TWC=iris.load_cube(os.path.join(savedir_tracking,'Mask_Segmentation_w_TWC.nc'),'segmentation_mask')

    logging.debug('data loaded from savefile')

    axis_extent=plotting_parameters['axis_extent']

    plot_tracks_mask_field_loop(track=track,field=W_max,mask=Mask_w,features=features,
                                plot_dir=os.path.join(plotdir,'w_max_Mask_w'),name='w_max_Mask_w',
                                 axis_extent=axis_extent,#figsize=figsize,orientation_colorbar='horizontal',pad_colorbar=0.2,
                                 vmin=0,vmax=30,
                                 plot_outline=True,plot_marker=True,marker_track=',',plot_number=True)

    make_animation(input_dir=os.path.join(plotdir,'w_max_Mask_w'),
                   output=os.path.join(plotdir,'Animation_w_max_Mask_w','w_max_Mask_w'),
                                          delay=200,make_gif=False)
    logging.debug('W_max with w mask plotted')
    
    plot_tracks_mask_field_loop(track=track,field=W_max,mask=Mask_w_TWC,features=features,
                                plot_dir=os.path.join(plotdir,'w_max_Mask_w_TWC'),name='w_max_Mask_w_TWC',
                                axis_extent=axis_extent,#figsize=figsize,orientation_colorbar='horizontal',pad_colorbar=0.2,
                                vmin=0,vmax=30,
                                plot_outline=True,plot_marker=True,marker_track='x',plot_number=True,plot_features=True)
    logging.debug('W_max with w_TWC mask plotted')
    
    make_animation(input_dir=os.path.join(plotdir,'w_max_Mask_w_TWC'),
                   output=os.path.join(plotdir,'Animation_w_max_Mask_w_TWC','w_max_Mask_w_TWC'),
                                          delay=200,make_gif=False)

    plot_tracks_mask_field_loop(track=track,field=TWP,mask=Mask_TWC,features=features,
                                plot_dir=os.path.join(plotdir,'TWP_Mask_TWC'),name='TWP_Mask_TWC',
                                axis_extent=axis_extent,#figsize=figsize,orientation_colorbar='horizontal',pad_colorbar=0.2,
                                vmin=0,vmax=30,
                                plot_outline=True,plot_marker=True,marker_track='x',plot_number=True,plot_features=True)
    logging.debug('TWP with w_TWC plotted')

    make_animation(input_dir=os.path.join(plotdir,'TWP_Mask_TWC'),
                   output=os.path.join(plotdir,'Animation_TWP_Mask_TWC','TWP_Mask_TWC'),
                                          delay=200,make_gif=False)

def plot_lifetime(top_savedir_tracking=None,top_plotdir=None):
    from html_animate import make_animation

    savedir_tracking=os.path.join(top_savedir_tracking)
    track_TWC=pd.read_hdf(os.path.join(savedir_tracking,'Track_TWC.h5'),'table')
    plot_dir=os.path.join(top_plotdir,'Lifetime')
    os.makedirs(plot_dir,exist_ok=True)
    plot_dir_cells=os.path.join(plot_dir,'Cells')
    os.makedirs(plot_dir_cells,exist_ok=True)
    fig4,ax4=plt.subplots(nrows=1,ncols=1,figsize=(15/2.54,15/2.54))
    Track_cell=track_TWC.groupby('cell')
    Lifetimes=(Track_cell['time_cell'].max()/pd.Timedelta(minutes=1)).rename("lifetime")
    Lifetimes.to_csv(os.path.join(plot_dir,'Lifetimes.txt'))
    
    for cell, Track in Track_cell:
        fig5,ax5=plt.subplots(nrows=1,ncols=1,figsize=(15/2.54,15/2.54))
        ax4.plot((Track['time_cell'].dt.total_seconds()/60.).values,Track['ncells'].values,ls='-',lw=0.3)
        ax5.plot((Track['time_cell'].dt.total_seconds()/ 60.).values,Track['ncells'].values,ls='-',lw=1,color='navy')
        ax5.set_xlabel('cell lifetime (min)')
        ax5.set_ylabel('ncells')
        ax5.set_xlim([0,2*1e9*3600])
        ax4.set_ylim([0,max(10,1.1*Track['ncells'].max())])
        ax5.set_xticks(1e9*3600*np.arange(0,2,0.25))
#                    ax5.xaxis.set_major_locator(hours)
#                    ax5.xaxis.set_major_formatter(formatter)
        ax5.set_title('cell: '+str(cell))
        filename=os.path.join(plot_dir_cells,'lifetime'+'_'+str(int(cell))+'.png')
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
                    output=os.path.join(plot_dir,'Animation_Lifetime_cells'),
                    delay=200,make_gif=False)

def plot_profiles_w(cell=None,savedir=None,plotdir=None):
    Track=pd.read_hdf(os.path.join(savedir,'Track.h5'),'table')
    for cell in np.unique(Track['cell']):
        fig6,ax6=plt.subplots(nrows=1,ncols=1)
        savefile_max=os.path.join(savedir,'W','maximum','W_maximum' + '_'+str(int(cell))+'.nc')
        savefile_min=os.path.join(savedir,'W','minimum','W_minimum' +'_'+str(int(cell))+'.nc')
        savefile_mean=os.path.join(savedir,'W','mean','W_mean' + '_'+str(int(cell))+'.nc')
        savefile_median=os.path.join(savedir,'W','median','W_median' + '_'+str(int(cell))+'.nc')
        os.makedirs(os.path.join(plotdir,'Profiles_w_individual'),exist_ok=True)
        os.makedirs(os.path.join(plotdir,'Profiles_w'),exist_ok=True)
        w_max=iris.load_cube(savefile_max)
        w_min=iris.load_cube(savefile_min)
        w_mean=iris.load_cube(savefile_mean)
        w_median=iris.load_cube(savefile_median)
        
        for i,time in enumerate(w_max.coord('time').points):
            logging.debug('plotting '+str(cell)+ ' ' +str(i))                        
            plot_max_ind=ax6.plot(w_max[i].data,w_max[i].coord('geopotential_height').points/1000,color='r',marker='.',markersize=2,linestyle='None',label='max')    
            plot_min_ind=ax6.plot(w_min[i].data,w_min[i].coord('geopotential_height').points/1000,color='b',marker='.',markersize=2,linestyle='None',label='min')    
            plot_mean_ind=ax6.plot(w_mean[i].data,w_mean[i].coord('geopotential_height').points/1000,color='k',marker='.',markersize=2,linestyle='None',label='mean')
            plot_median_ind=ax6.plot(w_median[i].data,w_median[i].coord('geopotential_height').points/1000,color='grey',marker='.',markersize=2,linestyle='None',label='median')
        
        ax6.set_xlim([0,30])
        ax6.set_ylim([0,15])
        ax6.set_xlabel('w (m/s)')
        ax6.set_ylabel('z (km)')
                
        ax6.legend(handles=[plot_min_ind,plot_max_ind,plot_mean_ind,plot_median_ind],loc=1)
        
        fig6.savefig(os.path.join(plotdir,'Profiles_w_individual','Profile_w_'+str(int(cell))+'.png'),dpi=600)
        plt.close(fig6)
        logging.debug('profiles plotted  for cell '+str(cell))

def plot_profile_CDNC(savedir=None,plotdir=None):
    Track=pd.read_hdf(os.path.join(savedir,'Track.h5'),'table')
    for cell in np.unique(Track['cell']):
        fig6,ax6=plt.subplots(nrows=1,ncols=1)
        
        savefile_max=os.path.join(savedir,'CDNC','maximum','CDNC_maximum' + '_'+str(int(cell))+'.nc')
        savefile_min=os.path.join(savedir,'CDNC','minimum','CDNC_minimum' +'_'+str(int(cell))+'.nc')
        savefile_mean=os.path.join(savedir,'CDNC','mean','CDNC_mean' + '_'+str(int(cell))+'.nc')
        savefile_median=os.path.join(savedir,'CDNC','median','CDNC_median' + '_'+str(int(cell))+'.nc')

        os.makedirs(os.path.join(plotdir,'Profiles_CDNC_individual'),exist_ok=True)
        os.makedirs(os.path.join(plotdir,'Profiles_w'),exist_ok=True)

        CDNC_max=iris.load_cube(savefile_max)
        CDNC_min=iris.load_cube(savefile_min)
        CDNC_mean=iris.load_cube(savefile_mean)
        CDNC_median=iris.load_cube(savefile_median)

        for i,time in enumerate(CDNC_max.coord('time').points):
            plot_max_ind=ax6.plot(CDNC_max[i].data/1e6,CDNC_max[i].coord('geopotential_height').points/1000,color='r',marker='.',markersize=2,linestyle='None',label='max')    
            plot_min_ind=ax6.plot(CDNC_min[i].data/1e6,CDNC_min[i].coord('geopotential_height').points/1000,color='b',marker='.',markersize=2,linestyle='None',label='min')    
            plot_mean_ind=ax6.plot(CDNC_mean[i].data/1e6,CDNC_mean[i].coord('geopotential_height').points/1000,color='k',marker='.',markersize=2,linestyle='None',label='mean')
            plot_median_ind=ax6.plot(CDNC_median[i].data/1e6,CDNC_median[i].coord('geopotential_height').points/1000,color='grey',marker='.',markersize=2,linestyle='None',label='median')

        ax6.set_xlim([0,5000])
        ax6.set_ylim([0,15])
        
        ax6.set_xlabel('CDNC (#/cm^2)')
        ax6.set_ylabel('z (km)')
        
        ax6.legend(handles=[plot_min_ind,plot_max_ind,plot_mean_ind,plot_median_ind],loc=1)


        fig6.savefig(os.path.join(plotdir,'Profiles_CDNC_individual','Profile_CDNC_'+str(int(cell))+'.png'),dpi=600)
        plt.close(fig6)
        logging.debug('profiles plotted for cell '+str(cell))

def color_TWP():
    cmap_TWP = colors.LinearSegmentedColormap.from_list("", ["white","lightblue",'lightsteelblue','cornflowerblue','royalblue','rebeccapurple','indigo'])
    bounds_TWP = [2, 3, 7, 9, 15]
    norm_TWP = colors.BoundaryNorm(bounds_TWP, cmap_TWP.N)
    levels_TWP=[0,0.5,1,1.5,2,3,4,5,10,15,20,30,40,50]
    return cmap_TWP,bounds_TWP,norm_TWP,levels_TWP

def color_w():
    cmap_w = colors.LinearSegmentedColormap.from_list("", ["forestgreen",'olivedrab','yellowgreen','khaki','gold', ])
    levels_w=[1,3,5,10,15]
    bounds_w=[0.5, 2, 4, 7, 12, 18]
    norm_w = colors.BoundaryNorm(bounds_w, cmap_w.N)
    return cmap_w,levels_w,bounds_w,norm_w

def angles_3D():

    ele=10.
    azim=300.
    return ele, azim


def plot_mask_mass_TWC(cell=None,savedir_tracking=None,savedir_data=None,plotdir=None,cell_selection=None,constraint_tracking_period=None,plotting_period=None):
    from html_animate import make_animation
    logging.debug('Start plotting tracks and masks for individual cells')
        
    features=pd.read_hdf(os.path.join(savedir_tracking,'Features.h5'),'table')
    track=pd.read_hdf(os.path.join(savedir_tracking,'Track_TWC.h5'),'table')

    W_mid_max=iris.load_cube(os.path.join(savedir_data,'Data_w_mid_max.nc')).extract(constraint_tracking_period)    
    TWP=iris.load_cube(os.path.join(savedir_data,'Data_WP.nc'),'TWP').extract(constraint_tracking_period)  

    mask=iris.load_cube(os.path.join(savedir_tracking,'Mask_Segmentation_TWC.nc'),'segmentation_mask')

    logging.debug('data loaded from savefile')
    cmap_TWP,levels_TWP,bounds_TWP,norm_TWP=color_TWP()
    cmap_w,levels_w,bounds_w,norm_w=color_w()
    if plotting_period:
        track_time=track.loc[(track['time']>plotting_period[0]) & (track['time']<plotting_period[1])]

    if cell_selection is None:
        cell_selection=track_time['cell'].dropna().unique()

    for cell in cell_selection:                        
        
        logging.debug('Start plotting mask for cell'+ str(cell))
                    
        savedir_cell=os.path.join(savedir_tracking,'cells',str(int(cell)))
        track_mass=pd.read_hdf(os.path.join(savedir_cell,'track_mass_integrated_TWC'+str(int(cell))+'.h5'),'table')
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
                                               plotdir=os.path.join(plotdir,'Masking_TWC_static_mass'))
    make_animation(input_dir=os.path.join(plotdir,'Masking_TWC_static_mass','Masking_TWC_static_mass_'+str(int(cell))),
                   output=os.path.join(plotdir,'Animations_Masking_TWC_static_mass','Masking_TWC_static_mass_'+str(int(cell))),
                   delay=200,make_gif=False)

def plot_mask_mass_w_TWC(cell=None,savedir_tracking=None,savedir_data=None,plotdir=None,cell_selection=None,constraint_tracking_period=None,plotting_period=None):
    from html_animate import make_animation
    logging.debug('Start plotting tracks and masks for individual cells')
        
    features=pd.read_hdf(os.path.join(savedir_tracking,'Features.h5'),'table')
    track=pd.read_hdf(os.path.join(savedir_tracking,'Track_TWC.h5'),'table')

    W_mid_max=iris.load_cube(os.path.join(savedir_data,'Data_w_mid_max.nc')).extract(constraint_tracking_period)    
    TWP=iris.load_cube(os.path.join(savedir_data,'Data_WP.nc'),'TWP').extract(constraint_tracking_period)  

    mask=iris.load_cube(os.path.join(savedir_tracking,'Mask_Segmentation_w_TWC.nc'),'segmentation_mask')

    logging.debug('data loaded from savefile')
    cmap_TWP,levels_TWP,bounds_TWP,norm_TWP=color_TWP()
    cmap_w,levels_w,bounds_w,norm_w=color_w()
    if plotting_period:
        track_time=track.loc[(track['time']>plotting_period[0]) & (track['time']<plotting_period[1])]

    if cell_selection is None:
        cell_selection=track_time['cell'].dropna().unique()

    for cell in cell_selection:                        
        
        logging.debug('Start plotting mask for cell'+ str(cell))
                    
        savedir_cell=os.path.join(savedir_tracking,'cells',str(int(cell)))
        track_mass=pd.read_hdf(os.path.join(savedir_cell,'track_mass_integrated_w_TWC'+str(int(cell))+'.h5'),'table')
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
                                               plotdir=os.path.join(plotdir,'Masking_TWC_static_mass'))
    make_animation(input_dir=os.path.join(plotdir,'Masking_w_TWC_static_mass','Masking_w_TWC_static_mass_'+str(int(cell))),
                   output=os.path.join(plotdir,'Animations_Masking_w_TWC_static_mass','Masking_w_TWC_static_mass_'+str(int(cell))),
                   delay=200,make_gif=False)

def plot_mask_w_TWC(savedir_tracking=None,savedir_data=None,plotdir=None,cell_selection=None,constraint_tracking_period=None,plotting_period=None):
    from html_animate import make_animation
    features=pd.read_hdf(os.path.join(savedir_tracking,'Features.h5'),'table')
    TWP=iris.load_cube(os.path.join(savedir_data,'Data_WP.nc'),'TWP').extract(constraint_tracking_period)
    W_mid_max=iris.load_cube(os.path.join(savedir_data,'Data_w_mid_max.nc')).extract(constraint_tracking_period)
    track=pd.read_hdf(os.path.join(savedir_tracking,'Track_w_TWC.h5'),'table')
    mask=iris.load_cube(os.path.join(savedir_tracking,'Mask_Segmentation_w_TWC.nc'),'segmentation_mask')
    if plotting_period:
        track_time=track.loc[(track['time']>plotting_period[0]) & (track['time']<plotting_period[1])]
    cmap_TWP,levels_TWP,bounds_TWP,norm_TWP=color_TWP()
    cmap_w,levels_w,bounds_w,norm_w=color_w()
    if cell_selection is None:
        cell_selection=track_time['cell'].dropna().unique()

    for cell in cell_selection:  
        logging.debug('Start plotting mask for cell'+ str(int(cell)))
        # Make plots centred around and moving with the tracked cell (10km each side)
        name='Masking_w_TWC_static_ncells'
        plot_mask_cell_track_static_timeseries(cell=cell,
                                       track=track, 
                                       cog=None,
                                       features=features,
                                       mask_total=mask,
                                       field_contour=W_mid_max, 
                                       label_field_contour='midlevel max w (m/s)',
                                       cmap_field_contour=cmap_w,norm_field_contour=norm_w,
                                       vmin_field_contour=0,vmax_field_contour=15,levels_field_contour=levels_w,
                                       field_filled=TWP, 
                                       label_field_filled='TWP (kg/m^2)',
                                       cmap_field_filled=cmap_TWP,
                                       vmin_field_filled=0,vmax_field_filled=50,levels_field_filled=levels_TWP,
                                       track_variable=track,variable='ncells',variable_label='number of cells (watershedding)',
                                       width=10000,n_extend=2,
                                       name= name+'_'+str(int(cell)), 
                                       plotdir=os.path.join(plotdir,name),
                                       )
        make_animation(input_dir=os.path.join(plotdir,name,name+'_'+str(int(cell))),
                       output=os.path.join(plotdir,'Animations'+'_'+name,name+'_'+str(int(cell))),
                       delay=200,make_gif=False)

def plot_mask_w(savedir_tracking=None,savedir_data=None,plotdir=None,cell_selection=None,constraint_tracking_period=None,plotting_period=None):
    from html_animate import make_animation
    features=pd.read_hdf(os.path.join(savedir_tracking,'Features.h5'),'table')
    W_max=iris.load_cube(os.path.join(savedir_data,'Data_w_max.nc')).extract(constraint_tracking_period)
    W_mid_max=iris.load_cube(os.path.join(savedir_data,'Data_w_mid_max.nc')).extract(constraint_tracking_period)
    track=pd.read_hdf(os.path.join(savedir_tracking,'Track_w.h5'),'table')
    mask=iris.load_cube(os.path.join(savedir_tracking,'Mask_Segmentation_w.nc'),'segmentation_mask')
    if plotting_period:
        track_time=track.loc[(track['time']>plotting_period[0]) & (track['time']<plotting_period[1])]
    cmap_TWP,levels_TWP,bounds_TWP,norm_TWP=color_TWP()
    cmap_w,levels_w,bounds_w,norm_w=color_w()

    if cell_selection is None:
        cell_selection=track_time['cell'].dropna().unique()
    else:
        cell_selection=cell_selection
    
    for cell in cell_selection:
        logging.debug('Start plotting mask for cell'+ str(int(cell)))
        name='Masking_w_static_ncells'

        plot_mask_cell_track_static_timeseries(cell=cell,
                                       track=track, 
                                       cog=None,
                                       features=features,
                                       mask_total=mask,
                                       field_contour=W_mid_max, 
                                       label_field_contour='midlevel max w (m/s)',
                                       cmap_field_contour=cmap_w,norm_field_contour=norm_w,
                                       vmin_field_contour=0,vmax_field_contour=15,levels_field_contour=levels_w,
                                       field_filled=W_max, 
                                       label_field_filled='max w (kg/m^2)',
                                       cmap_field_filled=cmap_TWP,
                                       vmin_field_filled=0,vmax_field_filled=15,levels_field_filled=levels_w,
                                       track_variable=track,variable='ncells',variable_label='number of cells (watershedding)',
                                       width=10000,n_extend=2,
                                       name= name+'_'+str(int(cell)), 
                                       plotdir=os.path.join(plotdir,name)
                                       )
        make_animation(input_dir=os.path.join(plotdir,name,name+'_'+str(int(cell))),
                       output=os.path.join(plotdir,'Animations_'+name,name+'_'+str(int(cell))),
                       delay=200,make_gif=False)

def plot_mask_TWC_3D(savedir_tracking=None,savedir_data=None,plotdir=None,cell_selection=None,constraint_tracking_period=None,plotting_period=None):
    from html_animate import make_animation
    features=pd.read_hdf(os.path.join(savedir_tracking,'Features.h5'),'table')
    track=pd.read_hdf(os.path.join(savedir_tracking,'Track_TWC.h5'),'table')
    TWP=iris.load_cube(os.path.join(savedir_data,'Data_WP.nc'),'TWP').extract(constraint_tracking_period)
    W_mid_max=iris.load_cube(os.path.join(savedir_data,'Data_w_mid_max.nc')).extract(constraint_tracking_period)
    mask=iris.load_cube(os.path.join(savedir_tracking,'Mask_Segmentation_TWC.nc'),'segmentation_mask')
    if plotting_period:
        track_time=track.loc[(track['time']>plotting_period[0]) & (track['time']<plotting_period[1])]
    cmap_TWP,levels_TWP,bounds_TWP,norm_TWP=color_TWP()
    cmap_w,levels_w,bounds_w,norm_w=color_w()
    ele, azim =angles_3D() 

    if cell_selection is None:
        cell_selection_TWC=track_time['cell'].dropna().unique()
    else:
        cell_selection_TWC=cell_selection
    
    for cell in cell_selection_TWC:
        logging.debug('Start plotting mask for cell'+ str(int(cell)))
        name='Masking_TWC_static_3D'

        plot_mask_cell_track_3Dstatic(cell=cell,
                                       track=track, 
                                       cog=None,
                                       features=features,
                                       mask_total=mask,
                                       field_contour=W_mid_max, 
                                       label_field_contour='midlevel max w (m/s)',
                                       cmap_field_contour=cmap_w,norm_field_contour=norm_w,
                                       vmin_field_contour=0,vmax_field_contour=15,levels_field_contour=levels_w,
                                       field_filled=TWP, 
                                       label_field_filled='TWP (kg/m^2)',
                                       cmap_field_filled=cmap_TWP,
                                       vmin_field_filled=0,vmax_field_filled=50,levels_field_filled=levels_TWP,
                                       width=10000,n_extend=2,
                                       name= name+'_'+str(int(cell)), 
                                       plotdir=os.path.join(plotdir,name)
                                       )
        
        make_animation(input_dir=os.path.join(plotdir,name,name+'_'+str(int(cell))),
                       output=os.path.join(plotdir,'Animations_'+name,name+'_'+str(int(cell))),
                       delay=200,make_gif=False)

def plot_mask_w_2D3D(savedir_tracking=None,savedir_data=None,plotdir=None,cell_selection=None,constraint_tracking_period=None,plotting_period=None):
    from html_animate import make_animation
    features=pd.read_hdf(os.path.join(savedir_tracking,'Features.h5'),'table')
    W_max=iris.load_cube(os.path.join(savedir_data,'Data_w_max.nc')).extract(constraint_tracking_period)
    W_mid_max=iris.load_cube(os.path.join(savedir_data,'Data_w_mid_max.nc')).extract(constraint_tracking_period)
    track=pd.read_hdf(os.path.join(savedir_tracking,'Track_w.h5'),'table')
    mask=iris.load_cube(os.path.join(savedir_tracking,'Mask_Segmentation_w.nc'),'segmentation_mask')
    if plotting_period:
        track_time=track.loc[(track['time']>plotting_period[0]) & (track['time']<plotting_period[1])]
    cmap_TWP,levels_TWP,bounds_TWP,norm_TWP=color_TWP()
    cmap_w,levels_w,bounds_w,norm_w=color_w()
    ele, azim =angles_3D() 

    if cell_selection is None:
        cell_selection_w=track_time['cell'].dropna().unique()
    else:
        cell_selection_w=cell_selection
    
    for cell in cell_selection_w:
        logging.debug('Start plotting mask for cell'+ str(int(cell)))
        name='Masking_w_static_2D3D'

        plot_mask_cell_track_2D3Dstatic(cell=cell,
                                       track=track, 
                                       cog=None,
                                       features=features,
                                       mask_total=mask,
                                       field_contour=W_mid_max, 
                                       label_field_contour='midlevel max w (m/s)',
                                       cmap_field_contour=cmap_w,norm_field_contour=norm_w,
                                       vmin_field_contour=0,vmax_field_contour=15,levels_field_contour=levels_w,
                                       field_filled=W_max, 
                                       label_field_filled='max w (kg/m^2)',
                                       cmap_field_filled=cmap_TWP,
                                       vmin_field_filled=0,vmax_field_filled=30,levels_field_filled=levels_w,
                                       width=10000,n_extend=2,
                                       name= name+'_'+str(int(cell)), 
                                       plotdir=os.path.join(plotdir,name),
                                       ele=ele,azim=azim
                                       )
        
        make_animation(input_dir=os.path.join(plotdir,name,name+'_'+str(int(cell))),
                       output=os.path.join(plotdir,'Animations_'+name,name+'_'+str(int(cell))),
                       delay=200,make_gif=False)

def plot_mask_w_TWC_3D(savedir_tracking=None,savedir_data=None,plotdir=None,cell_selection=None,constraint_tracking_period=None,plotting_period=None):
    from html_animate import make_animation
    features=pd.read_hdf(os.path.join(savedir_tracking,'Features.h5'),'table')
    TWP=iris.load_cube(os.path.join(savedir_data,'Data_WP.nc'),'TWP').extract(constraint_tracking_period)
    W_mid_max=iris.load_cube(os.path.join(savedir_data,'Data_w_mid_max.nc')).extract(constraint_tracking_period)
    track=pd.read_hdf(os.path.join(savedir_tracking,'Track_w.h5'),'table')
    mask=iris.load_cube(os.path.join(savedir_tracking,'Mask_Segmentation_w_TWC.nc'),'segmentation_mask')
    if plotting_period:
        track_time=track.loc[(track['time']>plotting_period[0]) & (track['time']<plotting_period[1])]
    cmap_TWP,levels_TWP,bounds_TWP,norm_TWP=color_TWP()
    cmap_w,levels_w,bounds_w,norm_w=color_w()
    ele, azim =angles_3D() 
    if cell_selection is None:
        cell_selection=track_time['cell'].dropna().unique()
    else:
        cell_selection=cell_selection
    
    for cell in cell_selection:
        logging.debug('Start plotting mask for cell'+ str(int(cell)))
        name='Masking_w_TWC_static_3D'

        plot_mask_cell_track_3Dstatic(cell=cell,
                                       track=track, 
                                       cog=None,
                                       features=features,
                                       mask_total=mask,
                                       field_contour=W_mid_max, 
                                       label_field_contour='midlevel max w (m/s)',
                                       cmap_field_contour=cmap_w,norm_field_contour=norm_w,
                                       vmin_field_contour=0,vmax_field_contour=15,levels_field_contour=levels_w,
                                       field_filled=TWP, 
                                       label_field_filled='TWP (kg/m^2)',
                                       cmap_field_filled=cmap_TWP,
                                       vmin_field_filled=0,vmax_field_filled=50,levels_field_filled=levels_TWP,
                                       width=10000,n_extend=2,
                                       name= name+'_'+str(int(cell)), 
                                       plotdir=os.path.join(plotdir,name)
                                       )
        
        make_animation(input_dir=os.path.join(plotdir,name,name+'_'+str(int(cell))),
                       output=os.path.join(plotdir,'Animations_'+name,name+'_'+str(int(cell))),
                       delay=200,make_gif=False)

def plot_mask_w_TWC_2D3D(savedir_tracking=None,savedir_data=None,plotdir=None,cell_selection=None,constraint_tracking_period=None,plotting_period=None):
    from html_animate import make_animation
    features=pd.read_hdf(os.path.join(savedir_tracking,'Features.h5'),'table')
    TWP=iris.load_cube(os.path.join(savedir_data,'Data_WP.nc'),'TWP').extract(constraint_tracking_period)
    W_mid_max=iris.load_cube(os.path.join(savedir_data,'Data_w_mid_max.nc')).extract(constraint_tracking_period)
    track=pd.read_hdf(os.path.join(savedir_tracking,'Track_w.h5'),'table')
    mask=iris.load_cube(os.path.join(savedir_tracking,'Mask_Segmentation_w_TWC.nc'),'segmentation_mask')
    if plotting_period:
        track_w_TWC_time=track.loc[(track['time']>plotting_period[0]) & (track['time']<plotting_period[1])]
    cmap_TWP,levels_TWP,bounds_TWP,norm_TWP=color_TWP()
    cmap_w,levels_w,bounds_w,norm_w=color_w()
    ele, azim =angles_3D() 

    if cell_selection is None:
        cell_selection_w_TWC=track_w_TWC_time['cell'].dropna().unique()
    else:
        cell_selection_w_TWC=cell_selection
    
    for cell in cell_selection_w_TWC:
        logging.debug('Start plotting mask for cell'+ str(int(cell)))
        name='Masking_w_TWC_static_2D3D'

        plot_mask_cell_track_2D3Dstatic(cell=cell,
                                       track=track, 
                                       cog=None,
                                       features=features,
                                       mask_total=mask,
                                       field_contour=W_mid_max, 
                                       label_field_contour='midlevel max w (m/s)',
                                       cmap_field_contour=cmap_w,norm_field_contour=norm_w,
                                       vmin_field_contour=0,vmax_field_contour=15,levels_field_contour=levels_w,
                                       field_filled=TWP, 
                                       label_field_filled='TWP (kg/m^2)',
                                       cmap_field_filled=cmap_TWP,
                                       vmin_field_filled=0,vmax_field_filled=50,levels_field_filled=levels_TWP,
                                       width=10000,n_extend=2,
                                       name= name+'_'+str(int(cell)), 
                                       plotdir=os.path.join(plotdir,name),
                                       ele=ele,azim=azim
                                       )
        
        make_animation(input_dir=os.path.join(plotdir,name,name+'_'+str(int(cell))),
                       output=os.path.join(plotdir,'Animations_'+name,name+'_'+str(int(cell))),
                       delay=200,make_gif=False)

def plot_mask_TWC_2D3D(savedir_tracking=None,savedir_data=None,plotdir=None,cell_selection=None,constraint_tracking_period=None,plotting_period=None):
    from html_animate import make_animation
    features=pd.read_hdf(os.path.join(savedir_tracking,'Features.h5'),'table')
    track_TWC=pd.read_hdf(os.path.join(savedir_tracking,'Track_TWC.h5'),'table')
    TWP=iris.load_cube(os.path.join(savedir_data,'Data_WP.nc'),'TWP').extract(constraint_tracking_period)
    W_mid_max=iris.load_cube(os.path.join(savedir_data,'Data_w_mid_max.nc')).extract(constraint_tracking_period)
    mask_TWC=iris.load_cube(os.path.join(savedir_tracking,'Mask_Segmentation_TWC.nc'),'segmentation_mask')
    if plotting_period:
        track_TWC_time=track_TWC.loc[(track_TWC['time']>plotting_period[0]) & (track_TWC['time']<plotting_period[1])]
    cmap_TWP,levels_TWP,bounds_TWP,norm_TWP=color_TWP()
    cmap_w,levels_w,bounds_w,norm_w=color_w()
    ele, azim =angles_3D() 
    if cell_selection is None:
        cell_selection_TWC=track_TWC_time['cell'].dropna().unique()
    else:
        cell_selection_TWC=cell_selection
    
    for cell in cell_selection_TWC:
        logging.debug('Start plotting mask for cell'+ str(int(cell)))
        name='Masking_TWC_static_2D3D'

        plot_mask_cell_track_2D3Dstatic(cell=cell,
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
                                       width=10000,n_extend=2,
                                       name= name+'_'+str(int(cell)), 
                                       plotdir=os.path.join(plotdir,name),
                                       ele=ele,azim=azim
                                       )
        
        make_animation(input_dir=os.path.join(plotdir,name,name+'_'+str(int(cell))),
                       output=os.path.join(plotdir,'Animations_'+name,name+'_'+str(int(cell))),
                       delay=200,make_gif=False)

def plot_mask_w_3D(savedir_tracking=None,savedir_data=None,plotdir=None,cell_selection=None,constraint_tracking_period=None,plotting_period=None):
    from html_animate import make_animation
    features=pd.read_hdf(os.path.join(savedir_tracking,'Features.h5'),'table')
    W_max=iris.load_cube(os.path.join(savedir_data,'Data_w_max.nc')).extract(constraint_tracking_period)
    W_mid_max=iris.load_cube(os.path.join(savedir_data,'Data_w_mid_max.nc')).extract(constraint_tracking_period)
    track=pd.read_hdf(os.path.join(savedir_tracking,'Track_w.h5'),'table')
    mask=iris.load_cube(os.path.join(savedir_tracking,'Mask_Segmentation_w.nc'),'segmentation_mask')
    if plotting_period:
        track_time=track.loc[(track['time']>plotting_period[0]) & (track['time']<plotting_period[1])]
    cmap_TWP,levels_TWP,bounds_TWP,norm_TWP=color_TWP()
    cmap_w,levels_w,bounds_w,norm_w=color_w()
    ele, azim =angles_3D() 

    if cell_selection is None:
        cell_selection_w=track_time['cell'].dropna().unique()
    else:
        cell_selection_w=cell_selection
    
    for cell in cell_selection_w:
        logging.debug('Start plotting mask for cell'+ str(int(cell)))
        name='Masking_w_static_3D'

        plot_mask_cell_track_3Dstatic(cell=cell,
                                       track=track, 
                                       cog=None,
                                       features=features,
                                       mask_total=mask,
                                       field_contour=W_mid_max, 
                                       label_field_contour='midlevel max w (m/s)',
                                       cmap_field_contour=cmap_w,norm_field_contour=norm_w,
                                       vmin_field_contour=0,vmax_field_contour=15,levels_field_contour=levels_w,
                                       field_filled=W_max, 
                                       label_field_filled='max w (kg/m^2)',
                                       cmap_field_filled=cmap_TWP,
                                       vmin_field_filled=0,vmax_field_filled=30,levels_field_filled=levels_w,
                                       width=10000,n_extend=2,
                                       name= name+'_'+str(int(cell)), 
                                       plotdir=os.path.join(plotdir,name)
                                       )
        
        make_animation(input_dir=os.path.join(plotdir,name,name+'_'+str(int(cell))),
                       output=os.path.join(plotdir,'Animations_'+name,name+'_'+str(int(cell))),
                       delay=200,make_gif=False)

def plot_mask_precipitation(savedir_tracking=None,savedir_data=None,plotdir=None,cell_selection=None,constraint_tracking_period=None,plotting_period=None,masking='w'):
    from html_animate import make_animation
    logging.debug('Start plotting tracks and masks for individual cells')
    features=pd.read_hdf(os.path.join(savedir_tracking,'Features.h5'),'table')
    track=pd.read_hdf(os.path.join(savedir_tracking,'Track.h5'),'table')
    W_mid_max=iris.load_cube(os.path.join(savedir_data,'Data_w_mid_max.nc')).extract(constraint_tracking_period)    
    TWP=iris.load_cube(os.path.join(savedir_data,'Data_WP.nc'),'TWP').extract(constraint_tracking_period)  
    mask=iris.load_cube(os.path.join(savedir_tracking,f'Mask_Segmentation_{masking}.nc'),'segmentation_mask')
    cmap_TWP,levels_TWP,bounds_TWP,norm_TWP=color_TWP()
    cmap_w,levels_w,bounds_w,norm_w=color_w()

    name=f'Masking_{masking}_static_precip'
    logging.debug('data loaded from savefile')    
    if plotting_period:
            track_time=track.loc[(track['time']>plotting_period[0]) & (track['time']<plotting_period[1])]
    
    if cell_selection is None:
        cell_selection=np.unique(track_time['cell'])
        
    logging.debug(f'cell_selection: {cell_selection}')
    
    # loop over individual cells for analysis
    for cell in cell_selection:
        name_cell=f'{name}_{int(cell)}'
        savedir_cell=os.path.join(savedir_tracking,'cells',str(int(cell)))
        precip=['surface_precipitation_instantaneous']
        track_precip=pd.read_hdf(os.path.join(savedir_cell,f'track_processes_lumped_integrated_{masking}_{int(cell)}.h5'),'table')
        plot_mask_cell_track_static_timeseries(cell=cell,
                                           track=track_precip, 
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
                                           track_variable=track_precip,variable=precip,
                                           variable_label=precip,variable_color='cornflowerblue',
                                           variable_ylabel='surface precipitation (kg/s/m^2)',variable_legend=True,
                                           width=10000,n_extend=2,
                                           name=name_cell, 
                                           plotdir=os.path.join(plotdir,name),
                                           )
        make_animation(input_dir=os.path.join(plotdir,name,name_cell),
                       output=os.path.join(plotdir,f'Animations_{name}',name_cell),
                       delay=200,make_gif=False)

def plot_mask_TWC(savedir_tracking=None,savedir_data=None,plotdir=None,cell_selection=None,constraint_tracking_period=None,plotting_period=None):
    from html_animate import make_animation
    features=pd.read_hdf(os.path.join(savedir_tracking,'Features.h5'),'table')
    track=pd.read_hdf(os.path.join(savedir_tracking,'Track_TWC.h5'),'table')
    TWP=iris.load_cube(os.path.join(savedir_data,'Data_WP.nc'),'TWP').extract(constraint_tracking_period)
    W_mid_max=iris.load_cube(os.path.join(savedir_data,'Data_w_mid_max.nc')).extract(constraint_tracking_period)
    mask=iris.load_cube(os.path.join(savedir_tracking,'Mask_Segmentation_TWC.nc'),'segmentation_mask')
    if plotting_period:
        track_time=track.loc[(track['time']>plotting_period[0]) & (track['time']<plotting_period[1])]
    cmap_TWP,levels_TWP,bounds_TWP,norm_TWP=color_TWP()
    cmap_w,levels_w,bounds_w,norm_w=color_w()
    if cell_selection is None:
        cell_selection=track_time['cell'].dropna().unique()
    for cell in cell_selection:  
        logging.debug('Start plotting mask for cell'+ str(int(cell)))
        # Make plots centred around and moving with the tracked cell (10km each side)
        name='Masking_TWC_static_ncells'
        plot_mask_cell_track_static_timeseries(cell=cell,
                                       track=track, 
                                       cog=None,
                                       features=features,
                                       mask_total=mask,
                                       field_contour=W_mid_max, 
                                       label_field_contour='midlevel max w (m/s)',
                                       cmap_field_contour=cmap_w,norm_field_contour=norm_w,
                                       vmin_field_contour=0,vmax_field_contour=15,levels_field_contour=levels_w,
                                       field_filled=TWP, 
                                       label_field_filled='TWP (kg/m^2)',
                                       cmap_field_filled=cmap_TWP,
                                       vmin_field_filled=0,vmax_field_filled=50,levels_field_filled=levels_TWP,
                                       track_variable=track,variable='ncells',variable_label='number of cells (watershedding)',
                                       width=10000,n_extend=2,
                                       name= name+'_'+str(int(cell)), 
                                       plotdir=os.path.join(plotdir,name),
                                       )
        make_animation(input_dir=os.path.join(plotdir,name,name+'_'+str(int(cell))),
                       output=os.path.join(plotdir,'Animations'+'_'+name,name+'_'+str(int(cell))),
                       delay=200,make_gif=False)
