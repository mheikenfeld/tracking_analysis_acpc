import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

import logging

from tobac import plot_tracks_mask_field_loop
from tobac import plot_mask_cell_track_static_timeseries
from tobac import plot_mask_cell_track_3Dstatic,plot_mask_cell_track_2D3Dstatic
from tobac import cell_statistics

from tobac import feature_detection_multithreshold
from tobac import segmentation_3D
from tobac import linking_trackpy
from tobac import get_spacings

from .cell_analysis import extract_cell_cubes_subset,extract_cell_cubes_subset_2D,interpolate_alongacross_mean
import iris
from iris.analysis import MIN,MAX,MEAN,MEDIAN,PERCENTILE
import os
import pandas as pd
import numpy as np

from matplotlib import colors
from multiprocessing.pool import ThreadPool
import dask

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
                      parameters_features=None,
                      parameters_linking=None,
                      parameters_tracking=None,
                      parameters_segmentation_TWC=None,
                      parameters_segmentation_w=None,
                      parameters_segmentation_w_2D=None,
                      plotting_parameters=None,
                      plotting_period=None,
                      n_core=1):
    
    print('n_core: '+str(n_core))
    logging.debug('ncore: '+str(n_core) )
    
    if n_core==1:
        dask.config.set(scheduler='single-threaded')
    elif n_core>1:
        dask.set_options(pool=ThreadPool(n_core-1))
    
    os.makedirs(top_savedir_tracking,exist_ok=True)
    f = open(os.path.join(top_savedir_tracking,'parameters_tracking.txt'), 'w')
    f.write('tracking_period: ' + repr(tracking_period) + '\n' )
    f.write('parameters_tracking: ' + repr(parameters_tracking) + '\n' )
    f.write('parameters_segmentation_TWC: ' + repr(parameters_segmentation_TWC) + '\n' )
    f.write('parameters_segmentation_w: ' + repr(parameters_segmentation_w) + '\n' )
    f.close()
    
    logging.debug('start tracking analysis')
    logging.debug('top_savedir_data: '+ top_savedir_data)
    logging.debug('top_savedir_tracking: '+ top_savedir_tracking)
    logging.debug('top_plotdir: '+ top_plotdir)

    logging.debug('number of processors used :'+ str(n_core))
    
    logging.debug(f'cell_selection : {cell_selection}')

    if loading_period:
        constraint_loading_period=iris.Constraint(time = lambda cell: loading_period[0] <= cell <  loading_period[1])
    elif loading_period is None:
        constraint_loading_period=iris.Constraint(None)

    if tracking_period:
        constraint_tracking_period=iris.Constraint(time = lambda cell: tracking_period[0] <= cell <  tracking_period[1])
    elif tracking_period is None:
        constraint_tracking_period=iris.Constraint(None)

    #Tracking setup:
    compression=True
    if 'load_data_w' in mode:
        load_w(filenames=filenames,
                  savedir=top_savedir_data,
                  constraint_loading_period=constraint_loading_period,
                  load_function=load_function,compression=compression)

    if 'load_data_WC' in mode:
        load_WC(filenames=filenames,
                  savedir=top_savedir_data,
                  constraint_loading_period=constraint_loading_period,
                  load_function=load_function,compression=compression)

    if 'load_data_WP' in mode:
        load_WP(filenames=filenames,
                  savedir=top_savedir_data,
                  constraint_loading_period=constraint_loading_period,
                  load_function=load_function,compression=compression)
        
    if 'load_data_Airmass' in mode:
        load_Airmass(filenames=filenames,
                  savedir=top_savedir_data,
                  constraint_loading_period=constraint_loading_period,
                  load_function=load_function,compression=compression)

    if 'load_data_temperature' in mode:
        load_temperature(filenames=filenames,
                  savedir=top_savedir_data,
                  constraint_loading_period=constraint_loading_period,
                  load_function=load_function,compression=compression)

    if 'load_data_winds' in mode:
        load_winds(filenames=filenames,
                  savedir=top_savedir_data,
                  constraint_loading_period=constraint_loading_period,
                  load_function=load_function,compression=compression)

    if 'load_data_precipitation' in mode:
        load_precipitation(filenames=filenames,
                  savedir=top_savedir_data,
                  constraint_loading_period=constraint_loading_period,
                  load_function=load_function,compression=compression)

    if 'load_data_hydrometeors' in mode:
        load_hydrometeors(filenames=filenames,
                  savedir=top_savedir_data,
                  constraint_loading_period=constraint_loading_period,
                  load_function=load_function,compression=compression)

    if 'load_processes' in mode:
        load_processes(filenames=filenames,
                       savedir=top_savedir_data,
                       constraint_loading_period=constraint_loading_period,
                       load_function=load_function,compression=compression)
        
    if 'load_processes_detailed' in mode:
        
        load_processes_detailed(filenames=filenames,
                       savedir=top_savedir_data,
                       constraint_loading_period=constraint_loading_period,
                       load_function=load_function,compression=compression)
        
    if 'load_processes_all' in mode:
        
        load_processes_all(filenames=filenames,
                       savedir=top_savedir_data,
                       constraint_loading_period=constraint_loading_period,
                       load_function=load_function,compression=compression)

    if 'detection' in mode:
        features=run_feature_detection(loaddir=os.path.join(top_savedir_data),savedir=os.path.join(top_savedir_tracking),
                                       constraint_tracking_period=constraint_tracking_period,
                                       parameters_features=parameters_features)
    
    if 'segmentation' in mode:
        run_segmentation_TWC(loaddir=os.path.join(top_savedir_data),savedir=os.path.join(top_savedir_tracking),
                             compression=True,
                             constraint_tracking_period=constraint_tracking_period,
                             parameters_segmentation_TWC=parameters_segmentation_TWC)
        
        run_segmentation_w(loaddir=os.path.join(top_savedir_data),savedir=os.path.join(top_savedir_tracking),
                           compression=True,
                           constraint_tracking_period=constraint_tracking_period,
                           parameters_segmentation_w=parameters_segmentation_w)
        
        run_segmentation_w_TWC(loaddir=os.path.join(top_savedir_data),savedir=os.path.join(top_savedir_tracking),
                               constraint_tracking_period=constraint_tracking_period,
                               compression=True,
                               parameters_segmentation_w=parameters_segmentation_w,
                               parameters_segmentation_TWC=parameters_segmentation_TWC)
        
    if 'segmentation_TWC' in mode:
        run_segmentation_TWC(loaddir=os.path.join(top_savedir_data),savedir=os.path.join(top_savedir_tracking),
                             compression=True,
                             constraint_tracking_period=constraint_tracking_period,
                             parameters_segmentation_TWC=parameters_segmentation_TWC)
        
    if 'segmentation_w' in mode:
        run_segmentation_w(loaddir=os.path.join(top_savedir_data),savedir=os.path.join(top_savedir_tracking),
                           compression=True,
                           constraint_tracking_period=constraint_tracking_period,
                           parameters_segmentation_w=parameters_segmentation_w)
        
    if 'segmentation_w_TWC' in mode:
        run_segmentation_w_TWC(loaddir=os.path.join(top_savedir_data),savedir=os.path.join(top_savedir_tracking),
                               constraint_tracking_period=constraint_tracking_period,
                               compression=True,
                               parameters_segmentation_w=parameters_segmentation_w,
                               parameters_segmentation_TWC=parameters_segmentation_TWC)

    if 'linking' in mode:
        run_linking(loaddir=os.path.join(top_savedir_data),savedir=os.path.join(top_savedir_tracking),
                    constraint_tracking_period=constraint_tracking_period,
                    parameters_linking=parameters_linking)
        
        #merge information from segmentation into track and save
        merge_features_track(savedir=os.path.join(top_savedir_tracking))

    if 'plot' in mode:
        plots_maps_loop(savedir_data=top_savedir_data,
                        savedir_tracking=top_savedir_tracking,
                        plotdir=top_plotdir,
                        plotting_parameters=plotting_parameters,
                        constraint_tracking_period=constraint_tracking_period)

    if ('plot_mask_TWC' in mode):
        plot_mask_TWC(savedir_tracking=top_savedir_tracking,savedir_data=top_savedir_data,plotdir=top_plotdir,
                      constraint_tracking_period=constraint_tracking_period,plotting_period=plotting_period)
        
    if ('plot_mask_w' in mode):
        plot_mask_w(savedir_tracking=top_savedir_tracking,savedir_data=top_savedir_data,plotdir=top_plotdir,
                                    cell_selection=cell_selection,constraint_tracking_period=constraint_tracking_period,plotting_period=plotting_period)
    if ('plot_mask_w_TWC' in mode):
        plot_mask_w_TWC(savedir_tracking=top_savedir_tracking,savedir_data=top_savedir_data,plotdir=top_plotdir,
                        cell_selection=cell_selection,constraint_tracking_period=constraint_tracking_period,plotting_period=plotting_period)

    if ('plot_mask_w_3D' in mode):
        plot_mask_w_3D(savedir_tracking=top_savedir_tracking,savedir_data=top_savedir_data,plotdir=top_plotdir,
                       cell_selection=cell_selection,constraint_tracking_period=constraint_tracking_period,plotting_period=plotting_period)

    if ('plot_mask_TWC_3D' in mode):
        plot_mask_TWC_3D(savedir_tracking=top_savedir_tracking,savedir_data=top_savedir_data,plotdir=top_plotdir,
                         cell_selection=cell_selection,constraint_tracking_period=constraint_tracking_period,plotting_period=plotting_period)

    if ('plot_mask_w_TWC_3D' in mode):
        plot_mask_w_TWC_3D(savedir_tracking=top_savedir_tracking,savedir_data=top_savedir_data,plotdir=top_plotdir,
                           cell_selection=cell_selection,constraint_tracking_period=constraint_tracking_period,plotting_period=plotting_period)

    if ('plot_mask_w_2D3D' in mode):
        plot_mask_w_2D3D(savedir_tracking=top_savedir_tracking,savedir_data=top_savedir_data,plotdir=top_plotdir,
                         cell_selection=cell_selection,constraint_tracking_period=constraint_tracking_period,plotting_period=plotting_period)

    if ('plot_mask_TWC_2D3D' in mode):
        plot_mask_TWC_2D3D(savedir_tracking=top_savedir_tracking,savedir_data=top_savedir_data,plotdir=top_plotdir,
                           cell_selection=cell_selection,constraint_tracking_period=constraint_tracking_period,plotting_period=plotting_period)

    if ('plot_mask_w_TWC_2D3D' in mode):
        plot_mask_w_TWC_2D3D(savedir_tracking=top_savedir_tracking,savedir_data=top_savedir_data,plotdir=top_plotdir,
                             cell_selection=cell_selection,constraint_tracking_period=constraint_tracking_period,plotting_period=plotting_period)

    if (('plot_mask_processes_TWC' not in mode) and ('interpolation_processes_TWC' in mode)):
        interpolation_processes_TWC(savedir_tracking=top_savedir_tracking,savedir_data=top_savedir_data,
                                    cell_selection=cell_selection,constraint_tracking_period=constraint_tracking_period,plotting_period=plotting_period)

    if (('plot_mask_processes_w_TWC'  not in mode) and ('interpolation_processes_w_TWC' in mode)):
        interpolation_processes_w_TWC(savedir_tracking=top_savedir_tracking,savedir_data=top_savedir_data,
                                    cell_selection=cell_selection,constraint_tracking_period=constraint_tracking_period,plotting_period=plotting_period)

    if (('plot_mask_processes_TWC' in mode) and ('interpolation_processes_TWC' not in mode)):
        plot_mask_processes_TWC(savedir_tracking=top_savedir_tracking,savedir_data=top_savedir_data,plotdir=top_plotdir,
                                cell_selection=cell_selection,constraint_tracking_period=constraint_tracking_period,plotting_period=plotting_period)

    if (('plot_mask_processes_w_TWC' in mode) and ('interpolation_processes_w_TWC' not in mode)):
        plot_mask_processes_w_TWC(savedir_tracking=top_savedir_tracking,savedir_data=top_savedir_data,plotdir=top_plotdir,
                                  cell_selection=cell_selection,constraint_tracking_period=constraint_tracking_period,plotting_period=plotting_period)

    if (('plot_mask_processes_TWC' in mode) and ('interpolation_processes_TWC' in mode)):
        interpolation_plot_processes_TWC(savedir_tracking=top_savedir_tracking,savedir_data=top_savedir_data,plotdir=top_plotdir,
                                  cell_selection=cell_selection,constraint_tracking_period=constraint_tracking_period,plotting_period=plotting_period)

    if (('plot_mask_processes_w_TWC' in mode) and ('plot_mask_processes_w_TWC' in mode)):
        interpolation_plot_processes_w_TWC(savedir_tracking=top_savedir_tracking,savedir_data=top_savedir_data,plotdir=top_plotdir,
                                  cell_selection=cell_selection,constraint_tracking_period=constraint_tracking_period,plotting_period=plotting_period)

    if (('plot_mask_hydrometeors_TWC' not in mode) and ('interpolation_hydrometeors_TWC' in mode)):
        interpolation_hydrometeors_TWC(savedir_tracking=top_savedir_tracking,savedir_data=top_savedir_data,
                                    cell_selection=cell_selection,constraint_tracking_period=constraint_tracking_period,plotting_period=plotting_period)

    if (('plot_mask_hydrometeors_w_TWC'  not in mode) and ('interpolation_hydrometeors_w_TWC' in mode)):
        interpolation_hydrometeors_w_TWC(savedir_tracking=top_savedir_tracking,savedir_data=top_savedir_data,
                                    cell_selection=cell_selection,constraint_tracking_period=constraint_tracking_period,plotting_period=plotting_period)

#     if (('plot_mask_hydrometeors_TWC' in mode) and ('interpolation_hydrometeors_TWC' not in mode)):
#         plot_mask_hydrometeors_TWC(savedir_tracking=top_savedir_tracking,savedir_data=top_savedir_data,plotdir=top_plotdir,
#                                 cell_selection=cell_selection,constraint_tracking_period=constraint_tracking_period,plotting_period=plotting_period)

    # if (('plot_mask_hydrometeors_w_TWC' in mode) and ('interpolation_hydrometeors_w_TWC' not in mode)):
    #     plot_mask_hydrometeors_w_TWC(savedir_tracking=top_savedir_tracking,savedir_data=top_savedir_data,plotdir=top_plotdir,
    #                               cell_selection=cell_selection,constraint_tracking_period=constraint_tracking_period,plotting_period=plotting_period)

    # if (('plot_mask_hydrometeors_TWC' in mode) and ('interpolation_hydrometeors_TWC' in mode)):
    #     interpolation_plot_hydrometeors_TWC(savedir_tracking=top_savedir_tracking,savedir_data=top_savedir_data,plotdir=top_plotdir,
    #                               cell_selection=cell_selection,constraint_tracking_period=constraint_tracking_period,plotting_period=plotting_period)

    # if (('plot_mask_hydrometeors_w_TWC' in mode) and ('plot_mask_hydrometeors_w_TWC' in mode)):
    #     interpolation_plot_hydrometeors_w_TWC(savedir_tracking=top_savedir_tracking,savedir_data=top_savedir_data,plotdir=top_plotdir,
    #                               cell_selection=cell_selection,constraint_tracking_period=constraint_tracking_period,plotting_period=plotting_period)



    if 'plot_lifetime' in mode:
        plot_lifetime(savedir=top_savedir_tracking,plotdir=top_plotdir)

    if ('profiles_CDNC' in mode):
        profiles_cdnc(savedir_tracking=top_savedir_tracking,savedir_data=top_savedir_data,filenames=filenames,load_function=load_function)

    if 'plot_profiles_CDNC' in mode:
        plot_profile_CDNC(savedir_tracking=top_savedir_tracking,plotdir=top_plotdir)

    if ('profiles_w' in mode):
        calculate_profiles_w(savedir_tracking=top_savedir_tracking,savedir_data=top_savedir_data,filenames=filenames,load_function=load_function)

    if 'plot_profiles_w' in mode:
        plot_profiles_w(savedir_tracking=top_savedir_tracking,plotdir=top_plotdir)

    if (('interpolation_mass_TWC' in mode) and ('plot_mask_mass_TWC' not in mode)):
            interpolation_mass_TWC(savedir_tracking=top_savedir_tracking,savedir_data=top_savedir_data,
                                  cell_selection=cell_selection,constraint_tracking_period=constraint_tracking_period)
            
    if (('plot_mask_mass_TWC' in mode) and ('interpolation_mass_TWC' not in mode)):
            plot_mask_mass_TWC(savedir_tracking=top_savedir_tracking,savedir_data=top_savedir_data,plotdir=top_plotdir,
                           cell_selection=cell_selection,constraint_tracking_period=constraint_tracking_period,plotting_period=plotting_period)
            
    if (('plot_mask_mass_TWC' in mode) and ('interpolation_mass_TWC' in mode)):
                interpolation_plot_mask_mass_TWC(savedir_tracking=top_savedir_tracking,savedir_data=top_savedir_data,plotdir=top_plotdir,
                                             cell_selection=cell_selection,constraint_tracking_period=constraint_tracking_period,plotting_period=plotting_period)

    if (('interpolation_mass_w_TWC' in mode) and ('plot_mask_mass_w_TWC' not in mode)):
            interpolation_mass_w_TWC(savedir_tracking=top_savedir_tracking,savedir_data=top_savedir_data,
                                  cell_selection=cell_selection,constraint_tracking_period=constraint_tracking_period)
            
    if (('plot_mask_mass_w_TWC' in mode) and ('interpolation_mass_w_TWC' not in mode)):
            plot_mask_mass_w_TWC(savedir_tracking=top_savedir_tracking,savedir_data=top_savedir_data,plotdir=top_plotdir,
                           cell_selection=cell_selection,constraint_tracking_period=constraint_tracking_period,plotting_period=plotting_period)
            
    if (('plot_mask_mass_w_TWC' in mode) and ('interpolation_mass_w_TWC' in mode)):
                interpolation_plot_mask_mass_w_TWC(savedir_tracking=top_savedir_tracking,savedir_data=top_savedir_data,plotdir=top_plotdir,
                                             cell_selection=cell_selection,constraint_tracking_period=constraint_tracking_period,plotting_period=plotting_period)

    # if (('interpolation_precipitation_TWC' in mode)):
    #         interpolation_precipitation_TWC(savedir_tracking=top_savedir_tracking,savedir_data=top_savedir_data,
    #                               cell_selection=cell_selection,constraint_tracking_period=constraint_tracking_period)

    # if (('interpolation_precipitation_w' in mode)):
    #         interpolation_precipitation_w(savedir_tracking=top_savedir_tracking,savedir_data=top_savedir_data,
    #                               cell_selection=cell_selection,constraint_tracking_period=constraint_tracking_period)

    # if (('interpolation_precipitation_w_TWC' in mode)):
    #         interpolation_precipitation_w_TWC(savedir_tracking=top_savedir_tracking,savedir_data=top_savedir_data,
    #                               cell_selection=cell_selection,constraint_tracking_period=constraint_tracking_period)

    if (('interpolation_precipitation_TWC' in mode)):
            interpolation_precipitation(savedir_tracking=top_savedir_tracking,savedir_data=top_savedir_data,
                                  cell_selection=cell_selection,constraint_tracking_period=constraint_tracking_period,masking='TWC')

    if (('interpolation_precipitation_w' in mode)):
            interpolation_precipitation(savedir_tracking=top_savedir_tracking,savedir_data=top_savedir_data,
                                  cell_selection=cell_selection,constraint_tracking_period=constraint_tracking_period,masking='w')

    if (('interpolation_precipitation_w_TWC' in mode)):
            interpolation_precipitation(savedir_tracking=top_savedir_tracking,savedir_data=top_savedir_data,
                                  cell_selection=cell_selection,constraint_tracking_period=constraint_tracking_period,masking='w_TWC')

    if ('plot_mask_precipitation_TWC' in mode):
        plot_mask_precipitation(savedir_tracking=top_savedir_tracking,savedir_data=top_savedir_data,plotdir=top_plotdir,
                                cell_selection=cell_selection,constraint_tracking_period=constraint_tracking_period,plotting_period=plotting_period,masking='TWC')

    if ('plot_mask_precipitation_w' in mode):
        plot_mask_precipitation(savedir_tracking=top_savedir_tracking,savedir_data=top_savedir_data,plotdir=top_plotdir,
                                cell_selection=cell_selection,constraint_tracking_period=constraint_tracking_period,plotting_period=plotting_period,masking='w')

    if ('plot_mask_precipitation_w_TWC' in mode):
        plot_mask_precipitation(savedir_tracking=top_savedir_tracking,savedir_data=top_savedir_data,plotdir=top_plotdir,
                                cell_selection=cell_selection,constraint_tracking_period=constraint_tracking_period,plotting_period=plotting_period,masking='w_TWC')
    if 'slices_processes_lumped' in mode:
        slices_processes_lumped(savedir_tracking=top_savedir_tracking,savedir_data=top_savedir_data,
                                  cell_selection=cell_selection,constraint_tracking_period=constraint_tracking_period)
    if 'slices_processes_all' in mode:
        slices_processes_all(savedir_tracking=top_savedir_tracking,savedir_data=top_savedir_data,
                                  cell_selection=cell_selection,constraint_tracking_period=constraint_tracking_period)
        
    if 'slices_hydrometeors_mass' in mode:
        slices_hydrometeors_mass(savedir_tracking=top_savedir_tracking,savedir_data=top_savedir_data,
                                  cell_selection=cell_selection,constraint_tracking_period=constraint_tracking_period)

    if 'slices_hydrometeors_number' in mode:
        slices_hydrometeors_number(savedir_tracking=top_savedir_tracking,savedir_data=top_savedir_data,
                                  cell_selection=cell_selection,constraint_tracking_period=constraint_tracking_period)
        
    if 'slices_aux' in mode:
        slices_aux(savedir_tracking=top_savedir_tracking,savedir_data=top_savedir_data,
                                  cell_selection=cell_selection,constraint_tracking_period=constraint_tracking_period)


logging.debug('done')

def add_geopotential_height(cube):
    from .make_geopotential_height_coord import geopotential_height_coord,geopotential_height_coord_stag
    coord_names=[coord.name() for coord in cube.coords()]
    if ('model_level_number' in coord_names) and ('geopotential_height' not in coord_names):
        if cube.coord('model_level_number').shape[0]==95:
            cube.add_aux_coord(geopotential_height_coord_stag,data_dims=cube.coord_dims('model_level_number'))
        if cube.coord('model_level_number').shape[0]==94:
            cube.add_aux_coord(geopotential_height_coord,data_dims=cube.coord_dims('model_level_number'))
    return cube

def load_w(filenames,savedir,constraint_loading_period,load_function,compression=True):
    logging.info('start loading data for tracking')
    os.makedirs(os.path.join(savedir),exist_ok=True)
    
    logging.debug(' start loading w')
    W=load_function(filenames,'w').extract(constraint_loading_period)
    if compression:
         zlib=True
         complevel=4
         packing=None
         chunksizes=list(W.shape)
         chunksizes[0]=1
    else:
         zlib=False
         complevel=None
         packing=None
         chunksizes=None

    iris.save([W],os.path.join(savedir,'Data_w.nc'),zlib=zlib,complevel=complevel,chunksizes=chunksizes,packing=packing)
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
    
def load_winds(filenames,savedir,constraint_loading_period,load_function,compression=True):
    logging.info('start loading data for tracking')
    os.makedirs(os.path.join(savedir),exist_ok=True)
    
    logging.debug(' start loading w')
    W=load_function(filenames,'w_unstaggered').extract(constraint_loading_period)
    logging.debug(' start loading u')
    U=load_function(filenames,'u_unstaggered').extract(constraint_loading_period)
    logging.debug(' start loading v')
    V=load_function(filenames,'v_unstaggered').extract(constraint_loading_period)

    if compression:
         zlib=True
         complevel=4
         packing=None
         chunksizes=list(W.shape)
         chunksizes[0]=1
    else:
         zlib=False
         complevel=None
         packing=None
         chunksizes=None

    iris.save([U,V,W],os.path.join(savedir,'Data_winds.nc'),zlib=zlib,complevel=complevel,chunksizes=chunksizes,packing=packing)
    logging.debug('winds loaded and saved')

def load_temperature(filenames,savedir,constraint_loading_period,load_function,compression=True):
    logging.info('start loading data for tracking')
    os.makedirs(os.path.join(savedir),exist_ok=True)
    
    logging.debug(' start loading theta')
    Theta=load_function(filenames,'potential_temperature').extract(constraint_loading_period)
    logging.debug(' start loading T')
    T=load_function(filenames,'air_temperature').extract(constraint_loading_period)

    if compression:
         zlib=True
         complevel=4
         packing=None
         chunksizes=list(Theta.shape)
         chunksizes[0]=1
    else:
         zlib=False
         complevel=None
         packing=None
         chunksizes=None

    iris.save([Theta,T],os.path.join(savedir,'Data_temperature.nc'),zlib=zlib,complevel=complevel,chunksizes=chunksizes,packing=packing)
    logging.debug('Temperature loaded and saved')

def load_WC(filenames,savedir,constraint_loading_period,load_function,compression=True):
    logging.debug(' start loading water content')
    os.makedirs(os.path.join(savedir),exist_ok=True)

    IWC=load_function(filenames,'IWC').extract(constraint_loading_period)  
    LWC=load_function(filenames,'LWC').extract(constraint_loading_period)  
    TWC=IWC+LWC
    TWC.rename('TWC')
    
    if compression:
         zlib=True
         complevel=4
         packing=None
         chunksizes=list(TWC.shape)
         chunksizes[0]=1
    else:
         zlib=False
         complevel=None
         packing=None
         chunksizes=None
         
    logging.debug('water content loaded start saving')

    iris.save([TWC,LWC,IWC],os.path.join(savedir,'Data_WC.nc'),zlib=zlib,complevel=complevel,chunksizes=chunksizes,packing=packing)
    logging.debug('water content loaded and saved')

def load_Airmass(filenames,savedir,constraint_loading_period,load_function,compression=True):
    logging.debug(' start loading airmass')
    os.makedirs(os.path.join(savedir),exist_ok=True)

    Airmass=load_function(filenames,'airmass').extract(constraint_loading_period) 
    logging.debug('airmass loaded')

    Airmass_path=load_function(filenames,'airmass_path').extract(constraint_loading_period) 
    logging.debug('airmass path loaded')

    rho=load_function(filenames,'air_density').extract(constraint_loading_period) 
    logging.debug('density loaded')

    
    if compression:
         zlib=True
         complevel=4
         packing=None
         chunksizes=list(Airmass.shape)
         chunksizes[0]=1
    else:
         zlib=False
         complevel=None
         packing=None
         chunksizes=None
         
    logging.debug('water content loaded start saving')

    iris.save([Airmass,Airmass_path,rho],os.path.join(savedir,'Data_Airmass.nc'),zlib=zlib,complevel=complevel,chunksizes=chunksizes,packing=packing)
    logging.debug('Airmass and density loaded and saved')



def load_WP(filenames,savedir,constraint_loading_period,load_function,compression=True):
    logging.debug(' start loading water path')
    os.makedirs(os.path.join(savedir),exist_ok=True)

    IWP=load_function(filenames,'IWP').extract(constraint_loading_period) 
    logging.debug('ice water path loaded')

    LWP=load_function(filenames,'LWP').extract(constraint_loading_period)      
    logging.debug('liquid water path loaded')

    TWP=IWP+LWP
    TWP.rename('TWP')
    
    logging.debug('water paths loaded, start saving')

    if compression:
         zlib=True
         complevel=4
         packing=None
         chunksizes=list(TWP.shape)
         chunksizes[0]=1
    else:
         zlib=False
         complevel=None
         packing=None
         chunksizes=None

    iris.save([TWP,LWP,IWP],os.path.join(savedir,'Data_WP.nc'),zlib=zlib,complevel=complevel,chunksizes=chunksizes,packing=packing)
    logging.debug('water path loaded and saved')

def load_precipitation(filenames,savedir,constraint_loading_period,load_function,compression=True):
    logging.debug(' start loading precipitation')
    os.makedirs(os.path.join(savedir),exist_ok=True)

    surface_precipitation=iris.cube.CubeList()
    surface_precipitation.append(load_function(filenames,'surface_precipitation_average').extract(constraint_loading_period))
    logging.debug('average from accumulated surface precipitation loaded')

    surface_precipitation.append(load_function(filenames,'surface_precipitation_accumulated').extract(constraint_loading_period))
    logging.debug('accumulated surface precipitation loaded')

    # surface_precipitation.append(load_function(filenames,'surface_precipitation_instantaneous').extract(constraint_loading_period))
    # logging.debug('instantaneous surface precipitation loaded')

    logging.debug('precipitation loaded, start saving')

    if compression:
         zlib=True
         complevel=4
         packing=None
         chunksizes=list(surface_precipitation[0].shape)
         chunksizes[0]=1
    else:
         zlib=False
         complevel=None
         packing=None
         chunksizes=None

    iris.save(surface_precipitation,os.path.join(savedir,'Data_Precip.nc'),zlib=zlib,complevel=complevel,chunksizes=chunksizes,packing=packing)
    logging.debug('precipitation loaded and saved')



def load_hydrometeors(filenames,savedir,constraint_loading_period,load_function,compression=True):
    logging.debug(' start loading hydrometeors')
    os.makedirs(os.path.join(savedir),exist_ok=True)

    hydrometeors_mass=load_function(filenames,'hydrometeor_mass').extract(constraint_loading_period) 
    logging.debug('hydrometeors_mass loaded')
    
    hydrometeors_number=load_function(filenames,'hydrometeor_number').extract(constraint_loading_period) 
    logging.debug('hydrometeors_number loaded')

    if compression:
         zlib=True
         complevel=4
         packing=None
         chunksizes=list(hydrometeors_mass[0].shape)
         chunksizes[0]=1
    else:
         zlib=False
         complevel=None
         packing=None
         chunksizes=None

    iris.save(hydrometeors_mass,os.path.join(savedir,'Data_hydrometeor_mass.nc'),zlib=zlib,complevel=complevel,chunksizes=chunksizes,packing=packing)
    iris.save(hydrometeors_number,os.path.join(savedir,'Data_hydrometeor_number.nc'),zlib=zlib,complevel=complevel,chunksizes=chunksizes,packing=packing)
    logging.debug('hydrometeors loaded and saved')


def load_processes(filenames,savedir,constraint_loading_period,load_function,compression=True,processes_type='processes_lumped'):
    logging.info('start loading data for tracking')
    os.makedirs(os.path.join(savedir),exist_ok=True)


    logging.debug(' start process rates')

    processes_lumped=load_function(filenames,'processes_lumped').extract(constraint_loading_period)    
    
    if compression:
         zlib=True
         complevel=4
         packing=None
         chunksizes=list(processes_lumped[0].shape)
         chunksizes[0]=1
    else:
         zlib=False
         complevel=None
         packing=None
         chunksizes=None
                 
    iris.save(processes_lumped,os.path.join(savedir,'Data_Processes_lumped.nc'),zlib=zlib,complevel=complevel,chunksizes=chunksizes,packing=packing)
    logging.debug('loaded and saved processes')
    del(processes_lumped)

    logging.debug('data loaded and saved')

def load_processes_detailed(filenames,savedir,constraint_loading_period,load_function,compression=True,processes_type='processes_lumped'):
    logging.info('start loading data for tracking')
    os.makedirs(os.path.join(savedir),exist_ok=True)


    logging.debug(' start process rates')

    processes_lumped=load_function(filenames,'processes_detailed').extract(constraint_loading_period)    
    
    if compression:
         zlib=True
         complevel=4
         packing=None
         chunksizes=list(processes_lumped[0].shape)
         chunksizes[0]=1
    else:
         zlib=False
         complevel=None
         packing=None
         chunksizes=None
                 
    iris.save(processes_lumped,os.path.join(savedir,'Data_Processes_lumped_detailed.nc'),zlib=zlib,complevel=complevel,chunksizes=chunksizes,packing=packing)
    logging.debug('loaded and saved processes')
    del(processes_lumped)

    logging.debug('data loaded and saved')

def load_processes_all(filenames,savedir,constraint_loading_period,load_function,compression=True):
    logging.info('start loading data for tracking')
    os.makedirs(os.path.join(savedir),exist_ok=True)

    logging.debug(' start loading process rates')
    processes=load_function(filenames,'processes_all').extract(constraint_loading_period)    
    
    if compression:
         zlib=True
         complevel=4
         packing=None
         chunksizes=list(processes[0].shape)
         chunksizes[0]=1
    else:
         zlib=False
         complevel=None
         packing=None
         chunksizes=None
                 
    iris.save(processes,os.path.join(savedir,'Data_Processes_all.nc'),zlib=zlib,complevel=complevel,chunksizes=chunksizes,packing=packing)
    logging.debug('loaded and saved processes')
    del(processes)

    logging.debug('data loaded and saved')


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

def run_feature_detection(loaddir=None,savedir=None,
                          constraint_tracking_period=None,
                          parameters_features=None):
    
    logging.debug('start loading data for tracking')
    os.makedirs(savedir,exist_ok=True)

    #W_max=iris.load_cube(os.path.join(loaddir,'Data_w_max.nc')).extract(constraint_tracking_period)  
    W_mid_max=iris.load_cube(os.path.join(loaddir,'Data_w_mid_max.nc')).extract(constraint_tracking_period)  

    dxy,dt=get_spacings(W_mid_max)
    logging.debug('data loaded for tracking')

    logging.info('start feature detection based on midlevel maximum vertical velocity:')
    
    method_detection=parameters_features.pop('method_detection')
    print(method_detection)
    if method_detection=='threshold_multi':
        features=feature_detection_multithreshold(W_mid_max,dxy,**parameters_features)
        features.to_hdf(os.path.join(savedir,'Features.h5'),'table')
    else:
        raise ValueError('method_detection unknown, only "threshold_multi" is implemented')
                        
    logging.info('start segmentation based in TWC')
    # Set to True for compression (decreases size of saved files to mb's but saving process takes really long)
    return features

def run_segmentation_TWC(loaddir=None, savedir=None,
                         compression=True,
                         constraint_tracking_period=None,
                         parameters_segmentation_TWC=None):
    # perform 3D segmentation based on updraft TWC
    logging.info('start watershedding TWC')

    # load results from previous steps:
    features=pd.read_hdf(os.path.join(savedir,'Features.h5'),'table')
    TWC=iris.load_cube(os.path.join(loaddir,'Data_WC.nc'),'TWC').extract(constraint_tracking_period)
    dxy,dt=get_spacings(TWC)
    mask,features_out=segmentation_3D(features,TWC,dxy,**parameters_segmentation_TWC)
     
    #Set compression for saving mask
    if compression:
         zlib=True
         complevel=4
         packing=None
         chunksizes=list(mask.shape)
         chunksizes[0]=1
         packing={'dtype': np.int32, 'scale_factor':1, 'add_offset':1}
    else:
         zlib=False
         complevel=None
         packing=None
         chunksizes=None
    logging.debug('segmentation TWC performed, start saving results to files')
    iris.save([mask],os.path.join(savedir,'Mask_Segmentation_TWC.nc'),zlib=zlib,complevel=complevel,packing=packing,chunksizes=chunksizes)                
    features_out.to_hdf(os.path.join(savedir,'Features_TWC.h5'),'table')
    logging.debug('segmentation TWC performed and saved')

def run_segmentation_w(loaddir=None, savedir=None,
                       compression=True,
                       constraint_tracking_period=None,
                       parameters_segmentation_w=None):
    # perform 3D segmentation based on updraft velocity
    logging.info('start watershedding w')
    # load results from previous steps:
    features=pd.read_hdf(os.path.join(savedir,'Features.h5'),'table')
    W=iris.load_cube(os.path.join(loaddir,'Data_w.nc')).extract(constraint_tracking_period)
    dxy,dt=get_spacings(W)
    mask,features_out=segmentation_3D(features,W,dxy,**parameters_segmentation_w)
    #Set compression for saving mask
    if compression:
         zlib=True
         complevel=4
         packing=None
         chunksizes=list(mask.shape)
         chunksizes[0]=1
         packing={'dtype': np.int32, 'scale_factor':1, 'add_offset':1}
    else:
         zlib=False
         complevel=None
         packing=None
         chunksizes=None

    logging.debug('segmentation w performed, start saving results to files')
    iris.save([mask],os.path.join(savedir,'Mask_Segmentation_w.nc'),zlib=zlib,complevel=complevel,packing=packing,chunksizes=chunksizes)
    features_out.to_hdf(os.path.join(savedir,'Features_w.h5'),'table')
    logging.debug('segmentation w performed and saved')

def run_segmentation_w_TWC(loaddir=None, savedir=None,
                           constraint_tracking_period=None,
                           compression=True,
                           parameters_segmentation_TWC=None,parameters_segmentation_w=None):
    # perform 3D segmentation based on updraft velocity and TWC
    logging.info('start watershedding w with TWC masking')
    # load results from previous steps:
    features=pd.read_hdf(os.path.join(savedir,'Features.h5'),'table')
    W=iris.load_cube(os.path.join(loaddir,'Data_w.nc')).extract(constraint_tracking_period)
    TWC=iris.load_cube(os.path.join(loaddir,'Data_WC.nc'),'TWC').extract(constraint_tracking_period)            
    dxy,dt=get_spacings(W)
    # Set field to zero for TWC below threshold
    if (W.shape==TWC.shape):
        W_TWC=W*(TWC.core_data()>parameters_segmentation_TWC['threshold'])
    # Workaround for WRF, where W is on a staggered grid...
    elif (W.shape[1]==(TWC.shape[1]+1)):
        W_TWC=W[:,:-1]*(TWC.core_data()>parameters_segmentation_TWC['threshold'])
    mask,features_out=segmentation_3D(features,W_TWC,dxy,**parameters_segmentation_w)
        #Set compression for saving mask
    if compression:
         zlib=True
         complevel=4
         packing=None
         chunksizes=list(mask.shape)
         chunksizes[0]=1
         packing={'dtype': np.int32, 'scale_factor':1, 'add_offset':1}
    else:
         zlib=False
         complevel=None
         packing=None
         chunksizes=None

    logging.debug('segmentation w_TWC performed, start saving results to files')
    iris.save([mask],os.path.join(savedir,'Mask_Segmentation_w_TWC.nc'),zlib=zlib,complevel=complevel,chunksizes=chunksizes,packing=packing)
    features_out.to_hdf(os.path.join(savedir,'Features_w_TWC.h5'),'table')
    logging.debug('segmentation w_TWC performed and saved')

def run_linking(loaddir=None,savedir=None,
                constraint_tracking_period=None,parameters_linking=None):
    # perform 3D segmentation based on updraft velocity and TWC
    logging.info('start trajectory linking')
    # load results from previous steps:
    features=pd.read_hdf(os.path.join(savedir,'Features.h5'),'table')
    W_mid_max=iris.load_cube(os.path.join(loaddir,'Data_w_mid_max.nc')).extract(constraint_tracking_period)  
    dxy,dt=get_spacings(W_mid_max)
    track=linking_trackpy(features,W_mid_max,dt=dt,dxy=dxy,**parameters_linking)
    #save trajectories:
    track.to_hdf(os.path.join(savedir,'Track.h5'),'table')
    return track

def merge_features_track(savedir):
    track=pd.read_hdf(os.path.join(savedir,'Track.h5'),'table')
    features_w=pd.read_hdf(os.path.join(savedir,'Features_w.h5'),'table')
    features_TWC=pd.read_hdf(os.path.join(savedir,'Features_TWC.h5'),'table')
    features_w_TWC=pd.read_hdf(os.path.join(savedir,'Features_w_TWC.h5'),'table')

    track_w=features_w.merge(track[['feature','cell','time_cell']],on='feature')
    track_w.to_hdf(os.path.join(savedir,'Track_w.h5'),'table')
    #merge information from segmentation into track and save
    track_TWC=features_TWC.merge(track[['feature','cell','time_cell']],on='feature')
    track_TWC.to_hdf(os.path.join(savedir,'Track_TWC.h5'),'table')
    #merge information from segmentation into track and save
    track_w_TWC=features_w_TWC.merge(track[['feature','cell','time_cell']],on='feature')
    track_w_TWC.to_hdf(os.path.join(savedir,'Track_w_TWC.h5'),'table')


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

def calculate_profiles_w(savedir_tracking=None,savedir_data=None,plotdir=None):
    track=pd.read_hdf(os.path.join(savedir_tracking,'Track.h5'),'table')
    mask_w=iris.load_cube(os.path.join(savedir_tracking,'Mask_Segmentation_w.nc'),'segmentation_mask')
    logging.debug('W loaded')

    W=iris.load_cube(os.path.join(savedir_data,'Data_w.nc'))
    W=add_geopotential_height(W)                

    dimensions=('x','y')
    cubelist_Aux=iris.cube.CubeList([W]) 
    for cell in np.unique(track['cell']):

        cell_statistics(cubelist_Aux,track,mask_w,dimensions=dimensions,cell=cell,aggregators=[MEAN,MAX,MIN,MEDIAN],output_path=savedir_tracking,output_name='W')
        cell_statistics(cubelist_Aux,track,mask_w,dimensions=dimensions,cell=cell,aggregators=[PERCENTILE],output_path=savedir_tracking,output_name='W',percent=[0,25,50,75,100])
        logging.debug('profile calculations done and saved for cell '+str(cell))

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

def profiles_cdnc(savedir_tracking=None,savedir_data=None,filenames=None,load_function=None):

    track=pd.read_hdf(os.path.join(savedir_tracking,'Track.h5'),'table')
    mask_w=iris.load_cube(os.path.join(savedir_tracking,'Mask_Segmentation_w.nc'),'segmentation_mask')
    logging.debug('start loading NCLOUD')

#    QCLOUD=load_function(filenames,'QCLOUD')
    NCLOUD=load_function(filenames,'NCLOUD')
    
    logging.debug('start loading air_density')
    rho=load_function(filenames,'air_density')
    
    logging.debug(str(NCLOUD))
    logging.debug(str(rho))
    CDNC=NCLOUD*rho
    CDNC.rename('cloud_droplet_number_concentration')
    CDNC=add_geopotential_height(CDNC)
    logging.debug(str(CDNC))

    if (CDNC.coord('model_level_number').shape[0] == mask_w.coord('model_level_number').shape[0]):
        mask_w_unstaggered=mask_w
    elif (CDNC.coord('model_level_number').shape[0] == mask_w.coord('model_level_number').shape[0]-1):
        mask_w_unstaggered=mask_w[:,1:,:,:]

    logging.debug('CDNC loaded')

    savedir=os.path.join(savedir_tracking)

    logging.debug('Start calculating profiles for individual cells (simple method)')


    dimensions=('x','y')
    cubelist_Aux=iris.cube.CubeList([CDNC])

    for cell in np.unique(track['cell']):
        cell_statistics(cubelist_Aux,track,mask_w_unstaggered,aggregators=[MEAN,MAX,MIN,MEDIAN],dimensions=dimensions,cell=cell,output_path=savedir,output_name='CDNC')                 
        cell_statistics(cubelist_Aux,track,mask_w_unstaggered,aggregators=[PERCENTILE],dimensions=dimensions,cell=cell,output_path=savedir,output_name='CDNC',percent=[0,25,50,75,100])                 
        logging.debug('profile calculations done and saved for cell '+str(cell))

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

def interpolation_mass_TWC(savedir_tracking=None,savedir_data=None,plotdir=None,cell_selection=None,constraint_tracking_period=None):
    logging.debug('loading tracks and Mask')

    track=pd.read_hdf(os.path.join(savedir_tracking,'Track.h5'),'table')
    mask=iris.load_cube(os.path.join(savedir_tracking,'Mask_Segmentation_TWC.nc'),'segmentation_mask')
    
    #load water contents and fix coordinate names:
    WC=iris.load(os.path.join(savedir_data,'Data_WC.nc')).extract(constraint_tracking_period)   
    logging.debug('loaded total water content from file')
    for cube in WC:
        for coord in cube.coords():
            coord.var_name=coord.name()

    # load airmass (used to calculate total values from mixing ratios) and fix coordinate names:
    airmass=iris.load_cube(os.path.join(savedir_data,'Data_Airmass.nc'),'airmass').extract(constraint_tracking_period) 
    for coord in airmass.coords():
        coord.var_name=coord.name()
    logging.debug('airmass loaded from file')

    if cell_selection is None:
        cell_selection=track['cell'].dropna().unique()

    for cell in cell_selection:                        
        savedir_cell=os.path.join(savedir_tracking,'cells',str(int(cell)))
        logging.debug( 'Start calculating Data for cell '+ str(int(cell)))
        savedir_cell=os.path.join(savedir_tracking,'cells',str(int(cell)))
        os.makedirs(savedir_cell,exist_ok=True)
        tracks=pd.read_hdf(os.path.join(savedir_tracking,'Track.h5'),'table')
        mask=iris.load_cube(os.path.join(savedir_tracking,'Mask_Segmentation_TWC.nc'),'segmentation_mask')


        height_levels=np.arange(0,20001,1000.)
        
        mass_sum=iris.cube.CubeList()

        cube=WC.extract_strict('TWC')

        for coord in airmass.coords():
            print(coord)
        for coord in cube.coords():
            print(coord)

        cube_sum=cube*airmass
        cube_sum.rename('total_water_mass')
        cube_sum=add_geopotential_height(cube_sum)                
        mass_sum.append(cube_sum)

        cube=WC.extract_strict('LWC')
        cube_sum=cube*airmass
        cube_sum.rename('liquid_water_mass')
        cube_sum=add_geopotential_height(cube_sum)                
        mass_sum.append(cube_sum)

        cube=WC.extract_strict('IWC')
        cube_sum=cube*airmass
        cube_sum.rename('frozen_water_mass')
        cube_sum=add_geopotential_height(cube_sum)                
        mass_sum.append(cube_sum)
            
        mass_cell_sum, mass_cell_integrated,track_mass_integrated=extract_cell_cubes_subset(cubelist_in=mass_sum,
                                                                                            mask=mask,
                                                                                            track=tracks,
                                                                                            cell=cell,
                                                                                            z_coord='model_level_number',
                                                                                            height_levels=height_levels)

        iris.save(mass_cell_sum,os.path.join(savedir_cell,'Mass_profile_TWC.nc'),zlib=True,complevel=4)
        iris.save(mass_cell_integrated,os.path.join(savedir_cell,'Mass_integrated_TWC.nc'),zlib=True,complevel=4)
        track_mass_integrated.to_hdf(os.path.join(savedir_cell,'track_mass_integrated_TWC.h5'),'table')
        logging.debug( 'mass calculated and saved to '+ savedir_cell)

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
        
def interpolation_mass_w_TWC(savedir_tracking=None,savedir_data=None,plotdir=None,cell_selection=None,constraint_tracking_period=None):
    logging.debug('loading tracks and Mask')

    track=pd.read_hdf(os.path.join(savedir_tracking,'Track.h5'),'table')
    mask=iris.load_cube(os.path.join(savedir_tracking,'Mask_Segmentation_TWC.nc'),'segmentation_mask')
    
    #load water contents and fix coordinate names:
    WC=iris.load(os.path.join(savedir_data,'Data_WC.nc')).extract(constraint_tracking_period)   
    logging.debug('loaded total water content from file')
    for cube in WC:
        for coord in cube.coords():
            coord.var_name=coord.name()

    # load airmass (used to calculate total values from mixing ratios) and fix coordinate names:
    airmass=iris.load_cube(os.path.join(savedir_data,'Data_Airmass.nc'),'airmass').extract(constraint_tracking_period) 
    for coord in airmass.coords():
        coord.var_name=coord.name()
    logging.debug('airmass loaded from file')

    if cell_selection is None:
        cell_selection=track['cell'].dropna().unique()

    for cell in cell_selection:                        
        savedir_cell=os.path.join(savedir_tracking,'cells',str(int(cell)))
        logging.debug( 'Start calculating Data for cell '+ str(int(cell)))
        savedir_cell=os.path.join(savedir_tracking,'cells',str(int(cell)))
        os.makedirs(savedir_cell,exist_ok=True)
        tracks=pd.read_hdf(os.path.join(savedir_tracking,'Track.h5'),'table')
        mask=iris.load_cube(os.path.join(savedir_tracking,'Mask_Segmentation_w_TWC.nc'),'segmentation_mask')


        height_levels=np.arange(0,20001,1000.)
        
        mass_sum=iris.cube.CubeList()

        cube=WC.extract_strict('TWC')

        for coord in airmass.coords():
            print(coord)
        for coord in cube.coords():
            print(coord)

        cube_sum=cube*airmass
        cube_sum.rename('total_water_mass')
        cube_sum=add_geopotential_height(cube_sum)                
        mass_sum.append(cube_sum)

        cube=WC.extract_strict('LWC')
        cube_sum=cube*airmass
        cube_sum.rename('liquid_water_mass')
        cube_sum=add_geopotential_height(cube_sum)                
        mass_sum.append(cube_sum)

        cube=WC.extract_strict('IWC')
        cube_sum=cube*airmass
        cube_sum.rename('frozen_water_mass')
        cube_sum=add_geopotential_height(cube_sum)                
        mass_sum.append(cube_sum)
            
        mass_cell_sum, mass_cell_integrated,track_mass_integrated=extract_cell_cubes_subset(cubelist_in=mass_sum,
                                                                                            mask=mask,
                                                                                            track=tracks,
                                                                                            cell=cell,
                                                                                            z_coord='model_level_number',
                                                                                            height_levels=height_levels)

        iris.save(mass_cell_sum,os.path.join(savedir_cell,'Mass_profile_w_TWC.nc'),zlib=True,complevel=4)
        iris.save(mass_cell_integrated,os.path.join(savedir_cell,'Mass_integrated_w_TWC.nc'),zlib=True,complevel=4)
        track_mass_integrated.to_hdf(os.path.join(savedir_cell,'track_mass_integrated_w_TWC'+str(int(cell))+'.h5'),'table')
        logging.debug( 'mass calculated and saved to '+ savedir_cell)

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

def interpolation_processes_w_TWC(savedir_tracking=None,savedir_data=None,cell_selection=None,constraint_tracking_period=None,plotting_period=None):
    logging.debug('Start calculations for individual cells')
        
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
        cell_selection=tracks['cell'].dropna().unique()

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

def interpolation_processes_TWC(savedir_tracking=None,savedir_data=None,cell_selection=None,constraint_tracking_period=None,plotting_period=None):
    logging.debug('Start calculations for individual cells')
        
    logging.debug('loading tracks and Mask')

    tracks=pd.read_hdf(os.path.join(savedir_tracking,'Track.h5'),'table')
    mask=iris.load_cube(os.path.join(savedir_tracking,'Mask_Segmentation_TWC.nc'),'segmentation_mask')

    processes_file=os.path.join(savedir_data,'Data_Processes_lumped.nc')
    if os.path.exists(processes_file):       
        logging.debug('start loading processes form preprocessed file:')
        processes_lumped=iris.load(processes_file).extract(constraint_tracking_period)        
#            else:#        
#                logging.debug(' start process rates from data files')
#                processes_lumped=load_function(filenames,'processes_lumped').extract(constraint_tracking_period)                 

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
        cell_selection=tracks['cell'].dropna().unique()

    logging.debug(f'cell_selection: {cell_selection}')


    # loop over individual cells for analysis
    for cell in cell_selection:                        
        logging.debug( 'Start calculating Data for cell '+ str(int(cell)))
        savedir_cell=os.path.join(savedir_tracking,'cells',str(int(cell)))
        os.makedirs(savedir_cell,exist_ok=True)

        processes_lumped_cell_sum, processes_lumped_cell_integrated,track_processes_lumped_integrated=extract_cell_cubes_subset(cubelist_in=processes_lumped_sum,
                                                                                                                             mask=mask,
                                                                                                                             track=tracks,
                                                                                                                             cell=cell,
                                                                                                                             z_coord='model_level_number',
                                                                                                                             height_levels=height_levels)
        
        iris.save(processes_lumped_cell_sum,os.path.join(savedir_cell,'Processes_lumped_profile_TWC.nc'),zlib=True,complevel=4)
        iris.save(processes_lumped_cell_integrated,os.path.join(savedir_cell,'Processes_lumped_integrated_TWC.nc'),zlib=True,complevel=4)
        track_processes_lumped_integrated.to_hdf(os.path.join(savedir_cell,'track_processes_lumped_integrated_TWC.h5'),'table')
        
        logging.debug( 'processes calculated and saved to '+ savedir_cell)

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

def interpolation_hydrometeors_TWC(savedir_tracking=None,savedir_data=None,cell_selection=None,constraint_tracking_period=None,plotting_period=None):
    logging.debug('Start calculations for individual cells')
        
    logging.debug('loading tracks and Mask')

    tracks=pd.read_hdf(os.path.join(savedir_tracking,'Track.h5'),'table')
    mask=iris.load_cube(os.path.join(savedir_tracking,'Mask_Segmentation_TWC.nc'),'segmentation_mask')

    hydrometeors_file=os.path.join(savedir_data,'Data_hydrometeor_mass.nc')
    if os.path.exists(hydrometeors_file):       
        logging.debug('start loading hydrometers form preprocessed file:')
        hydrometeors=iris.load(hydrometeors_file).extract(constraint_tracking_period)        
#            else:#        
#                logging.debug(' start process rates from data files')
#                processes_lumped=load_function(filenames,'processes_lumped').extract(constraint_tracking_period)                 

    for cube in hydrometeors:
        for coord in cube.coords():
            coord.var_name=coord.name()
            coord.coord_system=None

    logging.debug('loaded hydrometeors')

    logging.debug('start loading auxillary variables')
    
    
    # load airmass (used to calculate total values from mixing ratios)
    airmass=iris.load_cube(os.path.join(savedir_data,'Data_Airmass.nc'),'airmass').extract(constraint_tracking_period) 
    for coord in airmass.coords():
        coord.var_name=coord.name()
    logging.debug('airmass loaded from file')



    height_levels=np.arange(0,20001,1000.)
                    
    # Sum up process rates for individual cells
    hydrometeors_sum=iris.cube.CubeList()
    
    for cube in hydrometeors:
        cube_sum=cube*airmass
        cube_sum.rename(cube.name())
        cube_sum=add_geopotential_height(cube_sum)
        hydrometeors_sum.append(cube_sum)
#            
    if cell_selection is None:
        cell_selection=tracks['cell'].dropna().unique()

    logging.debug(f'cell_selection: {cell_selection}')
    # loop over individual cells for analysis
    for cell in cell_selection:                        
        logging.debug( 'Start calculating Data for cell '+ str(int(cell)))
        savedir_cell=os.path.join(savedir_tracking,'cells',str(int(cell)))
        os.makedirs(savedir_cell,exist_ok=True)

        hydrometeors_cell_sum, hydrometeors_cell_integrated,track_hydrometeors_integrated=extract_cell_cubes_subset(cubelist_in=hydrometeors_sum,
                                                                                                                             mask=mask,
                                                                                                                             track=tracks,
                                                                                                                             cell=cell,
                                                                                                                             z_coord='model_level_number',
                                                                                                                             height_levels=height_levels)
        
        iris.save(hydrometeors_cell_sum,os.path.join(savedir_cell,'Hydrometeors_profile_TWC.nc'),zlib=True,complevel=4)
        iris.save(hydrometeors_cell_integrated,os.path.join(savedir_cell,'Hydrometeors_integrated_TWC.nc'),zlib=True,complevel=4)
        track_hydrometeors_integrated.to_hdf(os.path.join(savedir_cell,'track_hydrometeors_integrated_TWC.h5'),'table')
        
        logging.debug( 'hydrometeors calculated and saved to '+ savedir_cell)





def interpolation_precipitation(savedir_tracking=None,savedir_data=None,cell_selection=None,constraint_tracking_period=None,plotting_period=None,masking='w'):
    logging.debug('Start calculations for individual cells')
        
    logging.debug('loading tracks and Mask')

    tracks=pd.read_hdf(os.path.join(savedir_tracking,'Track.h5'),'table')
    mask=iris.load_cube(os.path.join(savedir_tracking,f'Mask_Segmentation_{masking}.nc'),'segmentation_mask')

    precipitation=iris.load(os.path.join(savedir_data,'Data_Precip.nc')).extract(constraint_tracking_period)        

    logging.debug('loaded precipitation')
    
    if cell_selection is None:
        cell_selection=tracks['cell'].dropna().unique()

    logging.debug(f'cell_selection: {cell_selection}')

    # loop over individual cells for analysis
    for cell in cell_selection:                        
        logging.debug(f'Start calculating Data for cell {int(cell)}')
        savedir_cell=os.path.join(savedir_tracking,'cells',str(int(cell)))
        os.makedirs(savedir_cell,exist_ok=True)

#
        precipitation_sum,track_precipitation=extract_cell_cubes_subset_2D(cubelist_in=precipitation,
                                                         mask=mask,
                                                         track=tracks,
                                                         cell=cell,
                                                         z_coord='model_level_number')
        
        iris.save(precipitation_sum,os.path.join(savedir_cell,f'Precipitation_{masking}.nc'))

        track_precipitation.to_hdf(os.path.join(savedir_cell,f'track_precipitation_{masking}.h5'),'table')
        
        logging.debug(f'precipitation calculated and saved to {savedir_cell}')

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
        
        
def slices_processes_lumped(savedir_tracking=None,savedir_data=None,cell_selection=None,constraint_tracking_period=None,plotting_period=None):
    logging.debug('Start calculations for individual cells')
        
    logging.debug('loading tracks and Mask')

    tracks=pd.read_hdf(os.path.join(savedir_tracking,'Track.h5'),'table')

    processes_file=os.path.join(savedir_data,'Data_Processes_lumped.nc')
    if os.path.exists(processes_file):       
        logging.debug('start loading processes form preprocessed file:')
        processes=iris.load(processes_file).extract(constraint_tracking_period)
        
        for cube in processes:
            cube=add_geopotential_height(cube)
            
    logging.debug('loaded processes')

    height_levels=np.arange(0,20001,1000.)

    if cell_selection is None:
        cell_selection=tracks['cell'].dropna().unique()

    logging.debug(f'cell_selection: {cell_selection}')

    # loop over individual cells for analysis
    for cell in cell_selection:                        
        logging.debug( 'Start calculating Data for cell '+ str(int(cell)))
        savedir_cell=os.path.join(savedir_tracking,'cells',str(int(cell)))
        os.makedirs(savedir_cell,exist_ok=True)

        Processes_along,Processes_across=interpolate_alongacross_mean(processes,tracks,cell,dx=500,width=10000,z_coord='geopotential_height',height_level_borders=height_levels)
        iris.save(Processes_along,os.path.join(savedir_cell,f'Processes_lumped_along.nc'),zlib=True,complevel=4)
        iris.save(Processes_across,os.path.join(savedir_cell,f'Processes_lumped_across.nc'),zlib=True,complevel=4)
        
        
        logging.debug( 'Processes along and across calculated and saved to '+ savedir_cell)

def slices_processes_all(savedir_tracking=None,savedir_data=None,cell_selection=None,constraint_tracking_period=None,plotting_period=None):
    logging.debug('Start calculations for individual cells')
        
    logging.debug('loading tracks and Mask')

    tracks=pd.read_hdf(os.path.join(savedir_tracking,'Track.h5'),'table')

    processes_file=os.path.join(savedir_data,'Data_Processes_all.nc')
    if os.path.exists(processes_file):       
        logging.debug('start loading processes form preprocessed file:')
        processes=iris.load(processes_file).extract(constraint_tracking_period)
        
        for cube in processes:
            cube=add_geopotential_height(cube)
            
    print(processes[0])
    logging.debug('loaded processes')

    height_levels=np.arange(0,20001,1000.)

    if cell_selection is None:
        cell_selection=tracks['cell'].dropna().unique()

    logging.debug(f'cell_selection: {cell_selection}')

    # loop over individual cells for analysis
    for cell in cell_selection:                        
        logging.debug( 'Start calculating Data for cell '+ str(int(cell)))
        savedir_cell=os.path.join(savedir_tracking,'cells',str(int(cell)))
        os.makedirs(savedir_cell,exist_ok=True)

        Processes_along,Processes_across=interpolate_alongacross_mean(processes,tracks,cell,dx=500,width=10000,z_coord='geopotential_height',height_level_borders=height_levels)
        iris.save(Processes_along,os.path.join(savedir_cell,f'Processes_all_along.nc'),zlib=True,complevel=4)
        iris.save(Processes_across,os.path.join(savedir_cell,f'Processes_all_across.nc'),zlib=True,complevel=4)
        
        
        logging.debug( 'Processes along and across calculated and saved to '+ savedir_cell)

def slices_hydrometeors_mass(savedir_tracking=None,savedir_data=None,cell_selection=None,constraint_tracking_period=None,plotting_period=None):
    logging.debug('Start calculations for individual cells')
        
    logging.debug('loading tracks and Mask')

    tracks=pd.read_hdf(os.path.join(savedir_tracking,'Track.h5'),'table')

    hydrometeors_file=os.path.join(savedir_data,'Data_hydrometeor_mass.nc')
    if os.path.exists(hydrometeors_file):       
        logging.debug('start loading hydrometers form preprocessed file:')
        hydrometeors=iris.load(hydrometeors_file).extract(constraint_tracking_period)
        
        for cube in hydrometeors:
            cube=add_geopotential_height(cube)
            
    print(hydrometeors[0])

    logging.debug('loaded hydrometeors')

    height_levels=np.arange(0,20001,1000.)

    if cell_selection is None:
        cell_selection=tracks['cell'].dropna().unique()
        
    logging.debug(f'cell_selection: {cell_selection}')


    # loop over individual cells for analysis
    for cell in cell_selection:                        
        logging.debug( 'Start calculating Data for cell '+ str(int(cell)))
        savedir_cell=os.path.join(savedir_tracking,'cells',str(int(cell)))
        os.makedirs(savedir_cell,exist_ok=True)



        Hydrometeors_mass_along,Hydrometeors_mass_across=interpolate_alongacross_mean(hydrometeors,tracks,cell,dx=500,width=10000,z_coord='geopotential_height',height_level_borders=height_levels)
        iris.save(Hydrometeors_mass_along,os.path.join(savedir_cell,f'Hydrometeors_mass_along.nc'),zlib=True,complevel=4)
        iris.save(Hydrometeors_mass_across,os.path.join(savedir_cell,f'Hydrometeors_mass_across.nc'),zlib=True,complevel=4)
        
        logging.debug( 'hydrometeors mass along and across calculated and saved to '+ savedir_cell)

def slices_hydrometeors_number(savedir_tracking=None,savedir_data=None,cell_selection=None,constraint_tracking_period=None,plotting_period=None):
    logging.debug('Start calculations for individual cells')
        
    logging.debug('loading tracks and Mask')

    tracks=pd.read_hdf(os.path.join(savedir_tracking,'Track.h5'),'table')

    hydrometeors_file=os.path.join(savedir_data,'Data_hydrometeor_number.nc')
    if os.path.exists(hydrometeors_file):       
        logging.debug('start loading hydrometers form preprocessed file:')
        hydrometeors=iris.load(hydrometeors_file).extract(constraint_tracking_period)        
#            else:#        
#                logging.debug(' start process rates from data files')
#                processes_lumped=load_function(filenames,'processes_lumped').extract(constraint_tracking_period)                 

        for cube in hydrometeors:
            cube=add_geopotential_height(cube)

    logging.debug('loaded hydrometeors')

    height_levels=np.arange(0,20001,1000.)

    if cell_selection is None:
        cell_selection=tracks['cell'].dropna().unique()

    logging.debug(f'cell_selection: {cell_selection}')

    # loop over individual cells for analysis
    for cell in cell_selection:                        
        logging.debug( 'Start calculating Data for cell '+ str(int(cell)))
        savedir_cell=os.path.join(savedir_tracking,'cells',str(int(cell)))
        os.makedirs(savedir_cell,exist_ok=True)

        Hydrometeors_number_along,Hydrometeors_number_across=interpolate_alongacross_mean(hydrometeors,tracks,cell,dx=500,width=10000,z_coord='geopotential_height',height_level_borders=height_levels)
        iris.save(Hydrometeors_number_along,os.path.join(savedir_cell,f'Hydrometeors_number_along.nc'),zlib=True,complevel=4)
        iris.save(Hydrometeors_number_across,os.path.join(savedir_cell,f'Hydrometeors_number_across.nc'),zlib=True,complevel=4)
        
        logging.debug( 'hydrometeors number along and across calculated and saved to '+ savedir_cell)

def slices_aux(savedir_tracking=None,savedir_data=None,cell_selection=None,constraint_tracking_period=None,plotting_period=None):
    logging.debug('Start calculations for individual cells')
        
    logging.debug('loading tracks and Mask')

    tracks=pd.read_hdf(os.path.join(savedir_tracking,'Track.h5'),'table')

    WC_file=os.path.join(savedir_data,'Data_WC.nc')
    if os.path.exists(WC_file):       
        logging.debug('start loading water content form preprocessed file:')
        WC=iris.load(WC_file).extract(constraint_tracking_period)
        for cube in WC:
            cube=add_geopotential_height(cube)

        
    # winds_file=os.path.join(savedir_data,'Data_winds.nc')
    # if os.path.exists(winds_file):       
    #     logging.debug('start loading winds form preprocessed file:')
    #     winds=iris.load(winds_file).extract(constraint_tracking_period)
        # for cube in winds:
        #     cube=add_geopotential_height(cube)

    temperature_file=os.path.join(savedir_data,'Data_temperature.nc')
    if os.path.exists(temperature_file):       
        logging.debug('start loading temperature form preprocessed file:')
        temperature=iris.load(temperature_file).extract(constraint_tracking_period)
        for cube in temperature:
            cube=add_geopotential_height(cube)

    airmass_file=os.path.join(savedir_data,'Data_Airmass.nc')
    if os.path.exists(airmass_file):       
        logging.debug('start loading airmass form preprocessed file:')
        airmass=iris.load(airmass_file).extract(constraint_tracking_period)
        for cube in airmass:
            cube=add_geopotential_height(cube)

    Aux=iris.cube.CubeList()
    Aux.extend(WC)
    # Aux.extend(winds)
    Aux.extend(temperature)
    Aux.extend(airmass.extract('air_density'))

    logging.debug('loaded hydrometeors')

    height_levels=np.arange(0,20001,1000.)

    if cell_selection is None:
        cell_selection=tracks['cell'].dropna().unique()

    logging.debug(f'cell_selection: {cell_selection}')

    # loop over individual cells for analysis
    for cell in cell_selection:                        
        logging.debug( 'Start calculating Data for cell '+ str(int(cell)))
        savedir_cell=os.path.join(savedir_tracking,'cells',str(int(cell)))
        os.makedirs(savedir_cell,exist_ok=True)

        Aux_along,Aux_across=interpolate_alongacross_mean(Aux,tracks,cell,dx=500,width=10000,z_coord='geopotential_height',height_level_borders=height_levels)
        iris.save(Aux_along,os.path.join(savedir_cell,f'Aux_along.nc'),zlib=True,complevel=4)
        iris.save(Aux_across,os.path.join(savedir_cell,f'Aux_across.nc'),zlib=True,complevel=4)
        
        logging.debug( 'hydrometeors number along and across calculated and saved to '+ savedir_cell)
