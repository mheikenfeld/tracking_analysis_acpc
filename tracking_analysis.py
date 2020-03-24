import matplotlib.pyplot as plt
import logging

import iris
from iris.analysis import MIN,MAX,MEAN,MEDIAN,PERCENTILE
import os
import pandas as pd
import numpy as np

from matplotlib import colors
import dask
from multiprocessing.pool import ThreadPool

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
        run_feature_detection(loaddir=os.path.join(top_savedir_data),savedir=os.path.join(top_savedir_tracking),
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

#    if (('plot_mask_hydrometeors_w_TWC'  not in mode) and ('interpolation_hydrometeors_w_TWC' in mode)):
#        interpolation_hydrometeors_w_TWC(savedir_tracking=top_savedir_tracking,savedir_data=top_savedir_data,
#                                    cell_selection=cell_selection,constraint_tracking_period=constraint_tracking_period,plotting_period=plotting_period)
#
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