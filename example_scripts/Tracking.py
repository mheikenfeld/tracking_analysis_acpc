import logging
import warnings
import os
import pandas as pd
import sys
import glob
import datetime

sys.path.append('../')
from tracking_analysis_acpc import tracking_analysis
from tracking_analysis_acpc import set_load_function

#Disable some annoying warnings:
warnings.filterwarnings('ignore', category=UserWarning, append=True)
warnings.filterwarnings('ignore', category=RuntimeWarning, append=True)
warnings.filterwarnings('ignore',category=pd.io.pytables.PerformanceWarning)
warnings.filterwarnings('ignore', category=FutureWarning, append=True)

import sys
model=sys.argv[1]
case=sys.argv[2]


#Job_id=str(v_max)+'_'+str(w_max_threshold)

# Setup output log file:
if 'LSB_JOBID' in os.environ.keys():
    Job_id=str(os.environ['LSB_JOBID'])
else:
    Job_id=str(model)+str(case)
    

rootLogger = logging.getLogger()
rootLogger.setLevel(logging.DEBUG)
#rootLogger.setLevel(logging.INFO)
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
log_filename='log_'+Job_id+'.log'
fileHandler=logging.FileHandler(log_filename)
fileHandler.setFormatter(formatter)
rootLogger.addHandler(fileHandler)

if 'LSB_DJOB_NUMPROC'  in os.environ.keys():
    n_core=int(os.environ['LSB_DJOB_NUMPROC'])
else:
    n_core=4
 
logging.info('number of processors used:'+ str(n_core))
logging.info('Job Id:'+ str(Job_id))

        

acpc_workspace='/group_workspaces/jasmin2/acpc'

#%%
#Set path to data
directory={}
filename={}
directory['RAMS']={}
directory['RAMS']['CLN']=os.path.join(acpc_workspace,'houston_deep_convection/RAMS_CSU/CLN/x.out.1m')
directory['RAMS']['POL']=os.path.join(acpc_workspace,'houston_deep_convection/RAMS_CSU/POL/x.out.1m')

filename['RAMS']='a-A-2013-06-19*g3.h5'
#filename['RAMS']='a-A-2013-06-19-21*g3.h5'

directory['WRF']={}
directory['WRF']['CLN']=os.path.join(acpc_workspace,'houston_deep_convection/WRF_Oxford/CLN/1min/d03')
directory['WRF']['POL']=os.path.join(acpc_workspace,'houston_deep_convection/WRF_Oxford/POL/1min/d03')

filename['WRF']='wrfout_d03_2013-06-19*'
#filename['WRF']='wrfout_d03_2013-06-19_21*'


# Setup which parts of the processing to perform:  
mode=[]  
#mode.append('load_data')    
mode.append('tracking')    
mode.append('segmentation')    
mode.append('plot')    
#mode.append('statistics')    
#mode.append('profiles')
#mode.append('plot_profiles')
mode.append('plot_mask')
#mode.append('interpolation_mass')
#mode.append('plot_mask_mass')
mode.append('interpolation_processes')
mode.append('plot_mask_processes')
#mode.append('plot_lifetime')    


#Tracking setup:
parameters_tracking={}
parameters_tracking['method_detection']='threshold_multi'
parameters_tracking['method_linking']='predict'
parameters_tracking['position_threshold']='weighted_diff'
parameters_tracking['sigma_threshold']=0.5
parameters_tracking['adaptive_stop']=0.2
parameters_tracking['adaptive_step']=0.95
parameters_tracking['n_erosion_threshold']=4
parameters_tracking['stubs']=4
parameters_tracking['min_num']=4
parameters_tracking['extrapolate']=0
parameters_tracking['order']=1
parameters_tracking['v_max']=10
parameters_tracking['subnetwork_size']=100

parameters_tracking['threshold']=[5,10] #m/s

parameters_segmentation_TWC={}
parameters_segmentation_w={}
parameters_segmentation_w['method']='watershed'
parameters_segmentation_TWC['method']='watershed'


#Set output directories
top_savedir_data=os.path.join(acpc_workspace,'Analysis/mheiken/Tracking','Save_1min','data',model,case)
#
#top_savedir_tracking=os.path.join(acpc_workspace,'Analysis/mheiken/Tracking','Save_1min',model,case,'multithreshold')
#top_plotdir=os.path.join(acpc_workspace,'public/mheiken/Tracking','Plots_1min',model,case,'multithreshold')
#parameters_tracking['memory']=1

#top_savedir_tracking=os.path.join(acpc_workspace,'Analysis/mheiken/Tracking','Save_1min',model,case,'multithreshold_test')
#top_plotdir=os.path.join(acpc_workspace,'public/mheiken/Tracking','Plots_1min',model,case,'multithreshold_test')
#parameters_tracking['memory']=1

#top_savedir_tracking=os.path.join(acpc_workspace,'Analysis/mheiken/Tracking','Save_1min',model,case,'multithreshold_test_unfiltered')
#top_plotdir=os.path.join(acpc_workspace,'public/mheiken/Tracking','Plots_1min',model,case,'multithreshold_test_unfiltered')
#parameters_tracking['memory']=1
#parameters_tracking['min_num']=0



#top_savedir_tracking=os.path.join(acpc_workspace,'Analysis/mheiken/Tracking','Save_1min',model,case,'version_1_0')
#top_plotdir=os.path.join(acpc_workspace,'public/mheiken/Tracking','Plots_1min',model,case,'version_1_0')
#parameters_tracking['memory']=1
#parameters_tracking['min_num']=0
#parameters_tracking['min_distance']=2000
#
##
#top_savedir_tracking=os.path.join(acpc_workspace,'Analysis/mheiken/Tracking','Save_1min','version_1_1',model,case)
#top_plotdir=os.path.join(acpc_workspace,'public/mheiken/Tracking','Plots_1min','version_1_1',model,case)
#parameters_tracking['memory']=1
#parameters_tracking['min_num']=0
#parameters_tracking['min_distance']=5000
#parameters_segmentation_TWC['threshold']=1e-3  # kg/kg mixing ratio
#parameters_segmentation_w['threshold']=5 # m/s
#
####
#top_savedir_tracking=os.path.join(acpc_workspace,'Analysis/mheiken/Tracking','Save_1min','version_1_2',model,case)
#top_plotdir=os.path.join(acpc_workspace,'public/mheiken/Tracking','Plots_1min','version_1_2',model,case)
#parameters_tracking['memory']=1
#parameters_tracking['min_num']=0
#parameters_tracking['min_distance']=5000
#parameters_segmentation_TWC['threshold']=2e-3  # kg/kg mixing ratio
#parameters_segmentation_w['threshold']=5 # m/s
#parameters_tracking['n_erosion_threshold']=0

#top_savedir_tracking=os.path.join(acpc_workspace,'Analysis/mheiken/Tracking','Save_1min','version_1_3',model,case)
#top_plotdir=os.path.join(acpc_workspace,'public/mheiken/Tracking','Plots_1min','version_1_3',model,case)
#parameters_tracking['memory']=1
#parameters_tracking['min_num']=0
#parameters_tracking['min_distance']=2000
#parameters_segmentation_TWC['threshold']=2e-3  # kg/kg mixing ratio
#parameters_segmentation_w['threshold']=5 # m/s
#parameters_tracking['n_erosion_threshold']=0
####
#top_savedir_tracking=os.path.join(acpc_workspace,'Analysis/mheiken/Tracking','Save_1min','version_1_4',model,case)
#top_plotdir=os.path.join(acpc_workspace,'public/mheiken/Tracking','Plots_1min','version_1_4',model,case)
#parameters_tracking['memory']=1
#parameters_tracking['min_num']=0
#parameters_tracking['min_distance']=0
#parameters_segmentation_TWC['threshold']=2e-3  # kg/kg mixing ratio
#parameters_segmentation_w['threshold']=5 # m/s
#parameters_tracking['n_erosion_threshold']=0
#]
top_savedir_tracking=os.path.join(acpc_workspace,'Analysis/mheiken/Tracking','Save_1min','version_1_5',model,case)
top_plotdir=os.path.join(acpc_workspace,'public/mheiken/Tracking','Plots_1min','version_1_5',model,case)
parameters_tracking['memory']=1
parameters_tracking['min_num']=0
parameters_tracking['min_distance']=0
parameters_segmentation_TWC['threshold']=2e-3  # kg/kg mixing ratio
parameters_segmentation_w['threshold']=3 # m/s
parameters_tracking['n_erosion_threshold']=0
parameters_tracking['threshold']=[3,5,10] #m/s

#top_savedir_tracking=os.path.join(acpc_workspace,'Analysis/mheiken/Tracking','Save_1min',model,case,'multithreshold_30')
#top_plotdir=os.path.join(acpc_workspace,'public/mheiken/Tracking','Plots_1min',model,case,'multithreshold_30')
#parameters_tracking['memory']=0
#



plotting_parameters={}
plotting_parameters['axis_extent']=[-96.3,-93.9,28.4,30.5]

load_function=set_load_function(model)
filename=filename[model]
directory=directory[model][case]
imput_str=os.path.join(directory,filename)
logging.debug(imput_str)
filenames=glob.glob(imput_str)

loading_period=[datetime.datetime(2013,6,19,21,0),datetime.datetime(2013,6,20,0,0)]
tracking_period=[datetime.datetime(2013,6,19,21,0),datetime.datetime(2013,6,20,0,0)]
plotting_period=[datetime.datetime(2013,6,19,21,0),datetime.datetime(2013,6,20,0,0)]

tracking_analysis(filenames=filenames,
                  mode=mode,
                  top_savedir_data=top_savedir_data,
                  top_savedir_tracking=top_savedir_tracking,
                  top_plotdir=top_plotdir,
                  loading_period=loading_period,
                  tracking_period=tracking_period,
                  parameters_tracking=parameters_tracking,
                  parameters_segmentation_w=parameters_segmentation_w,
                  parameters_segmentation_TWC=parameters_segmentation_TWC,
                  load_function=load_function,
                  n_core=n_core,
                  plotting_parameters=plotting_parameters,
                  plotting_period=plotting_period
                  )
    
logging.debug('done')


