import matplotlib
matplotlib.use('Agg')

import os
import matplotlib.pyplot as plt
from matplotlib import dates
import pandas as pd 
import numpy as np
from cloudtrack import plot_lifetime_histogram
import logging
import sys
import datetime
sys.path.append('../')

from html_animate import make_animation

# Setup output log file:
if 'LSB_JOBID' in os.environ.keys():
    Job_id=str(os.environ['LSB_JOBID'])
else:
    Job_id='2' 

rootLogger = logging.getLogger()
rootLogger.setLevel(logging.DEBUG)
#rootLogger.setLevel(logging.INFO)
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
log_filename='log_'+Job_id+'.log'
fileHandler=logging.FileHandler(log_filename)
fileHandler.setFormatter(formatter)
rootLogger.addHandler(fileHandler)


models=[]
models.append('RAMS')
models.append('WRF')

cases=[]
cases.append('CLN')
cases.append('POL')

acpc_workspace='/group_workspaces/jasmin2/acpc'

savedir_top=os.path.join(acpc_workspace,'Analysis/mheiken/Tracking','Save_1min')
plotdir=os.path.join(acpc_workspace,'public/mheiken/Tracking','Plots_1min')


#version='multithreshold'
#version='multithreshold_test_unfiltered'
version='version_1_0'

Tracks={}
for model in models:
    Tracks[model]={}
    for case in cases:
                savedir_i=os.path.join(savedir_top,model,case,version)
                Track=pd.read_hdf(os.path.join(savedir_i,'Track_TWC.hdf'),'table') 
                #Restrict Tracks to those that are not affected/cut by time window or domain boundaries:
                if model is 'WRF':
                    x_min=5000
                    x_max=245000
                    y_min=5000
                    y_max=245000
                if model is 'RAMS':
                    x_min=-120000
                    x_max=120000
                    y_min=-120000
                    y_max=120000
                time_min=datetime.datetime(2013,6,19,21,10,0,0)
                time_max=datetime.datetime(2013,6,19,23,50,0,0)
                #Find all cells that are close to boundary either in time or space
                cells_remove=Track[ (Track['projection_x_coordinate']<x_min) 
                                 |(Track['projection_x_coordinate']>x_max)
                                 |(Track['projection_y_coordinate']<y_min) 
                                 |(Track['projection_y_coordinate']>y_max)
                                 |(Track['time']<time_min)
                                 |(Track['time']>time_max) 
                                 ]['particle'].unique()
                Track_filtered=Track[~Track['particle'].isin(cells_remove)]
                # Add filtered DataFrame to dict for each simulation
                Tracks[model][case]=Track_filtered
                

logging.info("data loaded from file")


#colors=['blue','green','gold','darkorange','red','darkred']
#linestyles=['-','--',':','-.']
    
color={}
color['WRF']='darkblue'
color['RAMS']='darkred'
ls={}
ls['CLN']='--'
ls['POL']='-'
pos={}
pos['WRF']={}
pos['RAMS']={}


plot_dir=os.path.join(plotdir,'Statistics',version)
os.makedirs(plot_dir,exist_ok=True)
width=3
bin_min=0
bin_max=120

bin_edges=np.arange(bin_min,bin_max,width)
#fig4,ax4=plt.subplots(nrows=1,ncols=1)        
fig3,ax3=plt.subplots(nrows=1,ncols=1,figsize=(15/2.54,15/2.54))
for model in models:
    for case in cases:
        Track_i=Tracks[model][case]
        os.makedirs(plot_dir,exist_ok=True)
        plot_lifetime_histogram(Track_i,axes=ax3,bin_edges=bin_edges,
                                color=color[model],ls=ls[case],
                                label=model+' '+case)
ax3.legend(fontsize=6)
ax3.set_xlabel('cell lifetime (min)')
ax3.set_ylabel('counts')
filename=os.path.join(plot_dir,'Histogram_lifetime_line.png')
fig3.savefig(filename,dpi=300)

logging.info("histogram lifetimes plotted")


#hours = dates.MinuteLocator(interval=15)
formatter = dates.DateFormatter('%M')

for model in models:
    for case in cases:
        Track_i=Tracks[model][case]
        plot_dir=os.path.join(plotdir,model,case,version,'Lifetime')
        os.makedirs(plot_dir,exist_ok=True)
        plot_dir_cells=os.path.join(plot_dir,'Cells')
        os.makedirs(plot_dir_cells,exist_ok=True)
        fig4,ax4=plt.subplots(nrows=1,ncols=1,figsize=(15/2.54,10/2.54))
        Track_particle=Track_i.groupby('particle')
        for particle, Track in Track_particle:
            fig5,ax5=plt.subplots(nrows=1,ncols=1,figsize=(15/2.54,15/2.54))
            Track.plot(ax=ax4,x='time_cell',y='ncells',ls='-',lw=0.3,legend=False)
            Track.plot(ax=ax5,x='time_cell',y='ncells',ls='-',lw=1,legend=False,color='navy')
            ax5.set_xlabel('cell lifetime (min)')
            ax5.set_ylabel('ncells')
            ax5.set_xlim([0,2*1e9*3600])
            ax4.set_ylim([0,max(10,1.1*Track['ncells'].max())])
            ax5.set_xticks(1e9*3600*np.arange(0,2,0.25))
#                    ax5.xaxis.set_major_locator(hours)
#                    ax5.xaxis.set_major_formatter(formatter)
            ax5.set_title('cell: '+str(particle))
            filename=os.path.join(plot_dir_cells,'lifetime_'+ model+'_'+case+'_'+str(int(particle))+'.png')
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
        filename=os.path.join(plot_dir,'lifetimes_'+ model+'_'+case+'_all'+'.png')
        fig4.savefig(filename,dpi=300)
        logging.info("lifetimes ncell plotted for "+ model+' ' + case)

        make_animation(input_dir=plot_dir_cells,
                       output=os.path.join(plot_dir,'Animation_Lifetime_cells.html'))



logging.info('done')