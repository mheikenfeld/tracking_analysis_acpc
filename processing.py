import iris
import logging
import os

from tobac import plot_tracks_mask_field_loop
from tobac import plot_mask_cell_track_static_timeseries
from tobac import plot_mask_cell_track_3Dstatic,plot_mask_cell_track_2D3Dstatic
from tobac import cell_statistics

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
