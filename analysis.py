import logging
import iris
import os
import pandas as pd

from tobac import feature_detection_multithreshold
from tobac import segmentation_3D
from tobac import linking_trackpy
from tobac import get_spacings

from .cell_analysis import extract_cell_cubes_subset,extract_cell_cubes_subset_2D,interpolate_alongacross_mean


#######################
# Interpolation profiles:

def interpolation_mask(cubelist_in,savedir_tracking=None,cell_selection=None,constraint_time=None,height_levels=None,type_mask=None,output_name=None)
    '''
    function description
    '''
    #load tracks and mass
    tracks=pd.read_hdf(os.path.join(savedir_tracking,'Track.h5'),'table')
    mask=iris.load_cube(os.path.join(savedir_tracking,f'Mask_Segmentation_{type_mask}.nc'),'segmentation_mask')

    if cell_selection is None:
        cell_selection=tracks['cell'].dropna().unique()

    logging.debug(f'cell_selection: {cell_selection}')

    # loop over individual cells for analysis
    for cell in cell_selection:                        
        logging.debug(f'Start calculating Data for cell {int(cell)}')
        savedir_cell=os.path.join(savedir_tracking,'cells',str(int(cell)))
        os.makedirs(savedir_cell,exist_ok=True)

        cubelist_cell_sum, cubelist_cell_integrated,track_integrated=extract_cell_cubes_subset(cubelist_in=cubelist_in,
                                                                                               mask=mask,
                                                                                               track=tracks,
                                                                                               cell=cell,
                                                                                               z_coord='model_level_number',
                                                                                               height_levels=height_levels)
        iris.save(cubelist_cell_sum,os.path.join(savedir_cell,f'{output_name}_profile_{type_mask}.nc'),zlib=True,complevel=4)
        iris.save(cubelist_cell_integrated,os.path.join(savedir_cell,f'{output_name}_integrated_{type_mask}.nc'),zlib=True,complevel=4)
        track_processes_lumped_integrated.to_hdf(os.path.join(savedir_cell,'track_{output_name}_integrated_{type_mask}.h5'),'table')
        
        logging.debug( f'{output_name} calculated and saved to {savedir_cell}'')

def interpolation_processes(savedir_tracking=None,savedir_data=None,cell_selection=None,constraint_time=None,height_levels=None,type_mask=None):
    '''
    function description
    '''
    output_name='Processes_lumped'
    data_file=os.path.join(savedir_data,f'Data_Processes_lumped.nc')
    if os.path.exists(data_file):       
        logging.debug(f'start loading Processes_lumped form preprocessed file:')
        processes_lumped=iris.load(data_file).extract(constraint_time)        

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

    # Sum up process rates for individual cells
    processes_lumped_sum=iris.cube.CubeList()
    
    for cube in processes_lumped:
        cube_sum=cube*airmass
        cube_sum.rename(cube.name())
        cube_sum=add_geopotential_height(cube_sum)
        processes_lumped_sum.append(cube_sum)
                                
    interpolation_mask(cubelist_in=processes_lumped_sum,
                       savedir_tracking=savedir_tracking,
                       output_name=output_name,
                       cell_selection=None,
                       constraint_time=None,
                       height_levels=None,
                       type_mask=None)               

def interpolation_hydrometeors_mass(savedir_tracking=None,savedir_data=None,cell_selection=None,constraint_tracking_period=None,plotting_period=None,type_mask=None,height_level=None):
    '''
    function description
    '''

    output_name='Hydrometeors_mass'
    hydrometeors_file=os.path.join(savedir_data,'Data_Hydrometeors_mass.nc')
    if os.path.exists(hydrometeors_file):       
        logging.debug('start loading hydrometers form preprocessed file:')
        hydrometeors=iris.load(data_file).extract(constraint_tracking_period)        

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

    # Sum up process rates for individual cells
    hydrometeors_sum=iris.cube.CubeList()
    
    for cube in hydrometeors:
        cube_sum=cube*airmass
        cube_sum.rename(cube.name())
        cube_sum=add_geopotential_height(cube_sum)
        hydrometeors_sum.append(cube_sum)

    interpolation_mask(cubelist_in=hydrometeors_sum,
                       savedir_tracking=savedir_tracking,
                       output_name=output_name,
                       cell_selection=None,
                       constraint_time=None,
                       height_levels=height_levels,
                       type_mask=None)

def interpolation_hydrometeors_number(savedir_tracking=None,savedir_data=None,cell_selection=None,constraint_tracking_period=None,plotting_period=None,type_mask=None,height_level=None):
    '''
    function description
    '''

    output_name='Hydrometeors_number'
    hydrometeors_file=os.path.join(savedir_data,'Data_Hydrometeors_number.nc')
    if os.path.exists(hydrometeors_file):       
        logging.debug('start loading hydrometers form preprocessed file:')
        hydrometeors=iris.load(data_file).extract(constraint_tracking_period)        

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

    # Sum up process rates for individual cells
    hydrometeors_sum=iris.cube.CubeList()
    
    for cube in hydrometeors:
        cube_sum=cube*airmass
        cube_sum.rename(cube.name())
        cube_sum=add_geopotential_height(cube_sum)
        hydrometeors_sum.append(cube_sum)

    interpolation_mask(cubelist_in=hydrometeors_sum,
                       savedir_tracking=savedir_tracking,
                       output_name=output_name,
                       cell_selection=None,
                       constraint_time=None,
                       height_levels=height_levels,
                       type_mask=None)

def interpolation_mass(savedir_tracking=None,savedir_data=None,plotdir=None,cell_selection=None,constraint_time=None,type_mask=None):
    '''
    function description
    '''
    output_name='Mass'
    #load water contents and fix coordinate names:
    WC=iris.load(os.path.join(savedir_data,'Data_WC.nc')).extract(constraint_time)   
    logging.debug('loaded total water content from file')
    for cube in WC:
        for coord in cube.coords():
            coord.var_name=coord.name()

    # load airmass (used to calculate total values from mixing ratios) and fix coordinate names:
    airmass=iris.load_cube(os.path.join(savedir_data,'Data_Airmass.nc'),'airmass').extract(constraint_time) 
    for coord in airmass.coords():
        coord.var_name=coord.name()
    logging.debug('airmass loaded from file')
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
            
    interpolation_mask(cubelist_in=mass_sum,
                       savedir_tracking=savedir_tracking,
                       output_name=output_name,
                       cell_selection=None,
                       constraint_time=None,
                       height_levels=height_levels,
                       type_mask=type_mask)
   
                      
                      
################################
###Slices:
                      
def slices(cubelist_in,savedir_tracking=None,savedir_data=None,cell_selection=None,constraint_tracking_period=None,plotting_period=None,height_levels=None):
    '''
    function description
    '''
                      
    tracks=pd.read_hdf(os.path.join(savedir_tracking,'Track.h5'),'table')
    
    # loop over individual cells for analysis
    for cell in cell_selection:                        
        logging.debug( 'Start calculating Data for cell '+ str(int(cell)))
        savedir_cell=os.path.join(savedir_tracking,'cells',str(int(cell)))
        os.makedirs(savedir_cell,exist_ok=True)

        cubelist_along,cubelist_across=interpolate_alongacross_mean(cubelist_in,tracks,cell,dx=500,width=10000,z_coord='geopotential_height',height_level_borders=height_levels)
        iris.save(cubelist_along,os.path.join(savedir_cell,f'{output_name}_along.nc'),zlib=True,complevel=4)
        iris.save(cubelist_across,os.path.join(savedir_cell,f'{output_name}_across.nc'),zlib=True,complevel=4)
        
        
        logging.debug(f'{output_name} along and across calculated and saved to {savedir_cell}')

def slices_processes_lumped(savedir_tracking=None,savedir_data=None,cell_selection=None,constraint_time=None,plotting_period=None,height_levels=None):
    '''
    function description
    '''
    output_name='Processes_lumped'

    processes_file=os.path.join(savedir_data,'Data_Processes_lumped.nc')
    if os.path.exists(processes_file):       
        logging.debug('start loading processes form preprocessed file:')
        processes=iris.load(processes_file).extract(constraint_time)
        for cube in processes:
            cube=add_geopotential_height(cube)
            
    logging.debug('loaded processes')

    slices(processes,savedir_tracking=None,savedir_data=None,cell_selection=None,constraint_time=None,height_levels=height_levels,output_name=output_name)
                      
def slices_processes_all(savedir_tracking=None,savedir_data=None,cell_selection=None,constraint_tracking_period=None,height_levels=None):
    '''
    function description
    '''
    output_name='Processes_all'
    tracks=pd.read_hdf(os.path.join(savedir_tracking,'Track.h5'),'table')

    processes_file=os.path.join(savedir_data,'Data_Processes_all.nc')
    if os.path.exists(processes_file):       
        logging.debug('start loading processes form preprocessed file:')
        processes=iris.load(processes_file).extract(constraint_tracking_period)
        
        for cube in processes:
            cube=add_geopotential_height(cube)
            
    logging.debug('loaded processes')

    slices(processes,savedir_tracking=None,savedir_data=None,cell_selection=None,constraint_time=None,height_levels=height_levels,output_name=output_name)

def slices_hydrometeors_mass(savedir_tracking=None,savedir_data=None,cell_selection=None,constraint_tracking_period=None,height_levels=None):
    '''
    function description
    '''
    output_name='Hydrometeors_mass'

    tracks=pd.read_hdf(os.path.join(savedir_tracking,'Track.h5'),'table')

    hydrometeors_file=os.path.join(savedir_data,'Data_hydrometeor_mass.nc')
    if os.path.exists(hydrometeors_file):       
        logging.debug('start loading hydrometers form preprocessed file:')
        hydrometeors=iris.load(hydrometeors_file).extract(constraint_tracking_period)
        
        for cube in hydrometeors:
            cube=add_geopotential_height(cube)
            
    logging.debug('loaded hydrometeors')

    slices(hydrometeors,savedir_tracking=None,savedir_data=None,cell_selection=None,constraint_time=None,height_levels=height_levels,output_name=output_name)

def slices_hydrometeors_number(savedir_tracking=None,savedir_data=None,cell_selection=None,constraint_tracking_period=None, height_levels=None):
    '''
    function description
    '''
    output_name='Hydrometeors_number'
    tracks=pd.read_hdf(os.path.join(savedir_tracking,'Track.h5'),'table')

    hydrometeors_file=os.path.join(savedir_data,'Data_hydrometeor_number.nc')
    if os.path.exists(hydrometeors_file):       
        logging.debug('start loading hydrometers form preprocessed file:')
        hydrometeors=iris.load(hydrometeors_file).extract(constraint_tracking_period)        

        for cube in hydrometeors:
            cube=add_geopotential_height(cube)

    logging.debug('loaded hydrometeors')

    slices(hydrometeors,savedir_tracking=None,savedir_data=None,cell_selection=None,constraint_time=None,height_levels=height_levels,output_name=output_name)

def slices_aux(savedir_tracking=None,savedir_data=None,cell_selection=None,constraint_tracking_period=None, height_levels=None):
    '''
    function description
    '''
    output_name='Aux'

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

    logging.debug(Aux,'loaded hydrometeors')

    slices(Aux, savedir_tracking=None, savedir_data=None,cell_selection=None,constraint_time=None, height_levels=height_levels)

def calculate_profiles_w(savedir_tracking=None,savedir_data=None,plotdir=None):
    '''
    function description
    '''
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


def interpolation_precipitation(savedir_tracking=None,savedir_data=None,cell_selection=None,constraint_tracking_period=None,plotting_period=None,masking='w'):
    '''
    function description
    '''
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

def profiles_cdnc(savedir_tracking=None,savedir_data=None,filenames=None,load_function=None):
    '''
    function description
    '''

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
        
def calculate_profiles_w(savedir_tracking=None,savedir_data=None,plotdir=None):
    '''
    function description
    '''
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
                      
                      
# def interpolation_processes_w_TWC(savedir_tracking=None,savedir_data=None,cell_selection=None,constraint_tracking_period=None,plotting_period=None):
#     logging.debug('Start calculations for individual cells')
        
#     logging.debug('loading tracks and Mask')

#     tracks=pd.read_hdf(os.path.join(savedir_tracking,'Track.h5'),'table')
#     mask=iris.load_cube(os.path.join(savedir_tracking,'Mask_Segmentation_w_TWC.nc'),'segmentation_mask')

#     processes_file=os.path.join(savedir_data,'Data_Processes_lumped.nc')
#     if os.path.exists(processes_file):       
#         logging.debug('start loading processes form preprocessed file:')
#         processes_lumped=iris.load(processes_file).extract(constraint_tracking_period)        
# #        else:#        
# #            logging.debug(' start process rates from data files')
# #            processes_lumped=load_function(filenames,'processes_lumped').extract(constraint_tracking_period)                 

#     for cube in processes_lumped:
#         for coord in cube.coords():
#             coord.var_name=coord.name()
#             coord.coord_system=None

#     logging.debug('loaded processes')

#     logging.debug('start loading auxillary variables')
    
    
#     # load airmass (used to calculate total values from mixing ratios)
#     airmass=iris.load_cube(os.path.join(savedir_data,'Data_Airmass.nc'),'airmass').extract(constraint_tracking_period) 
#     for coord in airmass.coords():
#         coord.var_name=coord.name()
#     logging.debug('airmass loaded from file')


#     height_levels=np.arange(0,20001,1000.)
                    
#     # Sum up process rates for individual cells
#     processes_lumped_sum=iris.cube.CubeList()
    
#     for cube in processes_lumped:
#         cube_sum=cube*airmass
#         cube_sum.rename(cube.name())
#         cube_sum=add_geopotential_height(cube_sum)
#         processes_lumped_sum.append(cube_sum)
# #            
#     if cell_selection is None:
#         cell_selection=tracks['cell'].dropna().unique()

#     logging.debug(f'cell_selection: {cell_selection}')

#     # loop over individual cells for analysis
#     for cell in cell_selection:                        
#         logging.debug( 'Start calculating Data for cell '+ str(int(cell)))
#         savedir_cell=os.path.join(savedir_tracking,'cells',str(int(cell)))
#         os.makedirs(savedir_cell,exist_ok=True)

# #
#         processes_lumped_cell_sum, processes_lumped_cell_integrated,track_processes_lumped_integrated=extract_cell_cubes_subset(cubelist_in=processes_lumped_sum,
#                                                                                                                              mask=mask,
#                                                                                                                              track=tracks,
#                                                                                                                              cell=cell,
#                                                                                                                              z_coord='model_level_number',
#                                                                                                                              height_levels=height_levels)
        
#         iris.save(processes_lumped_cell_sum,os.path.join(savedir_cell,'Processes_lumped_profile_w_TWC.nc'),zlib=True,complevel=4)
#         iris.save(processes_lumped_cell_integrated,os.path.join(savedir_cell,'Processes_lumped_integrated_w_TWC.nc'),zlib=True,complevel=4)
#         track_processes_lumped_integrated.to_hdf(os.path.join(savedir_cell,'track_processes_lumped_integrated_w_TWC.h5'),'table')
        
#         logging.debug( 'processes calculated and saved to '+ savedir_cell)



# def interpolation_processes_TWC(savedir_tracking=None,savedir_data=None,cell_selection=None,constraint_tracking_period=None):
#     logging.debug('Start calculations for individual cells')
        
#     logging.debug('loading tracks and Mask')



#     tracks=pd.read_hdf(os.path.join(savedir_tracking,'Track.h5'),'table')
#     mask=iris.load_cube(os.path.join(savedir_tracking,'Mask_Segmentation_TWC.nc'),'segmentation_mask')

#     processes_file=os.path.join(savedir_data,'Data_Processes_lumped.nc')
#     if os.path.exists(processes_file):       
#         logging.debug('start loading processes form preprocessed file:')
#         processes_lumped=iris.load(processes_file).extract(constraint_tracking_period)        
# #            else:#        
# #                logging.debug(' start process rates from data files')
# #                processes_lumped=load_function(filenames,'processes_lumped').extract(constraint_tracking_period)                 

#     for cube in processes_lumped:
#         for coord in cube.coords():
#             coord.var_name=coord.name()
#             coord.coord_system=None

#     logging.debug('loaded processes')

#     logging.debug('start loading auxillary variables')
    
    
#     # load airmass (used to calculate total values from mixing ratios)
#     airmass=iris.load_cube(os.path.join(savedir_data,'Data_Airmass.nc'),'airmass').extract(constraint_tracking_period) 
#     for coord in airmass.coords():
#         coord.var_name=coord.name()
#     logging.debug('airmass loaded from file')



#     height_levels=np.arange(0,20001,1000.)
                    
#     # Sum up process rates for individual cells
#     processes_lumped_sum=iris.cube.CubeList()

#     height_levels=np.arange(0,20001,1000.)
                    
#     # Sum up process rates for individual cells
#     processes_lumped_sum=iris.cube.CubeList()
    
#     for cube in processes_lumped:
#         cube_sum=cube*airmass
#         cube_sum.rename(cube.name())
#         cube_sum=add_geopotential_height(cube_sum)
#         processes_lumped_sum.append(cube_sum)
# #            
#     if cell_selection is None:
#         cell_selection=tracks['cell'].dropna().unique()

#     logging.debug(f'cell_selection: {cell_selection}')


#     # loop over individual cells for analysis
#     for cell in cell_selection:                        
#         logging.debug( 'Start calculating Data for cell '+ str(int(cell)))
#         savedir_cell=os.path.join(savedir_tracking,'cells',str(int(cell)))
#         os.makedirs(savedir_cell,exist_ok=True)

#         processes_lumped_cell_sum, processes_lumped_cell_integrated,track_processes_lumped_integrated=extract_cell_cubes_subset(cubelist_in=processes_lumped_sum,
#                                                                                                                              mask=mask,
#                                                                                                                              track=tracks,
#                                                                                                                              cell=cell,
#                                                                                                                              z_coord='model_level_number',
#                                                                                                                              height_levels=height_levels)
        
#         iris.save(processes_lumped_cell_sum,os.path.join(savedir_cell,'Processes_lumped_profile_TWC.nc'),zlib=True,complevel=4)
#         iris.save(processes_lumped_cell_integrated,os.path.join(savedir_cell,'Processes_lumped_integrated_TWC.nc'),zlib=True,complevel=4)
#         track_processes_lumped_integrated.to_hdf(os.path.join(savedir_cell,'track_processes_lumped_integrated_TWC.h5'),'table')
        
#         logging.debug( 'processes calculated and saved to '+ savedir_cell)

# def interpolation_mass_TWC(savedir_tracking=None,savedir_data=None,plotdir=None,cell_selection=None,constraint_tracking_period=None):
#     logging.debug('loading tracks and Mask')

#     track=pd.read_hdf(os.path.join(savedir_tracking,'Track.h5'),'table')
#     mask=iris.load_cube(os.path.join(savedir_tracking,'Mask_Segmentation_TWC.nc'),'segmentation_mask')
    
#     #load water contents and fix coordinate names:
#     WC=iris.load(os.path.join(savedir_data,'Data_WC.nc')).extract(constraint_tracking_period)   
#     logging.debug('loaded total water content from file')
#     for cube in WC:
#         for coord in cube.coords():
#             coord.var_name=coord.name()

#     # load airmass (used to calculate total values from mixing ratios) and fix coordinate names:
#     airmass=iris.load_cube(os.path.join(savedir_data,'Data_Airmass.nc'),'airmass').extract(constraint_tracking_period) 
#     for coord in airmass.coords():
#         coord.var_name=coord.name()
#     logging.debug('airmass loaded from file')

#     if cell_selection is None:
#         cell_selection=track['cell'].dropna().unique()

#     for cell in cell_selection:                        
#         savedir_cell=os.path.join(savedir_tracking,'cells',str(int(cell)))
#         logging.debug( 'Start calculating Data for cell '+ str(int(cell)))
#         savedir_cell=os.path.join(savedir_tracking,'cells',str(int(cell)))
#         os.makedirs(savedir_cell,exist_ok=True)
#         tracks=pd.read_hdf(os.path.join(savedir_tracking,'Track.h5'),'table')
#         mask=iris.load_cube(os.path.join(savedir_tracking,'Mask_Segmentation_TWC.nc'),'segmentation_mask')


#         height_levels=np.arange(0,20001,1000.)
        
#         mass_sum=iris.cube.CubeList()

#         cube=WC.extract_strict('TWC')

#         for coord in airmass.coords():
#             print(coord)
#         for coord in cube.coords():
#             print(coord)

#         cube_sum=cube*airmass
#         cube_sum.rename('total_water_mass')
#         cube_sum=add_geopotential_height(cube_sum)                
#         mass_sum.append(cube_sum)

#         cube=WC.extract_strict('LWC')
#         cube_sum=cube*airmass
#         cube_sum.rename('liquid_water_mass')
#         cube_sum=add_geopotential_height(cube_sum)                
#         mass_sum.append(cube_sum)

#         cube=WC.extract_strict('IWC')
#         cube_sum=cube*airmass
#         cube_sum.rename('frozen_water_mass')
#         cube_sum=add_geopotential_height(cube_sum)                
#         mass_sum.append(cube_sum)
            
#         mass_cell_sum, mass_cell_integrated,track_mass_integrated=extract_cell_cubes_subset(cubelist_in=mass_sum,
#                                                                                             mask=mask,
#                                                                                             track=tracks,
#                                                                                             cell=cell,
#                                                                                             z_coord='model_level_number',
#                                                                                             height_levels=height_levels)

#         iris.save(mass_cell_sum,os.path.join(savedir_cell,'Mass_profile_TWC.nc'),zlib=True,complevel=4)
#         iris.save(mass_cell_integrated,os.path.join(savedir_cell,'Mass_integrated_TWC.nc'),zlib=True,complevel=4)
#         track_mass_integrated.to_hdf(os.path.join(savedir_cell,'track_mass_integrated_TWC.h5'),'table')
#         logging.debug( 'mass calculated and saved to '+ savedir_cell)

# def interpolation_hydrometeors_TWC(savedir_tracking=None,savedir_data=None,cell_selection=None,constraint_tracking_period=None,plotting_period=None):
#     logging.debug('Start calculations for individual cells')
        
#     logging.debug('loading tracks and Mask')

#     tracks=pd.read_hdf(os.path.join(savedir_tracking,'Track.h5'),'table')
#     mask=iris.load_cube(os.path.join(savedir_tracking,'Mask_Segmentation_TWC.nc'),'segmentation_mask')

#     hydrometeors_file=os.path.join(savedir_data,'Data_hydrometeor_mass.nc')
#     if os.path.exists(hydrometeors_file):       
#         logging.debug('start loading hydrometers form preprocessed file:')
#         hydrometeors=iris.load(hydrometeors_file).extract(constraint_tracking_period)        
# #            else:#        
# #                logging.debug(' start process rates from data files')
# #                processes_lumped=load_function(filenames,'processes_lumped').extract(constraint_tracking_period)                 

#     for cube in hydrometeors:
#         for coord in cube.coords():
#             coord.var_name=coord.name()
#             coord.coord_system=None

#     logging.debug('loaded hydrometeors')

#     logging.debug('start loading auxillary variables')
    
    
#     # load airmass (used to calculate total values from mixing ratios)
#     airmass=iris.load_cube(os.path.join(savedir_data,'Data_Airmass.nc'),'airmass').extract(constraint_tracking_period) 
#     for coord in airmass.coords():
#         coord.var_name=coord.name()
#     logging.debug('airmass loaded from file')



#     height_levels=np.arange(0,20001,1000.)
                    
#     # Sum up process rates for individual cells
#     hydrometeors_sum=iris.cube.CubeList()
    
#     for cube in hydrometeors:
#         cube_sum=cube*airmass
#         cube_sum.rename(cube.name())
#         cube_sum=add_geopotential_height(cube_sum)
#         hydrometeors_sum.append(cube_sum)
# #            
#     if cell_selection is None:
#         cell_selection=tracks['cell'].dropna().unique()

#     logging.debug(f'cell_selection: {cell_selection}')
#     # loop over individual cells for analysis
#     for cell in cell_selection:                        
#         logging.debug( 'Start calculating Data for cell '+ str(int(cell)))
#         savedir_cell=os.path.join(savedir_tracking,'cells',str(int(cell)))
#         os.makedirs(savedir_cell,exist_ok=True)

#         hydrometeors_cell_sum, hydrometeors_cell_integrated,track_hydrometeors_integrated=extract_cell_cubes_subset(cubelist_in=hydrometeors_sum,
#                                                                                                                              mask=mask,
#                                                                                                                              track=tracks,
#                                                                                                                              cell=cell,
#                                                                                                                              z_coord='model_level_number',
#                                                                                                                              height_levels=height_levels)
        
#         iris.save(hydrometeors_cell_sum,os.path.join(savedir_cell,'Hydrometeors_profile_TWC.nc'),zlib=True,complevel=4)
#         iris.save(hydrometeors_cell_integrated,os.path.join(savedir_cell,'Hydrometeors_integrated_TWC.nc'),zlib=True,complevel=4)
#         track_hydrometeors_integrated.to_hdf(os.path.join(savedir_cell,'track_hydrometeors_integrated_TWC.h5'),'table')
        
#         logging.debug( 'hydrometeors calculated and saved to '+ savedir_cell)

