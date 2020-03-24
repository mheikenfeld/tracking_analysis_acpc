import iris
import logging
import os

def load_w(filenames,savedir,constraint_loading_period,load_function,compression=True):
    from iris.analysis import MAX
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
