def set_load_function(model):
    import ramscube
    import wrfcube
    import iris
    from mpdiag import calculate_wrf_mp_path,calculate_rams_mp_path,lump_processes
    #Define load functions for the two models and place them in dictionary
    #General load function    
    if model=='WRF':
        load_function=wrfcube.load
        def load_function_out(filename,variable_in):
            if variable_in is 'processes':
                processes=calculate_wrf_mp_path(filename,                                                    
                                          microphysics_scheme='morrison',processes='mass',
                                          signed=True,
                                          quantity='mixing ratio',constraint=None,
                                          debug_nproc=None)
                out=processes
            elif variable_in is 'processes_lumped':
                processes=calculate_wrf_mp_path(filename,                                                    
                                          microphysics_scheme='morrison',processes='mass',
                                          signed=True,
                                          quantity='mixing ratio',constraint=None,
                                          debug_nproc=None)
                processes_lumped=lump_processes(processes,microphysics_scheme='morrison',lumping='basic')
                out=processes_lumped
            elif variable_in is 'processes_lumped_detailed':
                processes=calculate_wrf_mp_path(filename,                                                    
                                          microphysics_scheme='morrison',processes='mass',
                                          signed=True,
                                          quantity='mixing ratio',constraint=None,
                                          debug_nproc=None)
                processes_lumped=lump_processes(processes,microphysics_scheme='morrison',lumping='detailed')
                out=processes_lumped
            elif variable_in is "hydrometeors":
                out=iris.cube.Cubelist()
                for hydrometeor in ['QCLOUD','QRAIN','QICE','QSNOW','QGRAUP']:
                    out.append(load_function(filename,hydrometeor))
            else:
                if variable_in=='w':
                    variable='W'
                elif variable_in=='NCLOUD':
                    variable='QNCLOUD'
                elif variable_in=='QCLOUD':
                    variable='QCLOUD'
                elif variable_in=='airmass':
                    variable='airmass'
                elif variable_in=='airmass_path':
                    variable='airmass_path'
                else: 
                    variable=variable_in
                out=load_function(filename,variable)
                out.rename(variable_in)
            return out
        return load_function_out
        
    if model=='RAMS':
        load_function=ramscube.load
        def load_function_out(filename,variable_in): 
            if variable_in is 'processes':
                processes=calculate_rams_mp_path(filename,                                                    
                                                 microphysics_scheme='RAMS',processes='mass',
                                                 signed=True,
                                                 quantity='mixing ratio',constraint=None,
                                                 debug_nproc=None)
                # if process rates are accumulated, divide by timestep (deemd ok for 1min data)
                time_coord=processes[0].coord('time')
                dt=(time_coord.units.num2date(time_coord.points[1])-time_coord.units.num2date(time_coord.points[0])).total_seconds()
                dt_coord=iris.coords.AuxCoord(dt,units='s')
                for i,process in enumerate(processes):
                    processes[i]=process/dt_coord
                    processes[i].rename(process.name())
                out=processes
            elif variable_in is 'processes_lumped':
                processes=calculate_rams_mp_path(filename,                                                    
                                                 microphysics_scheme='RAMS',processes='mass',
                                                 signed=True,
                                                 quantity='mixing ratio',constraint=None,
                                                 debug_nproc=None)
                processes_lumped=lump_processes(processes,microphysics_scheme='RAMS',lumping='basic')
                
                # if process rates are accumulated, divide by timestep (deemd ok for 1min data)
                time_coord=processes_lumped[0].coord('time')
                dt=(time_coord.units.num2date(time_coord.points[1])-time_coord.units.num2date(time_coord.points[0])).total_seconds()
                dt_coord=iris.coords.AuxCoord(dt,units='s')
                for i,process in enumerate(processes_lumped):
                    processes_lumped[i]=process/dt_coord
                    processes_lumped[i].rename(process.name())
                out=processes_lumped
                
            elif variable_in is 'processes_lumped_detailed':
                processes=calculate_rams_mp_path(filename,                                                    
                                                 microphysics_scheme='RAMS',processes='mass',
                                                 signed=True,
                                                 quantity='mixing ratio',constraint=None,
                                                 debug_nproc=None)
                processes_lumped=lump_processes(processes,microphysics_scheme='RAMS',lumping='detailed')
                
                # if process rates are accumulated, divide by timestep (deemd ok for 1min data)
                time_coord=processes_lumped[0].coord('time')
                dt=(time_coord.units.num2date(time_coord.points[1])-time_coord.units.num2date(time_coord.points[0])).total_seconds()
                dt_coord=iris.coords.AuxCoord(dt,units='s')
                for i,process in enumerate(processes_lumped):
                    processes_lumped[i]=process/dt_coord
                    processes_lumped[i].rename(process.name())
                out=processes_lumped
            elif variable_in is "hydrometeors":
                out=iris.cube.Cubelist()
                for hydrometeor in ['RCP','RDP','RRP','RPP','RSP','RAP','RGP','RHP']:
                    out.append(load_function(filename,hydrometeor))
            else:
                if variable_in=='w':
                    variable='WC'
                elif variable_in=='NCLOUD':
                    variable='CCP'
                elif variable_in=='QCLOUD':
                    variable='RCP'
                elif variable_in=='airmass':
                    variable='airmass'
                elif variable_in=='airmass_path':
                    variable='airmass_path'
                else: 
                    variable=variable_in
                out=load_function(filename,variable)
                out.rename(variable_in)

            return out

        return load_function_out
