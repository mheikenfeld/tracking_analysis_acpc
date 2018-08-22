def set_load_function(model):
    import ramscube
    import wrfcube
    from wrfmpdiag import calculate_wrf_mp_path,calculate_rams_mp_path,lump_processes
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
                out=processes

            elif variable_in is 'processes_lumped':
                processes=calculate_rams_mp_path(filename,                                                    
                                                 microphysics_scheme='RAMS',processes='mass',
                                                 signed=True,
                                                 quantity='mixing ratio',constraint=None,
                                                 debug_nproc=None)
                processes_lumped=lump_processes(processes,microphysics_scheme='RAMS',lumping='basic')
                out=processes_lumped

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