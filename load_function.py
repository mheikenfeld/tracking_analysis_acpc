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
            elif variable_in is 'processes_volume':
                processes=calculate_wrf_mp_path(filename,                                                    
                                          microphysics_scheme='morrison',processes='mass',
                                          signed=True,
                                          quantity='volume',constraint=None,
                                          debug_nproc=None)
                out=processes
            elif variable_in is 'processes_all':
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
            elif variable_in is "hydrometeor_mass":
                out=iris.cube.CubeList()
                for hydrometeor in ['QCLOUD','QRAIN','QICE','QSNOW','QGRAUP']:
                    out.append(load_function(filename,hydrometeor))
            elif variable_in is "hydrometeor_number":
                out=iris.cube.CubeList()
                for hydrometeor in ['QNCLOUD','QNRAIN','QNICE','QNSNOW','QNGRAUPEL']:
                    out.append(load_function(filename,hydrometeor))

            else:
                if variable_in=='w':
                    variable='W'
                elif variable_in=='u':
                    variable='U'
                elif variable_in=='v':
                    variable='V'
                elif variable_in=='w_unstaggered':
                    variable='w_unstaggered'
                elif variable_in=='u_unstaggered':
                    variable='u_unstaggered'
                elif variable_in=='v_unstaggered':
                    variable='v_unstaggered'
                elif variable_in=='NCLOUD':
                    variable='QNCLOUD'
                elif variable_in=='QCLOUD':
                    variable='QCLOUD'
                elif variable_in=='airmass':
                    variable='airmass'
                elif variable_in=='airmass_path':
                    variable='airmass_path'
                elif variable_in=='air_density':
                    variable='air_density'
                elif variable_in=='OLR':
                    variable='OLR'
                elif variable_in=='air_temperature':
                    variable='air_temperature'
                elif variable_in=='potential_temperature':
                    variable='potential_temperature'
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
                                                 microphysics_scheme='rams',processes='mass_grouped',
                                                 signed=True,
                                                 quantity='mixing ratio',constraint=None,
                                                 debug_nproc=None)
                out=processes
            elif variable_in is 'processes_volume':
                processes=calculate_rams_mp_path(filename,                                                    
                                                 microphysics_scheme='rams',processes='mass_grouped',
                                                 signed=True,
                                                 quantity='volume',constraint=None,
                                                 debug_nproc=None)
                out=processes
            elif variable_in is 'processes_all':
                processes=calculate_rams_mp_path(filename,                                                    
                                                 microphysics_scheme='rams',processes='mass',
                                                 signed=True,
                                                 quantity='mixing ratio',constraint=None,
                                                 debug_nproc=None)
                out=processes
            elif variable_in is 'processes_lumped':
                processes=calculate_rams_mp_path(filename,                                                    
                                                 microphysics_scheme='rams',processes='mass_grouped',
                                                 signed=True,
                                                 quantity='mixing ratio',constraint=None,
                                                 debug_nproc=None)
                processes_lumped=lump_processes(processes,microphysics_scheme='rams',lumping='basic')
                out=processes_lumped
                
            elif variable_in is 'processes_lumped_detailed':
                processes=calculate_rams_mp_path(filename,                                                    
                                                 microphysics_scheme='rams',processes='mass_grouped',
                                                 signed=True,
                                                 quantity='mixing ratio',constraint=None,
                                                 debug_nproc=None)
                processes_lumped=lump_processes(processes,microphysics_scheme='rams',lumping='detailed')
                out=processes_lumped
                
            elif variable_in is "hydrometeor_mass":
                out=iris.cube.CubeList()
                for hydrometeor in ['RCP','RDP','RRP','RPP','RSP','RAP','RGP','RHP']:
                    out.append(load_function(filename,hydrometeor))
            elif variable_in is "hydrometeor_number":
                out=iris.cube.CubeList()
                for hydrometeor in ['CCP','CDP','CRP','CPP','CSP','CAP','CGP','CHP']:
                    out.append(load_function(filename,hydrometeor))

            else:
                if variable_in=='w':
                    variable='WC'
                elif variable_in=='u':
                    variable='UC'
                elif variable_in=='v':
                    variable='VC'
                elif variable_in=='w_unstaggered':
                    variable='WC'
                elif variable_in=='u_unstaggered':
                    variable='UC'
                elif variable_in=='v_unstaggered':
                    variable='VC'
                elif variable_in=='NCLOUD':
                    variable='CCP'
                elif variable_in=='QCLOUD':
                    variable='RCP'
                elif variable_in=='airmass':
                    variable='airmass'
                elif variable_in=='airmass_path':
                    variable='airmass_path'                
                elif variable_in=='air_density':
                    variable='air_density'
                elif variable_in=='OLR':
                    variable='OLR'
                elif variable_in=='air_temperature':
                    variable='air_temperature'
                elif variable_in=='potential_temperature':
                    variable='THETA'
                else: 
                    variable=variable_in
                out=load_function(filename,variable)
                out.rename(variable_in)

            return out

        return load_function_out
