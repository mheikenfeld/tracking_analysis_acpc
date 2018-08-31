
#    if 'cog' in mode:                    
#
#        logging.debug('Start COG calculations for individual cells')
#        
#        savedir=os.path.join(top_savedir_tracking)
#        input_str=os.path.join(directory,filename)
#        logging.debug(input_str)
#    
#
#        logging.debug('Start calculating COG')
#
#        Tracks=pd.read_hdf(os.path.join(savedir,'Track.hdf'),'table')
#        Mask=iris.load_cube(os.path.join(savedir,'Watershed_TWC.nc'),'watershedding_output_mask')
#
#        TWC=iris.load_cube(os.path.join(savedir,'Data_WC.nc'),'total water content')    
#        LWC=iris.load_cube(os.path.join(savedir,'Data_WC.nc'),'liquid water content')    
#        IWC=iris.load_cube(os.path.join(savedir,'Data_WC.nc'),'ice water content')    
##       Airmass=iris.load_cube(os.path.join(savedir,'Data_WC.nc'),'mass_of_air')    
#
#        logging.debug('water contents loaded')
#
#        
#        logging.debug('Data loaded')
#
#        M_total=TWC*Airmass
#        M_total.rename('Mass_total')
#        M_liquid=LWC*Airmass
#        M_liquid.rename('Mass_total')
#        M_frozen=IWC*Airmass
#        M_frozen.rename('Mass_frozen')
        
        
        
#        # serial calculation of all COG in one go
#        Tracks_COG_total=calculate_cog(Tracks,M_total,Mask)
#        logging.debug('COG total loaded')
#        Tracks_COG_liquid=calculate_cog(Tracks,M_liquid,Mask)
#        logging.debug('COG liquid loaded')
#        Tracks_COG_frozen=calculate_cog(Tracks,M_frozen,Mask)
#        logging.debug('COG frozen loaded')
#
#        # loop over particles for possible parallelisation using dask
#        Tracks_COG_total_list=[]
#        Tracks_COG_liquid_list=[]
#        Tracks_COG_frozen_list=[]
#
#        cells=list(np.unique(Tracks['particle']))
#        if cell_selection is not None:
#            cells=list(set(cells) & set(cell_selection))
#            
#        
#        for particle in cells:
#            cog_cell(Tracks,M_total,M_liquid,M_frozen,Mask,savedir,particle)
#        
#        #Set savedir for COG calculations
#        savedir_cog=os.path.join(savedir,'COG')
#        os.makedirs(savedir_cog,exist_ok=True)
#
#        # Merge individualc cells into one DataFrame:
#        Tracks_COG_total=pd.concat(Tracks_COG_total_list)
#        Tracks_COG_liquid=pd.concat(Tracks_COG_liquid_list)
#        Tracks_COG_frozen=pd.concat(Tracks_COG_frozen_list)
#        # set files to save data to:
#        savefile_COG_total=os.path.join(savedir_cog,'COG_total.h5')
#        savefile_COG_liquid=os.path.join(savedir_cog,'COG_liquid.h5')
#        savefile_COG_frozen=os.path.join(savedir_cog,'COG_frozen.h5')
#        # Save data to hdf files:
#        Tracks_COG_total.to_hdf(savefile_COG_total,'table')
#        Tracks_COG_liquid.to_hdf(savefile_COG_liquid,'table')
#        Tracks_COG_frozen.to_hdf(savefile_COG_frozen,'table')
#        logging.debug('COG merged and saved to '+ savedir_cog)
#

#        COG_untracked_total=calculate_cog_untracked(M_total,Mask)
#        COG_untracked_liquid=calculate_cog_untracked(M_liquid,Mask)
#        COG_untracked_frozen=calculate_cog_untracked(M_frozen,Mask)   
#        # set files to save data to:
#        savefile_COG_untracked_total=os.path.join(savedir_cog,'COG_untracked_total.h5')
#        savefile_COG_untracked_liquid=os.path.join(savedir_cog,'COG_untracked_liquid.h5')
#        savefile_COG_untracked_frozen=os.path.join(savedir_cog,'COG_untracked_frozen.h5')
#        # Save data to hdf files:
#        COG_untracked_total.to_hdf(savefile_COG_untracked_total,'table')
#        COG_untracked_liquid.to_hdf(savefile_COG_untracked_liquid,'table')
#        COG_untracked_frozen.to_hdf(savefile_COG_untracked_frozen,'table')
#        logging.debug('untracked COG calculated and saved to '+ savedir_cog)
#    
#        # Calculate domainwide COG data
#        COG_domain_total=calculate_cog_domain(M_total)
#        COG_domain_liquid=calculate_cog_domain(M_liquid)
#        COG_domain_frozen=calculate_cog_domain(M_frozen)
#        # set files to save data to:
#        savefile_COG_domain_total=os.path.join(savedir_cog,'COG_domain_total.h5')
#        savefile_COG_domain_liquid=os.path.join(savedir_cog,'COG_domain_liquid.h5')
#        savefile_COG_domain_frozen=os.path.join(savedir_cog,'COG_domain_frozen.h5')
#        # Save data to hdf files:
#        COG_domain_total.to_hdf(savefile_COG_domain_total,'table')
#        COG_domain_liquid.to_hdf(savefile_COG_domain_liquid,'table')
#        COG_domain_frozen.to_hdf(savefile_COG_domain_frozen,'table')
#        logging.debug('domainwide COG calculated and saved to '+ savedir_cog)
#    
#        logging.debug('COG calculated and saved to '+ savedir_cog)
    
#        #%%
#    if 'profiles' in mode:
#        savedir=os.path.join(top_savedir_tracking)
#        plotdir=os.path.join(top_plotdir)
#    
#
#        logging.debug('Start calculating profiles for individual cells (simple method)')
#
#        Track=pd.read_hdf(os.path.join(savedir,'Track.hdf'),'table')
#        Mask_w=iris.load_cube(os.path.join(savedir,'Watershed_w.nc'),'watershedding_output_mask')
#        Mask_T=Mask_w[:,1:,:,:]
#
#        logging.debug('Tracks and mask loaded')
#
##                W=iris.load_cube(os.path.join(savedir,'Data_w.nc'))
##                          
##                logging.debug('vertical velocity loaded')
#        
#        QCLOUD=load_function(filenames,'QCLOUD')
#        logging.debug('CDNC loaded')
#
#        CDNC=load_function(filenames,'NCLOUD')
#        
#        logging.debug('CDNC loaded')
##    ##         
#        logging.debug('start loading processes')
#
#        Processes=load_function(filenames,'processes')
#        Processes_lumped=lump_processes(Processes,microphysics_scheme=mp[model],others=False)
#        
#        logging.debug('microphysics process rates loaded')
##    
#        Airmass=load_function(filenames,'airmass')
##
#        Processes_sum=iris.cube.CubeList()
#        for process in Processes_lumped:
#            cube=Airmass*process
#            cube.rename(process.name())
#            Processes_sum.append(cube)
#
#        dimensions=('x','y')
#        cubelist_Aux=iris.cube.CubeList([CDNC])        
#        cell_statistics(cubelist_Aux,Track,Mask_T,dimensions,aggregators=[MEAN,MAX,MIN,MEDIAN],output_path=savedir,output_name='CDNC')                 
#        cell_statistics(cubelist_Aux,Track,Mask_T,dimensions,aggregators=[PERCENTILE],output_path=savedir,output_name='CDNC',percent=[0,25,50,75,100])                 
#
#        cell_statistics(Processes_lumped,Track,Mask_T,dimensions,aggregators=[MEAN,MAX,MIN,MEDIAN],output_path=savedir,output_name='Processes')                
#        
#        dimensions=('x','y','model_level_number')
#        cell_statistics(Processes_sum,Track,Mask_w,dimensions,aggregators=[SUM],output_path=savedir,output_name='Processes_integrated')                
##
#        logging.debug('loaded and saved profiles')
#                
#        logging.debug('profile calculations done and saved')
#
#    
