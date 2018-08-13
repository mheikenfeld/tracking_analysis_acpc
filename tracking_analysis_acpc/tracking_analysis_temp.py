
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
    
#    if 'interpolation' in mode:     
                            
#        logging.debug('Start calculations for individual cells')
#        savedir=os.path.join(top_savedir_tracking)
#        input_str=os.path.join(directory,filename)
#        logging.debug(input_str)
#    
#        logging.debug('Start calculating profiles for individual cells (simple method)')
#
#        Tracks=pd.read_hdf(os.path.join(savedir,'Track.hdf'),'table')
##                    Mask_w=iris.load_cube(os.path.join(savedir,'Watershed_w.nc'),'watershedding_output_mask')
##                    Mask=Mask_w[:,1:,:,:]
#        Mask=iris.load_cube(os.path.join(savedir,'Watershed_TWC.nc'),'watershedding_output_mask')
#
#
#
#        logging.debug('start loading processes')
#                
#        add_coordinates=['xy','z']
#        Processes=calculate_wrf_mp_path(filenames,                                                    
#                                                microphysics_scheme=mp['WRF'],processes='mass',
#                                                signed=True,
#                                                quantity='volume',constraint=None,add_coordinates=add_coordinates,
#                                                debug_nproc=None,verbose=False)
#        logging.debug('loaded processes')
#
#        logging.debug('start loading auxillary variables')
#
#        #Load auxillary variables for current time:
#        Aux, Aux_2D,Aux_sum=load_AUX(filenames,microphysics_scheme=mp['WRF'],
#                            constraint=None,add_coordinates=add_coordinates,
#                            essential_only=False)
##    
#
#        logging.debug('loaded auxillary variables')
#        
#        cells=list(np.unique(Tracks['particle']))
#        if cell_selection is not None:
#            cells=list(set(cells) & set(cell_selection))
#
#        for particle in cells:                        
#            
#            logging.debug( 'Start calculating Data for cell '+ str(int(particle)))
#            savedir_cell=os.path.join(savedir,'cells',str(int(particle)))
#            os.makedirs(savedir_cell,exist_ok=True)
#
#            savefile_along=os.path.join(savedir_cell,'Processes_along'+'_'+str(int(particle))+'.nc')
#            savefile_across=os.path.join(savedir_cell,'Processes_across'+'_'+str(int(particle))+'.nc')
#            savefile_Aux_along=os.path.join(savedir_cell,'Aux_along'+'_'+str(int(particle))+'.nc')
#            savefile_Aux_across=os.path.join(savedir_cell,'Aux_across'+'_'+str(int(particle))+'.nc')
##                                savefile_Aux_2D=os.path.join(savedir,subpath,'Aux_2D'+'_' +subpath +'.nc')
#            savefile_Aux_2D_cell_averaged=os.path.join(savedir_cell,'Aux_2D_cell_averaged'+'_'+str(int(particle))+'.nc')
#            savefile_Aux_2D_cell_sum=os.path.join(savedir_cell,'Aux_2D_cell_sum'+'_'+str(int(particle))+'.nc')
#            savefile_Processes_cell_sum=os.path.join(savedir_cell,'Processes_cell_sum'+'_'+str(int(particle))+'.nc')
#            savefile_Aux_cell_averaged=os.path.join(savedir_cell,'Aux_cell_averaged'+'_'+str(int(particle))+'.nc')
#            savefile_Aux_cell_sum=os.path.join(savedir_cell,'Aux_cell_sum'+'_'+str(int(particle))+'.nc')
#
#            Mask_particle=mask_particle(Mask,particle,masked=False)
#            Mask_particle_surface=mask_particle_surface(Mask,particle,masked=False,z_coord='model_level_number')
#            Track=Tracks[Tracks['particle']==particle]
#
#            width=20
#            dx=1000
#            height_levels=np.arange(0,20001,1000.)
#                                        
#            Interpolation_out=track_cell_orientation_subset(Processes,Aux,Aux_2D,Aux_sum,                                                     
#                                                            track=Track,mask_particle=Mask_particle,mask_particle_surface=Mask_particle_surface,
#                                                            width=width,dx=dx, z_coord='geopotential_height', height_levels=height_levels)
#
#            iris.save(Interpolation_out['Processes_along'],savefile_along)
#            iris.save(Interpolation_out['Processes_across'],savefile_across)
#            iris.save(Interpolation_out['Aux_along'],savefile_Aux_along)
#            iris.save(Interpolation_out['Aux_across'],savefile_Aux_across)
#            iris.save(Interpolation_out['Aux_2D_cell_averaged'],savefile_Aux_2D_cell_averaged)
#            iris.save(Interpolation_out['Aux_2D_cell_averaged'],savefile_Aux_2D_cell_sum)
#            iris.save(Interpolation_out['Processes_cell_sum'],savefile_Processes_cell_sum)
#            iris.save(Interpolation_out['Aux_cell_averaged'],savefile_Aux_cell_averaged)
#            iris.save(Interpolation_out['Aux_cell_sum'],savefile_Aux_cell_sum)
#
##                        iris.save(Aux_domain_averaged,savefile_Aux_domain_averaged)
#            logging.debug( 'Data calculated and saved to '+ savedir_cell)
#        logging.debug( 'Data calculated and saved to file for all cells')

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
#    if 'plot_profiles' in mode:
#        savedir=os.path.join(top_savedir_tracking)
#        plotdir=os.path.join(top_plotdir)
#        
#        Track=pd.read_hdf(os.path.join(savedir,'Track.hdf'),'table')
#        fig5,ax5=plt.subplots(nrows=1,ncols=1)
#
#        for particle in np.unique(Track['particle']):            
#            fig6,ax6=plt.subplots(nrows=1,ncols=1)
#
##                savefile_max=os.path.join(savedir,'W_CDNC','Max','W_CDNC'+'_Max'+'_'+str(int(particle))+'.nc')
##                savefile_min=os.path.join(savedir,'W_CDNC','Min','W_CDNC'+'_Min'+'_'+str(int(particle))+'.nc')
##                savefile_mean=os.path.join(savedir,'W_CDNC','Mean','W_CDNC'+'_Mean'+'_'+str(int(particle))+'.nc')
#            
#            savefile_max=os.path.join(savedir,'Max','Max'+'_'+str(int(particle))+'.nc')
#            savefile_min=os.path.join(savedir,'Min','Min'+'_'+str(int(particle))+'.nc')
#            savefile_mean=os.path.join(savedir,'Mean','Mean'+'_'+str(int(particle))+'.nc')
#
#            os.makedirs(os.path.join(plotdir,'Profiles_w_individual'),exist_ok=True)
#            os.makedirs(os.path.join(plotdir,'Profiles_w'),exist_ok=True)
#
#            w_max=iris.load_cube(savefile_max,'w')
#            w_min=iris.load_cube(savefile_min,'w')
#            w_mean=iris.load_cube(savefile_mean,'w')
#            
#            for i,time in enumerate(w_max.coord('time').points):
#                logging.debug('plotting '+str(particle)+ ' ' +str(i))                        
#                logging.debug(str(w_max[i]))
#                plot_max=ax5.plot(w_max[i].data,w_max[i].coord('model_level_number').points,color='r',marker='.',markersize=2,linestyle='None',label='max')    
#                plot_min=ax5.plot(w_min[i].data,w_min[i].coord('model_level_number').points,color='b',marker='.',markersize=2,linestyle='None',label='min')    
#                plot_mean=ax5.plot(w_mean[i].data,w_mean[i].coord('model_level_number').points,color='k',marker='.',markersize=2,linestyle='None',label='mean')    
#                plot_max_ind=ax6.plot(w_max[i].data,w_max[i].coord('model_level_number').points,color='r',marker='.',markersize=2,linestyle='None',label='max')    
#                plot_min_ind=ax6.plot(w_min[i].data,w_min[i].coord('model_level_number').points,color='b',marker='.',markersize=2,linestyle='None',label='min')    
#                plot_mean_ind=ax6.plot(w_mean[i].data,w_mean[i].coord('model_level_number').points,color='k',marker='.',markersize=2,linestyle='None',label='mean')
#                ax6.set_xlim([0,30])
#                ax6.set_ylim([0,100])
#
##                    iplt.plot(w_max[i].coord('geopotential_height').points,'rx',ax=ax5)    
##                    iplt.plot(w_min[i].coord('geopotential_height').points,'bx',ax=ax5)    
##                    iplt.plot(w_mean[i].coord('geopotential_height').points,'kx',ax=ax5)    
##                    iplt.plot(w_max[i].coord('geopotential_height').points,'rx',ax=ax6)    
##                    iplt.plot(w_min[i].coord('geopotential_height').points,'bx',ax=ax6)    
##                    iplt.plot(w_mean[i].coord('geopotential_height').points,'kx',ax=ax6)    
##                        ax6.legend()
##                        ax5.legend()
#            fig6.savefig(os.path.join(plotdir,'Profiles_w_individual','Profile_w_'+str(int(particle))+'.png'),dpi=600)
#            plt.close(fig6)
#        ax5.set_xlim([0,30])
#        ax5.set_ylim([0,15])
#
#        fig5.savefig(os.path.join(plotdir,'Profiles_w','Profile_w'+'.png'),dpi=600)
#        plt.close(fig5)
#    
##
#    