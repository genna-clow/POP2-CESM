! module satellite_sim

! use io_types,                     only : nml_in, nml_filename
! use kinds_mod,                    only : char_len, int_kind 
! use named_field_mod,              only : named_field_get !! get chlorophyll 
! use prognostic,                   only : TRACER, curtime !! get temp and salt fields
! use tavg,                         only : define_tavg_field
! use tavg,                         only : accumulate_tavg_field

! logical       :: run_satellite_sim = .false.  
! logical       :: chl = .false. 
! logical       :: temp = .false. 
! real (r8)     :: sample_time = .false. 
! integer       :: swath_width = .false. 

! ! Namelist definition
! namelist /satellite_sim_nml/ chl, chl_sample_time, chl_swath_width, temp, temp_sample_time, temp_swath_width


! contains


! subroutine init_satellite_sim()
    
!     #ifdef SAT_SIM
!         integer (int_kind) :: tavg_sat_Chl
!         integer (int_kind) :: tavg_sat_Chl_wgt
!         integer (int_kind) :: tavg_sat_SST
!         integer (int_kind) :: tavg_sat_SST_wgt
!         integer (int_kind) :: tavg_sat_PIC
!         integer (int_kind) :: tavg_sat_PIC_wgt

!         ! Open and read Namelist file.
!         open (nml_in, file=nml_filename, iostat=nml_error)
!         read (unit=nml_in, nml=satellite_sim_nml, iostat=nml_error) ! satellite_sim

!         if (temp) then
!             call define_tavg_field(tavg_SST_sat,'sst_satellite', 2,              &
!                                 long_name='Satellite-observed SST',   &
!                                 units='K', grid_loc='2110',      &
!                                 coordinates='TLONG TLAT time')
!             call define_tavg_field(tavg_SST_sat_wgt,'sst_satellite_wgt',2,              &
!                                 long_name='Weight for Satellite-observed SST',   &
!                                 units='K', grid_loc='2110',      &
!                                 coordinates='TLONG TLAT time')
!             if (temp_sample_time) then 
!                 call define_tavg_field(tavg_SST_sat_swath,'sst_satellite_swath', 2,              &
!                                     long_name='Satellite-observed SST with swath sampling',   &
!                                     units='degC', grid_loc='2110',      &
!                                     coordinates='TLONG TLAT time')
!                 call define_tavg_field(tavg_SST_sat_swath_wgt,'sst_satellite_swath_wgt',2,              &
!                                     long_name='Weight for Satellite-observed SST with swath sampling',   &
!                                     units='degC', grid_loc='2110',      &
!                                     coordinates='TLONG TLAT time')
!             endif
!         endif


!         if (chl) then
!             call define_tavg_field(tavg_chl_sat,'chl_satellite', 2,              &
!                                 long_name='Satellite-observed chlorophyll',   &
!                                 units='mg m3', grid_loc='2110',      &
!                                 coordinates='TLONG TLAT time')

!             call define_tavg_field(tavg_chl_sat_wgt,'chl_satellite_wgt',2,              &
!                                 long_name='Weight for satellite-observed chlorophyll',   &
!                                 units='mg m3', grid_loc='2110',      &
!                                 coordinates='TLONG TLAT time')
!             if (temp_sample_time) then 
!                 call define_tavg_field(tavg_chl_sat,'chl_satellite_swath', 2,              &
!                                     long_name='Satellite-observed chlorophyll with swath sampling',   &
!                                     units='degC', grid_loc='2110',      &
!                                     coordinates='TLONG TLAT time')
!                 call define_tavg_field(tavg_chl_sat_swath_wgt,'sst_satellite_swath_wgt',2,              &
!                                     long_name='Weight for satellite-observed chlorophyll with swath sampling',   &
!                                     units='degC', grid_loc='2110',      &
!                                     coordinates='TLONG TLAT time')
!             endif
            
!         endif
!     #endif
! end subroutine init_satellite_sim


! subroutine accumulate_weighted_field() 


!     if (sample_time) then 
!         call calculate_swath(sample_time, swath_width)
!     endif


!     if (chl) then
!         call named_field_get(totChl_surf_nf_ind, bid, CHL(:,:)) 
!         call accumulate_tavg_field(CHL(:,:), tavg_Chl, bid, 1) 
!     endif
   
!     if (temp) then
!         TCUR = TRACER(:,:,:,:,curtime,iblock)
!         call accumulate_tavg_field(TCUR(:,:,k,1), tavg_SST_sat, iblock,1) 
!         ! TCUR = tracers at current timestep
!         ! TCUR(:,:,k,2) = salt
!         ! TCUR(:,:,k,1) = temp
!         ! k = levels, 1 = surface
!         ! parameter 1 = surface
!     endif
        

! end subroutine accumulate_weighted_field 


! subroutine calculate_swath(sample_time, swath_width)

!     real (r8)                               :: loc_time(nx_block, ny_block, max_blocks_clinic) ! local time of day 
!     real (r8)                               :: sat_loc_time ! approximate time of satellite fly-over 
!     real (r8)                               :: sat_lon ! longitude of satellite 
!     real (r8)                               :: swath(nx_block, ny_block, max_blocks_clinic) ! swath mask 
!     real (r8)                               :: dlon(nx_block, ny_block, max_blocks_clinic) ! distance to satellite longitude in deg 
!     real (r8)                               :: dx(nx_block, ny_block, max_blocks_clinic) ! distance to satellite longitude in km 

!     sat_lon = 15._r8*(sample_time-(frac_day*24._r8))
!     swath = c0
!     dlon = mod((TLOND - sat_lon + 180._r8), 360.0_r8) - 180._r8
!     dx = dlon * (c1/radian) * COS(TLAT) * radius
!     where (abs(dx)<(swath_width*p5))
!        swath = c1
!     endwhere

! end subroutine calculate_swath 


! module tavg_sat_vars_type 
!     type tavg_sat_vars                          
!         real (r8) :: data 
!         integer :: tavg_ind 
!         real (r8) :: data_wgts 
!         integer :: tavg_wgt_ind               
!     end type tavg_sat_vars 

! end module tavg_sat_vars_type 



! subroutine calculate_satellite_chlor(isccp, sample_swath) 

 
!     use tavg_sat_vars_type 


!     real (r8)                               :: CHL(nx_block,ny_block) ! total surface chlorophyll conc. 
!     real (r8)                               :: Chl_sat_weight(nx_block, ny_block) ! weight to mask areas not viewable by satellite 
!     real (r8)                               :: cloud_weight_modis(nx_block, ny_block) ! cloud weight for MODIS 
!     real (r8)                               :: Chl_sat_weight_modis(nx_block, ny_block) ! chl weight for MODIS 
!     real (r8)                               :: cloud_weight_isccp(nx_block, ny_block) ! cloud weight for ISCCP 
!     real (r8)                               :: Chl_sat_weight_isccp(nx_block, ny_block) ! chl weight for ISCCP 
!     type (tavg_sat_vars_type)               :: baseline, isccp_chlor, modis_chlor, modis_chlor_swath 

   

!     call named_field_get(totChl_surf_nf_ind, bid, CHL(:,:)) 
!     ! Accumulate total chlorophyll field 
!     call accumulate_tavg_field(CHL(:,:), tavg_Chl, bid, 1) 

 

!     no_swath = c1 
!     ! Set weight to 1 
!     Chl_sat_weight = c1 
!     swath = call calculate_swath 
!     !! Accumulate daylight chlorophyll 
!     ! Calculate daylight weight 
!     where (COSZEN(:,:,bid) .le. 0.342) 
!        Chl_sat_weight = c0 
!     endwhere 

 

!     ! Calculate ice weight 
!     Chl_sat_weight =  Chl_sat_weight*(c1-IFRAC(:,:,bid)) 
!     ! Calculate cloud/ice weight 
!     ! Assume that sea ice and clouds have random overlap 
!     Chl_sat_weight_modis = Chl_sat_weight*(c1-(CLOUDFRAC_MODIS(:,:,bid)/100._r8)) 
!     Chl_sat_weight_isccp = Chl_sat_weight*(c1-(CLOUDFRAC_ISCCP(:,:,bid)/100._r8)) 

 

!     ! Define variables: 
!     !!! baseline chlor 
!     baseline = tavg_sat_vars(CHL(:,:,bid), tavg_SatChl_nocld, Chl_sat_weight(:,:), tavg_SatChl_weight_nocld) 
!     call accumulate_weighted_field(baseline, no_swath) 


!     !!! cloud fractions 
!     ! Accumulate weighted cloud fraction (daylight-only) 
!     call accumulate_tavg_field(CLOUDFRAC_MODIS(:,:,bid)*Chl_sat_weight(:,:), tavg_cloudfrac_modis, bid, 1) 
!     call accumulate_tavg_field(CLOUDFRAC_ISCCP(:,:,bid)*Chl_sat_weight(:,:), tavg_cloudfrac_isccp, bid, 1) 


!     !!! chlorophyll 
!     isccp_chlor = tavg_sat_vars(CHL(:,:,bid), tavg_isccp_Chl, Chl_sat_weight_isccp(:,:), tavg_isccp_Chl_weight) 
!     call accumulate_weighted_field(isccp, no_swath) 

!     modis_chlor = tavg_sat_vars(CHL(:,:,bid), tavg_modis_Chl, Chl_sat_weight_modis(:,:), tavg_modis_Chl_weight) 
!     call accumulate_weighted_field(modis_chlor, no_swath) 

!     modis_chlor_swath = tavg_sat_vars(CHL(:,:,bid), tavg_modis_Chl, Chl_sat_weight_modis(:,:), tavg_modis_Chl_weight) 
!     call accumulate_weighted_field(modis_chlor_swath, swath) 

! end 