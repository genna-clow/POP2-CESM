!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

module forcing_fields

!BOP
! !MODULE: forcing_fields

! !DESCRIPTION:
!  Contains the forcing fields necessary for supporting high-level coupling.
!  These fields originally resided in modules forcing and forcing_coupled.

! !REVISION HISTORY:
!  SVN:$Id$
!
! !USES:

   use kinds_mod
   use blocks,      only: nx_block, ny_block
   use constants,   only: c0
   use domain_size, only: max_blocks_clinic,nt
   
   ! new imports: 
   ! use timers
   use time_management ! defines tday00 and nsteps_per_interval

   implicit none
   save

!EOP
!BOC
! !PUBLIC DATA MEMBERS:

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic),public ::  &
      EVAP_F = c0,       &! evaporation   flux    from cpl (kg/m2/s)
      PREC_F = c0,       &! precipitation flux    from cpl (kg/m2/s)
                          ! (rain + snow)
      SNOW_F = c0,       &! snow          flux    from cpl (kg/m2/s)
      MELT_F = c0,       &! melt          flux    from cpl (kg/m2/s)
      ROFF_F = c0,       &! river runoff  flux    from cpl (kg/m2/s)
      IOFF_F = c0,       &! ice   runoff  flux    from cpl (kg/m2/s)
      SALT_F = c0,       &! salt          flux    from cpl (kg(salt)/m2/s)
      SENH_F = c0,       &! sensible heat flux    from cpl (W/m2   )
      LWUP_F = c0,       &! longwave heat flux up from cpl (W/m2   )
      LWDN_F = c0,       &! longwave heat flux dn from cpl (W/m2   )
      MELTH_F= c0         ! melt     heat flux    from cpl (W/m2   )


   integer(kind=int_kind), public :: &
      ATM_CO2_PROG_nf_ind = 0, & ! bottom atm level prognostic co2
      ATM_CO2_DIAG_nf_ind = 0    ! bottom atm level diagnostic co2

  integer(kind=int_kind), public :: &
       ATM_NHx_nf_ind = 0, & ! bottom atm level NHx flux
       ATM_NOy_nf_ind = 0    ! bottom atm level NOy flux

   real (r8), dimension(nx_block,ny_block,2,max_blocks_clinic), &
      public, target :: &
      SMF,  &!  surface momentum fluxes (wind stress)
      SMFT   !  surface momentum fluxes on T points if avail

   real (r8), dimension(nx_block,ny_block,nt,max_blocks_clinic), &
      public, target :: &
      STF,      &! surface tracer fluxes
      STF_RIV,  &! riverine tracer fluxes
      TFW        ! tracer content in freshwater flux

   logical (log_kind), dimension(nt), public :: &
      lhas_vflux,   & ! true if a tracer uses virtual fluxes
      lhas_riv_flux   ! true if a tracer has a riverine flux

   integer(kind=int_kind), public :: &
      vflux_tracer_cnt ! number of tracers for which lhas_vflux is .true.

   logical (log_kind), public :: &
      lsmft_avail   ! true if SMFT is an available field

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic), &
      public, target ::  &
      IFRAC,             &! ice fraction; not initialized in this routine
      U10_SQR,           &! 10m wind speed squared; not initialized in this routine
      ATM_PRESS,         &! atmospheric pressure forcing
      FW,FW_OLD,         &! freshwater flux at T points (cm/s)
                          ! FW_OLD is at time n-1
      LAMULT,            &! Langmuir multiplier
      LASL,              &! surface layer averaged Langmuir number
      USTOKES,           &! surface Stokes drift x component
      VSTOKES             ! surface Stokes drift y component

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic), &
      public, target ::  &
      ATM_FINE_DUST_FLUX,       &! fine dust flux from atm from cpl (g/cm2/s)
      ATM_COARSE_DUST_FLUX,     &! coarse dust flux from atm from cpl (g/cm2/s)
      SEAICE_DUST_FLUX,         &! coarse dust flux from seaice from cpl (g/cm2/s)
      ATM_BLACK_CARBON_FLUX,    &! black carbon flux from atm from cpl (g/cm2/s)
      SEAICE_BLACK_CARBON_FLUX   ! black carbon flux from seaice from cpl (g/cm2/s)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! new section copied from forcing_coupled.F90
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! copied nn from overflows.F90
   integer   (int_kind)  ::  &
      nn,                    & ! ovf tracer index
      timer_compute_cosz,    & ! 
      iblock

   real (r8) ::  &
      tday00_interval_beg,    & ! model time at beginning of coupling interval
      orb_eccen,              & ! Earth eccentricity
      orb_obliqr,             & ! Earth Obliquity
      orb_lambm0,             & ! longitude of perihelion at v-equinox
      orb_mvelpp                ! Earths Moving vernal equinox of orbit +pi

   real (r8) :: cosz_day        ! time where cosz is computed

   real (r8), dimension(:,:,:), public, allocatable :: & ! added public
      QSW_COSZ_WGHT,      & ! weights
      QSW_COSZ_WGHT_NORM    ! normalization for QSW_COSZ_WGHT

   integer (int_kind), private ::   &
      cpl_ts                ! flag id for coupled_ts flag

 contains

!-----------------------------------------------------------------------
!  Compute QSW_COSZ_WGHT_NORM.
!-----------------------------------------------------------------------
subroutine calc_QSW_COSZ_WGHT
!       ! delete if statement because we always want to calculate cosz
      ! if ( qsw_distrb_iopt == qsw_distrb_iopt_cosz ) then
   tday00_interval_beg = tday00

   !$OMP PARALLEL DO PRIVATE(iblock,nn,cosz_day)
   do iblock = 1, nblocks_clinic

      QSW_COSZ_WGHT_NORM(:,:,iblock) = c0

      do nn = 1, nsteps_per_interval
         cosz_day = tday00_interval_beg + interval_cum_dayfrac(nn-1) &
            - interval_cum_dayfrac(nsteps_per_interval)
   
         call compute_cosz(cosz_day, iblock, QSW_COSZ_WGHT(:,:,iblock))

         if (interval_avg_ts(nn)) then
            QSW_COSZ_WGHT_NORM(:,:,iblock) = &
               QSW_COSZ_WGHT_NORM(:,:,iblock) &
               + p5 * QSW_COSZ_WGHT(:,:,iblock)
         else
            QSW_COSZ_WGHT_NORM(:,:,iblock) = &
               QSW_COSZ_WGHT_NORM(:,:,iblock) &
               + QSW_COSZ_WGHT(:,:,iblock)
         endif

      enddo

      where (QSW_COSZ_WGHT_NORM(:,:,iblock) > c0) &
         QSW_COSZ_WGHT_NORM(:,:,iblock) = &
            (fullsteps_per_interval + p5 * halfsteps_per_interval) &
            / QSW_COSZ_WGHT_NORM(:,:,iblock)

   enddo
   !$OMP END PARALLEL DO

end subroutine calc_QSW_COSZ_WGHT
!***********************************************************************
!BOP
! !IROUTINE: compute_cosz
! !INTERFACE:

subroutine compute_cosz(tday, iblock, COSZ)

   ! !DESCRIPTION:
   !  This subroutine computes cos of the solar zenith angle.
   !  Negative values are set to zero.
   !
   ! !REVISION HISTORY:
   !  same as module
   !
   ! !USES:
   
      use shr_orb_mod, only: shr_orb_decl, shr_orb_cosz
   
   ! !INPUT PARAMETERS:
   
      real (r8), intent(in) :: tday
      integer (int_kind), intent(in) :: iblock
   
   ! !OUTPUT PARAMETERS:
   
      real (r8), dimension(:,:), intent(out) :: COSZ
   
   !EOP
   !BOC
   !-----------------------------------------------------------------------
   !
   !  local variables
   !
   !-----------------------------------------------------------------------
   
      integer (int_kind) ::   &
         i, j            ! loop indices
   
      real (r8) :: &
         calday,       & ! Calendar day, including fraction
         delta,        & ! Solar declination angle in rad
         eccf            ! Earth-sun distance factor (ie. (1/r)**2)
   
   !-----------------------------------------------------------------------
   
      ! call timer_start(timer_compute_cosz, block_id=iblock)
   
   !  shr_orb code assumes Jan 1 = calday 1, unlike Jan 1 = tday 0
      calday = tday + c1
   
      call shr_orb_decl(calday, orb_eccen, orb_mvelpp, orb_lambm0, &
                        orb_obliqr, delta, eccf)
   
      do j = 1, ny_block
         do i = 1, nx_block
            COSZ(i,j) = shr_orb_cosz(calday, TLAT(i,j,iblock), &
                                       TLON(i,j,iblock), delta)
            COSZ(i,j) = max(c0, COSZ(i,j))
         enddo
      enddo
   
      ! call timer_stop(timer_compute_cosz, block_id=iblock)
   
   !-----------------------------------------------------------------------
   !EOC
   
end subroutine compute_cosz

!***********************************************************************

end module forcing_fields

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
