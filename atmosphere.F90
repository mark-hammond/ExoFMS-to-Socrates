!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                   !!
!!                   GNU General Public License                      !!
!!                                                                   !!
!! This file is part of the Flexible Modeling System (FMS).          !!
!!                                                                   !!
!! FMS is free software; you can redistribute it and/or modify       !!
!! it and are expected to follow the terms of the GNU General Public !!
!! License as published by the Free Software Foundation.             !!
!!                                                                   !!
!! FMS is distributed in the hope that it will be useful,            !!
!! but WITHOUT ANY WARRANTY; without even the implied warranty of    !!
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the     !!
!! GNU General Public License for more details.                      !!
!!                                                                   !!
!! You should have received a copy of the GNU General Public License !!
!! along with FMS; if not, write to:                                 !!
!!          Free Software Foundation, Inc.                           !!
!!          59 Temple Place, Suite 330                               !!
!!          Boston, MA  02111-1307  USA                              !!
!! or see:                                                           !!
!!          http://www.gnu.org/licenses/gpl.txt                      !!
!!                                                                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Modified by MDH to interface with SOCRATES (UK Met Office) radiative transfer routine


MODULE atmosphere_mod

  !-----------------------------------------------------------------------
  !
  !    interface for FV dynamical core with Held-Suarez forcing
  !
  !-----------------------------------------------------------------------


  USE time_manager_mod, ONLY: time_type, get_time, set_time, OPERATOR(+)

  USE fms_mod,          ONLY: file_exist, open_namelist_file,   &
       error_mesg, FATAL,                &
       check_nml_error, stdlog,          &
       write_version_number,             &
       close_file, set_domain

  USE hs_forcing_mod,   ONLY: hs_forcing_init
  USE constants_mod, ONLY: omega, cp_air, rdgas, kappa, radius, grav, rvgas, &
       cp_vapor, cp_water, latent_heat, rho_cp, rho0, r_vir, cp_vir

  !------------------
  ! FV specific codes:
  !------------------
  !  use fv_pack
  !use            fv_pack, only: nlon, mlat, nlev, beglat, endlat, beglon, endlon, &
  !                              u, v, pt, q, ua, va, delp, phis, ps, pe, peln,    &
  !                              pk, pkz, omga, u_srf, v_srf, rlonb, rlatb,        &
  !                              cold_start, ncnst, pnats, consv_te, ptop,         &
  !                              fv_init, fv_domain, fv_end, change_time, map_dt,  &
  !                              adiabatic
  USE fv_pack, ONLY: nlon, mlat, nlev, beglat, endlat, beglon, endlon, &
       rlonb, rlatb, rlat, rlon,       &
       cold_start, ncnst, pnats, consv_te, ptop,         &
       fv_init, fv_domain, fv_end, change_time, map_dt,  &
       adiabatic, restart_format, &
       get_eta_level, p_var, ak, bk, ng_d, master, &
       nt_phys
  USE update_fv_phys_mod, ONLY: update_fv_phys
  USE fv_arrays_mod, ONLY: fv_array_sync

  USE fv_diagnostics, ONLY: fv_diag_init, fv_diag, fv_time
  USE timingModule,   ONLY: timing_on, timing_off
  USE fv_restart_mod, ONLY: fv_restart, write_fv_rst
  USE fv_dynamics_mod, ONLY: fv_dynamics
  USE fv_phys_mod, ONLY: fv_phys

  USE  diag_manager_mod, ONLY: register_diag_field, send_data
  !  use  radiation_mod, only: radiation_init, radiation_down, radiation_up, &
  !                                    radiation_end
  USE  surface_flux_mod, ONLY: surface_flux
  USE dry_convection_mod, ONLY: dry_convection
  USE  qe_moist_convection_mod, ONLY: moist_convection, compute_k, d622
  USE  simple_sat_vapor_pres_mod, ONLY: escomp

  !-----------------------------------------------------------------------

  !-------
  !SOCRATES specific
  !-------
  USE realtype_rd
  USE soc_constants_mod
  USE read_control_mod
  USE socrates_calc_mod
  USE socrates_interface_mod
  USE compress_spectrum_mod
  USE def_spectrum
  USE def_dimen,   ONLY: StrDim
  USE def_control, ONLY: StrCtrl,  allocate_control,   deallocate_control
  USE def_atm,     ONLY: StrAtm,   allocate_atm,       deallocate_atm
  USE def_cld,     ONLY: StrCld,   allocate_cld,       deallocate_cld, &
       allocate_cld_prsc,  deallocate_cld_prsc, &
       allocate_cld_mcica, deallocate_cld_mcica
  USE def_aer,     ONLY: StrAer,   allocate_aer,       deallocate_aer, &
       allocate_aer_prsc,  deallocate_aer_prsc
  USE def_bound,   ONLY: StrBound, allocate_bound,     deallocate_bound
  USE def_out,     ONLY: StrOut,                       deallocate_out
  !-----------------------------------------------------------------------

  IMPLICIT NONE
  PRIVATE

  PUBLIC   atmosphere_init, atmosphere,  atmosphere_end

  !-----------------------------------------------------------------------

  CHARACTER(len=128) :: version = '$Id: atmosphere.F90,v 13.0.2.2 2006/05/19 16:44:39 wfc Exp $'
  CHARACTER(len=128) :: tag = '$Name:  $'
  CHARACTER(len=10), PARAMETER :: mod_name='atmosphere'

  !-----------------------------------------------------------------------
  !---- private data ----

  TYPE        (time_type) :: Time_step_atmos
  REAL                    :: dt_atmos
  INTEGER :: sec
  INTEGER days, seconds

  !-----------------------------------------------------------------------
  !df 1401010
  INTEGER :: i, j, k
  INTEGER :: ij, nx, tsiz, isiz
  INTEGER :: is, ie
  !real :: eff_heat_capacity = rho_cp*1.
  REAL :: mld = 1.
  INTEGER :: id_t_surf, id_conv_rain, id_cond_rain, id_pme
  INTEGER :: id_conv_rain_profile, id_cond_rain_profile
  INTEGER :: id_flux_t, id_flux_q, id_flux_r, id_rh
  LOGICAL :: used, soc_mode, hires_mode
  REAL, ALLOCATABLE, DIMENSION(:,:)   :: t_surf, fms_stellar_flux, lats, lons
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: p_half, p_full, rh
  REAL, ALLOCATABLE, DIMENSION(:,:)   :: flux_t, flux_q, flux_r
  REAL, ALLOCATABLE, DIMENSION(:,:)   :: conv_rain, cond_rain, pme
  REAL, ALLOCATABLE, DIMENSION(:,:)   :: net_surf_sw_down, surf_lw_down
  REAL, ALLOCATABLE, DIMENSION(:,:)   :: output_net_surf_sw_down, output_surf_lw_down
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: conv_rain_profile, cond_rain_profile, output_heating_rate

  !-------------------------------

  ! Dimensions:
  TYPE(StrDim) :: dimen


  ! Atmospheric input:
  TYPE(StrAtm) :: atm_input

  !----------------------------------

  !-----------------------------------------------------------------------
CONTAINS

  !#######################################################################

  SUBROUTINE atmosphere_init ( Time_init, Time, Time_step )

#include <fv_arrays.h>

    TYPE (time_type), INTENT(in) :: Time_step
    TYPE (time_type), INTENT(in) :: Time_init
    TYPE (time_type), INTENT(in) :: Time
    LOGICAL cold_start

    ! local:
    INTEGER axes(4)
    INTEGER ss, ds

#include <fv_point.inc>
    !----- write version and namelist to log file -----

    CALL write_version_number ( version, tag )

    !---- compute physics/atmos time step in seconds ----

    Time_step_atmos = Time_step
    CALL get_time (Time_step_atmos, sec)
    dt_atmos = REAL(sec)

    !----- initialize FV dynamical core -----

    CALL fv_init( sec )
    CALL fv_restart( days, seconds )

    IF ( .NOT. cold_start ) THEN
       !
       ! Check consistency in Time
       !
       fv_time = set_time (seconds, days)
       CALL get_time (Time, ss,  ds)

       IF( seconds /= ss .OR. days /= ds ) THEN
          !       write(6,*) 'FMS:', ds, ss
          !       write(6,*) 'FV:', days, seconds
          CALL error_mesg('FV_init:', &
               'Time inconsistent between fv_rst and INPUT/atmos_model.res', &
               FATAL)
       ENDIF
    ELSE
       fv_time = time
    ENDIF

    CALL fv_diag_init( axes, Time )

    IF( nlev > 1 ) CALL hs_forcing_init ( axes, Time )

    !-----------------------------------------------------------------------
    !df 1401010
    ALLOCATE(t_surf      (1:nlon, beglat:endlat))
    ALLOCATE(fms_stellar_flux      (1:nlon, beglat:endlat))
    ALLOCATE(lats      (1:nlon, beglat:endlat))
    ALLOCATE(lons      (1:nlon, beglat:endlat))
    !t_surf(:,:) = tsini
    !----------------------------------------------------------------------
    isiz = nlon/4
    nx = nlon/isiz
    tsiz = nx * ( endlat - beglat + 1 )         ! Total loop length
    ! Calling phys one latitude at a time using a big/fat OpenMP loop
    ! For cache performance, this could be changed to finer decomposition if needed
    !$omp parallel do private(ij,i,j,k,m,is,ie,p_full,p_half)
    IF (.NOT. file_exist('INPUT/atmos_tracers.res.nc')) THEN
       DO ij=1,tsiz
          j  = beglat + (ij-1) / nx
          is = 1 + isiz * MOD(ij-1, nx)
          ie = is + isiz - 1
          t_surf(is:ie,j) = pt(is:ie,j,nlev) + 1.
       END DO
    ELSE
       DO ij=1,tsiz
          j  = beglat + (ij-1) / nx
          is = 1 + isiz * MOD(ij-1, nx)
          ie = is + isiz - 1
          t_surf(is:ie,j) = q(is:ie,j,nlev,2) !-20.
          IF (master) WRITE (6,*) "restart Ts"
       END DO
    END IF

    q(:,:,nlev,2) = 0. !1e-4

    ALLOCATE(p_half      (1:nlon, beglat:endlat, nlev+1))
    ALLOCATE(p_full      (1:nlon, beglat:endlat, nlev))
    ALLOCATE(flux_t      (1:nlon, beglat:endlat))
    ALLOCATE(flux_q      (1:nlon, beglat:endlat))
    ALLOCATE(flux_r      (1:nlon, beglat:endlat))
    ALLOCATE(conv_rain   (1:nlon, beglat:endlat))
    ALLOCATE(cond_rain   (1:nlon, beglat:endlat))
    ALLOCATE(pme         (1:nlon, beglat:endlat))
    ALLOCATE(net_surf_sw_down   (1:nlon, beglat:endlat))
    ALLOCATE(output_net_surf_sw_down   (1:nlon, beglat:endlat))
    ALLOCATE(surf_lw_down       (1:nlon, beglat:endlat))
    ALLOCATE(output_surf_lw_down       (1:nlon, beglat:endlat))
    ALLOCATE(conv_rain_profile  (1:nlon, beglat:endlat, nlev))
    ALLOCATE(cond_rain_profile  (1:nlon, beglat:endlat, nlev))
    ALLOCATE(rh                 (1:nlon, beglat:endlat, nlev))
    ALLOCATE(output_heating_rate   (1:nlon, beglat:endlat, nlev))

    p_half = 0.; p_full =0.
    !     ----- register diagnostic fields -----
    id_t_surf = register_diag_field(mod_name, 't_surf',        &
         axes(1:2), Time, 'Surface temperature','Kelvin')

    id_flux_t = register_diag_field(mod_name, 'flux_t',        &
         axes(1:2), Time, 'Sensible heat flux','W/m**2')

    id_flux_q = register_diag_field(mod_name, 'flux_q',        &
         axes(1:2), Time, 'Evaporation flux','W/m**2')

    id_flux_r = register_diag_field(mod_name, 'flux_r',        &
         axes(1:2), Time, 'Upward IR at surface','W/m**2')

    id_conv_rain = register_diag_field(mod_name, 'conv_rain',        &
         axes(1:2), Time, 'Convection rain','kg/s/m**2')

    id_cond_rain = register_diag_field(mod_name, 'cond_rain',        &
         axes(1:2), Time, 'Condensation rain','kg/s/m**2')

    id_pme = register_diag_field(mod_name, 'pme',        &
         axes(1:2), Time, 'P-E','kg/s/m**2')

    id_conv_rain_profile = register_diag_field(mod_name, 'conv_profile',  &
         axes(1:3), Time, 'Convection rain profile','kg/s/m**2')

    id_cond_rain_profile = register_diag_field(mod_name, 'cond_profile',  &
         axes(1:3), Time, 'Condensation rain profile','kg/s/m**2')

    id_rh = register_diag_field(mod_name, 'rh',        &
         axes(1:3), Time, 'Relative humidity','%')

    !-----------------------------------------------------------------------

    IF (beglat == 1) THEN
                                                                           
PRINT*, ' '
       PRINT*, '-----------------------------------'

PRINT*, ' '

       PRINT*,' ________                 ________  ____    ____   ______   '
       PRINT*,'|_   __  |               |_   __  ||_   \  /   _|.` ____ \  '
       PRINT*,'  | |_ \_| _   __   .--.   | |_ \_|  |   \/   |  | (___ \_| '
       PRINT*,'  |  _| _ [ \ [  ]/ .``\ \ |  _|     | |\  /| |   _.____`.  '
       PRINT*,' _| |__/ | > `  < | \__. |_| |_     _| |_\/_| |_ | \____) | '
       PRINT*,'|________|[__]`\_] `.__.`|_____|   |_____||_____| \______.` '

       PRINT*, ' '

       PRINT*, 'Initisalised ExoFMS v0.1'
       PRINT*, nlon, 'longitude points'
       PRINT*, nlev, 'vertical levels'

       PRINT*, ' '
       PRINT*, '-----------------------------------'
       PRINT*, ' '
    END IF

    CALL socrates_init(1, nlon, beglat, endlat, nlev, axes, Time,rlat(:,:))


  END SUBROUTINE atmosphere_init


  !#######################################################################

  SUBROUTINE atmosphere (Time)
#include <fv_arrays.h>
    TYPE(time_type), INTENT(in) :: Time

    ! local:
    REAL    zvir
    LOGICAL p_map                      ! Perform only partial remapping
#include <fv_point.inc>

    !df 141011 used in my physics schemes---------------------------------
    REAL, DIMENSION(1:nlon, beglat:endlat, nlev) :: tg_tmp, qg_tmp, cp, &
         u_tmp, v_tmp, rh_tmp, esat, &
         dt_ug, dt_vg, dt_tg, dt_qg, &
         diff_u, diff_v, atm_mass, atm_rho
    REAL, DIMENSION(1:nlon, beglat:endlat) :: q1, p1, dp1, dmass, lmass, rl, rt, ct, kappa1, t11, &
         zh, Ee, sigma, psg_tmp, Ek, dt_psg

    REAL, DIMENSION(1:nlon, beglat:endlat) ::      &
         gust,                 &   !gustiness constant
         flux_u,               &   ! surface flux of zonal mom.
         flux_v,               &   ! surface flux of meridional mom.
         drag_m,               &   ! momentum drag coefficient
         drag_t,               &   ! heat drag coefficient
         drag_q,               &   ! moisture drag coefficient
         w_atm,                &   ! wind speed
         ustar,                &   ! friction velocity
         bstar,                &   ! buoyancy scale
         qstar,                &   ! moisture scale
         dhdt_surf,            &   ! d(sensible heat flux)/d(surface temp)
         dedt_surf,            &   ! d(latent heat flux)/d(surface temp)???
         dedq_surf,            &   ! d(latent heat flux)/d(surface moisture)???
         drdt_surf,            &   ! d(upward longwave)/d(surface temp)
         dhdt_atm,             &   ! d(sensible heat flux)/d(atmos.temp)
         dedq_atm,             &   ! d(latent heat flux)/d(atmospheric mixing rat.)
         dtaudv_atm                ! d(stress component)/d(atmos wind)
    LOGICAL, DIMENSION(1:nlon, beglat:endlat) ::    land, avail

    REAL, DIMENSION(nlev) :: t_prev, q_prev, pfull_prev, t_after, q_after, &
         pfull_after, rain_profile, lnpfull_prev, &
         lnp_full, water, u_after, v_after, &
         tmean, qmean

    REAL, DIMENSION(nlev+1) :: phalf_prev, phalf_after, lnp_half

    REAL, DIMENSION(1:nlon, nlev) :: tm, qm, um, vm
    REAL, DIMENSION(1:nlon)       :: psm, tsm
    REAL :: psmean, tsmean

    INTEGER :: k,n
    INTEGER :: ngroup, nn, l
    !integer :: i, j, k, n
    !integer ij, nx, tsiz, isiz
    !integer is, ie

    !SOCRATES
    LOGICAL :: input_l_planet_grey_surface = .FALSE.
    REAL(r_def) :: input_planet_albedo = 0.06
    REAL(r_def) :: input_planet_emissivity = 0.94
    INTEGER(i_def), PARAMETER :: n_profile = 1!144!432
    INTEGER(i_def), PARAMETER :: n_layer = 40
    INTEGER(i_def) :: input_n_cloud_layer = 40
    INTEGER(i_def) :: input_n_aer_mode = 40
    INTEGER(i_def) :: input_cld_subcol_gen = 40
    INTEGER(i_def) :: input_cld_subcol_req = 40
    REAL(r_def) :: soc_stellar_constant = 3500000.0
    !INTEGER, PARAMETER :: n_layer = 40

    ! Fluxes and radiances calculated
    !  REAL  (RealK), ALLOCATABLE :: flux_diffuse_diown(:,:,:)
    !       Diffuse downward flux
    !  REAL  (RealK), ALLOCATABLE :: flux_net(:,:,:)
    !       Net flux
    !  REAL  (RealK), ALLOCATABLE :: heating_rate(:,:,:)
    !       Heating rates

    REAL, DIMENSION(n_profile,n_layer) :: input_p, input_t, input_mixing_ratio, &
         input_d_mass, input_density, input_layer_heat_capacity, &
         soc_heating_rate, input_o3_mixing_ratio
    REAL, DIMENSION(n_profile,0:n_layer) :: input_p_level, input_t_level, soc_flux_direct, &
         soc_flux_down, soc_flux_up, output_flux_net
    REAL, DIMENSION(n_profile) :: input_t_surf, input_cos_zenith_angle, input_solar_irrad, &
         input_orog_corr


    REAL :: Ep, delta_t_surf, Emass
    REAL :: lh, rain1
    REAL :: k3, k4, Ep_2, water_1, water_2
    REAL :: vcoeff, delta_t
    REAL :: ftop
    REAL gmean

    land = .FALSE.; avail = .TRUE.; gust = 1.
    dt_ug = 0.; dt_vg = 0.; dt_tg = 0.; dt_qg = 0.; dt_psg = 0.
    diff_u = 0.; diff_v = 0.
    delta_t = dt_atmos
    !-----------------------------------------------------------------------

    fv_time = Time + Time_step_atmos
    CALL get_time (fv_time, seconds,  days)

    !---- call fv dynamics -----
    zvir = 0.         ! no virtual effect if not full physics
    !  zvir = d608

    IF ( MOD(seconds, map_dt) == 0 ) THEN
       p_map = .FALSE.
    ELSE
       p_map = .TRUE.
    ENDIF
    CALL timing_on('fv_dynamics')

#ifndef USE_LIMA
    CALL fv_dynamics (nlon,    mlat,   nlev,    beglat,   endlat,    &
         ncnst,   pnats,  p_map,   consv_te,            &
         u,       v,      delp,    pt,       q,         &
         ps,      pe,     pk,      pkz,      phis,      &
         omga,    peln,   ptop,    omega,    sec,       &
         zvir,    cp_air, rdgas,   kappa,  radius, ua, va, fv_time )
#else
    CALL fv_dynamics (nlon,    mlat,   nlev,    beglat,   endlat,    &
         ncnst,   pnats,  p_map,   consv_te,            &
         u,       v,      delp,    pt,       q,         &
         ps,      pe,     pk,      pkz,      phis,      &
         omga,    peln,   ptop,    omega,    sec,       &
         zvir,    cp_air, rdgas,   kappa,  radius, ua, va )
#endif

    CALL timing_off('fv_dynamics')

    IF( nlev /=1 .AND. .NOT. adiabatic ) THEN
       CALL timing_on('FV_PHYS')
       !                call fv_phys ( fv_time, sec )
       !df physics 141011
       !----------------------------------------------------------------------
       isiz = nlon/4
       nx = nlon/isiz
       tsiz = nx * ( endlat - beglat + 1 )         ! Total loop length
       ! Calling phys one latitude at a time using a big/fat OpenMP loop
       ! For cache performance, this could be changed to finer decomposition if needed
       !$omp parallel do private(ij,i,j,k,m,is,ie,p_full,p_half)
       DO ij=1,tsiz
          j  = beglat + (ij-1) / nx
          is = 1 + isiz * MOD(ij-1, nx)
          ie = is + isiz - 1
          !       do k=1,nlev
          DO i=is,ie
             CALL get_eta_level(nlev, ps(i,j), p_full(i,j,:), p_half(i,j,:))
             u_tmp(i,j,:) = ua(i,j,:)
             v_tmp(i,j,:) = va(i,j,:)
             tg_tmp(i,j,:) = pt(i,j,:)
             qg_tmp(i,j,:) = q(i,j,:,1)
             psg_tmp(i,j) = ps(i,j)
          END DO
       END DO
       n = nlev
       cp = cp_air * (1 - qg_tmp) + cp_vapor * qg_tmp


       !-----------
       ! Socrates interface - qwe

       !Set tide-locked flux - should be set by namelist eventually!
       soc_stellar_constant = 3500000.0
       fms_stellar_flux = soc_stellar_constant*COS(rlat)*COS(rlon)
       WHERE (fms_stellar_flux < 0.0) fms_stellar_flux = 0.0

       ! Retrieve output_heating_rate, and downward surface SW and LW fluxes
       soc_mode = .TRUE.
       hires_mode = .FALSE.
       CALL socrates_interface(Time, rlat, rlon, soc_mode, hires_mode,    &
            tg_tmp, t_surf, p_full, p_half, n_profile, n_layer,     &
            output_heating_rate, output_net_surf_sw_down, output_surf_lw_down, fms_stellar_flux )

       ! NB - output_heating_rate wrong!
       !output_heating_rate(:,:,39:) = 0.0
       output_heating_rate(:,:,:5) = 0.1*output_heating_rate(:,:,:5)
       tg_tmp = tg_tmp + output_heating_rate*delta_t!, (/144,  40/)) * delta_

       surf_lw_down(:,:) = output_surf_lw_down(:,:)

       IF (1==2) THEN
          ! Retrieve output_heating_rate, and downward surface SW and LW fluxes
          soc_mode = .FALSE.
          CALL socrates_interface(Time, rlat, rlon, soc_mode, hires_mode,    &
               tg_tmp, t_surf, p_full, p_half, n_profile, n_layer,     &
               output_heating_rate, output_net_surf_sw_down, output_surf_lw_down, fms_stellar_flux )

          output_heating_rate(:,:,:5) = 0.1*output_heating_rate(:,:,:5)
          !output_heating_rate(:,:,39:) = 0.0
          tg_tmp = tg_tmp + output_heating_rate*delta_t!, (/144,  40/)) * delta_t


          net_surf_sw_down(:,:) = output_net_surf_sw_down(:,:)

          ! NB net_surf_sw_down and surf_lw_down have now been set
          ! THey are used below
          !-----------
       END IF

       !net_surf_sw_down = fms_stellar_flux

       omga(:,:,:5)=omga(:,:,:)*0.99
       u(:,:,:5)=u(:,:,:)*0.99
       v(:,:,:5)=v(:,:,:)*0.99
       !WHERE (tg_tmp < 800.0) tg_tmp = 800.0


       !------------------------------------------------------


       CALL surface_flux(       &
            tg_tmp(:,:,nlev),       &
            qg_tmp(:,:,nlev),       &
            u_tmp(:,:,nlev),            &
            v_tmp(:,:,nlev),            &
            p_full(:,:,nlev),       &
            dt_psg(:,:),            & !zatm
            p_half(:,:,nlev+1),     &
            t_surf(:,:),            &
            t_surf(:,:),            &
            dt_psg(:,:),            & !q_surf
            dt_psg(:,:),            & !u_surf
            dt_psg(:,:),            & !v_surf
            dt_psg(:,:),            & !rough_mom
            dt_psg(:,:),            & !rough_heat
            dt_psg(:,:),            & !rough_moist
            gust(:,:),              &
            flux_t(:,:),            &
            flux_q(:,:),            &
            flux_r(:,:),            &
            flux_u(:,:),                              &
            flux_v(:,:),                              &
            drag_m(:,:),                              &
            drag_t(:,:),                              &
            drag_q(:,:),                              &
            w_atm(:,:),                              &
            ustar(:,:),                              &
            bstar(:,:),                              &
            qstar(:,:),                              &
            dhdt_surf(:,:),                              &
            dedt_surf(:,:),                              &
            dedq_surf(:,:),                              &
            drdt_surf(:,:),                              &
            dhdt_atm(:,:),                              &
            dedq_atm(:,:),                              &
            dtaudv_atm(:,:),                              &
            delta_t,                              &
            land(:,:),                              &
            avail(:,:)  )

       ! boundary layer scheme
       tg_tmp(:,:,n) = tg_tmp(:,:,n) + flux_t * delta_t &
            * grav / (p_half(:,:,n+1) - p_half(:,:,n)) &
            / cp(:,:,n)


       !dt_ug(:, :, n) = - ug(:,:,n,previous) / 86400.
       !dt_vg(:, :, n) = - vg(:,:,n,previous) / 86400.  !* 2.
       !Ek = (p_half(:,:,n+1) - p_half(:,:,n)) / grav/ 86400.  &
       !     * (ug(:,:,n,previous)**2 + vg(:,:,n,previous)**2)  !t_surf
       !do i = is, ie
       !   do j = js, je
       !      if (flux_q(i,j) .le. 0) flux_q(i,j) = 0.
       !   end do
       !end do
       !Rayleigh damping
       vcoeff = -1./(1. - 0.7) / 86400.
       DO k = 1, nlev
          sigma(:,:) = p_full(:,:,k) / psg_tmp(:,:)
          WHERE (sigma(:,:) > 0.7)
             dt_ug(:,:,k) = u_tmp(:,:,k) *vcoeff *(sigma(:,:) - 0.7)
             dt_vg(:,:,k) = v_tmp(:,:,k) *vcoeff *(sigma(:,:) - 0.7)
             !   elsewhere (p_full(:,:,k) < 1e2)
             !      dt_ug(:,:,k) = -u_tmp(:,:,k) / 86400.
             !      dt_vg(:,:,k) = -v_tmp(:,:,k) / 86400.
          ELSEWHERE
             dt_ug(:,:,k) = 0.
             dt_vg(:,:,k) = 0.
          endwhere
       END DO
       dt_ug(:,:,1:5) = -u_tmp(:,:,1:5) / 86400. /10.
       dt_vg(:,:,1:5) = -v_tmp(:,:,1:5) / 86400. /10.



       tg_tmp = tg_tmp - delta_t / cp * &
            (dt_ug*(u_tmp(:,:,:)+dt_ug*delta_t/2.)+ &
            dt_vg*(v_tmp(:,:,:)+dt_vg*delta_t/2.) )
       u_tmp = u_tmp + dt_ug * delta_t
       v_tmp = v_tmp + dt_vg * delta_t
       WHERE (flux_q(:,:) < 0.0)
          flux_q(:,:) = 0.
       endwhere


       !evaporation, convection
       !$omp parallel do private(ij,i,j,k,m,is,ie,p_full,p_half)
       DO ij=1,tsiz
          j  = beglat + (ij-1) / nx
          is = 1 + isiz * MOD(ij-1, nx)
          ie = is + isiz - 1

          DO i=is,ie
             ! moist convection------------------------------------------------------
             t_prev = tg_tmp(i,j,:) *1.
             q_prev = qg_tmp(i,j,:) *1.
             phalf_prev = p_half(i,j,:) *1.
             pfull_prev = p_full(i,j,:) *1.
             !dry convection
             CALL dry_convection(t_prev,  pfull_prev, &
                  t_after)

             !large-scale condensation
             !prev is actually after
             ! update surface temperature and pressure

             ! TEST - 0*SWD, 0*LWD
             delta_t_surf = (surf_lw_down(i,j) + net_surf_sw_down(i,j)  &
                  - flux_t(i,j) - flux_r(i,j) ) &
                  * delta_t / rho_cp / mld !eff_heat_capacity
             !PRINT*, '-----------------------'
             !PRINT*, 't_surf'
             !PRINT*, t_surf(i,j)
             !PRINT*, 'lw'
             !    delta_t_surf = (surf_lw_down(i,j) + net_surf_sw_down(i,j))  &
             !                    * delta_t / rho_cp / mld !eff_heat_capacity

             !    delta_t_surf = (surf_lw_down(i,j) + net_surf_sw_down(i,j)-0.5*5.67e-8*t_surf(i,j)**4) &
             !                   * delta_t / 200000.0
             t_surf(i,j) = t_surf(i,j) + delta_t_surf
             !correct the energy imbalance due to vertical interpolation

             tg_tmp(i,j,:) = t_after
             !qg_tmp(i,j,:) = q_after
             !u_tmp(i,j,:) = u_after
             !v_tmp(i,j,:) = v_after
             !    p_half(i,j,:) = phalf_after
             !    p_full(i,j,:) = pfull_after
             CALL escomp(tg_tmp(i,j,:), esat(i,j,:))
             rh_tmp(i,j,:) = qg_tmp(i,j,:)*pfull_after /&
                  (0.622 + 0.378 * qg_tmp(i,j,:)) /esat(i,j,:) * 100

          END DO

          ps(is:ie,j)    = psg_tmp(is:ie,j)
          pt(is:ie,j,:)  = tg_tmp(is:ie,j,:)
          q(is:ie,j,:,1) = qg_tmp(is:ie,j,:)
          !u(is:ie,j,:)   = u_tmp(is:ie,j,:)
          !v(is:ie,j,:)   = v_tmp(is:ie,j,:)
          rh(is:ie,j,:)  = rh_tmp(is:ie,j,:)
          u_dt(is:ie,j,:) = (u_tmp(is:ie,j,:)-ua(is:ie,j,:))/delta_t
          v_dt(is:ie,j,:) = (v_tmp(is:ie,j,:)-va(is:ie,j,:))/delta_t
          t_dt(is:ie,j,:) = 0.
          q_dt(is:ie,j,:,:) = 0.
       END DO


       DO j=beglat,endlat
          IF (j==1) THEN
             !            psmean   = sum(ps(:,j+1), dim=1)/nlon
             tmean(:) = SUM(pt(:,j+1,:), dim=1)/nlon
             qmean(:) = SUM(q(:,j+1,:,1),dim=1)/nlon
             DO i=1,nlon
                !              ps(i,j)   = psmean
                pt(i,j,:) = tmean(:)*1.
                q(i,j,:,1)= qmean(:)*1.
             END DO
          ELSE IF (j==mlat) THEN
             !            psmean   = sum(ps(:,j-1), dim=1)/nlon
             tmean(:) = SUM(pt(:,j-1,:), dim=1) / nlon
             qmean(:) = SUM(q(:,j-1,:,1),dim=1)/nlon
             DO i=1,nlon
                !                ps(i,j)   = psmean
                pt(i,j,:) = tmean(:)*1.
                q(i,j,:,1)= qmean(:)*1.
             END DO
          END IF


          DO i=1,nlon
             DO k=1,nlev
                delp(i,j,k) = ak(k+1) - ak(k) + ps(i,j) * (bk(k+1) - bk(k))
             ENDDO
          ENDDO
       ENDDO
       CALL p_var(nlon, mlat, nlev, beglat, endlat, ptop, delp, ps,   &
            pe, peln, pk, pkz, kappa, q, ng_d, ncnst, .FALSE. )
       !----------------------------------------------------------------------
       CALL update_fv_phys ( delta_t, nt_phys, .FALSE., .FALSE., Time )
       !--------------------------------------------------------------
       CALL timing_off('FV_PHYS')
    ENDIF

    !---- diagnostics for FV dynamics -----

    CALL timing_on('FV_DIAG')

    CALL fv_diag(fv_time, nlon, mlat, nlev, beglat, endlat, ncnst, zvir,   &
         dt_atmos, .TRUE.)

    CALL timing_off('FV_DIAG')

    !--------------------------------------------------------
    !df 141010
    IF(id_t_surf > 0) used = send_data(id_t_surf, t_surf, Time)
    IF(id_flux_t > 0) used = send_data(id_flux_t, flux_t, Time)
    IF(id_flux_q > 0) used = send_data(id_flux_q, flux_q, Time)
    IF(id_flux_r > 0) used = send_data(id_flux_r, flux_r, Time)
    IF(id_conv_rain > 0) used = send_data(id_conv_rain, conv_rain, Time)
    IF(id_cond_rain > 0) used = send_data(id_cond_rain, cond_rain, Time)
    IF(id_pme       > 0) used = send_data(id_pme      , pme      , Time)
    IF(id_conv_rain_profile > 0) used = send_data(id_conv_rain_profile, conv_rain_profile, Time)
    IF(id_cond_rain_profile > 0) used = send_data(id_cond_rain_profile, cond_rain_profile, Time)

    IF(id_rh > 0)     used = send_data(id_rh, rh, Time)

    !--------------------------------------------------------
  END SUBROUTINE atmosphere


  SUBROUTINE atmosphere_end

#include <fv_arrays.h>
    isiz = nlon/4
    nx = nlon/isiz
    tsiz = nx * ( endlat - beglat + 1 )         ! Total loop length
    ! Calling phys one latitude at a time using a big/fat OpenMP loop
    ! For cache performance, this could be changed to finer decomposition if needed
    !$omp parallel do private(ij,i,j,k,m,is,ie,p_full,p_half)
    DO ij=1,tsiz
       j  = beglat + (ij-1) / nx
       is = 1 + isiz * MOD(ij-1, nx)
       ie = is + isiz - 1
       q(is:ie,j,nlev,2) = t_surf(is:ie,j)

       !!call socrates_hires_init(1, nlon, beglat, endlat, nlev, axes, Time,rlat(:,:))

       !!call socrates_hires_interface(Time, rlat, rlon,     &
       !!     tg_tmp, t_surf, p_full, p_half, n_profile, n_layer,     &
       !!     output_heating_rate, net_surf_sw_down, surf_lw_down, fms_stellar_flux )

    END DO

    !----- initialize domains for writing global physics data -----

    CALL set_domain ( fv_domain )
    CALL get_time (fv_time, seconds,  days)
    CALL write_fv_rst( 'RESTART/fv_rst.res', days, seconds, grav, &
         restart_format )


    CALL fv_end(days, seconds)

    DEALLOCATE(t_surf)
    DEALLOCATE(p_half, p_full)
    DEALLOCATE(flux_t, flux_q, flux_r)
    DEALLOCATE(conv_rain, cond_rain)
    DEALLOCATE(net_surf_sw_down, surf_lw_down)
    DEALLOCATE(conv_rain_profile, cond_rain_profile, rh)

    !    call radiation_end

  END SUBROUTINE atmosphere_end

END MODULE atmosphere_mod
