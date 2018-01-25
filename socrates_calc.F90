MODULE socrates_calc_mod

  ! Socrates calculation interface modules
  ! MDH added FMS diagnostics

  !----------
  !DIAG ExoFMS diagnostics
  USE    diag_manager_mod,   ONLY: register_diag_field, send_data

  ! ExoFMS time
  USE    time_manager_mod,   ONLY: time_type, &
       OPERATOR(+), OPERATOR(-), OPERATOR(/=)
  !----------
  IMPLICIT NONE

  !----------

CONTAINS

  ! ==================================================================================


  ! Set up the call to the Socrates radiation scheme
  ! -----------------------------------------------------------------------------
  SUBROUTINE socrates_calc(Time_diag,control, spectrum,                                    &
       n_profile, n_layer, n_cloud_layer, n_aer_mode,                               &
       cld_subcol_gen, cld_subcol_req,                                              &
       p_layer, t_layer, t_layer_boundaries, d_mass, density,                       &
       h2o, o3,                                                                     &
       t_rad_surf, cos_zenith_angle, solar_irrad, orog_corr,                        &
       l_planet_grey_surface, planet_albedo, planet_emissivity,                     &
       layer_heat_capacity,                                                         &
       flux_direct, flux_down, flux_up, heating_rate, spectral_olr)

    USE rad_pcf
    USE def_control,  ONLY: StrCtrl
    USE def_spectrum, ONLY: StrSpecData
    USE def_dimen,    ONLY: StrDim
    USE def_atm,      ONLY: StrAtm,   deallocate_atm
    USE def_bound,    ONLY: StrBound, deallocate_bound
    USE def_cld,      ONLY: StrCld,   deallocate_cld, deallocate_cld_prsc
    USE def_aer,      ONLY: StrAer,   deallocate_aer, deallocate_aer_prsc
    USE def_out,      ONLY: StrOut,   deallocate_out

    USE set_control_mod, ONLY: set_control
    USE set_dimen_mod,   ONLY: set_dimen
    USE set_atm_mod,     ONLY: set_atm
    USE set_bound_mod,   ONLY: set_bound
    USE set_cld_mod,     ONLY: set_cld
    USE set_aer_mod,     ONLY: set_aer

    USE soc_constants_mod,   ONLY: i_def, r_def

    IMPLICIT NONE

    !DIAG ExoFMS diagnostic time
    TYPE(time_type), INTENT(in)         :: Time_diag

    ! Spectral data:
    TYPE (StrSpecData), INTENT(in) :: spectrum

    ! Controlling options:
    TYPE (StrCtrl), INTENT(inout) :: control

    INTEGER(i_def), INTENT(in) :: n_profile
    !   Number of columns to operate on
    INTEGER(i_def), INTENT(in) :: n_layer
    !   Number of layers for radiation
    INTEGER(i_def), INTENT(in) :: n_cloud_layer
    !   Number of potentially cloudy layers
    INTEGER(i_def), INTENT(in) :: n_aer_mode
    !   Number of aerosol modes
    INTEGER(i_def), INTENT(in) :: cld_subcol_gen
    !   Number of sub-grid cloud columns generated
    INTEGER(i_def), INTENT(in) :: cld_subcol_req
    !   Number of sub-grid cloud columns required

    REAL(r_def), INTENT(in) :: p_layer(n_profile, n_layer)
    !   Pressure at layer centres
    REAL(r_def), INTENT(in) :: t_layer(n_profile, n_layer)
    !   Temperature at layer centres
    REAL(r_def), INTENT(in) :: t_layer_boundaries(n_profile, 0:n_layer)
    !   Temperature at layer boundaries
    REAL(r_def), INTENT(in) :: d_mass(n_profile, n_layer)
    !   Mass of layer (kg m-2)
    REAL(r_def), INTENT(in) :: density(n_profile, n_layer)
    !   Density of layer (kg m-3)
    REAL(r_def), INTENT(in) :: h2o(n_profile, n_layer)
    !   Mass mixing ratio of water vapour
    REAL(r_def), INTENT(in) :: o3(n_profile, n_layer)
    !   Mass mixing ratio of ozone

    REAL(r_def), INTENT(in) :: t_rad_surf(n_profile)
    !   Effective radiative temperature over whole grid-box
    REAL(r_def), INTENT(in) :: cos_zenith_angle(n_profile)
    !   Cosine of solar zenith angle
    REAL(r_def), INTENT(in) :: solar_irrad(n_profile)
    !   Solar irradiance at top-of-atmosphere (mean over timestep)
    REAL(r_def), INTENT(in) :: orog_corr(n_profile)
    !   Orographic correction factor

    LOGICAL, INTENT(in) :: l_planet_grey_surface
    !   Set a single grey albedo / emissivity for the surface
    REAL(r_def), INTENT(in) :: planet_albedo
    !   Surface albedo used for SW calculations
    REAL(r_def), INTENT(in) :: planet_emissivity
    !   Surface emissivity used for LW calculations

    REAL(r_def), INTENT(in) :: layer_heat_capacity(n_profile, n_layer)
    !   Heat capacity of layer

    REAL(r_def), INTENT(out) :: flux_direct(n_profile, 0:n_layer)
    !   Direct (unscattered) downwards flux (Wm-2)
    REAL(r_def), INTENT(out) :: flux_down(n_profile, 0:n_layer)
    !   Downwards flux (Wm-2)
    REAL(r_def), INTENT(out) :: flux_up(n_profile, 0:n_layer)
    !   Upwards flux (Wm-2)
    REAL(r_def), INTENT(out) :: heating_rate(n_profile, n_layer)
    !   Heating rate (Ks-1)

    REAL(r_def), INTENT(out) :: spectral_olr(n_profile, 135)
    !   Spectral OLR


    ! Dimensions:
    TYPE (StrDim) :: dimen

    ! Atmospheric properties:
    TYPE(StrAtm) :: atm

    ! Boundary conditions:
    TYPE(StrBound) :: bound

    ! Cloud properties:
    TYPE(StrCld) :: cld

    ! Aerosol properties:
    TYPE(StrAer) :: aer

    ! Output fields from core radiation code:
    TYPE(StrOut) :: radout

    INTEGER(i_def) :: l, i
    !   Loop variablesi

    !DIAG Diagnostic
    LOGICAL :: used



    CALL set_control(control)

    CALL set_dimen(control, dimen, spectrum, n_profile, n_layer,                   &
         n_cloud_layer, n_aer_mode, cld_subcol_gen, cld_subcol_req)

    CALL set_atm(control, dimen, spectrum, atm, n_profile, n_layer,                &
         p_layer, t_layer, t_layer_boundaries, d_mass, density, h2o, o3)

    CALL set_bound(control, dimen, spectrum, bound, n_profile,                     &
         t_rad_surf, cos_zenith_angle, solar_irrad, orog_corr,                        &
         l_planet_grey_surface, planet_albedo, planet_emissivity)

    CALL set_cld(control, dimen, spectrum, cld, n_profile)

    CALL set_aer(control, dimen, spectrum, aer, n_profile)


    ! DEPENDS ON: radiance_calc
    CALL radiance_calc(control, dimen, spectrum, atm, cld, aer, bound, radout)


    ! set heating rates and diagnostics
    DO l=1, n_profile
       DO i=1, n_layer
          heating_rate(l, i) = (radout%flux_down(l,i-1,1)-radout%flux_down(l,i,1)    &
               + radout%flux_up(l,i,1)-radout%flux_up(l,i-1,1))       &
               / layer_heat_capacity(l, i)
       END DO
    END DO



    DO l=1, n_profile
       DO i=0, n_layer
          flux_direct(l, i) = radout%flux_direct(l, i, 1)
          flux_down(l, i)   = radout%flux_down(l, i, 1)
          flux_up(l, i)     = radout%flux_up(l, i, 1)
       END DO
       spectral_olr(l,:) = radout%flux_up_clear_band(l,0,:)
    END DO


    CALL deallocate_out(radout)
    CALL deallocate_aer_prsc(aer)
    CALL deallocate_aer(aer)
    CALL deallocate_cld_prsc(cld)
    CALL deallocate_cld(cld)
    CALL deallocate_bound(bound)
    CALL deallocate_atm(atm)

  END SUBROUTINE socrates_calc
END MODULE socrates_calc_mod
