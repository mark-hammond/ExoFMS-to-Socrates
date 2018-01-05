MODULE socrates_interface_mod

  ! Socrates calculation interface module
  ! Takes FMS time, spectra, temperature, and pressure
  ! Outputs FMS heating rate, and downwards surface LW and SW
  ! MDH 07/12/17

  !----------

  ! ExoFMS diagnostics
  USE  diag_manager_mod, ONLY: register_diag_field, send_data

  ! ExoFMS time
  USE time_manager_mod, ONLY: time_type, OPERATOR(+), OPERATOR(-), OPERATOR(/=)

  ! Socrates modules
  USE read_control_mod
  USE def_control, ONLY: StrCtrl,  allocate_control,   deallocate_control
  USE def_spectrum

  IMPLICIT NONE

  ! Input spectra
  TYPE (StrSpecData) :: spectrum_lw
  TYPE (StrSpecData) :: spectrum_sw

  ! Control options:
  TYPE(StrCtrl) :: control
  TYPE(StrCtrl) :: control_sw
  TYPE(StrCtrl) :: control_lw

  ! Diagnostic IDs, name, and missing value
  INTEGER :: id_soc_olr, id_soc_olr_spectrum_lw, id_soc_surf_spectrum_sw
  INTEGER :: id_soc_heating_sw, id_soc_heating_lw, id_soc_heating_rate
  INTEGER :: id_soc_flux_up_lw, id_soc_flux_down_sw
  CHARACTER(len=10), PARAMETER :: soc_mod_name = 'socrates'
  REAL :: missing_value = -999


CONTAINS

  SUBROUTINE socrates_init(is, ie, js, je, num_levels, axes, Time, lat)
    !! Initialises Socrates spectra, arrays, and constants

    ! Arguments
    INTEGER, INTENT(in), DIMENSION(4) :: axes
    !! NB axes refers to the handles of the axes defined in fv_diagnostics
    TYPE(time_type), INTENT(in)       :: Time
    INTEGER, INTENT(in)               :: is, ie, js, je, num_levels
    REAL, INTENT(in) , DIMENSION(:,:)   :: lat


    ! Socrates spectral files -- should be set by namelist
    control_lw%spectral_file = '~/spec_file_co2_co_lowres'
    control_sw%spectral_file = '~/spec_file_co2_co_lowres'

    ! Read in spectral files
    CALL read_spectrum(control_lw%spectral_file,spectrum_lw)
    CALL read_spectrum(control_sw%spectral_file,spectrum_sw)

    ! Set Socrates configuration
    CALL read_control(control_lw,spectrum_lw)
    CALL read_control(control_sw,spectrum_sw)

    ! Specify LW and SW setups
    control_sw%isolir=1
    control_lw%isolir=2


    ! Register diagostic fields
    id_soc_olr = &
         register_diag_field ( soc_mod_name, 'soc_olr', axes(1:2), Time, &
         'outgoing longwave radiation', &
         'watts/m2', missing_value=missing_value               )

    id_soc_olr_spectrum_lw = &
         register_diag_field ( soc_mod_name, 'soc_olr_spectrum_lw',(/ axes(1:2), axes(5)/) , Time, &
         'socrates substellar LW OLR spectrum', &
         'watts/m2', missing_value=missing_value               )

    id_soc_surf_spectrum_sw = &
         register_diag_field ( soc_mod_name, 'soc_surf_spectrum_sw',(/ axes(1:2), axes(5)/) , Time, &
         'socrates substellar SW surface spectrum', &
         'watts/m2', missing_value=missing_value               )

    id_soc_heating_lw = &
         register_diag_field ( soc_mod_name, 'soc_heating_lw', axes(1:3), Time, &
         'socrates LW heating rate', &
         'J/s', missing_value=missing_value               )

    id_soc_heating_sw = &
         register_diag_field ( soc_mod_name, 'soc_heating_sw', axes(1:3), Time, &
         'socrates SW heating rate', &
         'J/s', missing_value=missing_value               )

    id_soc_heating_rate = &
         register_diag_field ( soc_mod_name, 'soc_heating_rate', axes(1:3), Time, &
         'socrates total heating rate', &
         'J/s', missing_value=missing_value               )

    id_soc_flux_up_lw = &
         register_diag_field ( soc_mod_name, 'soc_flux_up_lw', axes(1:3), Time, &
         'socrates LW flux up', &
         'watts/m2', missing_value=missing_value               )

    id_soc_flux_down_sw = &
         register_diag_field ( soc_mod_name, 'soc_flux_down_sw', axes(1:3), Time, &
         'socrates SW flux down', &
         'watts/m2', missing_value=missing_value               )

    RETURN
  END SUBROUTINE socrates_init
  ! ==================================================================================


  ! Set up the call to the Socrates radiation scheme
  ! -----------------------------------------------------------------------------
  SUBROUTINE socrates_interface(Time_diag, rlat, rlon, soc_mode,       &
       fms_temp, fms_t_surf, fms_p_full, fms_p_half, n_profile, n_layer,        &
       output_heating_rate, fms_net_surf_sw_down, fms_surf_lw_down, fms_stellar_flux )

    USE realtype_rd
    USE soc_constants_mod
    USE read_control_mod
    USE socrates_calc_mod
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

    ! Input time
    TYPE(time_type), INTENT(in)         :: Time_diag

    INTEGER(i_def), INTENT(in) :: n_profile, n_layer
    LOGICAL, INTENT(in) :: soc_mode
    INTEGER(i_def) :: nlat

    ! Input arrays
    REAL(r_def), INTENT(in) :: fms_temp(:,:,:)
    REAL(r_def), INTENT(in) :: fms_p_full(:,:,:)
    REAL(r_def), INTENT(in) :: fms_p_half(:,:,:)
    REAL(r_def), INTENT(in) :: fms_t_surf(:,:)
    REAL(r_def), INTENT(in) :: fms_stellar_flux(:,:)
    REAL(r_def), INTENT(in) :: rlon(:,:)
    REAL(r_def), INTENT(in) :: rlat(:,:)

    ! Output arrays
    REAL(r_def), INTENT(out) :: fms_net_surf_sw_down(:,:)
    REAL(r_def), INTENT(out) :: fms_surf_lw_down(:,:)
    REAL(r_def), INTENT(out) :: output_heating_rate(:,:,:)
    REAL(r_def) :: output_heating_rate_lw(144,3,40)
    REAL(r_def) :: output_heating_rate_sw(144,3,40)
    REAL(r_def) :: output_soc_flux_up_lw(144,3,40)
    REAL(r_def) :: output_soc_flux_down_sw(144,3,40)


    ! Arrays to send to Socrates
    REAL, DIMENSION(n_layer) :: input_p, input_t, input_mixing_ratio, &
         input_d_mass, input_density, input_layer_heat_capacity, &
         soc_heating_rate, input_o3_mixing_ratio, &
         soc_heating_rate_lw, soc_heating_rate_sw
    REAL, DIMENSION(0:n_layer) :: input_p_level, input_t_level, soc_flux_direct, &
         soc_flux_down_sw, soc_flux_up_sw, output_flux_net, &
         soc_flux_down_lw, soc_flux_up_lw
    REAL, DIMENSION(n_profile) :: input_t_surf, input_cos_zenith_angle, input_solar_irrad, &
         input_orog_corr


    ! Socrates options
    LOGICAL :: input_l_planet_grey_surface = .TRUE.
    REAL(r_def) :: input_planet_albedo = 0.06
    REAL(r_def) :: input_planet_emissivity = 0.5
    INTEGER(i_def) :: input_n_cloud_layer
    INTEGER(i_def) :: input_n_aer_mode
    INTEGER(i_def) :: input_cld_subcol_gen
    INTEGER(i_def) :: input_cld_subcol_req

    !-------------------------------
    ! Socrates input options -- to become namelist
    REAL :: soc_stellar_constant
    LOGICAL :: soc_tide_locked
    NAMELIST/socrates_nml/ soc_tide_locked, soc_stellar_constant

    ! Dimensions:
    TYPE(StrDim) :: dimen
    TYPE(StrAtm) :: atm_input

    ! Loop variables
    INTEGER(i_def) :: i, j, k, l, n, lon

    !DIAG Diagnostic
    LOGICAL :: used

    !----------------------------

    ! Set array sizes
    input_n_cloud_layer = n_layer
    input_n_aer_mode = n_layer
    input_cld_subcol_gen = n_layer
    input_cld_subcol_req = n_layer


    DO lon = 1,144
       DO n = 1, 3

          nlat = n

          !Set input T, p, p_level, and mixing ratio profiles
          input_t = fms_temp(lon,nlat,:)
          input_p = fms_p_full(lon,nlat,:)
          input_p_level = fms_p_half(lon,nlat,:)
          input_mixing_ratio = 1.E-3
          input_o3_mixing_ratio = 1.E-6



          !-------------

          !Default parameters
          input_cos_zenith_angle = 0.7
          input_orog_corr = 0.0
          input_layer_heat_capacity = 29.07

          !Set tide-locked flux - should be set by namelist eventually!
          input_solar_irrad = fms_stellar_flux(lon,nlat)
          input_t_surf = fms_t_surf(lon,nlat)


          !--------------

          ! Set input t_level by scaling t - NEEDS TO CHANGE!
          DO i = nlat, nlat
             DO k = 0,n_layer
                input_t_level(k) = 0.5*(input_t(k+1)+input_t(k))
             END DO
             input_t_level(40) = input_t(40) + input_t(40) - input_t_level(39)
             input_t_level(0) = input_t(1) - (input_t_level(1) - input_t(1))
          END DO



          ! Set input dry mass, density, and heat capacity profiles
          DO i=n_layer, 1, -1
             input_d_mass(i) = (input_p_level(i)-input_p_level(i-1))/23.0
             input_density(i) = input_p(i)/(8.31*input_t(i))
             input_layer_heat_capacity(i) = input_d_mass(i)*1005.0
          END DO


          ! Zero heating rate
          soc_heating_rate = 0.0
          soc_heating_rate_lw = 0.0
          soc_heating_rate_sw = 0.0

          ! Check if longwave or shortwave mode
          IF (soc_mode == .TRUE.) THEN
             control_lw%isolir = 2
             CALL read_control(control_lw, spectrum_lw)
             CALL socrates_calc(Time_diag, control_lw, spectrum_lw,                            &
                  n_profile, n_layer, input_n_cloud_layer, input_n_aer_mode,                   &
                  input_cld_subcol_gen, input_cld_subcol_req,                                  &
                  input_p, input_t, input_t_level, input_d_mass, input_density,                &
                  input_mixing_ratio, input_o3_mixing_ratio,                                   &
                  input_t_surf, input_cos_zenith_angle, input_solar_irrad, input_orog_corr,    &
                  input_l_planet_grey_surface, input_planet_albedo, input_planet_emissivity,   &
                  input_layer_heat_capacity,                                                   &
                  soc_flux_direct, soc_flux_down_lw, soc_flux_up_lw, soc_heating_rate_lw)

             ! Set output arrays
             fms_surf_lw_down(lon,nlat) = soc_flux_down_lw(40)
             output_heating_rate_lw(lon,nlat,:) = soc_heating_rate_lw
             output_soc_flux_up_lw(lon,nlat,:) = soc_flux_up_lw
             output_heating_rate(lon,nlat,:) = soc_heating_rate_lw

          ENDIF
          !--------------

! Shortwave mode
          IF (soc_mode == .FALSE.) THEN
             control_sw%isolir = 1
             CALL read_control(control_sw, spectrum_sw)

             CALL socrates_calc(Time_diag, control_sw, spectrum_sw,                                          &
                  n_profile, n_layer, input_n_cloud_layer, input_n_aer_mode,                   &
                  input_cld_subcol_gen, input_cld_subcol_req,                                  &
                  input_p, input_t, input_t_level, input_d_mass, input_density,                &
                  input_mixing_ratio, input_o3_mixing_ratio,                                      &
                  input_t_surf, input_cos_zenith_angle, input_solar_irrad, input_orog_corr,    &
                  input_l_planet_grey_surface, input_planet_albedo, input_planet_emissivity,   &
                  input_layer_heat_capacity,                                                   &
                  soc_flux_direct, soc_flux_down_sw, soc_flux_up_sw, soc_heating_rate_sw)

             ! Set output arrays
             output_heating_rate_sw(lon,nlat,:) = soc_heating_rate_sw
             fms_net_surf_sw_down(lon,nlat) = soc_flux_down_sw(40)
             output_heating_rate(lon,nlat,:) = soc_heating_rate_sw
             output_soc_flux_down_sw(lon,nlat,:) = soc_flux_down_sw

          ENDIF

       END DO

    END DO

    output_heating_rate(:,:,1) = 0.0

    ! Send diagnostics
    IF (soc_mode == .TRUE.) THEN
       used = send_data ( id_soc_heating_lw, output_heating_rate_lw, Time_diag)
       used = send_data ( id_soc_flux_up_lw, output_soc_flux_up_lw, Time_diag)
    ELSE
       used = send_data ( id_soc_heating_sw, output_heating_rate_sw, Time_diag)
       used = send_data ( id_soc_flux_down_sw, output_soc_flux_down_sw, Time_diag)
    ENDIF


  END SUBROUTINE socrates_interface
END MODULE socrates_interface_mod
