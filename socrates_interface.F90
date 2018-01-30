MODULE socrates_interface_mod

  ! Socrates calculation interface module
  ! Takes FMS time, spectra, temperature, and pressure
  ! Outputs FMS heating rate, and downwards surface LW and SW
  ! MDH 30/01/18

  !----------

  ! ExoFMS diagnostics
  USE  diag_manager_mod, ONLY: register_diag_field, send_data

  ! ExoFMS time
  USE time_manager_mod, ONLY: time_type, OPERATOR(+), OPERATOR(-), OPERATOR(/=)

  ! Socrates modules
  USE read_control_mod
  USE def_control, ONLY: StrCtrl,  allocate_control,   deallocate_control
  USE def_spectrum

  ! MDH Socrates modules
  USE tau_output_mod
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

  ! Socrates inputs from namelist
  REAL :: stellar_constant = 1370.0
  LOGICAL :: tidally_locked = .TRUE.
  NAMELIST/socrates_nml/ stellar_constant, tidally_locked



CONTAINS

  SUBROUTINE socrates_init(is, ie, js, je, num_levels, axes, Time, lat)
    !! Initialises Socrates spectra, arrays, and constants

    ! Arguments
    INTEGER, INTENT(in), DIMENSION(4) :: axes
    !! NB axes refers to the handles of the axes defined in fv_diagnostics
    TYPE(time_type), INTENT(in)       :: Time
    INTEGER, INTENT(in)               :: is, ie, js, je, num_levels
    REAL, INTENT(in) , DIMENSION(:,:)   :: lat
    !-------------------------------------------------------------------------------------

    ! Socrates spectral files -- should be set by namelist
    control_lw%spectral_file = '~/Work/spec_file_co_gcm'
    control_sw%spectral_file = '~/Work/spec_file_co_gcm'

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



    ! Print Socrates init data from one processor
    IF (js == 1) THEN
       PRINT*, ' '
       PRINT*, '-----------------------------------'

       PRINT*,'  ____                       _            '
       PRINT*,' / ___|  ___   ___ _ __ __ _| |_ ___  ___ '
       PRINT*,' \___ \ / _ \ / __|  __/ _` | __/ _ \/ __|'
       PRINT*,'  ___) | (_) | (__| | | (_| | ||  __/\__ \'
       PRINT*,' |____/ \___/ \___|_|  \__,_|\__\___||___/'

       PRINT*, ' '

       PRINT*, 'Initialised Socrates vX.x'
       PRINT*, 'Stellar constant = ', stellar_constant
       PRINT*, 'Longwave spectral file = ', TRIM(control_lw%spectral_file)
       PRINT*, 'Shortwave spectral file = ', TRIM(control_sw%spectral_file)
       PRINT*, ' '
       PRINT*, '-----------------------------------'
       PRINT*, ' '
    end if

    return
  end subroutine socrates_init
  ! ==================================================================================


  ! Set up the call to the Socrates radiation scheme
  ! -----------------------------------------------------------------------------
  subroutine socrates_interface(Time_diag, rlat, rlon, soc_mode,       &
       fms_temp, fms_t_surf, fms_p_full, fms_p_half, n_profile, n_layer,        &
       output_heating_rate, fms_net_surf_sw_down, fms_surf_lw_down, fms_stellar_flux )

    use realtype_rd
    use soc_constants_mod
    use read_control_mod
    use socrates_calc_mod
    use compress_spectrum_mod
    use def_spectrum
    use def_dimen,   only: StrDim
    use def_control, only: StrCtrl,  allocate_control,   deallocate_control
    use def_atm,     only: StrAtm,   allocate_atm,       deallocate_atm
    use def_cld,     only: StrCld,   allocate_cld,       deallocate_cld, &
         allocate_cld_prsc,  deallocate_cld_prsc, &
         allocate_cld_mcica, deallocate_cld_mcica
    use def_aer,     only: StrAer,   allocate_aer,       deallocate_aer, &
         allocate_aer_prsc,  deallocate_aer_prsc
    use def_bound,   only: StrBound, allocate_bound,     deallocate_bound
    use def_out,     only: StrOut,                       deallocate_out
    !-----------------------------------------------------------------------
    implicit none

    ! Input time
    type(time_type), intent(in)         :: Time_diag

    INTEGER(i_def), intent(in) :: n_profile, n_layer
    logical, intent(in) :: soc_mode
    INTEGER(i_def) :: nlat

    ! Input arrays
    real(r_def), intent(in) :: fms_temp(:,:,:)
    real(r_def), intent(in) :: fms_p_full(:,:,:)
    real(r_def), intent(in) :: fms_p_half(:,:,:)
    real(r_def), intent(in) :: fms_t_surf(:,:)
    real(r_def), intent(in) :: fms_stellar_flux(:,:)
    real(r_def), intent(in) :: rlon(:,:)
    real(r_def), intent(in) :: rlat(:,:)

    ! Output arrays
    real(r_def), intent(out) :: fms_net_surf_sw_down(:,:)
    real(r_def), intent(out) :: fms_surf_lw_down(:,:)
    real(r_def), intent(out) :: output_heating_rate(:,:,:)
    real(r_def) :: output_heating_rate_lw(144,3,40)
    real(r_def) :: output_heating_rate_sw(144,3,40)
    real(r_def) :: output_soc_flux_up_lw(144,3,40)
    real(r_def) :: output_soc_flux_down_sw(144,3,40)
    real(r_def) :: output_soc_spectral_olr(144,3,621)

    ! Hi-res output
    INTEGER, PARAMETER :: out_unit=20
    CHARACTER(len=200) :: file_name

    ! Arrays to send to Socrates
    real, dimension(n_layer) :: input_p, input_t, input_mixing_ratio, &
         input_d_mass, input_density, input_layer_heat_capacity, &
         soc_heating_rate, input_o3_mixing_ratio, &
         soc_heating_rate_lw, soc_heating_rate_sw
    real, dimension(0:n_layer) :: input_p_level, input_t_level, soc_flux_direct, &
         soc_flux_down_sw, soc_flux_up_sw, output_flux_net, &
         soc_flux_down_lw, soc_flux_up_lw
    real, dimension(n_profile) :: input_t_surf, input_cos_zenith_angle, input_solar_irrad, &
         input_orog_corr


    ! Socrates options
    logical :: input_l_planet_grey_surface = .true.
    real(r_def) :: input_planet_albedo = 0.06
    real(r_def) :: input_planet_emissivity = 0.5
    integer(i_def) :: input_n_cloud_layer
    integer(i_def) :: input_n_aer_mode
    integer(i_def) :: input_cld_subcol_gen
    integer(i_def) :: input_cld_subcol_req


    ! Dimensions:
    type(StrDim) :: dimen
    type(StrAtm) :: atm_input

    ! Loop variables
    integer(i_def) :: i, j, k, l, n, lon

    !DIAG Diagnostic
    logical :: used

    !----------------------------

    ! Set array sizes
    input_n_cloud_layer = n_layer
    input_n_aer_mode = n_layer
    input_cld_subcol_gen = n_layer
    input_cld_subcol_req = n_layer



    do lon = 1,144
       do n = 1, 3

          nlat = n

          !Set input T, p, p_level, and mixing ratio profiles
          input_t = fms_temp(lon,nlat,:)
          input_p = fms_p_full(lon,nlat,:)
          input_p_level = fms_p_half(lon,nlat,:)

          ! TODO: Remove or edit
          input_mixing_ratio = 0.0!
          input_o3_mixing_ratio = 0.0!


          !-------------

          !Default parameters
          input_cos_zenith_angle = 0.7
          input_orog_corr = 0.0
          input_layer_heat_capacity = 29.07

          !Set tide-locked flux - should be set by namelist eventually!
          input_solar_irrad = fms_stellar_flux(lon,nlat)!RESHAPE(fms_stellar_flux, (/n_profile/))
          input_t_surf = fms_t_surf(lon,nlat)!RESHAPE(fms_t_surf, (/n_profile/))


          !--------------

          !Set input t_level by scaling t - NEEDS TO CHANGE!
          DO i = nlat, nlat
             DO k = 0,n_layer
                input_t_level(k) = 0.5*(input_t(k+1)+input_t(k))
             END DO
             input_t_level(40) = input_t(40) + input_t(40) - input_t_level(39)
             input_t_level(0) = input_t(1) - (input_t_level(1) - input_t(1))
          END DO



          !Set input dry mass, density, and heat capacity profiles
          DO i=n_layer, 1, -1
             input_d_mass(i) = (input_p_level(i)-input_p_level(i-1))/23.0
             input_density(i) = input_p(i)/(8.31*input_t(i))!1000.!atm%p(l ,i) / 1000
             !KLUDGE
             input_layer_heat_capacity(i) = input_d_mass(i)*1303.1!17.0*1005.0
          END DO


          ! Zero heating rate
          soc_heating_rate = 0.0
          soc_heating_rate_lw = 0.0
          soc_heating_rate_sw = 0.0

          if (soc_mode == .TRUE.) then
             control_lw%isolir = 2
             CALL read_control(control_lw, spectrum_lw)
             !CALL compress_spectrum(control_lw, spectrum_lw)

             CALL socrates_calc(Time_diag, control_lw, spectrum_lw,                                          &
                  n_profile, n_layer, input_n_cloud_layer, input_n_aer_mode,                   &
                  input_cld_subcol_gen, input_cld_subcol_req,                                  &
                  input_p, input_t, input_t_level, input_d_mass, input_density,                &
                  input_mixing_ratio, input_o3_mixing_ratio,                                      &
                  input_t_surf, input_cos_zenith_angle, input_solar_irrad, input_orog_corr,    &
                  input_l_planet_grey_surface, input_planet_albedo, input_planet_emissivity,   &
                  input_layer_heat_capacity,                                                   &
                  soc_flux_direct, soc_flux_down_lw, soc_flux_up_lw, soc_heating_rate_lw)

             ! Set output arrays
             fms_surf_lw_down(lon,nlat) = soc_flux_down_lw(40)


             output_heating_rate_lw(lon,nlat,:) = soc_heating_rate_lw
             output_soc_flux_up_lw(lon,nlat,:) = soc_flux_up_lw
             output_heating_rate(lon,nlat,:) = soc_heating_rate_lw


          endif
          !--------------


          if (soc_mode == .FALSE.) then
             ! SW calculation
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

          endif

       end do


       !--------------

    end do

    ! Call routine to write out stored tau values
    IF (1==2) THEN
       CALL write_tau_output(INT(100.0*(1.8+rlat(1,1))))
    END IF

    ! Write spectral OLR column-by-column
    ! TODO Improve structure
    IF (1==2) THEN

       WRITE(file_name,'(a,i4.4,a)') "/network/group/aopp/testvol2/plan/fms-scratch-mdh/SPEC_OLR_PURE_CO_",INT(100.0*(1.8+rlat(1,1))),".TXT"
       ! Open file
       OPEN (unit=out_unit,file=TRIM(file_name),action="write",status="replace")

       ! Iterate through lons and lats, writing OLR spectrum
       DO j = 1, 144
          DO i = 1, 3
             WRITE (out_unit,*) output_soc_spectral_olr(j,i,:)
          END DO
       END DO
       ! Close file
       CLOSE (out_unit)
       !--------------
    END IF


    ! Send diagnostics
    if (soc_mode == .TRUE.) then
       used = send_data ( id_soc_heating_lw, output_heating_rate_lw, Time_diag)
       used = send_data ( id_soc_flux_up_lw, output_soc_flux_up_lw, Time_diag)
    else
       used = send_data ( id_soc_heating_sw, output_heating_rate_sw, Time_diag)
       used = send_data ( id_soc_flux_down_sw, output_soc_flux_down_sw, Time_diag)
    endif


  end subroutine socrates_interface
end module socrates_interface_mod
