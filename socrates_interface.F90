module socrates_interface_mod

! Socrates calculation interface module
! Takes FMS time, spectra, temperature, and pressure
! Outputs FMS heating rate, and downwards surface LW and SW
! MDH 07/12/17

!----------

! ExoFMS diagnostics
use  diag_manager_mod, only: register_diag_field, send_data

! ExoFMS time
use time_manager_mod, only: time_type, operator(+), operator(-), operator(/=)

! Socrates modules
use read_control_mod
use def_control, only: StrCtrl,  allocate_control,   deallocate_control
use def_spectrum

implicit none

! Input spectra
type (StrSpecData) :: spectrum_lw
type (StrSpecData) :: spectrum_sw

! Control options:
type(StrCtrl) :: control
type(StrCtrl) :: control_sw
type(StrCtrl) :: control_lw

! Diagnostic IDs, name, and missing value
integer :: id_soc_olr, id_soc_olr_spectrum_lw, id_soc_surf_spectrum_sw
integer :: id_soc_heating_sw, id_soc_heating_lw, id_soc_heating_rate
integer :: id_soc_flux_up_lw, id_soc_flux_down_sw
character(len=10), parameter :: soc_mod_name = 'socrates'
real :: missing_value = -999

! Socrates inputs from namelist
!real :: stellar_flux = 1370.0
!logical :: tidally_locked = .true.
!namelist/socrates_nml/ stellar_flux, tide_locked

contains

subroutine socrates_init(is, ie, js, je, num_levels, axes, Time, lat)
      !! Initialises Socrates spectra, arrays, and constants

! Arguments
integer, intent(in), dimension(4) :: axes
      !! NB axes refers to the handles of the axes defined in fv_diagnostics
type(time_type), intent(in)       :: Time
integer, intent(in)               :: is, ie, js, je, num_levels
real, intent(in) , dimension(:,:)   :: lat
!-------------------------------------------------------------------------------------

! Read in namelist
!unit = open_file ('input.nml', action='read')
!ierr=1
!do while (ierr /= 0)
!   read  (unit, nml=radiation_nml, iostat=io, end=10)
!   ierr = check_nml_error (io, 'socrates_nml')
!enddo
!10 call close_file (unit)

!-----------------------------------------------------------------------

! Socrates spectral files -- should be set by namelist
control_lw%spectral_file = '~/Work/spec_file_co_gcm'
control_sw%spectral_file = '~/Work/spec_file_co_gcm'

! Read in spectral files
call read_spectrum(control_lw%spectral_file,spectrum_lw)
call read_spectrum(control_sw%spectral_file,spectrum_sw)

! Set Socrates configuration
call read_control(control_lw,spectrum_lw)
call read_control(control_sw,spectrum_sw)

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

!-------------------------------
! Socrates input options -- to become namelist
real :: soc_stellar_constant
logical :: soc_tide_locked
namelist/socrates_nml/ soc_tide_locked, soc_stellar_constant

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
input_t = fms_temp(lon,nlat,:)!RESHAPE(fms_temp, (/432, 40/))
input_p = fms_p_full(lon,nlat,:)!RESHAPE(fms_p_full, (/432, 40/))
input_p_level = fms_p_half(lon,nlat,:)!RESHAPE(fms_p_half, (/432, 41/))
input_mixing_ratio = 0.0!00000001!10000.0!0.01
input_o3_mixing_ratio = 0.0!10000.0!0.01


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
fms_surf_lw_down(lon,nlat) = soc_flux_down_lw(40)!RESHAPE(soc_flux_down(:,40) , (/144,3/))
!fms_surf_lw_down = fms_surf_lw_down * 0.0


output_heating_rate_lw(lon,nlat,:) = soc_heating_rate_lw
output_soc_flux_up_lw(lon,nlat,:) = soc_flux_up_lw

output_heating_rate(lon,nlat,:) = soc_heating_rate_lw


!used = send_data ( id_soc_heating_lw, output_heating_rate_lw, Time_diag)
!used = send_data ( id_soc_flux_up_lw, RESHAPE(soc_flux_up_lw, (/144,3,40/)), Time_diag)
!used = send_data ( id_soc_flux_down_lw, RESHAPE(soc_flux_down_lw, (/144,3,40/)), Time_diag)
 !   used = send_data ( id_soc_heating_sw, RESHAPE(soc_flux_up, (/144,3,40/)), Time_diag)
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

!     used = send_data ( id_soc_olr_spectrum_lw, RESHAPE(radout%flux_up_band(:,0,:), (/144,3,20/)), Time_diag)
!used = send_data ( id_soc_heating_sw, output_heating_rate_sw, Time_diag)
endif

end do


!--------------

end do

!output_heating_rate(:,:,1) = 0.0
!output_heating_rate(:,:,:10) = 0.2*output_heating_rate(:,:,:10)

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
