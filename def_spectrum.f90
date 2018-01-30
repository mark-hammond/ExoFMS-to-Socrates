! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module to declare a structure of spectral data.
!
! Description:
!   This module contains the heirarchical declaration of structures
!   of spectral data.
!
!------------------------------------------------------------------------------
! CAUTION - Any changes made to this routine need to be mirrored in the
!           setup_spectra_mod module in the UM.
!------------------------------------------------------------------------------
MODULE def_spectrum

USE realtype_rd

IMPLICIT NONE


INTEGER, PARAMETER :: n_dim = 22
!   Number of dimensions in StrSpecDim

TYPE StrSpecDim
  INTEGER :: nd_type
!   Size allocated for spectral blocks
  INTEGER :: nd_band
!   Size allocated for spectral bands
  INTEGER :: nd_exclude
!   Size allocated for excluded bands
  INTEGER :: nd_k_term
!   Size allocated for k-terms
  INTEGER :: nd_species
!   Size allocated for gaseous species
  INTEGER :: nd_scale_variable
!   Size allocated for scaling variables
  INTEGER :: nd_continuum
!   Size allocated for continua
  INTEGER :: nd_drop_type
!   Size allocated for drop types
  INTEGER :: nd_ice_type
!   Size allocated for ice crystal types
  INTEGER :: nd_aerosol_species
!   Size allocated for aerosol species
  INTEGER :: nd_aerosol_mr
!   Size allocated for aerosol mixing ratios
  INTEGER :: nd_thermal_coeff
!   Size allocated for thermal coefficients
  INTEGER :: nd_cloud_parameter
!   Size allocated for cloud parameters
  INTEGER :: nd_humidity
!   Size allocated for humidities
  INTEGER :: nd_aod_wavel
!   Number of wavelengths for aerosol optical depths
  INTEGER :: nd_phase_term
!   Size allocated for terms in the phase function
  INTEGER :: nd_tmp
!   Number of reference temperature for k-terms
  INTEGER :: nd_pre
!   Number of reference pressures for k-terms
  INTEGER :: nd_mix
!   Number of eta for mixture absorbing species
  INTEGER :: nd_band_mix_gas
!   Number of bands where mixed species exist
  INTEGER :: nd_sub_band
!   Size allocated for spectral sub-bands (for spectral variability)
  INTEGER :: nd_times
!   Size allocated for times (for spectral variability)
END TYPE StrSPecDim


TYPE StrSpecBasic
  LOGICAL, ALLOCATABLE      :: l_present(:)
!   Blocks of spectral data in the file
  INTEGER                   :: n_band
!   Number of spectral bands used
  REAL (RealK), ALLOCATABLE :: wavelength_long(:)
!   Lower wavelength limits for the band
  REAL (RealK), ALLOCATABLE :: wavelength_short(:)
!   Higher wavelengths limits for the band
  INTEGER, ALLOCATABLE      :: n_band_exclude(:)
!   Number of exclusions from each band
  INTEGER, ALLOCATABLE      :: index_exclude(:, :)
!   List of excluded bands within each region
END TYPE StrSpecBasic


TYPE StrSpecSolar
  REAL (RealK), ALLOCATABLE :: solar_flux_band(:)
!   Fraction of the solar spectrum in each band
  REAL (RealK), ALLOCATABLE :: solar_flux_band_ses(:, :)
!   Fraction of the solar spectrum for each k-term
  REAL (RealK), ALLOCATABLE :: weight_blue(:)
!   Fraction of the surface flux designated as "blue" in each band
END TYPE StrSpecSolar


TYPE StrSpecRayleigh
  INTEGER :: i_rayleigh_scheme
!   Type of Rayleigh scattering
  REAL (RealK), ALLOCATABLE :: rayleigh_coeff(:)
!   Rayleigh scattering coefficients in each band for total gas
  INTEGER :: n_gas_rayleigh
!   Total number of Rayleigh scattering gases
  INTEGER, ALLOCATABLE      :: index_rayleigh(:)
!   Index of gases for which Rayleigh scattering coefficients are tabulated
  REAL (RealK), ALLOCATABLE :: rayleigh_coeff_gas(:,:)
!   Rayleigh scattering coefficients for each gas in each band
END TYPE StrSpecRayleigh


TYPE StrSpecGas
  INTEGER  :: n_absorb
!   Total number of gaseous absorbers
  INTEGER, ALLOCATABLE      :: n_band_absorb(:)
!   Number of gaseous absorbers in each band
  INTEGER, ALLOCATABLE      :: index_absorb(:, :)
!   Number of gaseous absorbers
  INTEGER, ALLOCATABLE      :: type_absorb(:)
!   Actual types of each gas in the spectral file
  INTEGER, ALLOCATABLE      :: n_mix_gas(:)
!   Number of mixed gases in a band
  INTEGER, ALLOCATABLE      :: index_mix_gas(:, :)
!   Index of mixed absorbers in each band
  INTEGER, ALLOCATABLE      :: num_mix(:)
!   Number of binary parameter for interpolation of absorption
!   coefficient for mixture of two species
  INTEGER, ALLOCATABLE      :: mix_gas_band(:)
!   Sequence band number (not real band number) of mixed species
  INTEGER, ALLOCATABLE      :: num_ref_p(:, :)
!   Number of reference pressures
  INTEGER, ALLOCATABLE      :: num_ref_t(:, :)
!   Number of reference temperatures
  INTEGER, ALLOCATABLE      :: i_band_k(:, :)
!   Number of k-terms in each band for each gas
  INTEGER, ALLOCATABLE      :: i_band_k_ses(:)
!   Number of k-terms in band for each gas
  INTEGER, ALLOCATABLE      :: i_scale_k(:, :)
!   Type of scaling applied to each k-term
  INTEGER, ALLOCATABLE      :: i_scale_fnc(:, :)
!   Type of scaling function
  INTEGER, ALLOCATABLE      :: i_scat(:, :, :)
!   Method of scattering treatment for each k-term

  REAL (RealK), ALLOCATABLE :: k(:, :, :)
!   Absorption coefficients of k-terms
  REAL (RealK), ALLOCATABLE :: w(:, :, :)
!   Weights for k-terms
  REAL (RealK), ALLOCATABLE :: scale(:, :, :, :)
!   Scaling parameters for each absorber and term
  REAL (RealK), ALLOCATABLE :: p_ref(:, :)
!   Reference pressures for scaling functions
  REAL (RealK), ALLOCATABLE :: t_ref(:, :)
!   Reference temperatures for scaling functions

  REAL (RealK), ALLOCATABLE :: p_lookup(:)
  REAL (RealK), ALLOCATABLE :: t_lookup(:, :)
  REAL (RealK), ALLOCATABLE :: k_lookup(:, :, :, :, :)
  REAL (RealK), ALLOCATABLE :: w_ses(:, :)
  REAL (RealK), ALLOCATABLE :: k_mix_gas(:, :, :, :, :)
!   Absorption coefficients for mixture species
  REAL (RealK), ALLOCATABLE :: f_mix(:)
!   Mixing ratio of mixed absorber amount

  LOGICAL, ALLOCATABLE      :: l_doppler(:)
!   Flag for Doppler broadening for each species
  REAL (RealK), ALLOCATABLE :: doppler_cor(:)
!   Doppler correction terms
END TYPE StrSpecGas


TYPE StrSpecPlanck
  INTEGER                   :: n_deg_fit
!   Degree of the fit to the Planckian function
  REAL (RealK), ALLOCATABLE :: thermal_coeff(:, :)
!   Coefficients in polynomial fit to source function
  REAL (RealK), ALLOCATABLE :: theta_planck_tbl(:)
!   Temperatures at which the band-integrated Planck function
!   has been evaluated.
  REAL (RealK)              :: t_ref_planck
!   Reference temperature for the Plackian function
  LOGICAL                   :: l_planck_tbl
!   Flag for using a look-up table instead of a polynomial
END TYPE StrSpecPlanck


TYPE StrSpecCont
  INTEGER, ALLOCATABLE      :: n_band_continuum(:)
!   Number of continua in each band
  INTEGER, ALLOCATABLE      :: index_continuum(:, :)
!   List of continua in each band
  INTEGER                   :: index_water
!   Index of water vapour of continua in each band
  INTEGER, ALLOCATABLE      :: i_scale_fnc_cont(:, :)
!   Types of scaling functions for continua

  REAL (RealK), ALLOCATABLE :: k_cont(:, :)
!   ABsorption coefficients for continuum absorption
  REAL (RealK), ALLOCATABLE :: scale_cont(:, :, :)
!   Reference temperature for the Plackian function
  REAL (RealK), ALLOCATABLE :: p_ref_cont(:, :)
!   Reference pressures for continuum scaling functions
  REAL (RealK), ALLOCATABLE :: t_ref_cont(:, :)
!   Reference temperatures for continuum scaling functions
  REAL (RealK), ALLOCATABLE :: k_cont_ses(:, :, :, :)
  REAL (RealK), ALLOCATABLE :: k_h2oc(:, :, :, :)
!   Absorption coefficient for water vapour continuum
END TYPE StrSpecCont


TYPE StrSpecDrop
  LOGICAL, ALLOCATABLE      :: l_drop_type(:)
!   Flags for types of droplets present
  INTEGER, ALLOCATABLE      :: i_drop_parm(:)
!   Form of parametrization for each type of droplet
  INTEGER, ALLOCATABLE      :: n_phf(:)
!   Number of moments of the phase fuction fitted (N. B. This
!   array is not set for parametrizations which are implicitly
!   restricted to the asymmetry.)

  REAL (RealK), ALLOCATABLE :: parm_list(:, :, :)
!   Parameters used to fit the optical properties of droplets
  REAL (RealK), ALLOCATABLE :: parm_min_dim(:)
!   Minimum dimension permissible in the parametrization
  REAL (RealK), ALLOCATABLE :: parm_max_dim(:)
!   Maximum dimension permissible in the parametrization
END TYPE StrSpecDrop


TYPE StrSpecAerosol
  LOGICAL, ALLOCATABLE      :: l_aero_spec(:)
!   Flags for species of aerosol present

  INTEGER                   :: n_aerosol
!   Number of aerosol species present in spectral file
  INTEGER                   :: n_aerosol_mr
!   Number of aerosol species present in mixing ratio array
  INTEGER, ALLOCATABLE      :: type_aerosol(:)
!   Actual types of aerosols in the spectral file
  INTEGER, ALLOCATABLE      :: i_aerosol_parm(:)
!   Parametrization scheme used for each aerosol
  INTEGER, ALLOCATABLE      :: n_aerosol_phf_term(:)
!   Number of terms in the phase function
  INTEGER, ALLOCATABLE      :: nhumidity(:)
!   Number of values of humidity

  REAL (RealK), ALLOCATABLE :: abs(:, :, :)
!   Absortption by aerosols
  REAL (RealK), ALLOCATABLE :: scat(:, :, :)
!   Scattering by aerosols
  REAL (RealK), ALLOCATABLE :: phf_fnc(:, :, :, :)
!   Phase functions of aerosols
  REAL (RealK), ALLOCATABLE :: humidities(:, :)
!   Humdities of each component

! Fields for aerosol optical depth:
  INTEGER                   :: n_aod_wavel
!   Number of wavelengths
  INTEGER, ALLOCATABLE      :: i_aod_type(:)
!   Relationship between aerosol component and type
  REAL (RealK), ALLOCATABLE :: aod_wavel(:)
!   Wavelengths for the aod
  REAL (RealK), ALLOCATABLE :: aod_abs(:, :, :)
!   Monochromatic specific absorption coefficient
  REAL (RealK), ALLOCATABLE :: aod_scat(:, :, :)
!   Monochromatic specific scattering coefficient
END TYPE StrSpecAerosol


TYPE StrSpecIce
  LOGICAL, ALLOCATABLE      :: l_ice_type(:)
!   Flags for types of ice crystals present
  INTEGER, ALLOCATABLE      :: i_ice_parm(:)
!   Form of parametrization for each type of ice crystal
  INTEGER, ALLOCATABLE      :: n_phf(:)
!   Number of moments of the phase fuction fitted

  REAL (RealK), ALLOCATABLE :: parm_list(:, :, :)
!   Parameters used to fit the optical properties of ice crystals
  REAL (RealK), ALLOCATABLE :: parm_min_dim(:)
!   Minimum dimension permissible in the parametrization
  REAL (RealK), ALLOCATABLE :: parm_max_dim(:)
!   Maximum dimension permissible in the parametrization
END TYPE StrSpecIce


TYPE StrSpecVar
  INTEGER                   :: n_sub_band
!   Number of sub-bands used
  INTEGER                   :: n_times
!   Number of times at which the solar spectrum is given
  INTEGER                   :: n_repeat_times
!   Number of times over which to periodically repeat data into the future
  INTEGER                   :: n_rayleigh_coeff
!   Number of Rayleigh coefficients that vary
  INTEGER, ALLOCATABLE      :: index_sub_band(:, :)
!   Index of k-terms associated with each sub-band
  REAL (RealK), ALLOCATABLE :: wavelength_sub_band(:, :)
!   Wavelength limits for the sub-band

  INTEGER, ALLOCATABLE      :: time(:, :)
!   Times: year, month, day of month, seconds in day
  REAL (RealK), ALLOCATABLE :: total_solar_flux(:)
!   Total solar flux in Wm-2 at 1 AU for each time
  REAL (RealK), ALLOCATABLE :: solar_flux_sub_band(:, :)
!   Fraction of the solar spectrum in each sub-band for each time
  REAL (RealK), ALLOCATABLE :: rayleigh_coeff(:, :)
!   Rayleigh scattering coefficients in each sub-band for each time
END TYPE StrSpecVar


TYPE StrSpecData
  TYPE (StrSpecDim)               :: Dim
  TYPE (StrSpecBasic)             :: Basic
  TYPE (StrSpecSolar)             :: Solar
  TYPE (StrSpecRayleigh)          :: Rayleigh
  TYPE (StrSpecGas)               :: Gas
  TYPE (StrSpecPlanck)            :: Planck
  TYPE (StrSpecCont)              :: Cont
  TYPE (StrSpecDrop)              :: Drop
  TYPE (StrSpecAerosol)           :: Aerosol
  TYPE (StrSpecIce)               :: Ice
  TYPE (StrSpecVar)               :: Var
END TYPE StrSpecData


CONTAINS
!------------------------------------------------------------------------------
SUBROUTINE allocate_spectrum(Sp)

USE missing_data_mod, ONLY: rmdi

IMPLICIT NONE

TYPE (StrSpecData), INTENT(INOUT) :: Sp


! Basic
IF (.NOT. ALLOCATED(Sp%Basic%l_present)) THEN
  ALLOCATE(Sp%Basic%l_present(0:Sp%Dim%nd_type))
  Sp%Basic%l_present = .FALSE.
END IF

IF (.NOT. ALLOCATED(Sp%Basic%wavelength_long)) &
  ALLOCATE(Sp%Basic%wavelength_long( Sp%Dim%nd_band ))

IF (.NOT. ALLOCATED(Sp%Basic%wavelength_short)) &
  ALLOCATE(Sp%Basic%wavelength_short( Sp%Dim%nd_band ))

IF (.NOT. ALLOCATED(Sp%Basic%n_band_exclude)) &
  ALLOCATE(Sp%Basic%n_band_exclude( Sp%Dim%nd_band ))

IF (.NOT. ALLOCATED(Sp%Basic%index_exclude)) &
  ALLOCATE(Sp%Basic%index_exclude( Sp%Dim%nd_exclude, Sp%Dim%nd_band ))

! Solar
IF (.NOT. ALLOCATED(Sp%Solar%solar_flux_band)) &
  ALLOCATE(Sp%Solar%solar_flux_band( Sp%Dim%nd_band ))

IF (.NOT. ALLOCATED(Sp%Solar%solar_flux_band_ses)) &
  ALLOCATE(Sp%Solar%solar_flux_band_ses( Sp%Dim%nd_k_term, Sp%Dim%nd_band ))

IF (.NOT. ALLOCATED(Sp%Solar%weight_blue)) THEN
  ALLOCATE(Sp%Solar%weight_blue( Sp%Dim%nd_band ))
  Sp%Solar%weight_blue = rmdi
END IF

! Rayleigh
IF (.NOT. ALLOCATED(Sp%Rayleigh%rayleigh_coeff)) &
  ALLOCATE(Sp%Rayleigh%rayleigh_coeff( Sp%Dim%nd_band ))
IF (.NOT. ALLOCATED(Sp%Rayleigh%index_rayleigh)) &
  ALLOCATE(Sp%Rayleigh%index_rayleigh( Sp%Dim%nd_species ))
IF (.NOT. ALLOCATED(Sp%Rayleigh%rayleigh_coeff_gas)) &
  ALLOCATE(Sp%Rayleigh%rayleigh_coeff_gas( Sp%Dim%nd_species, Sp%Dim%nd_band ))

! Gas
IF (.NOT. ALLOCATED(Sp%Gas%n_band_absorb)) &
  ALLOCATE(Sp%Gas%n_band_absorb( Sp%Dim%nd_band ))

IF (.NOT. ALLOCATED(Sp%Gas%index_absorb)) &
  ALLOCATE(Sp%Gas%index_absorb( Sp%Dim%nd_species, Sp%Dim%nd_band ))

IF (.NOT. ALLOCATED(Sp%Gas%type_absorb)) &
  ALLOCATE(Sp%Gas%type_absorb( Sp%Dim%nd_species ))

IF (.NOT. ALLOCATED(Sp%Gas%n_mix_gas)) &
  ALLOCATE(Sp%Gas%n_mix_gas( Sp%Dim%nd_band ))

IF (.NOT. ALLOCATED(Sp%Gas%index_mix_gas)) &
  ALLOCATE(Sp%Gas%index_mix_gas( 2, Sp%Dim%nd_band_mix_gas ))

IF (.NOT. ALLOCATED(Sp%Gas%num_mix)) &
  ALLOCATE(Sp%Gas%num_mix( Sp%Dim%nd_band ))

IF (.NOT. ALLOCATED(Sp%Gas%mix_gas_band)) &
  ALLOCATE(Sp%Gas%mix_gas_band( Sp%Dim%nd_band ))

IF (.NOT. ALLOCATED(Sp%Gas%num_ref_p)) &
  ALLOCATE(Sp%Gas%num_ref_p( Sp%Dim%nd_species, Sp%Dim%nd_band ))

IF (.NOT. ALLOCATED(Sp%Gas%num_ref_t)) &
  ALLOCATE(Sp%Gas%num_ref_t( Sp%Dim%nd_species, Sp%Dim%nd_band ))

IF (.NOT. ALLOCATED(Sp%Gas%i_band_k)) THEN
  ALLOCATE(Sp%Gas%i_band_k( Sp%Dim%nd_band, Sp%Dim%nd_species ))
  Sp%Gas%i_band_k=0
END IF

IF (.NOT. ALLOCATED(Sp%Gas%i_band_k_ses)) THEN
  ALLOCATE(Sp%Gas%i_band_k_ses( Sp%Dim%nd_band ))
  Sp%Gas%i_band_k_ses=0
END IF

IF (.NOT. ALLOCATED(Sp%Gas%i_scale_k)) &
  ALLOCATE(Sp%Gas%i_scale_k( Sp%Dim%nd_band, Sp%Dim%nd_species ))

IF (.NOT. ALLOCATED(Sp%Gas%i_scale_fnc)) &
  ALLOCATE(Sp%Gas%i_scale_fnc( Sp%Dim%nd_band, Sp%Dim%nd_species ))

IF (.NOT. ALLOCATED(Sp%Gas%i_scat)) &
  ALLOCATE(Sp%Gas%i_scat( Sp%Dim%nd_k_term, Sp%Dim%nd_band, Sp%Dim%nd_species ))

IF (.NOT. ALLOCATED(Sp%Gas%k)) &
  ALLOCATE(Sp%Gas%k( Sp%Dim%nd_k_term, Sp%Dim%nd_band, Sp%Dim%nd_species ))

IF (.NOT. ALLOCATED(Sp%Gas%w)) &
  ALLOCATE(Sp%Gas%w( Sp%Dim%nd_k_term, Sp%Dim%nd_band, Sp%Dim%nd_species ))

IF (.NOT. ALLOCATED(Sp%Gas%scale)) &
  ALLOCATE(Sp%Gas%scale( Sp%Dim%nd_scale_variable, Sp%Dim%nd_k_term, &
                         Sp%Dim%nd_band, Sp%Dim%nd_species ))

IF (.NOT. ALLOCATED(Sp%Gas%p_ref)) &
  ALLOCATE(Sp%Gas%p_ref( Sp%Dim%nd_species, Sp%Dim%nd_band ))

IF (.NOT. ALLOCATED(Sp%Gas%t_ref)) &
  ALLOCATE(Sp%Gas%t_ref( Sp%Dim%nd_species, Sp%Dim%nd_band ))

IF (.NOT. ALLOCATED(Sp%Gas%p_lookup)) &
  ALLOCATE(Sp%Gas%p_lookup( Sp%Dim%nd_pre ))

IF (.NOT. ALLOCATED(Sp%Gas%t_lookup)) &
  ALLOCATE(Sp%Gas%t_lookup( Sp%Dim%nd_tmp, Sp%Dim%nd_pre ))

IF (.NOT. ALLOCATED(Sp%Gas%k_lookup)) &
  ALLOCATE(Sp%Gas%k_lookup( Sp%Dim%nd_tmp, Sp%Dim%nd_pre, Sp%Dim%nd_k_term, &
                            Sp%Dim%nd_species, Sp%Dim%nd_band ))

IF (.NOT. ALLOCATED(Sp%Gas%w_ses)) &
  ALLOCATE(Sp%Gas%w_ses( Sp%Dim%nd_k_term, Sp%Dim%nd_band ))

IF (.NOT. ALLOCATED(Sp%Gas%k_mix_gas)) &
  ALLOCATE(Sp%Gas%k_mix_gas( Sp%Dim%nd_pre, Sp%Dim%nd_tmp, Sp%Dim%nd_mix, &
                             Sp%Dim%nd_k_term, Sp%Dim%nd_band_mix_gas ))

IF (.NOT. ALLOCATED(Sp%Gas%f_mix)) &
  ALLOCATE(Sp%Gas%f_mix( Sp%Dim%nd_band ))

IF (.NOT. ALLOCATED(Sp%Gas%l_doppler)) THEN
  ALLOCATE(Sp%Gas%l_doppler( Sp%Dim%nd_species ))
  Sp%Gas%l_doppler = .FALSE.
END IF

IF (.NOT. ALLOCATED(Sp%Gas%doppler_cor)) &
  ALLOCATE(Sp%Gas%doppler_cor( Sp%Dim%nd_species ))

! Planck
IF (.NOT. ALLOCATED(Sp%Planck%thermal_coeff)) &
  ALLOCATE(Sp%Planck%thermal_coeff( 0:Sp%Dim%nd_thermal_coeff-1, &
                                    Sp%Dim%nd_band ))
IF (.NOT. ALLOCATED(Sp%Planck%theta_planck_tbl)) &
  ALLOCATE(Sp%Planck%theta_planck_tbl( 0:Sp%Dim%nd_thermal_coeff-1 ))

! Cont
IF (.NOT. ALLOCATED(Sp%Cont%n_band_continuum)) &
  ALLOCATE(Sp%Cont%n_band_continuum( Sp%Dim%nd_band ))

IF (.NOT. ALLOCATED(Sp%Cont%index_continuum)) &
  ALLOCATE(Sp%Cont%index_continuum( Sp%Dim%nd_band, Sp%Dim%nd_continuum ))

IF (.NOT. ALLOCATED(Sp%Cont%i_scale_fnc_cont)) &
  ALLOCATE(Sp%Cont%i_scale_fnc_cont( Sp%Dim%nd_band, Sp%Dim%nd_continuum ))

IF (.NOT. ALLOCATED(Sp%Cont%k_cont)) &
  ALLOCATE(Sp%Cont%k_cont( Sp%Dim%nd_band, Sp%Dim%nd_continuum ))

IF (.NOT. ALLOCATED(Sp%Cont%scale_cont)) &
  ALLOCATE(Sp%Cont%scale_cont( Sp%Dim%nd_scale_variable, &
                               Sp%Dim%nd_band, Sp%Dim%nd_continuum ))

IF (.NOT. ALLOCATED(Sp%Cont%p_ref_cont)) &
  ALLOCATE(Sp%Cont%p_ref_cont( Sp%Dim%nd_continuum, Sp%Dim%nd_band ))

IF (.NOT. ALLOCATED(Sp%Cont%t_ref_cont)) &
  ALLOCATE(Sp%Cont%t_ref_cont( Sp%Dim%nd_continuum, Sp%Dim%nd_band ))

IF (.NOT. ALLOCATED(Sp%Cont%k_cont_ses)) &
  ALLOCATE(Sp%Cont%k_cont_ses( Sp%Dim%nd_k_term, Sp%Dim%nd_tmp, &
                               Sp%Dim%nd_band, Sp%Dim%nd_continuum ))

IF (.NOT. ALLOCATED(Sp%Cont%k_h2oc)) &
  ALLOCATE(Sp%Cont%k_h2oc( Sp%Dim%nd_pre, Sp%Dim%nd_tmp, &
                           Sp%Dim%nd_k_term, Sp%Dim%nd_band ))

! Drop
IF (.NOT. ALLOCATED(Sp%Drop%l_drop_type)) THEN
  ALLOCATE(Sp%Drop%l_drop_type( Sp%Dim%nd_drop_type ))
  Sp%Drop%l_drop_type = .FALSE.
END IF

IF (.NOT. ALLOCATED(Sp%Drop%i_drop_parm)) &
  ALLOCATE(Sp%Drop%i_drop_parm( Sp%Dim%nd_drop_type ))

IF (.NOT. ALLOCATED(Sp%Drop%n_phf)) &
  ALLOCATE(Sp%Drop%n_phf( Sp%Dim%nd_drop_type ))

IF (.NOT. ALLOCATED(Sp%Drop%parm_list)) &
  ALLOCATE(Sp%Drop%parm_list( Sp%Dim%nd_cloud_parameter, Sp%Dim%nd_band, &
                              Sp%Dim%nd_drop_type ))

IF (.NOT. ALLOCATED(Sp%Drop%parm_min_dim)) &
  ALLOCATE(Sp%Drop%parm_min_dim( Sp%Dim%nd_drop_type ))

IF (.NOT. ALLOCATED(Sp%Drop%parm_max_dim)) &
  ALLOCATE(Sp%Drop%parm_max_dim( Sp%Dim%nd_drop_type ))

! Aerosol
IF (.NOT. ALLOCATED(Sp%Aerosol%l_aero_spec)) THEN
  ALLOCATE(Sp%Aerosol%l_aero_spec( Sp%Dim%nd_aerosol_species ))
  Sp%Aerosol%l_aero_spec = .FALSE.
END IF

IF (.NOT. ALLOCATED(Sp%Aerosol%type_aerosol)) &
  ALLOCATE(Sp%Aerosol%type_aerosol( Sp%Dim%nd_aerosol_species ))

IF (.NOT. ALLOCATED(Sp%Aerosol%i_aerosol_parm)) &
  ALLOCATE(Sp%Aerosol%i_aerosol_parm( Sp%Dim%nd_aerosol_species ))

IF (.NOT. ALLOCATED(Sp%Aerosol%n_aerosol_phf_term)) &
  ALLOCATE(Sp%Aerosol%n_aerosol_phf_term( Sp%Dim%nd_aerosol_species ))

IF (.NOT. ALLOCATED(Sp%Aerosol%nhumidity)) &
  ALLOCATE(Sp%Aerosol%nhumidity( Sp%Dim%nd_aerosol_species ))

IF (.NOT. ALLOCATED(Sp%Aerosol%abs)) &
  ALLOCATE(Sp%Aerosol%abs( Sp%Dim%nd_humidity, Sp%Dim%nd_aerosol_species, &
                           Sp%Dim%nd_band ))

IF (.NOT. ALLOCATED(Sp%Aerosol%scat)) &
  ALLOCATE(Sp%Aerosol%scat( Sp%Dim%nd_humidity, Sp%Dim%nd_aerosol_species, &
                            Sp%Dim%nd_band ))

IF (.NOT. ALLOCATED(Sp%Aerosol%phf_fnc)) &
  ALLOCATE(Sp%Aerosol%phf_fnc( Sp%Dim%nd_humidity, Sp%Dim%nd_phase_term, &
                               Sp%Dim%nd_aerosol_species, Sp%Dim%nd_band ))

IF (.NOT. ALLOCATED(Sp%Aerosol%humidities)) &
  ALLOCATE(Sp%Aerosol%humidities( Sp%Dim%nd_humidity, &
                                  Sp%Dim%nd_aerosol_species ))

IF (.NOT. ALLOCATED(Sp%Aerosol%i_aod_type)) &
  ALLOCATE(Sp%Aerosol%i_aod_type( Sp%Dim%nd_aerosol_species ))

IF (.NOT. ALLOCATED(Sp%Aerosol%aod_wavel)) &
  ALLOCATE(Sp%Aerosol%aod_wavel( Sp%Dim%nd_aod_wavel ))

IF (.NOT. ALLOCATED(Sp%Aerosol%aod_abs)) &
  ALLOCATE(Sp%Aerosol%aod_abs( Sp%Dim%nd_humidity, Sp%Dim%nd_aerosol_species, &
                               Sp%Dim%nd_aod_wavel ))

IF (.NOT. ALLOCATED(Sp%Aerosol%aod_scat)) &
  ALLOCATE(Sp%Aerosol%aod_scat( Sp%Dim%nd_humidity, Sp%Dim%nd_aerosol_species, &
                                Sp%Dim%nd_aod_wavel ))

! Ice
IF (.NOT. ALLOCATED(Sp%Ice%l_ice_type)) THEN
  ALLOCATE(Sp%Ice%l_ice_type( Sp%Dim%nd_ice_type ))
  Sp%Ice%l_ice_type = .FALSE.
END IF

IF (.NOT. ALLOCATED(Sp%Ice%i_ice_parm)) &
  ALLOCATE(Sp%Ice%i_ice_parm( Sp%Dim%nd_ice_type ))

IF (.NOT. ALLOCATED(Sp%Ice%n_phf)) &
  ALLOCATE(Sp%Ice%n_phf( Sp%Dim%nd_ice_type ))

IF (.NOT. ALLOCATED(Sp%Ice%parm_list)) &
  ALLOCATE(Sp%Ice%parm_list( Sp%Dim%nd_cloud_parameter, Sp%Dim%nd_band, &
                             Sp%Dim%nd_ice_type ))

IF (.NOT. ALLOCATED(Sp%Ice%parm_min_dim)) &
  ALLOCATE(Sp%Ice%parm_min_dim( Sp%Dim%nd_ice_type ))

IF (.NOT. ALLOCATED(Sp%Ice%parm_max_dim)) &
  ALLOCATE(Sp%Ice%parm_max_dim( Sp%Dim%nd_ice_type ))

! Spectral variability
IF (.NOT. ALLOCATED(Sp%Var%index_sub_band)) &
  ALLOCATE(Sp%Var%index_sub_band( 2, Sp%Dim%nd_sub_band ))

IF (.NOT. ALLOCATED(Sp%Var%wavelength_sub_band)) &
  ALLOCATE(Sp%Var%wavelength_sub_band( 2, Sp%Dim%nd_sub_band ))

IF (.NOT. ALLOCATED(Sp%Var%time)) &
  ALLOCATE(Sp%Var%time( 4, Sp%Dim%nd_times ))

IF (.NOT. ALLOCATED(Sp%Var%total_solar_flux)) &
  ALLOCATE(Sp%Var%total_solar_flux( Sp%Dim%nd_times ))

IF (.NOT. ALLOCATED(Sp%Var%solar_flux_sub_band)) &
  ALLOCATE(Sp%Var%solar_flux_sub_band( Sp%Dim%nd_sub_band, Sp%Dim%nd_times ))

IF (.NOT. ALLOCATED(Sp%Var%rayleigh_coeff)) &
  ALLOCATE(Sp%Var%rayleigh_coeff( Sp%Dim%nd_sub_band, 0:Sp%Dim%nd_times ))

END SUBROUTINE allocate_spectrum
!------------------------------------------------------------------------------
SUBROUTINE deallocate_spectrum(Sp)

IMPLICIT NONE

TYPE (StrSpecData), INTENT(INOUT) :: Sp

! Spectral variability
IF (ALLOCATED(Sp%Var%rayleigh_coeff)) &
   DEALLOCATE(Sp%Var%rayleigh_coeff)
IF (ALLOCATED(Sp%Var%solar_flux_sub_band)) &
   DEALLOCATE(Sp%Var%solar_flux_sub_band)
IF (ALLOCATED(Sp%Var%total_solar_flux)) &
   DEALLOCATE(Sp%Var%total_solar_flux)
IF (ALLOCATED(Sp%Var%time)) &
   DEALLOCATE(Sp%Var%time)
IF (ALLOCATED(Sp%Var%wavelength_sub_band)) &
   DEALLOCATE(Sp%Var%wavelength_sub_band)
IF (ALLOCATED(Sp%Var%index_sub_band)) &
   DEALLOCATE(Sp%Var%index_sub_band)

! Ice
IF (ALLOCATED(Sp%Ice%parm_max_dim)) &
   DEALLOCATE(Sp%Ice%parm_max_dim)
IF (ALLOCATED(Sp%Ice%parm_min_dim)) &
   DEALLOCATE(Sp%Ice%parm_min_dim)
IF (ALLOCATED(Sp%Ice%parm_list)) &
   DEALLOCATE(Sp%Ice%parm_list)
IF (ALLOCATED(Sp%Ice%n_phf)) &
   DEALLOCATE(Sp%Ice%n_phf)
IF (ALLOCATED(Sp%Ice%i_ice_parm)) &
   DEALLOCATE(Sp%Ice%i_ice_parm)
IF (ALLOCATED(Sp%Ice%l_ice_type)) &
   DEALLOCATE(Sp%Ice%l_ice_type)

! Aerosol
IF (ALLOCATED(Sp%Aerosol%aod_scat)) &
   DEALLOCATE(Sp%Aerosol%aod_scat)
IF (ALLOCATED(Sp%Aerosol%aod_abs)) &
   DEALLOCATE(Sp%Aerosol%aod_abs)
IF (ALLOCATED(Sp%Aerosol%aod_wavel)) &
   DEALLOCATE(Sp%Aerosol%aod_wavel)
IF (ALLOCATED(Sp%Aerosol%i_aod_type)) &
   DEALLOCATE(Sp%Aerosol%i_aod_type)
IF (ALLOCATED(Sp%Aerosol%humidities)) &
   DEALLOCATE(Sp%Aerosol%humidities)
IF (ALLOCATED(Sp%Aerosol%phf_fnc)) &
   DEALLOCATE(Sp%Aerosol%phf_fnc)
IF (ALLOCATED(Sp%Aerosol%scat)) &
   DEALLOCATE(Sp%Aerosol%scat)
IF (ALLOCATED(Sp%Aerosol%abs)) &
   DEALLOCATE(Sp%Aerosol%abs)
IF (ALLOCATED(Sp%Aerosol%nhumidity)) &
   DEALLOCATE(Sp%Aerosol%nhumidity)
IF (ALLOCATED(Sp%Aerosol%n_aerosol_phf_term)) &
   DEALLOCATE(Sp%Aerosol%n_aerosol_phf_term)
IF (ALLOCATED(Sp%Aerosol%i_aerosol_parm)) &
   DEALLOCATE(Sp%Aerosol%i_aerosol_parm)
IF (ALLOCATED(Sp%Aerosol%type_aerosol)) &
   DEALLOCATE(Sp%Aerosol%type_aerosol)
IF (ALLOCATED(Sp%Aerosol%l_aero_spec)) &
   DEALLOCATE(Sp%Aerosol%l_aero_spec)

! Drop
IF (ALLOCATED(Sp%Drop%parm_max_dim)) &
   DEALLOCATE(Sp%Drop%parm_max_dim)
IF (ALLOCATED(Sp%Drop%parm_min_dim)) &
   DEALLOCATE(Sp%Drop%parm_min_dim)
IF (ALLOCATED(Sp%Drop%parm_list)) &
   DEALLOCATE(Sp%Drop%parm_list)
IF (ALLOCATED(Sp%Drop%n_phf)) &
   DEALLOCATE(Sp%Drop%n_phf)
IF (ALLOCATED(Sp%Drop%i_drop_parm)) &
   DEALLOCATE(Sp%Drop%i_drop_parm)
IF (ALLOCATED(Sp%Drop%l_drop_type)) &
   DEALLOCATE(Sp%Drop%l_drop_type)

! Cont
IF (ALLOCATED(Sp%Cont%k_h2oc)) &
   DEALLOCATE(Sp%Cont%k_h2oc)
IF (ALLOCATED(Sp%Cont%k_cont_ses)) &
   DEALLOCATE(Sp%Cont%k_cont_ses)
IF (ALLOCATED(Sp%Cont%t_ref_cont)) &
   DEALLOCATE(Sp%Cont%t_ref_cont)
IF (ALLOCATED(Sp%Cont%p_ref_cont)) &
   DEALLOCATE(Sp%Cont%p_ref_cont)
IF (ALLOCATED(Sp%Cont%scale_cont)) &
   DEALLOCATE(Sp%Cont%scale_cont)
IF (ALLOCATED(Sp%Cont%k_cont)) &
   DEALLOCATE(Sp%Cont%k_cont)
IF (ALLOCATED(Sp%Cont%i_scale_fnc_cont)) &
   DEALLOCATE(Sp%Cont%i_scale_fnc_cont)
IF (ALLOCATED(Sp%Cont%index_continuum)) &
   DEALLOCATE(Sp%Cont%index_continuum)
IF (ALLOCATED(Sp%Cont%n_band_continuum)) &
   DEALLOCATE(Sp%Cont%n_band_continuum)

! Planck
IF (ALLOCATED(Sp%Planck%theta_planck_tbl)) &
   DEALLOCATE(Sp%Planck%theta_planck_tbl)
IF (ALLOCATED(Sp%Planck%thermal_coeff)) &
   DEALLOCATE(Sp%Planck%thermal_coeff)

! Gas
IF (ALLOCATED(Sp%Gas%doppler_cor)) &
   DEALLOCATE(Sp%Gas%doppler_cor)
IF (ALLOCATED(Sp%Gas%l_doppler)) &
   DEALLOCATE(Sp%Gas%l_doppler)
IF (ALLOCATED(Sp%Gas%f_mix)) &
   DEALLOCATE(Sp%Gas%f_mix)
IF (ALLOCATED(Sp%Gas%k_mix_gas)) &
   DEALLOCATE(Sp%Gas%k_mix_gas)
IF (ALLOCATED(Sp%Gas%w_ses)) &
   DEALLOCATE(Sp%Gas%w_ses)
IF (ALLOCATED(Sp%Gas%k_lookup)) &
   DEALLOCATE(Sp%Gas%k_lookup)
IF (ALLOCATED(Sp%Gas%t_lookup)) &
   DEALLOCATE(Sp%Gas%t_lookup)
IF (ALLOCATED(Sp%Gas%p_lookup)) &
   DEALLOCATE(Sp%Gas%p_lookup)
IF (ALLOCATED(Sp%Gas%t_ref)) &
   DEALLOCATE(Sp%Gas%t_ref)
IF (ALLOCATED(Sp%Gas%p_ref)) &
   DEALLOCATE(Sp%Gas%p_ref)
IF (ALLOCATED(Sp%Gas%scale)) &
   DEALLOCATE(Sp%Gas%scale)
IF (ALLOCATED(Sp%Gas%w)) &
   DEALLOCATE(Sp%Gas%w)
IF (ALLOCATED(Sp%Gas%k)) &
   DEALLOCATE(Sp%Gas%k)
IF (ALLOCATED(Sp%Gas%i_scat)) &
   DEALLOCATE(Sp%Gas%i_scat)
IF (ALLOCATED(Sp%Gas%i_scale_fnc)) &
   DEALLOCATE(Sp%Gas%i_scale_fnc)
IF (ALLOCATED(Sp%Gas%i_scale_k)) &
   DEALLOCATE(Sp%Gas%i_scale_k)
IF (ALLOCATED(Sp%Gas%i_band_k_ses)) &
   DEALLOCATE(Sp%Gas%i_band_k_ses)
IF (ALLOCATED(Sp%Gas%i_band_k)) &
   DEALLOCATE(Sp%Gas%i_band_k)
IF (ALLOCATED(Sp%Gas%num_ref_t)) &
   DEALLOCATE(Sp%Gas%num_ref_t)
IF (ALLOCATED(Sp%Gas%num_ref_p)) &
   DEALLOCATE(Sp%Gas%num_ref_p)
IF (ALLOCATED(Sp%Gas%mix_gas_band)) &
   DEALLOCATE(Sp%Gas%mix_gas_band)
IF (ALLOCATED(Sp%Gas%num_mix)) &
   DEALLOCATE(Sp%Gas%num_mix)
IF (ALLOCATED(Sp%Gas%index_mix_gas)) &
   DEALLOCATE(Sp%Gas%index_mix_gas)
IF (ALLOCATED(Sp%Gas%n_mix_gas)) &
   DEALLOCATE(Sp%Gas%n_mix_gas)
IF (ALLOCATED(Sp%Gas%type_absorb)) &
   DEALLOCATE(Sp%Gas%type_absorb)
IF (ALLOCATED(Sp%Gas%index_absorb)) &
   DEALLOCATE(Sp%Gas%index_absorb)
IF (ALLOCATED(Sp%Gas%n_band_absorb)) &
   DEALLOCATE(Sp%Gas%n_band_absorb)

! Rayleigh
IF (ALLOCATED(Sp%Rayleigh%rayleigh_coeff)) &
   DEALLOCATE(Sp%Rayleigh%rayleigh_coeff)
IF (ALLOCATED(Sp%Rayleigh%index_rayleigh)) &
   DEALLOCATE(Sp%Rayleigh%index_rayleigh)
IF (ALLOCATED(Sp%Rayleigh%rayleigh_coeff_gas)) &
   DEALLOCATE(Sp%Rayleigh%rayleigh_coeff_gas)

! Solar
IF (ALLOCATED(Sp%Solar%weight_blue)) &
   DEALLOCATE(Sp%Solar%weight_blue)
IF (ALLOCATED(Sp%Solar%solar_flux_band_ses)) &
   DEALLOCATE(Sp%Solar%solar_flux_band_ses)
IF (ALLOCATED(Sp%Solar%solar_flux_band)) &
   DEALLOCATE(Sp%Solar%solar_flux_band)

! Basic
IF (ALLOCATED(Sp%Basic%index_exclude)) &
   DEALLOCATE(Sp%Basic%index_exclude)
IF (ALLOCATED(Sp%Basic%n_band_exclude)) &
   DEALLOCATE(Sp%Basic%n_band_exclude)
IF (ALLOCATED(Sp%Basic%wavelength_short)) &
   DEALLOCATE(Sp%Basic%wavelength_short)
IF (ALLOCATED(Sp%Basic%wavelength_long)) &
   DEALLOCATE(Sp%Basic%wavelength_long)
IF (ALLOCATED(Sp%Basic%l_present)) &
   DEALLOCATE(Sp%Basic%l_present)

END SUBROUTINE deallocate_spectrum
!------------------------------------------------------------------------------

END MODULE def_spectrum
