MODULE tau_output_mod

  ! Optical thickness output module
  ! MDH 15/01/18

  !----------

  REAL, SAVE :: spectral_tau(144,3,40)
  INTEGER, SAVE :: total_lat_n
  INTEGER, SAVE :: total_lon_n
  INTEGER, SAVE :: total_n

CONTAINS

  !-------------------------------------

  SUBROUTINE send_tau_output(received_tau, lon_number, lat_number)


    REAL, INTENT(in) :: received_tau(40)

    INTEGER, INTENT(in) :: lon_number
    INTEGER, INTENT(in) :: lat_number


    !----------

!    PRINT*, 'ook'
!    PRINT*, received_tau

    total_n = total_n + 1


IF (total_n < 433) THEN
    total_lat_n = MOD((total_n-1),3)+1
    !total_lon_n = CEILING(REAL(total_n),3) + total_lat_n
    total_lon_n = ((total_n - total_lat_n) / 3) + 1


    spectral_tau(total_lon_n, total_lat_n, :) =  received_tau(:)

  END IF

!    PRINT*, 'ook'
!    PRINT*, received_tau(:)


  END SUBROUTINE send_tau_output

  !-------------------------------------


  SUBROUTINE write_tau_output(processor_lat)

    INTEGER, INTENT(in) :: processor_lat
    INTEGER :: i, j, k

    ! FIle handle and name
    INTEGER, PARAMETER :: out_unit=20
    CHARACTER(len=200) :: file_name

PRINT*, 'ook!'
PRINT*, spectral_tau(1,1,:)

    ! Generate file name
    WRITE(file_name,'(a,i4.4,a)') "/network/group/aopp/testvol2/plan/fms-scratch-mdh/tau_output/TAU_OUT_50_", processor_lat,".TXT"

    ! Open file
    OPEN (unit=out_unit,file=TRIM(file_name),action="write",status="replace")

    ! Iterate through lons and lats, writing spectral tau
    DO j = 1, 144
       DO i = 1, 3
          DO k = 1,40
             WRITE (out_unit,*) spectral_tau(j,i,k)
          END DO
       END DO
    END DO

    ! Close file
    CLOSE (out_unit)

  END SUBROUTINE write_tau_output

  !-------------------------------------

END MODULE tau_output_mod
