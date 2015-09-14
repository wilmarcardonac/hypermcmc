Module fiducial

Implicit none

save 
!#############################
! Parameters of Fiducial model
!#############################

!Real*8,parameter :: omega_b = 2.225d-2
!Real*8,parameter :: omega_cdm = 1.198d-1
!Real*8,parameter :: n_s = 9.645d-1
!Real*8,parameter :: A_s = 2.20652d-9
!Real*8,parameter :: H0 = 6.727d1
!Real*8,parameter :: m_ncdm = 6.0d-2
!Real*8,parameter :: MG_beta2 = 1.00d0
!Real*8,parameter :: N_ur = 0.d0
!Real*8,parameter :: N_ncdm = 1.d0
!Real*8,parameter :: deg_ncdm = 3.d0
!Real*8,parameter :: tau = 0.089d0

Real*8,parameter :: prior_A = 12.5d0
Real*8,parameter :: prior_bw = -3.d0
Real*8,parameter :: prior_sigma_int = 1.d-1

!################################################
! 1-sigma values for parameters in fiducial model
!################################################

!Real*8,parameter :: sigma_omega_b = 1.6d-4
!Real*8,parameter :: sigma_omega_cdm = 1.5d-3
!Real*8,parameter :: sigma_n_s = 4.9d-3
!Real*8,parameter :: sigma_A_s = 1.541d-11
!Real*8,parameter :: sigma_H0 = 6.6d-1
!Real*8,parameter :: sigma_m_ncdm = 5.d-3
!Real*8,parameter :: sigma_MG_beta2 = 5.d-1

Real*8,parameter :: sigma_A = 1.8d-2
Real*8,parameter :: sigma_bw = 6.d-2
!Real*8,parameter :: sigma_sigma_int = 2.d2

!#####################
! Other specifications 
!#####################

!Integer*4,parameter :: nbins = 10
!Real*8,parameter :: zmin = 0.1d0
!Real*8,parameter :: zmax = 2.0d0
!Real*8,parameter :: dz = 1.0d-3
!Integer*4,parameter :: lmax = 400
!Integer*4,parameter :: lmin = 2
!Real*8,parameter :: gal_per_sqarcmn = 30.d0
Real*8,parameter :: Pi = 3.141592653589793d0
!Integer*4,parameter :: n_points = 5 ! Number of points per cosmological parameter
!Real*8,parameter :: fsky = 1.5d4/4.1253d4
!Real*8,parameter :: theoreticalerror = 0.d0 !5.d-2

!################
! MCMC parameters
!################

Integer*4,parameter :: number_iterations = 11d6        ! total number of iterations in MCMC run
Integer*4,parameter :: number_of_parameters = 2       ! number of cosmological parameters 
Integer*4,parameter :: jumping_factor_update = 5d2    ! number of steps before updating jumping factor (if needed)
Integer*4,parameter :: covariance_matrix_update = 1d4 ! number of steps before updating covariance matrix (if needed)
Real*8,parameter :: step_size_changes = 1.d-1        ! It helps to change step size 
Integer*4,parameter :: steps_taken_before_definite_run = 1d7
logical,parameter :: separate_dataA = .true.
logical,parameter :: separate_dataB = .true.
logical,parameter :: separate_dataC = .true.
logical,parameter :: include_dataA = .true.
logical,parameter :: include_dataB = .true.
logical,parameter :: include_dataC = .true.

End Module fiducial
