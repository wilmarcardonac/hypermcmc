Module fiducial

    Implicit none

    save 

    !#############################
    ! PARAMETERS OF FIDUCIAL MODEL
    !#############################

    Real*8,parameter    :: prior_A = 12.5d0
    Real*8,parameter    :: prior_bw = -3.d0
    Real*8,parameter    :: prior_sigma_int = 1.d-1

    !################################################
    ! 1-SIGMA VALUES FOR PARAMETERS IN FIDUCIAL MODEL
    !################################################

    Real*8,parameter    :: sigma_A = 1.8d-2
    Real*8,parameter    :: sigma_bw = 6.d-2
    !Real*8,parameter :: sigma_sigma_int = 2.d2

    !#####################
    ! OTHER SPECIFICATIONS
    !#####################

    Real*8,parameter    :: Pi = 3.141592653589793d0

    !################
    ! MCMC PARAMETERS
    !################

    Integer*4,parameter :: number_iterations = 11d5              ! TOTAL NUMBER OF ITERATIONS IN MCMC RUN
    Integer*4,parameter :: number_of_parameters = 2              ! NUMBER OF PARAMETERS
    Integer*4,parameter :: jumping_factor_update = 2d2           ! NUMBER OF TAKEN STEPS BEFORE UPDATING JUMPING FACTOR (IF NEEDED)
    Integer*4,parameter :: covariance_matrix_update = 2d3        ! STEPS TAKEN BEFORE UPDATING COVARIANCE MATRIX (IF NEEDED)
    Integer*4,parameter :: steps_taken_before_definite_run = 3d5 ! STEPS TAKEN BEFORE DEFINITE RUN

    Real*8,parameter    :: step_size_changes = 1.d-1             ! CHANGES IN STEP SIZE

    Logical,parameter   :: separate_dataA = .true.               ! INCLUDE DATA SET A AS SINGLE POINTS IF SET IT TRUE
    Logical,parameter   :: separate_dataB = .true.               ! INCLUDE DATA SET B AS SINGLE POINTS IF SET IT TRUE
    Logical,parameter   :: separate_dataC = .true.               ! INCLUDE DATA SET C AS SINGLE POINTS IF SET IT TRUE
    Logical,parameter   :: include_dataA = .true.                ! INCLUDE DATA SET A IF SET IT TRUE
    Logical,parameter   :: include_dataB = .true.                ! INCLUDE DATA SET B IF SET IT TRUE
    Logical,parameter   :: include_dataC = .true.                ! INCLUDE DATA SET C IF SET IT TRUE
    Logical,parameter   :: start_from_fiducial = .true.         ! START MCMC ANALYSIS FROM FIDUCIAL POINT IF SET IT TRUE 
    Logical,parameter   :: testing_Gaussian_likelihood = .false. ! TEST GAUSSIAN LIKELIHOOD IF SET IT TRUE
    Logical,parameter   :: using_hyperparameters = .true.        ! USE HYPERPARAMETERS IF SET IT TRUE
    Logical,parameter   :: using_jeffreys_prior = .true. ! USE JEFFREYS PRIOR IF SET IT TRUE, OTHERWISE USE UNIFORM PRIOR [0,1] 

    Character(len=*),parameter :: path_to_datafileA = './data/dataA.txt'    ! PATH TO DATA SET A
    Character(len=*),parameter :: path_to_datafileB = './data/dataB.txt'    ! PATH TO DATA SET B 
    Character(len=*),parameter :: path_to_datafileC = './data/dataC.txt'    ! PATH TO DATA SET C
    Character(len=*),parameter :: path_to_datafileAB = './data/dataAB.txt'    ! PATH TO JOINT DATA SET AB
    !Character(len=*),parameter :: path_to_datafileABC = './data/dataABC.txt'    ! PATH TO JOINT DATA SET ABC
    Character*16,parameter :: phrase = 'randomizer'    ! PHRASE NEEDED BY RANDOM NUMBER GENERATOR 

End Module fiducial
