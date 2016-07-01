Module fiducial

    Implicit none

    save 

    !#############################
    ! PARAMETERS OF FIDUCIAL MODEL
    !#############################

    Real*8,parameter    :: prior_A = 12.5d0
    Real*8,parameter    :: prior_bw = -3.26d0
    Real*8,parameter    :: prior_sigma_int = 0.1d0!0.20d0 SN Ia hosts. MUST BE ABOUT 0.1
    Real*8,parameter    :: prior_sigma_int_LMC = 0.113d0 ! SAME VALUE AS IN EQUATION (4a) OF EFSTATHIOU'S PAPER
    Real*8,parameter    :: prior_sigma_int_MW = 0.10d0 ! SAME VALUE AS IN SUBSECTION 4.2 OF EFSTATHIOU'S PAPER
    Real*8,parameter    :: prior_alpha_j = 5.d-1
    Real*8,parameter    :: R = 0.410d0                  ! TAKEN FROM PAGE 7 IN R11
    Real*8,parameter    :: a_v = 0.71273d0               ! TAKEN FROM PAGE 16 IN R16. R16 USED B BAND (R11 V BAND), HENCE 'a_v' HERE STANDS FOR a_B 
    Real*8,parameter    :: a_cal = 0.d0                 ! AUXILIAR PARAMETER TO HAVE DIAGNONAL COVARIANCE MATRIX FOR LMC CEPHEID VARIABLES
    Real*8,parameter    :: NGC4258_distance = 7.54d0    ! TAKEN FROM PAGE 16 IN R16. UNITS : MPC
    Real*8,parameter    :: NGC4258_distance_2015 = 7.08d0    ! TAKEN FROM PAGE 1 IN J. POLSHAW. UNITS : MPC
    Real*8,parameter    :: mu_0_NGC4258 = 5.d0*log10(NGC4258_distance) + 25.d0 ! DEFINITION OF DISTANCE MODULUS
    Real*8,parameter    :: mu_0_NGC4258_2015 = 5.d0*log10(NGC4258_distance_2015) + 25.d0 ! DEFINITION OF DISTANCE MODULUS
    Real*8,parameter    :: LMC_distance = 49.97d-3       ! TAKEN FROM PAGE 76 IN PIETRZYNSKI. UNITS : MPC
    Real*8,parameter    :: mu_0_LMC = 5.d0*log10(LMC_distance) + 25.d0 
    Real*8,parameter    :: LMC_distance_2015 = 51.82d-3  ! FROM ABSTRACT IN M. M. FAUSNAUGH (2015). UNITS : MPC
    Real*8,parameter    :: mu_0_LMC_2015 = 5.d0*log10(LMC_distance_2015) + 25.d0
    Real*8,parameter    :: mu_0_M31 = 24.36d0   ! TAKEN FROM R16
    Real*8,parameter    :: meanOH_LMC = 8.5d0           ! MEAN METALLICITY FOR LMC CEPHEID VARIABLES ASSUMED BY EFSTATHIOU
    Real*8,parameter    :: meanOH_MW = 8.9d0           ! MEAN METALLICITY FOR MW CEPHEID VARIABLES ASSUMED BY EFSTATHIOU
    Real*8,parameter    :: prior_zpH = 28.d0
    Real*8,parameter    :: prior_bH = -2.7d0
    Real*8,parameter    :: prior_mu1 = 29.135d0 ! M101  FROM TABLE 5 IN R16 
    Real*8,parameter    :: prior_mu2 = 32.497d0 ! N1015 FROM TABLE 5 IN R16 
    Real*8,parameter    :: prior_mu3 = 32.523d0 ! N1309 FROM TABLE 5 IN R16 
    Real*8,parameter    :: prior_mu4 = 31.307d0 ! N1365 FROM TABLE 5 IN R16 
    Real*8,parameter    :: prior_mu5 = 31.311d0 ! N1448 FROM TABLE 5 IN R16
    Real*8,parameter    :: prior_mu6 = 31.511d0 ! N2442 FROM TABLE 5 IN R16 
    Real*8,parameter    :: prior_mu7 = 32.498d0 ! N3021 FROM TABLE 5 IN R16 
    Real*8,parameter    :: prior_mu8 = 32.072d0 ! N3370 FROM TABLE 5 IN R16 
    Real*8,parameter    :: prior_mu9 = 31.908d0 ! N3447 FROM TABLE 5 IN R16
    Real*8,parameter    :: prior_mu10 = 31.587d0! N3972 FROM TABLE 5 IN R16
    Real*8,parameter    :: prior_mu11 = 31.737d0! N3982 FROM TABLE 5 IN R16
    Real*8,parameter    :: prior_mu12 = 31.290d0! N4038 FROM TABLE 5 IN R16
    Real*8,parameter    :: prior_mu13 = 31.080d0! N4424 FROM TABLE 5 IN R16
    Real*8,parameter    :: prior_mu14 = 30.906d0! N4536 FROM TABLE 5 IN R16
    Real*8,parameter    :: prior_mu15 = 31.532d0! N4639 FROM TABLE 5 IN R16
    Real*8,parameter    :: prior_mu16 = 31.786d0! N5584 FROM TABLE 5 IN R16
    Real*8,parameter    :: prior_mu17 = 32.263d0! N5917 FROM TABLE 5 IN R16
    Real*8,parameter    :: prior_mu18 = 31.499d0! N7250 FROM TABLE 5 IN R16
    Real*8,parameter    :: prior_mu19 = 32.919d0! U9391 FROM TABLE 5 IN R16
    Real*8,parameter    :: prior_mu20 = 29.387d0! N4258
    Real*8,parameter    :: prior_mu21 = 24.36d0 ! M31
    Real*8,parameter    :: prior_mu22 = 18.5d0  ! LMC
    Real*8,parameter    :: prior_zpw = 29.d0
    Real*8,parameter    :: prior_zpw4258 = 30.5d0  ! CENTRAL VALUE OF PRIOR ON zp_{w,4258} 
    Real*8,parameter    :: prior_zpwLMC = 20.98d0
    Real*8,parameter    :: prior_Mw = -5.88d0
    Real*8,parameter    :: prior_Zw = 0.d0         ! CENTRAL VALUE FOR PRIOR ON Zw 
    Real*8,parameter    :: prior_H0 = 70.0d0
    Real*8,parameter    :: prior_bw_from_LMC = -3.3d0 ! CENTRAL VALUE FOR PRIOR ON bw

    !################################################
    ! 1-SIGMA VALUES FOR PARAMETERS IN FIDUCIAL MODEL
    !################################################

    Real*8,parameter    :: sigma_A = 1.8d-2
    Real*8,parameter    :: sigma_bw = 2.d-1
    Real*8,parameter    :: sigma_sigma_int = 1.d-1 !log10(sigma_int)
    Real*8,parameter    :: sigma_alpha_j = 1.d-3
    Real*8,parameter    :: sigma_a_v = 0.00176d0        ! TAKEN FROM PAGE 15 IN R16 (AGAIN THIS STANDS FOR 'sigma_a_v')
    Real*8,parameter    :: sigma_a_cal = 0.04d0         ! TAKEN FROM PAGE 10 IN EFSTATHIOU'S PAPER
    Real*8,parameter    :: sigma_NGC4258_quadrature =  0.1972d0    ! TAKEN FROM PAGE 16 IN R16. UNITS : MPC
    Real*8,parameter    :: sigma_NGC4258_2015 =  0.86d0    ! TAKEN FROM PAGE 1 IN J. POLSHAW. UNITS : MPC
    Real*8,parameter    :: sigma_mu_0_NGC4258 = 5.d0/log(10.d0)/NGC4258_distance*sigma_NGC4258_quadrature ! ERROR ON DISTANCE MODULUS
    Real*8,parameter    :: sigma_mu_0_NGC4258_2015 = 5.d0/log(10.d0)/NGC4258_distance_2015*sigma_NGC4258_2015 ! ERROR ON DISTANCE MODULUS
    Real*8,parameter    :: sigma_LMC_quadrature = 1.13d-3 ! TAKEN FROM PAGE 76 IN PIETRZYNSKI. UNITS : MPC
    Real*8,parameter    :: sigma_LMC_2015 = 3.23d-3  ! FROM ABSTRACT IN M. M. FAUSNAUGH (2015). UNITS : MPC
    Real*8,parameter    :: sigma_mu_0_LMC_2015 = 5.d0/log(10.d0)/LMC_distance_2015*sigma_LMC_2015
    Real*8,parameter    :: sigma_mu_0_LMC = 5.d0/log(10.d0)/LMC_distance*sigma_LMC_quadrature ! ERROR ON DISTANCE MODULUS
    Real*8,parameter    :: sigma_mu_0_M31 = 0.08d0        ! ERROR ON M31 DISTANCE MODULUS
    Real*8,parameter    :: sigma_zpH = 0.2d0
    Real*8,parameter    :: sigma_bH = 0.1d0
    Real*8,parameter    :: sigma_mu1 = sigma_mu_0_NGC4258 ! 0.045d0 
    Real*8,parameter    :: sigma_mu2 = sigma_mu_0_NGC4258 ! 0.081d0 
    Real*8,parameter    :: sigma_mu3 = sigma_mu_0_NGC4258 ! 0.055d0 
    Real*8,parameter    :: sigma_mu4 = sigma_mu_0_NGC4258 ! 0.057d0 
    Real*8,parameter    :: sigma_mu5 = sigma_mu_0_NGC4258 ! 0.045d0 
    Real*8,parameter    :: sigma_mu6 = sigma_mu_0_NGC4258 ! 0.053d0
    Real*8,parameter    :: sigma_mu7 = sigma_mu_0_NGC4258 ! 0.090d0
    Real*8,parameter    :: sigma_mu8 = sigma_mu_0_NGC4258 ! 0.049d0 
    Real*8,parameter    :: sigma_mu9 = sigma_mu_0_NGC4258 ! 0.043d0
    Real*8,parameter    :: sigma_mu10 = sigma_mu_0_NGC4258 ! 0.070d0 
    Real*8,parameter    :: sigma_mu11 = sigma_mu_0_NGC4258 ! 0.069d0
    Real*8,parameter    :: sigma_mu12 = sigma_mu_0_NGC4258 ! 0.112d0
    Real*8,parameter    :: sigma_mu13 = sigma_mu_0_NGC4258 ! 0.292d0
    Real*8,parameter    :: sigma_mu14 = sigma_mu_0_NGC4258 ! 0.053d0
    Real*8,parameter    :: sigma_mu15 = sigma_mu_0_NGC4258 ! 0.071d0
    Real*8,parameter    :: sigma_mu16 = sigma_mu_0_NGC4258 ! 0.046d0
    Real*8,parameter    :: sigma_mu17 = sigma_mu_0_NGC4258 ! 0.102d0
    Real*8,parameter    :: sigma_mu18 = sigma_mu_0_NGC4258 ! 0.078d0
    Real*8,parameter    :: sigma_mu19 = sigma_mu_0_NGC4258 ! 0.063d0 
    Real*8,parameter    :: sigma_mu20 = sigma_mu_0_NGC4258
    Real*8,parameter    :: sigma_mu21 = sigma_mu_0_M31
    Real*8,parameter    :: sigma_mu22 = sigma_mu_0_LMC
    Real*8,parameter    :: sigma_zpw = 1.d-1
    Real*8,parameter    :: sigma_zpw4258 = 0.1d0 
    Real*8,parameter    :: sigma_zpwLMC = 0.78d0
    Real*8,parameter    :: sigma_Mw = 0.05d0
    Real*8,parameter    :: sigma_Zw = 0.25d0
    Real*8,parameter    :: sigma_Zw_prior = 0.25d0 ! 0.02d0 STRONG PRIOR; 0.25d0 WEAK PRIOR
    Real*8,parameter    :: sigma_H0 = 2.4d0
    Real*8,parameter    :: sigma_bw_prior = 0.05d0

    !#####################
    ! OTHER SPECIFICATIONS
    !#####################

    Real*8,parameter    :: Pi = 3.141592653589793d0

    !################
    ! MCMC PARAMETERS
    !################

    Integer*4,parameter :: number_iterations = 11000000              ! TOTAL NUMBER OF ITERATIONS IN MCMC RUN
    Integer*4,parameter :: number_model_parameters = 28 ! NUMBER OF PARAMETERS IN MODEL : 2 FOR LMC ALONE, 10 FOR R11 DATA WITHOUT METALLICITY,
    ! 3 FOR CEPHEIDS ALONE (INCLUDING METALLICITY DEPENDENCE, AND FIXED sigma_int), 12 FOR ALL R11 CEPHEIDS, 14 FOR R11 DATA USING NGC4258 AS AN ANCHOR 
    ! INCLUDING METALLICITY AND REDDENING-FREE MAGNITUDE, 16 FOR ALL R11 CEPHEIDS + LMC CEPHEIDS AND USING LMC AS ANCHOR, 15 FOR ALL R11 CEPHEIDS +
    ! MW CEPHEIDS ANS USING MW AS ANCHOR, 16 FOR ALL R11 CEPHEIDS + NGC4258 AND LMC AS ANCHORS, 15 FOR ALL R11 CEPHEIDS + NGC4258 AND MW AS ANCHORS,
    ! 16 FOR ALL R11 CEPHEIDS + MW AND LMC AS ANCHORS, 16 FOR ALL R11 CEPHEIDS + NGC4258, LMC AND MW AS ANCHORS, 
    ! 28 FOR ALL R16 CEPHEIDS + NGC4258, LMC, MW, AND M31 AS ANCHORS, 26 FOR ALL R11 CEPHEIDS + NGC4258, LMC AND MW AS ANCHORS, SIGMA INT PER HOST GALAXY 
    ! BUT NOT INCLUDING METALLICITY, 4 FOR CEPHEIDS ALONE (INCLUDING METALLICITY DEPENDENCE, AND VARYING sigma_int)
    Integer*4,parameter :: number_hyperparameters = 0           ! NUMBER OF HYPER-PARAMETERS (MUST MATCH TOTAL NUMBER OF POINTS) 
    Integer*4,parameter :: number_of_parameters = number_model_parameters + number_hyperparameters ! TOTAL NUMBER OF PARAMETERS IN MODEL
    Integer*4,parameter :: jumping_factor_update = 100           ! NUMBER OF TAKEN STEPS BEFORE UPDATING JUMPING FACTOR (IF NEEDED)
    Integer*4,parameter :: covariance_matrix_update = 20000        ! STEPS TAKEN BEFORE UPDATING COVARIANCE MATRIX (IF NEEDED)
    Integer*4,parameter :: steps_taken_before_definite_run = 1000000 ! STEPS TAKEN BEFORE DEFINITE RUN
    Integer*4,parameter :: number_of_hosts_galaxies = 20 ! TOTAL NUMBER OF HOSTS GALAXIES AS IN R16 (NUMBER INCLUDES NGC4258)
    Integer*4,parameter :: UNIT_EXE_FILE = 90           ! UNIT NUMBER FOR EXECUTION INFORMATION FILE
    Integer*4,parameter :: UNIT_RANGES_FILE = 91           ! UNIT NUMBER FOR RANGES FILE
    Integer*4,parameter :: UNIT_PARAMNAMES_FILE = 92           ! UNIT NUMBER FOR PARAMNAMES FILE
    Integer*4,parameter :: UNIT_MCMC_FILE = 93           ! UNIT NUMBER FOR MCMC OUTPUT FILE (CALIBRATING PHASE)
    Integer*4,parameter :: UNIT_MCMC_FINAL_FILE = 94           ! UNIT NUMBER FOR MCMC FINAL OUTPUT FILE
    Integer*4,parameter :: UNIT_HP_FILE = 95           ! UNIT EFFECTIVE HPS FILE

    Real*8,parameter    :: step_size_changes = 1.d-2             ! CHANGES IN STEP SIZE
    Real*8,parameter    :: cepheid_Period_limit = 205.d0 !205.d0 !60.d0           ! DISREGARD CEPHEID VARIABLES WITH PERIOD GREATER THAN cepheid_Period_limit
    Real*8,parameter    :: cepheid_lower_Period_limit = 0.d0                    ! DISREGARD CEPHEID VARIABLES WITH PERIOD SHORTER THAN cepheid_lower_Period_limit

    Logical,parameter   :: separate_dataA = .false.!.true.               ! INCLUDE DATA SET A AS SINGLE POINTS IF SET IT TRUE
    Logical,parameter   :: separate_dataB = .false.!.true.               ! INCLUDE DATA SET B AS SINGLE POINTS IF SET IT TRUE
    Logical,parameter   :: separate_dataC = .false.!.true.               ! INCLUDE DATA SET C AS SINGLE POINTS IF SET IT TRUE
    Logical,parameter   :: include_dataA = .false.!.true.                ! INCLUDE DATA SET A IF SET IT TRUE
    Logical,parameter   :: include_dataB = .false.!.true.                ! INCLUDE DATA SET B IF SET IT TRUE
    Logical,parameter   :: include_dataC = .false.!.true.                ! INCLUDE DATA SET C IF SET IT TRUE
    Logical,parameter   :: include_table2_R11 = .true.            ! INCLUDE TABLE 2 IN R11 IF SET IT TRUE
    Logical,parameter   :: start_from_fiducial = .true.          ! START MCMC ANALYSIS FROM FIDUCIAL POINT IF SET IT TRUE 
    Logical,parameter   :: testing_Gaussian_likelihood = .false. ! TEST GAUSSIAN LIKELIHOOD IF SET IT TRUE
    Logical,parameter   :: using_hyperparameters = .true.        ! USE HYPER-PARAMETERS IF SET IT TRUE
    Logical,parameter   :: using_jeffreys_prior = .false.        ! USE JEFFREYS PRIOR IF SET IT TRUE, OTHERWISE USE UNIFORM PRIOR [0,1] 
    Logical,parameter   :: hyperparameters_as_mcmc = .false.      ! SET HYPER-PARAMETERS AS MCMC PARAMETERS IF SET IT TRUE
    Logical,parameter   :: use_NGC4258_as_anchor = .true.       ! USE NFC4258 AS ANCHOR IF SET IT TRUE
    Logical,parameter   :: use_LMC_as_anchor = .true.           ! USE LMC AS ANCHOR IF SET IT TRUE
    Logical,parameter   :: use_M31_as_anchor = .true.          ! USE M31 AS ANCHOR IF SET IT TRUE
    Logical,parameter   :: use_MW_as_anchor = .true.            ! USE MW AS ANCHOR IF SET IT TRUE
    Logical,parameter   :: use_metallicity = .true.             ! USE METALLICITY DEPENDENCE IF SET IT TRUE
    Logical,parameter   :: use_H_band = .false.!.true.                   ! USE H BAND IF SET IT TRUE, OTHERWISE USE W BAND
    Logical,parameter   :: use_HP_in_SNIa = .true.               ! USE HPs WHEN COMPUTING SNIa CHI2
    Logical,parameter   :: use_HP_in_av = .false.                ! USE HPs WHEN COMPUTINNG av CHI2
    Logical,parameter   :: use_HP_in_anchor = .true.            ! USE HPs WHEN COMPUTING ANCHOR CHI2
    Logical,parameter   :: use_HP_per_host = .false.              ! USE HPs FOR EACH HOST IN R11 IF SET IT TRUE
    Logical,parameter   :: use_HP_per_cepheid = .true.           ! USE HPs FOR EACH CEPHEID IN R11 IF SET IT TRUE
    Logical,parameter   :: use_HP_per_MW_cepheid = .true.       ! USE HPs FOR EACH CEPHEID IN MW IF SET IT TRUE
    Logical,parameter   :: use_HP_for_MW_dataset = .false.       ! USE HP FOR MW DATASET IF SET IT TRUE (JEFFREY'S PRIOR)
    Logical,parameter   :: doing_R11_analysis = .true.           ! DO R16 ANALYSIS IF SET IT TRUE, OTHERWISE DO EFSTATHIOU'S SECTION 2 (LMC CEPHEIDS ALONE)
    Logical,parameter   :: include_only_cepheids = .false.       ! INCLUDE ONLY CEPHEIDS DATA IF SET IT TRUE
    Logical,parameter   :: all_R11_hosts = .false.             ! INCLUDE ALL CEPHEIDS IN R11 SAMPLE SIMULTANEOUSLY IF SET IT TRUE
    Logical,parameter   :: use_prior_on_zpw4258 = .false. !.true.       ! USE PRIOR ON zp_{w,4258} IS SET IT TRUE
    Logical,parameter   :: use_prior_on_Zw = .true.              ! USE PRIOR ON Zw IF SET IT TRUE 
    Logical,parameter   :: use_prior_on_bw = .false.              ! USE PRIOR ON bw IF SET IT TRUE
    Logical,parameter   :: use_HP_in_Zw = .false.                 ! USE HPs WHEN USING PRIOR ON THE METALLICITY IF SET IT TRUE 
    Logical,parameter   :: varying_sigma_int = .true.            ! TRUE IF VARYING sigma_int IN MCMC WHEN NO sigma_int_per_R11_host, SET TO FALSE OTHERWISE
    Logical,parameter   :: sigma_int_per_R11_host = .true.        ! TRUE FOR MAIN ANALYSIS: IT INCLUDES SIGMA INT PER R11 HOST 
    Logical,parameter   :: include_mu_0_NGC4258_2015 = .true.    ! TRUE TO INCLUDE DISTANCE TO NGC4258 FROM 2015, FALSE TO INCLUDE ONLY MEASUREMENT FROM 2016

    Character(len=*),parameter :: path_to_datafileA = './data/dataA.txt'    ! PATH TO DATA SET A
    Character(len=*),parameter :: path_to_datafileB = './data/dataB.txt'    ! PATH TO DATA SET B 
    Character(len=*),parameter :: path_to_datafileC = './data/dataC.txt'    ! PATH TO DATA SET C
    Character(len=*),parameter :: path_to_datafileAB = './data/dataAB.txt'    ! PATH TO JOINT DATA SET AB
    Character(len=*),parameter :: path_to_datafileABC = './data/dataABC.txt'    ! PATH TO JOINT DATA SET ABC
    Character(len=*),parameter :: path_to_datafileM31 = './data/M31.txt'    ! PATH TO M31 CEPHEID DATA 
    Character*16,parameter     :: phrase = 'randomizer'    ! PHRASE NEEDED BY RANDOM NUMBER GENERATOR 
    character(len=*),parameter :: path_to_table2_R11 = './data/table4_R16.txt' ! PATH TO DATA OF TABLE 4 IN R16
    character(len=*),parameter :: path_to_table3_R11 = './data/table5_R16.txt' ! PATH TO DATA OF TABLE 5 IN R16
    character(len=*),parameter :: path_to_table2_LEEUWEN = './data/table2_Leeuwen.txt' ! PATH TO DATA OF TABLE 2 IN LEEUWEN
    Character(len=5),dimension(number_of_hosts_galaxies), parameter :: host = ['m101 ','n1015','n1309','n1365',&
         'n1448','n2442','n3021','n3370','n3447','n3972','n3982','n4038','n4424','n4536','n4639','n5584',&
         'n5917','n7250','u9391','n4258'] ! HOST GALAXIES IN SAME ORDER LISTED IN TABLE 5 OF R16
    Character(len=*),parameter :: EXECUTION_INFORMATION = './output/chains/execution_information.txt' ! PATH TO EXECUTION INFORMATION FILE

End Module fiducial
