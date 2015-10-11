Program mcmc 

!####################
! LOAD NEEDED MODULES 
!####################

    use fiducial
    use arrays
    use functions 

!#################################
! DECLARE VARIABLES AND PARAMETERS
!#################################

    Implicit none

    Integer*4 :: m,n,i,q                                        ! INTERGER FOR SHORT LOOPS 
    Integer*4 :: seed1,seed2                                      ! SEEDS FOR RANDOM NUMBER GENERATOR 
    Integer*4 :: number_accepted_points,number_rejected_points    ! MCMC PARAMETERS
    Integer*4 :: weight                                           ! IT COUNTS THE NUMBER OF STEPS TAKEN BEFORE MOVING TO A NEW POINT IN MCMC 

    Real*4 :: average_acceptance_probability           ! SAVES ACCEPTANCE PROBABILITY 
    Real*4 :: genunf                            ! RANDOM UNIFOR DEVIATES 
    Real*8 :: jumping_factor                           ! SAVES JUMPING FACTOR FOR MCMC (INCREASE IF WANT BIGGER STEP SIZE, DECREASE OTHERWISE) 
    Real*8 :: random_uniform                           ! SAVES RANDOM UNIFORM DEVIATE BETWEEN 0 AND 1 
    Real*8 :: old_loglikelihood,current_loglikelihood  ! TEMPORALY SAVES LIKELIHOOD VALUES 
    Real*4,dimension(number_of_parameters*(number_of_parameters+3)/2 + 1) :: parm     ! ARRAY NEEDED BY RANDON NUMBER GENERATOR 
    Real*8,dimension(number_of_parameters,number_of_parameters)           :: Covgauss ! COVARIANCE MATRIX OF GAUSSIAN LIKELIHOOD 
    Real*8,dimension(number_of_parameters,number_of_parameters)           :: Covguess ! COVARIANCE MATRIX 
    Real*4,dimension(number_of_parameters) :: work,x_old,x_new                        ! ARRAYS NEEDED BY RANDOM NUMBER GENERATOR 
    Real*8,dimension(number_of_parameters) :: bestfit,means                           ! SAVES BESTFIT AND MEANS OF PARAMETERS 

    Logical :: not_good_aap,non_plausible_parameters   ! IT CONTROLS PLAUSIBLE VALUES OF COSMOLOGICAL PARAMETERS
    Logical,dimension(number_of_parameters) :: plausibility  

    Character(len=10) :: string ! STORES STRINGS FOR INTEGERS
    Character(len=12),dimension(number_hyperparameters) :: alpha_string

!##########################################################
! ASSIGNMENTS AND INITIALIZATION OF RANDOM NUMBER GENERATOR
!##########################################################

    weight = 1

    number_rejected_points = 0

    number_accepted_points = 0

    call initialize()                               ! INITIALIZE RANDOM NUMBER GENERATOR

    call phrtsd(phrase,seed1,seed2)                 ! GENERATE SEEDS FOR RANDOM NUMBERS FROM PHRASE

    call set_initial_seed(seed1,seed2)              ! SET INITIAL SEEDS FOR RANDOM NUMBER GENERATOR

    ! ALLOCATING MEMORY FOR POINTS IN PARAMETER SPACE AND ACCEPTANCE PROBABILITY
    allocate (old_point(1:number_of_parameters),current_point(1:number_of_parameters),&
    acceptance_probability(number_iterations),stat = status1)

    If (testing_Gaussian_likelihood) then

        ! SETTING COVARIANCE MATRIX
        Do m=1,number_of_parameters

            Do n=1,number_of_parameters

                If (m .eq. n) then     

                    Covgauss(m,n) = 1.d0

                Else 

                    Covgauss(m,n) = 0.d0

                End If

            End Do

        End Do
        ! COVARIANCE MATRIX SET
        jumping_factor = 2.38d0/sqrt(dble(number_of_parameters))    ! MODIFY ACCORDING TO WANTED INITIAL ACCEPTANCE PROBABILITY
        ! COVARIANCE MATRIX ADJUSTED
        Covgauss = jumping_factor*Covgauss

        open(15,file='./output/execution_information.txt')    ! OPEN FILE FOR EXECUTION INFORMATION 

    Else

        ! SETTING COVARIANCE MATRIX
        Covguess = 0.d0
  
        If (doing_R11_analysis) then

            If ( ( use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. use_MW_as_anchor) then
    
               print *,'USE OF THREE ANCHORS SIMULTANEOUSLY NOT IMPLEMENTED YET'

               stop

            Else If ( ( use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then

               print *,'NGC4258+LMC NOT IMPLEMENTED YET'

               stop

            Else If ( ( use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. use_MW_as_anchor ) then

               print *,'NGC4258+MW NOT IMPLEMENTED YET'

               stop

            Else If ( ( .not.use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. use_MW_as_anchor ) then

               print *,'MW+LMC NOT IMPLEMENTED YET'

               stop

            Else If ( ( .not.use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. use_MW_as_anchor ) then

               print *,'MW NOT IMPLEMENTED YET'

               stop

            Else If ( ( .not.use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then

               print *,'LMC NOT IMPLEMENTED YET'

               stop

            Else If ( ( use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then

               If (use_metallicity) then 

                  If (use_H_band) then
                     
                     print *,'H BAND NOT IMPLEMENTED INCLUDING METALLICITY DEPENDENCE'
                     
                     stop
                     
                  Else

                     If (number_model_parameters .eq. 13) then

                        Covguess(1,1) = sigma_mu1**2 

                        Covguess(2,2) = sigma_mu2**2 

                        Covguess(3,3) = sigma_mu3**2 

                        Covguess(4,4) = sigma_mu4**2 

                        Covguess(5,5) = sigma_mu5**2 

                        Covguess(6,6) = sigma_mu6**2 

                        Covguess(7,7) = sigma_mu7**2 

                        Covguess(8,8) = sigma_mu8**2 

                        Covguess(9,9) = sigma_mu9**2 

                        Covguess(10,10) = sigma_zpw**2

                        Covguess(11,11) = sigma_bw**2 

                        Covguess(12,12) = sigma_m0v_ref**2 

                        Covguess(13,13) = sigma_Zw**2 

                     Else
                        
                        print *,'WRONG NUMBER OF MODEL PARAMETERS. CHECK FIDUCIAL MODULE AND COMPARE WITH EQUATION (18) IN R09'

                        stop

                     End If

                  End If

               Else

                  If (use_H_band) then
                     
                     If (number_model_parameters .eq. 10) then

                        Covguess(1,1) = sigma_zpH**2 

                        Covguess(2,2) = sigma_zpH**2 

                        Covguess(3,3) = sigma_zpH**2 

                        Covguess(4,4) = sigma_zpH**2 

                        Covguess(5,5) = sigma_zpH**2 

                        Covguess(6,6) = sigma_zpH**2 

                        Covguess(7,7) = sigma_zpH**2 

                        Covguess(8,8) = sigma_zpH**2 

                        Covguess(9,9) = sigma_zpH**2 

                        Covguess(10,10) = sigma_bH**2

                     Else
                        
                        print *,'WRONG NUMBER OF MODEL PARAMETERS. CHECK FIDUCIAL MODULE AND COMPARE WITH EQUATION (3) IN R09'

                        stop

                     End If

                  Else

                     print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE'
                     
                     stop

                  End If

               End If

            Else

               print *, 'USER MUST SET TRUE AT LEAST ONE ANCHOR DISTANCE IN FIDUCIAL MODULE'

               stop

            End If

    
        Else

            Covguess(1,1) = sigma_A**2 

            Covguess(2,2) = sigma_bw**2

            !Covguess(3,3) = sigma_sigma_int**2
 
        End If

        If (hyperparameters_as_mcmc) then
            ! SETTING PIECE OF COVARIANCE MATRIX FOR HYPER-PARAMETERS
            Do m=number_model_parameters+1,number_of_parameters

                Covguess(m,m) = sigma_alpha_j**2

            End Do
            
        End If
        
        ! COVARIANCE MATRIX SET

        jumping_factor = 2.38d0/sqrt(dble(number_of_parameters))!*1.d-3    !    MODIFY ACCORDING TO WANTED ACCEPTANCE PROBABILITY

        ! COVARIANCE MATRIX ADJUSTED 
        Covguess = jumping_factor*Covguess

        If (using_hyperparameters) then

            If (include_dataA) then

                call read_data_EfstathiouA(path_to_datafileA)

            End If

            If (include_dataB) then

                call read_data_EfstathiouB(path_to_datafileB)

            End If

            If (include_dataC) then

                call read_data_EfstathiouC(path_to_datafileC)

            End If

            If (include_table2_R11) then

                call read_table2_R11(path_to_table2_R11)

                call read_table3_R11(path_to_table3_R11)

            End If
            
            If (hyperparameters_as_mcmc) then

                open(15,file='./output/execution_information_HP_as_MCMC.txt')    ! OPEN FILE FOR EXECUTION INFORMATION
            
                write(15,*) 'WORKING WITH HYPER-PARAMETERS'     

                write(15,*) 'HYPER-PARAMETERS USED AS MCMC PARAMETERS'

                call read_data_Efstathiou(path_to_datafileABC)

                If (number_hyperparameters .ne. (size(NameA)+size(NameB)+size(NameC))) then
                
                    write(15,*) 'NUMBER OF HYPER-PARAMETERS MUST MATCH TOTAL NUMBER OF DATA POINTS ',&
                    size(NameA)+size(NameB)+size(NameC)

                    write(15,*) 'CHECK FIDUCIAL MODULE'
       
                    stop

                End If

            Else

                If (doing_R11_analysis) then

                   open(15,file='./output/execution_information_HP_R11.txt')    ! OPEN FILE FOR EXECUTION INFORMATION
            
                   write(15,*) 'WORKING WITH HYPER-PARAMETERS AND DOING R11 ANALYSIS'
                   
                Else

                   open(15,file='./output/execution_information_HP.txt')    ! OPEN FILE FOR EXECUTION INFORMATION
            
                   write(15,*) 'WORKING WITH HYPER-PARAMETERS'
   
                End If

                If (number_hyperparameters .ne. 0) then
                
                    write(15,*) 'NUMBER OF HYPER-PARAMETERS MUST BE ZERO WHEN SETTING FALSE "hyperparameters_as_mcmc" '
 
                    write(15,*) 'CHECK FIDUCIAL MODULE'
       
                    stop

                End If

            End If


        Else

            call read_data_Efstathiou(path_to_datafileAB)

            If (include_table2_R11) then

                call read_table2_R11(path_to_table2_R11)

            End If

            open(15,file='./output/execution_information.txt')    ! OPEN FILE FOR EXECUTION INFORMATION 

        End If

    End If

!##################################
! MARKOV CHAIN MONTE CARLO ANALYSIS
!##################################

    write(15,*) 'STARTING MCMC ANALYSIS'

    If (testing_Gaussian_likelihood) then

        write(15,*) 'TESTING CODE WITH GAUSSIAN LIKELIHOOD'

        open(13,file='./output/mcmc_final_output.txt')

        open(17,file='./output/mcmc_final_output.ranges')    !    OPEN FILE WITH HARD BOUNDS NEEDED BY GETDIST

        write(17,*) 'A    N    N '

        write(17,*) 'bw    N    N '

    !    write(17,*) 'sigma_int    N    N '

        close(17)

        Do i=1,number_of_parameters

            x_old(i) = genunf(-1.,1.)

        End Do

        Do m=1,number_of_parameters
    
            If (m .eq. 4) then

                old_point(m) = dble(x_old(m)) !exp(dble(x_old(m)))/(1.d1**1.d1)

            else

                old_point(m) = dble(x_old(m))

            End If

        End Do

        old_loglikelihood = log_Gaussian_likelihood(old_point)

    Else 
    !###########################################
    ! GENERATE A RANDOM POINT IN PARAMETER SPACE
    ! RANDOM NUMBER GENERATOR WORKS WITH SINGLE PRECISION WHEREAS THIS CODES USES DOUBLE PRECISION; CHANGES ARE CORRESPONDINGLY MADE.
    !################################################################################################################################

        If (start_from_fiducial) then

            write(15,*) 'STARTING FROM FIDUCIAL POINT'

            If (doing_R11_analysis) then

               If ( ( use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. use_MW_as_anchor) then
    
                  print *,'USE OF THREE ANCHORS SIMULTANEOUSLY NOT IMPLEMENTED YET'

                  stop

               Else If ( ( use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then

                  print *,'NGC4258+LMC NOT IMPLEMENTED YET'

                  stop

               Else If ( ( use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. use_MW_as_anchor ) then

                  print *,'NGC4258+MW NOT IMPLEMENTED YET'

                  stop

               Else If ( ( .not.use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. use_MW_as_anchor ) then

                  print *,'MW+LMC NOT IMPLEMENTED YET'

                  stop

               Else If ( ( .not.use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. use_MW_as_anchor ) then

                  print *,'MW NOT IMPLEMENTED YET'

                  stop

               Else If ( ( .not.use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then

                  print *,'LMC NOT IMPLEMENTED YET'

                  stop

               Else If ( ( use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then

                  If (use_metallicity) then 

                     If (use_H_band) then
                     
                        print *,'H BAND NOT IMPLEMENTED INCLUDING METALLICITY DEPENDENCE'
                     
                        stop

                     Else

                        old_point(1) = prior_mu1

                        old_point(2) = prior_mu2

                        old_point(3) = prior_mu3 

                        old_point(4) = prior_mu4

                        old_point(5) = prior_mu5

                        old_point(6) = prior_mu6

                        old_point(7) = prior_mu7

                        old_point(8) = prior_mu8

                        old_point(9) = prior_mu9

                        old_point(10) = prior_zpw

                        old_point(11) = prior_bw

                        old_point(12) = prior_m0v_ref

                        old_point(13) = prior_Zw
                        
                     End If

                  Else

                     If (use_H_band) then
                     
                        old_point(1) = prior_zpH

                        old_point(2) = prior_zpH

                        old_point(3) = prior_zpH 

                        old_point(4) = prior_zpH

                        old_point(5) = prior_zpH

                        old_point(6) = prior_zpH

                        old_point(7) = prior_zpH

                        old_point(8) = prior_zpH

                        old_point(9) = prior_zpH

                        old_point(10) = prior_bH

                     Else

                        print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE'
                     
                        stop
                        
                     End If

                  End If

               Else

                  print *, 'USER MUST SET TRUE AT LEAST ONE ANCHOR DISTANCE IN FIDUCIAL MODULE'

                  stop

               End If
    
            Else

               old_point(1) = prior_A         ! A 
               old_point(2) = prior_bw        ! bw 
               !old_point(3) = prior_sigma_int ! sigma_int 
 
            End If

            If (hyperparameters_as_mcmc) then
                ! SETTING INITIAL POINT FOR HYPER-PARAMETERS
                Do m=number_model_parameters+1,number_of_parameters

                    old_point(m) = prior_alpha_j

                End Do
            
            End If
    
            Do m=1,number_of_parameters

                If (m .gt. number_model_parameters) then
                    
                    If (using_jeffreys_prior) then

                        x_old(m) = real(log10(old_point(m)))    ! CONVERT TO log10(\alpha_j)

                    Else

                        x_old(m) = real(old_point(m))

                    End If

                else

                    x_old(m) = real(old_point(m))

                End If

            End Do

        Else

            If (doing_R11_analysis) then

               If ( ( use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. use_MW_as_anchor) then
    
                  print *,'USE OF THREE ANCHORS SIMULTANEOUSLY NOT IMPLEMENTED YET'

                  stop

               Else If ( ( use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then

                  print *,'NGC4258+LMC NOT IMPLEMENTED YET'

                  stop

               Else If ( ( use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. use_MW_as_anchor ) then

                  print *,'NGC4258+MW NOT IMPLEMENTED YET'

                  stop

               Else If ( ( .not.use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. use_MW_as_anchor ) then

                  print *,'MW+LMC NOT IMPLEMENTED YET'

                  stop

               Else If ( ( .not.use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. use_MW_as_anchor ) then

                  print *,'MW NOT IMPLEMENTED YET'

                  stop

               Else If ( ( .not.use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then

                  print *,'LMC NOT IMPLEMENTED YET'

                  stop

               Else If ( ( use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then

                  If (use_metallicity) then 

                     If (use_H_band) then

                        print *,'H BAND NOT IMPLEMENTED INCLUDING METALLICITY DEPENDENCE'
                     
                        stop

                     Else
                     
                        x_old(1) = genunf(real(prior_mu1 - sigma_mu1),real(prior_mu1 + sigma_mu1))

                        x_old(2) = genunf(real(prior_mu2 - sigma_mu2),real(prior_mu2 + sigma_mu2))

                        x_old(3) = genunf(real(prior_mu3 - sigma_mu3),real(prior_mu3 + sigma_mu3))

                        x_old(4) = genunf(real(prior_mu4 - sigma_mu4),real(prior_mu4 + sigma_mu4))

                        x_old(5) = genunf(real(prior_mu5 - sigma_mu5),real(prior_mu5 + sigma_mu5))

                        x_old(6) = genunf(real(prior_mu6 - sigma_mu6),real(prior_mu6 + sigma_mu6))

                        x_old(7) = genunf(real(prior_mu7 - sigma_mu7),real(prior_mu7 + sigma_mu7))

                        x_old(8) = genunf(real(prior_mu8 - sigma_mu8),real(prior_mu8 + sigma_mu8))

                        x_old(9) = genunf(real(prior_mu9 - sigma_mu9),real(prior_mu9 + sigma_mu9))

                        x_old(10) = genunf(real(prior_zpw - sigma_zpw),real(prior_zpw + sigma_zpw))

                        x_old(11) = genunf(real(prior_bw - sigma_bw),real(prior_bw + sigma_bw))

                        x_old(12) = genunf(real(prior_m0v_ref - sigma_m0v_ref),real(prior_m0v_ref + sigma_m0v_ref))

                        x_old(13) = genunf(real(prior_Zw - sigma_Zw),real(prior_Zw + sigma_Zw))

                     End If

                  Else

                     If (use_H_band) then
                     
                        x_old(1) = genunf(real(prior_zpH - sigma_zpH),real(prior_zpH + sigma_zpH))

                        x_old(2) = genunf(real(prior_zpH - sigma_zpH),real(prior_zpH + sigma_zpH))

                        x_old(3) = genunf(real(prior_zpH - sigma_zpH),real(prior_zpH + sigma_zpH))

                        x_old(4) = genunf(real(prior_zpH - sigma_zpH),real(prior_zpH + sigma_zpH))

                        x_old(5) = genunf(real(prior_zpH - sigma_zpH),real(prior_zpH + sigma_zpH))

                        x_old(6) = genunf(real(prior_zpH - sigma_zpH),real(prior_zpH + sigma_zpH))

                        x_old(7) = genunf(real(prior_zpH - sigma_zpH),real(prior_zpH + sigma_zpH))

                        x_old(8) = genunf(real(prior_zpH - sigma_zpH),real(prior_zpH + sigma_zpH))

                        x_old(9) = genunf(real(prior_zpH - sigma_zpH),real(prior_zpH + sigma_zpH))

                        x_old(10) = genunf(real(prior_bH - sigma_bH),real(prior_bH + sigma_bH))

                     Else

                        print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE'
                     
                        stop
                        
                     End If

                  End If

               Else

                  print *, 'USER MUST SET TRUE AT LEAST ONE ANCHOR DISTANCE IN FIDUCIAL MODULE'

                  stop

               End If
    
            Else

               x_old(1) = genunf(real(prior_A - sigma_A),real(prior_A + sigma_A))         ! A
               x_old(2) = genunf(real(prior_bw - sigma_bw),real(prior_bw + sigma_bw)) ! bw
               !x_old(3) = genunf(real(-10.d0),real(0.d0)) ! log10(sigma_int)

            End If

            If (hyperparameters_as_mcmc) then
                ! SETTING INITIAL POINT FOR HYPER-PARAMETERS
                Do m=number_model_parameters+1,number_of_parameters
                    
                    If (using_jeffreys_prior) then

                        x_old(m) = genunf(real(-1.d0),real(1.d0)) ! log10(\alpha_j)

                    Else

                        x_old(m) = genunf(real(0.d0),real(1.d0))

                    End If

                End Do
            
            End If

            Do m=1,number_of_parameters

                If (m .gt. number_model_parameters) then
                  
                    If (using_jeffreys_prior) then

                        old_point(m) = 10**(dble(x_old(m)))    !    CONVERT TO \alpha_j

                    Else

                        old_point(m) = dble(x_old(m)) 

                    End If

                else

                    old_point(m) = dble(x_old(m))

                End If

            End Do

        End If

        If (using_hyperparameters) then
           ! OPEN FILE TO STORE MCMC COMPUTATION
           open(13,file='./output/mcmc_final_output_HP.txt')

           write(15,*) 'WORKING WITH HYPER-PARAMETERS'

           write(15,*) 'COMPUTING LOG_LIKELIHOOD FOR INITIAL POINT'
            
           ! COMPUTE INITIAL LIKELIHOOD
           If (doing_R11_analysis) then
              
              If (use_H_band) then

                 old_loglikelihood = log_R11_likelihood_H(old_point(1:number_model_parameters-1),&
                      old_point(number_model_parameters),prior_sigma_int)

              Else

                 old_loglikelihood = log_R11_likelihood_W(old_point(1:number_model_parameters-4),&
                      old_point(number_model_parameters-3),old_point(number_model_parameters-2),&
                      old_point(number_model_parameters-1),old_point(number_model_parameters),prior_sigma_int)

              End If

           Else

              old_loglikelihood = log_Efstathiou_likelihood_hyperparameters(old_point(1),old_point(2),prior_sigma_int)

           End If

           open(16,file='./output/mcmc_final_output_HP.paramnames')    !    OPEN FILE WITH PARAMETER NAMES NEEDED BY GETDIST

           open(17,file='./output/mcmc_final_output_HP.ranges')    !    OPEN FILE WITH HARD BOUNDS NEEDED BY GETDIST

        Else
            ! OPEN FILE TO STORE MCMC COMPUTATION
            open(13,file='./output/mcmc_final_output.txt')

            write(15,*) 'NOT WORKING WITH HYPER-PARAMETERS' 

            write(15,*) 'COMPUTING LOG_LIKELIHOOD FOR INITIAL POINT'

            ! COMPUTE INITIAL LIKELIHOOD
            old_loglikelihood = log_Efstathiou_likelihood(old_point(1),old_point(2),prior_sigma_int) 

            open(16,file='./output/mcmc_final_output.paramnames')    !    OPEN FILE WITH PARAMETER NAMES NEEDED BY GETDIST

            open(17,file='./output/mcmc_final_output.ranges')    !    OPEN FILE WITH HARD BOUNDS NEEDED BY GETDIST 

        End If

        If (doing_R11_analysis) then

           If ( ( use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. use_MW_as_anchor) then
    
              print *,'USE OF THREE ANCHORS SIMULTANEOUSLY NOT IMPLEMENTED YET'

              stop

           Else If ( ( use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then

              print *,'NGC4258+LMC NOT IMPLEMENTED YET'

              stop

           Else If ( ( use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. use_MW_as_anchor ) then

              print *,'NGC4258+MW NOT IMPLEMENTED YET'

              stop

           Else If ( ( .not.use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. use_MW_as_anchor ) then

              print *,'MW+LMC NOT IMPLEMENTED YET'

              stop

           Else If ( ( .not.use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. use_MW_as_anchor ) then

              print *,'MW NOT IMPLEMENTED YET'

              stop

           Else If ( ( .not.use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then

              print *,'LMC NOT IMPLEMENTED YET'

              stop

           Else If ( ( use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then

              If (use_metallicity) then 

                 If (use_H_band) then

                    print *,'H BAND NOT IMPLEMENTED INCLUDING METALLICITY DEPENDENCE'
                     
                    stop

                 Else

                    write(16,*) 'mu01    \mu_{0,1}'

                    write(16,*) 'mu02    \mu_{0,2}'

                    write(16,*) 'mu03    \mu_{0,3}'

                    write(16,*) 'mu04    \mu_{0,4}'

                    write(16,*) 'mu05    \mu_{0,5}'

                    write(16,*) 'mu06    \mu_{0,6}'

                    write(16,*) 'mu07    \mu_{0,7}'

                    write(16,*) 'mu08    \mu_{0,8}'

                    write(16,*) 'mu04258    \mu_{0,4258}'

                    write(16,*) 'zpw4258    zp_{w,4258}'

                    write(16,*) 'bw    b_w'

                    write(16,*) 'm0v4258    m^0_{v,4258}'

                    write(16,*) 'Zw    Z_w'

                 End If

              Else

                 If (use_H_band) then

                    write(16,*) 'zpH1    zp_{H1}'

                    write(16,*) 'zpH2    zp_{H2}'

                    write(16,*) 'zpH3    zp_{H3}'

                    write(16,*) 'zpH4    zp_{H4}'

                    write(16,*) 'zpH5    zp_{H5}'

                    write(16,*) 'zpH6    zp_{H6}'

                    write(16,*) 'zpH7    zp_{H7}'

                    write(16,*) 'zpH8    zp_{H8}'

                    write(16,*) 'zpH4258    zp_{H4258}'

                    write(16,*) 'bH    b_H'

                 Else

                    print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE'
                     
                    stop
                        
                 End If

              End If

           Else

              print *, 'USER MUST SET TRUE AT LEAST ONE ANCHOR DISTANCE IN FIDUCIAL MODULE'

              stop

           End If
    
        Else

           write(16,*) 'A    A'

           write(16,*) 'bw    b_{w}'

           !write(16,*) 'sigma_int    \sigma_{int}'

        End If

        If (hyperparameters_as_mcmc) then
                ! WRITING PARAMNAMES FOR HYPER-PARAMETERS
                Do m=number_model_parameters+1,number_of_parameters

                    write(string,'(i2.2)') m-number_model_parameters

                    write(16,*) 'alpha_'//trim(string)//'    \alpha_{'//trim(string)//'}'

                    alpha_string(m-number_model_parameters) = 'alpha_'//trim(string)//'    '

                End Do
            
        End If

        close(16)

        If (doing_R11_analysis) then

           If ( ( use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. use_MW_as_anchor) then
    
              print *,'USE OF THREE ANCHORS SIMULTANEOUSLY NOT IMPLEMENTED YET'

              stop

           Else If ( ( use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then

              print *,'NGC4258+LMC NOT IMPLEMENTED YET'

              stop

           Else If ( ( use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. use_MW_as_anchor ) then

              print *,'NGC4258+MW NOT IMPLEMENTED YET'

              stop

           Else If ( ( .not.use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. use_MW_as_anchor ) then

              print *,'MW+LMC NOT IMPLEMENTED YET'

              stop

           Else If ( ( .not.use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. use_MW_as_anchor ) then

              print *,'MW NOT IMPLEMENTED YET'

              stop

           Else If ( ( .not.use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then

              print *,'LMC NOT IMPLEMENTED YET'

              stop

           Else If ( ( use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then

              If (use_metallicity) then 

                 If (use_H_band) then

                    print *,'H BAND NOT IMPLEMENTED INCLUDING METALLICITY DEPENDENCE'
                     
                    stop
                    
                 Else

                    write(17,*) 'mu01    0.    50.'

                    write(17,*) 'mu02    0.    50.'

                    write(17,*) 'mu03    0.    50.'

                    write(17,*) 'mu04    0.    50.'

                    write(17,*) 'mu05    0.    50.'

                    write(17,*) 'mu06    0.    50.'

                    write(17,*) 'mu07    0.    50.'

                    write(17,*) 'mu08    0.    50.'

                    write(17,*) 'mu04258    0.    50.'

                    write(17,*) 'zpw4258    0.    50.'

                    write(17,*) 'bw    -20.    0.'

                    write(17,*) 'm0v4258    0.    14.'

                    write(17,*) 'Zw    -2.    2.'

                 End If

              Else

                 If (use_H_band) then

                    write(17,*) 'zpH1    0.    50.'

                    write(17,*) 'zpH2    0.    50.'

                    write(17,*) 'zpH3    0.    50.'

                    write(17,*) 'zpH4    0.    50.'

                    write(17,*) 'zpH5    0.    50.'

                    write(17,*) 'zpH6    0.    50.'

                    write(17,*) 'zpH7    0.    50.'

                    write(17,*) 'zpH8    0.    50.'

                    write(17,*) 'zpH4258    0.    50.'

                    write(17,*) 'bH    -20.    0.'

                 Else

                    print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE'
                     
                    stop
                        
                 End If

              End If

           Else

              print *, 'USER MUST SET TRUE AT LEAST ONE ANCHOR DISTANCE IN FIDUCIAL MODULE'

              stop

           End If
    
        Else

           write(17,*) 'A    0.    50. '

           write(17,*) 'bw    -20.    0. '

           !    write(17,*) 'sigma_int    1.e-10    1 '

        End If

        If (hyperparameters_as_mcmc) then
                ! WRITING HARD BOUNDS FOR HYPER-PARAMETERS
                Do m=number_model_parameters+1,number_of_parameters

                    write(string,'(i2.2)') m-number_model_parameters
                    
                    If (using_jeffreys_prior) then

                        print *,'MODIFY THIS PART ACCORDING TO PRIOR'

                        stop

                    Else

                       write(17,*) 'alpha_'//trim(string)//'    0.    1.'

                    End If

                End Do
            
        End If

        close(17)

    End If

    ! OPEN TEMPORARY FILE TO SAVE CHAIN
    open(14,file='./output/mcmc_output.txt')     

    write(13,*) '# NUMBER OF ITERATIONS IN MCMC : ', number_iterations - steps_taken_before_definite_run

    If (start_from_fiducial .and. (.not.testing_Gaussian_likelihood .and. .not.doing_R11_analysis)) then

        write(15,*) '# FIDUCIAL MODEL IS (PARAMETERS ARE ORDERED AS IN CHAINS FILES) :', prior_A, prior_bw, prior_sigma_int

        write(15,'(a37,es18.10)') '# ln(L/L_max) AT THE FIDUCIAL MODEL :', old_loglikelihood

    End If

    If (hyperparameters_as_mcmc) then

        write(13,*) '# WEIGHT   -ln(L/L_{max})    A    bw   ', alpha_string(1:number_hyperparameters)
 
    Else

       If (doing_R11_analysis) then

          If ( ( use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. use_MW_as_anchor) then
    
             print *,'USE OF THREE ANCHORS SIMULTANEOUSLY NOT IMPLEMENTED YET'

             stop

          Else If ( ( use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then

             print *,'NGC4258+LMC NOT IMPLEMENTED YET'

             stop

          Else If ( ( use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. use_MW_as_anchor ) then

             print *,'NGC4258+MW NOT IMPLEMENTED YET'

             stop

          Else If ( ( .not.use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. use_MW_as_anchor ) then

             print *,'MW+LMC NOT IMPLEMENTED YET'

             stop

          Else If ( ( .not.use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. use_MW_as_anchor ) then

             print *,'MW NOT IMPLEMENTED YET'

             stop

          Else If ( ( .not.use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then

             print *,'LMC NOT IMPLEMENTED YET'

             stop

          Else If ( ( use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then

             If (use_metallicity) then 

                If (use_H_band) then

                   print *,'H BAND NOT IMPLEMENTED INCLUDING METALLICITY DEPENDENCE'
                     
                   stop
                   
                Else

                   write(13,*) '# WEIGHT   -ln(L/L_{max})    mu01    mu02    mu03'//trim(&
                   '    mu04    mu05    mu06    mu07    mu08    mu04258    zpw4258    bw    m0v4258    Zw')//'' 

                End If

             Else

                If (use_H_band) then

                   write(13,*) '# WEIGHT   -ln(L/L_{max})    zpH1    zpH2    zpH3'//trim(&
                   '    zpH4    zpH5    zpH6    zpH7    zpH8    zpH4258    bH')//'' 

                Else

                   print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE'
                     
                   stop
                        
                End If

             End If

          Else

             print *, 'USER MUST SET TRUE AT LEAST ONE ANCHOR DISTANCE IN FIDUCIAL MODULE'

             stop

          End If
    
       Else

          !write(13,*) '# Weight   -ln(L/L_{max})    A    bw    sigma_int ' 
          write(13,*) '# WEIGHT   -ln(L/L_{max})    A    bw ' 

       End If

    End If

    write(15,*)'STARTING SAMPLING OF PARAMETER SPACE'

    ! LOOP TO EXPLORE PARAMETER SPACE STARTS HERE
    Do m=1,number_iterations

        ! GENERATE NEW POINT IN PARAMETER SPACE FROM MULTIVARIATE DISTRIBUTION; CODE USES RANLIB LIBRARY (BE CAREFUL WITH X_OLD AND OLD_POINT DEFINITIONS)
        Do q=1,number_iterations 

            If (testing_Gaussian_likelihood) then

                call setgmn(x_old,real(Covgauss),number_of_parameters,parm) 
 
                call genmn(parm,x_new,work)

                exit

            Else

                call setgmn(x_old,real(Covguess),number_of_parameters,parm) 

                call genmn(parm,x_new,work)

            End If

            If (doing_R11_analysis) then

               If ( ( use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. use_MW_as_anchor) then
    
                  print *,'USE OF THREE ANCHORS SIMULTANEOUSLY NOT IMPLEMENTED YET'

                  stop

               Else If ( ( use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then

                  print *,'NGC4258+LMC NOT IMPLEMENTED YET'

                  stop

               Else If ( ( use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. use_MW_as_anchor ) then

                  print *,'NGC4258+MW NOT IMPLEMENTED YET'

                  stop

               Else If ( ( .not.use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. use_MW_as_anchor ) then

                  print *,'MW+LMC NOT IMPLEMENTED YET'

                  stop

               Else If ( ( .not.use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. use_MW_as_anchor ) then

                  print *,'MW NOT IMPLEMENTED YET'

                  stop

               Else If ( ( .not.use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then

                  print *,'LMC NOT IMPLEMENTED YET'

                  stop

               Else If ( ( use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then

                  If (use_metallicity) then 

                     If (use_H_band) then

                        print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE'
                     
                        stop

                     Else

                        plausibility(1) = (x_new(1) .le. real(0.d0)) .or. (x_new(1) .ge. real(5.d1))

                        plausibility(2) = (x_new(2) .le. real(0.d0)) .or. (x_new(2) .ge. real(5.d1))

                        plausibility(3) = (x_new(3) .le. real(0.d0)) .or. (x_new(3) .ge. real(5.d1))

                        plausibility(4) = (x_new(4) .le. real(0.d0)) .or. (x_new(4) .ge. real(5.d1))

                        plausibility(5) = (x_new(5) .le. real(0.d0)) .or. (x_new(5) .ge. real(5.d1))

                        plausibility(6) = (x_new(6) .le. real(0.d0)) .or. (x_new(6) .ge. real(5.d1))

                        plausibility(7) = (x_new(7) .le. real(0.d0)) .or. (x_new(7) .ge. real(5.d1))

                        plausibility(8) = (x_new(8) .le. real(0.d0)) .or. (x_new(8) .ge. real(5.d1))

                        plausibility(9) = (x_new(9) .le. real(0.d0)) .or. (x_new(9) .ge. real(50.d0))

                        plausibility(10) =  (x_new(10) .le. real(0.d0)) .or. (x_new(10) .ge. real(50.d0)) 

                        plausibility(11) =  (x_new(11) .le. real(-20.d0)) .or. (x_new(11) .ge. real(0.d0)) 

                        plausibility(12) =  (x_new(12) .le. real(0.d0)) .or. (x_new(12) .ge. real(14.d0)) 

                        plausibility(13) =  (x_new(13) .le. real(-2.d0)) .or. (x_new(13) .ge. real(2.d0)) 

                     End If

                  Else

                     If (use_H_band) then

                        plausibility(1) = (x_new(1) .le. real(0.d0)) .or. (x_new(1) .ge. real(5.d1))

                        plausibility(2) = (x_new(2) .le. real(0.d0)) .or. (x_new(2) .ge. real(5.d1))

                        plausibility(3) = (x_new(3) .le. real(0.d0)) .or. (x_new(3) .ge. real(5.d1))

                        plausibility(4) = (x_new(4) .le. real(0.d0)) .or. (x_new(4) .ge. real(5.d1))

                        plausibility(5) = (x_new(5) .le. real(0.d0)) .or. (x_new(5) .ge. real(5.d1))

                        plausibility(6) = (x_new(6) .le. real(0.d0)) .or. (x_new(6) .ge. real(5.d1))

                        plausibility(7) = (x_new(7) .le. real(0.d0)) .or. (x_new(7) .ge. real(5.d1))

                        plausibility(8) = (x_new(8) .le. real(0.d0)) .or. (x_new(8) .ge. real(5.d1))

                        plausibility(9) = (x_new(9) .le. real(0.d0)) .or. (x_new(9) .ge. real(50.d0))

                        plausibility(10) =  (x_new(10) .le. real(-20.d0)) .or. (x_new(10) .ge. real(0.d0)) 

                     Else

                        print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE'
                     
                        stop
                        
                     End If

                  End If

               Else

                  print *, 'USER MUST SET TRUE AT LEAST ONE ANCHOR DISTANCE IN FIDUCIAL MODULE'

                  stop

               End If
    
            Else

               plausibility(1) = (x_new(1) .le. real(0.d0)) .or. (x_new(1) .ge. real(5.d1))
               plausibility(2) = (x_new(2) .le. real(-2.d1)) .or. (x_new(2) .ge. real(0.d0))
               !plausibility(3) =  (x_new(3) .gt. real(0.d0)) .or. (x_new(3) .lt. real(-10.d0))    ! limit log10(sigma_int)

            End If

            If (hyperparameters_as_mcmc) then
                ! CHECKING PLAUSIBILITY FOR HYPER-PARAMETERS
                Do n=number_model_parameters+1,number_of_parameters
                   
                    If (x_new(n) .gt. real(1.d0)) then
                     
                        x_new(n) = x_old(n)

                    Else If (x_new(n) .lt. real(0.d0)) then

                        x_new(n) = x_old(n)

                    End If

                    plausibility(n) =  (x_new(n) .gt. real(1.d0)) .or. (x_new(n) .lt. real(0.d0))    ! limit alpha_j

                End Do
            
            End If

            Do n=1,number_of_parameters

                If (plausibility(n)) then

                   x_new(n) = x_old(n)
!                    non_plausible_parameters = .true.

!                    exit

!                Else if (n .eq. number_of_parameters) then

                    non_plausible_parameters = .false.

                End If

            End Do

            If (non_plausible_parameters .and. (q .ne. number_iterations)) then

                call genmn(parm,x_new,work)

            Else If (q .eq. number_iterations) then

                write(15,*) 'LOOP TO GENERATE MULTIVARIATE GAUSSIAN DEVIATE HIT MAXIMUM NUMBER OF'
                write(15,*) 'ITERATIONS WITHOUT FINDING AN ALLOWED POINT'

                stop

            Else 

                exit

            End If

        End Do
        ! NEW POINT GENERATED 

        Do n=1,number_of_parameters

            If (n .gt. number_model_parameters) then

                If (using_jeffreys_prior) then

                    current_point(n) = 10**(dble(x_new(n))) ! CONVERTING LOG10(alpha_j) to alpha_j 

                Else

                    current_point(n) = dble(x_new(n))
              
                End If

            Else

                current_point(n) = dble(x_new(n))

            End If

        End Do

        ! EVALUATE LOG_LIKELIHOOD FOR CURRENT POINT IN PARAMETER SPACE
        If (testing_Gaussian_likelihood) then

            current_loglikelihood = log_Gaussian_likelihood(current_point)

        Else

            If (using_hyperparameters) then    

               If (doing_R11_analysis) then

                  If (use_H_band) then

                     current_loglikelihood = log_R11_likelihood_H(current_point(1:number_model_parameters-1),&
                          current_point(number_model_parameters),prior_sigma_int)

                  Else

                     current_loglikelihood = log_R11_likelihood_W(current_point(1:number_model_parameters-4),&
                          current_point(number_model_parameters-3),current_point(number_model_parameters-2),&
                          current_point(number_model_parameters-1),current_point(number_model_parameters),prior_sigma_int)

                  End If

               Else

                  current_loglikelihood = log_Efstathiou_likelihood_hyperparameters(current_point(1),&
                       current_point(2),prior_sigma_int)

               End If

            Else

                current_loglikelihood = log_Efstathiou_likelihood(current_point(1),current_point(2),prior_sigma_int)

            End If

        End If
        ! LOG_LIKELIHOOD FOR CURRENT POINT COMPUTED

        !MAKE DECISION ABOUT CURRENT POINT : ACCEPT OR REJECT IT
        If (current_loglikelihood .ge. old_loglikelihood) then ! ACCEPT CURRENT POINT

            If (m .gt. steps_taken_before_definite_run) then

                number_accepted_points = number_accepted_points + 1         

            End If

            ! COMPUTING ACCEPTANCE PROBABILITY FOR CURRENT POINT
            acceptance_probability(m) = min(1.d0,exp(current_loglikelihood - old_loglikelihood))    
        
            If (m .le. steps_taken_before_definite_run) then ! WRITE OUT INFORMATION IN TEMPORARY FILE
               
               If (doing_R11_analysis) then

                  write(14,*) weight,-old_loglikelihood,old_point(1:number_of_parameters)

               Else

                  write(14,*) weight,-old_loglikelihood,old_point(1:number_of_parameters)

               End If

            Else ! WRITE OUT INFORMATION IN DEFINITE FILE

               If (doing_R11_analysis) then

                  write(13,*) weight,-old_loglikelihood,old_point(1:number_of_parameters)


               Else

                  write(13,*) weight,-old_loglikelihood,old_point(1:number_of_parameters)

               End If

            End If
       
            weight = 1    

            old_loglikelihood = current_loglikelihood
        
            Do i=1,number_of_parameters 

                old_point(i) = current_point(i)

                If (i .gt. number_model_parameters) then

                    If (using_jeffreys_prior) then

                        x_old(i) = log10(real(old_point(i))) ! CONVERTING alpha_j TO  log10(alpha_j)

                    Else

                        x_old(i) = real(old_point(i))

                    End If

                Else

                    x_old(i) = real(old_point(i))

                End If

            End Do
   
        Else ! ACCEPT OR REJECT THE CURRENT POINT ACCORDING TO :

            random_uniform = dble(genunf(real(0.),real(1.)))

            If ( random_uniform .le. exp(current_loglikelihood-old_loglikelihood)) then ! ACCEPT CURRENT POINT

                If (m .gt. steps_taken_before_definite_run) then

                    number_accepted_points = number_accepted_points + 1         

                End If
                
                acceptance_probability(m) = min(1.d0,dexp(current_loglikelihood - old_loglikelihood))    

                If (m .le. steps_taken_before_definite_run) then ! WRITE OUT INFORMATION TO TEMPORARY FILE

                   If (doing_R11_analysis) then

                      write(14,*) weight,-old_loglikelihood,old_point(1:number_of_parameters)

                   Else

                      write(14,*) weight,-old_loglikelihood,old_point(1:number_of_parameters)

                   End If

                Else ! WRITE OUT INFORMATION TO DEFINITE FILE
                   
                   If (doing_R11_analysis) then

                      write(13,*) weight,-old_loglikelihood,old_point(1:number_of_parameters)

                   Else

                      write(13,*) weight,-old_loglikelihood,old_point(1:number_of_parameters)

                   End If

                End If

                weight = 1

                old_loglikelihood = current_loglikelihood

                Do i=1,number_of_parameters 

                    old_point(i) = current_point(i)

                    If (i .gt. number_model_parameters) then
                        
                        If (using_jeffreys_prior) then

                            x_old(i) = real(log10(old_point(i))) ! CONVERTING alpha_j TO log10(alpha_j)

                        Else

                            x_old(i) = real(old_point(i))

                        End If

                    Else

                        x_old(i) = real(old_point(i))

                    End If

                End Do

            Else   ! REJECT CURRENT POINT 

                If (m .gt. steps_taken_before_definite_run) then

                    number_rejected_points = number_rejected_points + 1            

                End If

                acceptance_probability(m) = min(1.d0,exp(current_loglikelihood - old_loglikelihood))    

                weight = weight + 1

                Do i=1,number_of_parameters 

                    If (i .gt. number_model_parameters) then

                        If (using_jeffreys_prior) then

                            x_old(i) = real(log10(old_point(i))) ! CONVERT alpha_j TO log10(alpha_j)

                        Else

                            x_old(i) = real(old_point(i))

                        End If

                    Else

                        x_old(i) = real(old_point(i))

                    End If

                End Do

            End If

        End If
        ! DECISION ABOUT CURRENT POINT MADE

        ! COMPUTING AVERAGE ACCEPTANCE PROBABILITY AND UPDATING BOTH COVARIANCE MATRIX AND JUMPING FACTOR (IF NEEDED)
        If ((mod(m,jumping_factor_update) .eq. 0) .and. (m .le. steps_taken_before_definite_run) ) then

            average_acceptance_probability = sum(acceptance_probability(m-jumping_factor_update+1:m))&
            /real(jumping_factor_update)

!            write(15,*) 'CURRENT AVERAGE ACCEPTANCE PROBABILITY = ',average_acceptance_probability
            
            ! UPDATE JUMPING FACTOR IF NEEDED        
            If (average_acceptance_probability .lt. 0.1) then 
               
                jumping_factor = (1.d0 - step_size_changes)    !    Decreasing step size

                If (testing_Gaussian_likelihood) then

                    Covgauss = jumping_factor*Covgauss

                Else

                    Covguess = jumping_factor*Covguess

                End If

            Else if (average_acceptance_probability .gt. 0.4) then

                jumping_factor = (1.d0 + step_size_changes)    !    Increasing step size 

                If (testing_Gaussian_likelihood) then

                    Covgauss = jumping_factor*Covgauss

                Else

                    Covguess = jumping_factor*Covguess

                End If

            End If
            ! JUMPING FACTOR UPDATED (IF IT WAS NEEDED)
             
            not_good_aap = (average_acceptance_probability .lt. 0.1) .or. (average_acceptance_probability .gt. 0.4)

            If ( (mod(m,covariance_matrix_update) .eq. 0) .and. not_good_aap) then
                
                call stat('./output/mcmc_output.txt',buff,status1)

                If ((status1 .eq. 0) .and. (buff(8) .gt. 0)) then
              
                    If (testing_Gaussian_likelihood) then

                        call system('cd output; python compute_covariance_matrix_Gaussian.py')

                        call read_covariance_matrix_mcmc(Covgauss)

                        close(14)

                        call system('rm ./output/mcmc_output.txt')

                        open(14,file='./output/mcmc_output.txt')

                    Else

                        If (using_hyperparameters) then

                            If (hyperparameters_as_mcmc) then
                                
                                call system('cd output; python compute_covariance_matrix_HP.py')
                                
                            Else

                                call system('cd output; python compute_covariance_matrix.py')

                            End If

                        Else

                            call system('cd output; python compute_covariance_matrix.py')

                        End If

!                        write(15,*) 'CURRENT COVARIANCE MATRIX BEFORE',Covguess

                        call read_covariance_matrix_mcmc(Covguess)

!                        write(15,*) 'CURRENT COVARIANCE MATRIX AFTER',Covguess

                        close(14)

                        call system('rm ./output/mcmc_output.txt')

                        open(14,file='./output/mcmc_output.txt')

                    End If

                End If

            End If

        End If

    End Do
    ! LOOP TO EXPLORE PARAMETER SPACE ENDED

    write(15,*) 'NUMBER OF REJECTED POINTS = ', number_rejected_points

    write(15,*) 'ACCEPTANCE RATIO = ', dble(number_iterations - steps_taken_before_definite_run - number_rejected_points)/&
    dble(number_iterations - steps_taken_before_definite_run)
 
    ! CLOSE FILE STORING CHAIN
    close(13)
    ! CLOSE TEMPORARY FILE FOR CHAINS
    close(14)

    !ANALYZE SAMPLES, MAKE FIGURES, COMPUTE BESTFIT AND HYPER-PARAMETERS (IF NEEDED)
    If (testing_Gaussian_likelihood) then

        call system('cd analyzer; python analyze.py')

    Else

        If (using_hyperparameters) then

            If (hyperparameters_as_mcmc) then

                call system('cd analyzer; python analyze_HP_as_MCMC.py')

            Else
               
               If (doing_R11_analysis) then  !MUST IMPLEMENT OTHER OPTIONS LATER!!!!!!!!!!!!!!!!!

                  If (use_H_band) then

                     call system('cd analyzer; python analyze_HP_R11_H.py')

                  Else

                     call system('cd analyzer; python analyze_HP_R11_W.py')

                  End If

               Else

                  call system('cd analyzer; python analyze_HP.py')

               End If

            End If

        Else

            call system('cd analyzer; python analyze.py')

        End If    

    End If

    call read_bestfit_mcmc(bestfit)

    call read_means_mcmc(means)

    If (doing_R11_analysis) then

       If ( ( use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. use_MW_as_anchor) then
    
          print *,'USE OF THREE ANCHORS SIMULTANEOUSLY NOT IMPLEMENTED YET'

          stop

       Else If ( ( use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then

          print *,'NGC4258+LMC NOT IMPLEMENTED YET'

          stop

       Else If ( ( use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. use_MW_as_anchor ) then

          print *,'NGC4258+MW NOT IMPLEMENTED YET'

          stop
          
       Else If ( ( .not.use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. use_MW_as_anchor ) then

          print *,'MW+LMC NOT IMPLEMENTED YET'

          stop

       Else If ( ( .not.use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. use_MW_as_anchor ) then

          print *,'MW NOT IMPLEMENTED YET'

          stop

       Else If ( ( .not.use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then

          print *,'LMC NOT IMPLEMENTED YET'

          stop

       Else If ( ( use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then

          If (use_metallicity) then 

             If (use_H_band) then

                print *,'H BAND NOT IMPLEMENTED INCLUDING METALLICITY DEPENDENCE'
                     
                stop

             Else

                write(15,*) 'BESTFIT IS : '

                write(15,*) 'mu01 = ', bestfit(1)

                write(15,*) 'mu02 = ', bestfit(2)

                write(15,*) 'mu03 = ', bestfit(3)

                write(15,*) 'mu04 = ', bestfit(4)

                write(15,*) 'mu05 = ', bestfit(5)

                write(15,*) 'mu06 = ', bestfit(6)

                write(15,*) 'mu07 = ', bestfit(7)

                write(15,*) 'mu08 = ', bestfit(8)

                write(15,*) 'mu04258 = ', bestfit(9)

                write(15,*) 'zpw4258 = ', bestfit(10)

                write(15,*) 'bw = ', bestfit(11)

                write(15,*) 'm0v4258 = ', bestfit(12)

                write(15,*) 'Zw = ', bestfit(13)

             End If

          Else

             If (use_H_band) then

                write(15,*) 'BESTFIT IS : '

                write(15,*) 'bH = ', bestfit(number_of_parameters)

             Else

                print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE'
                     
                stop
                        
             End If

          End If

       Else

          print *, 'USER MUST SET TRUE AT LEAST ONE ANCHOR DISTANCE IN FIDUCIAL MODULE'

          stop

       End If
    
    Else

       write(15,*) 'BESTFIT IS : '

       write(15,*) 'A = ', bestfit(1)

       write(15,*) 'bw = ', bestfit(2)

       !write(15,*) 'sigma_int = ', bestfit(3)

    End If


    If (hyperparameters_as_mcmc .and. using_hyperparameters) then
    ! WRITING BESTFIT FOR HYPER-PARAMETERS
        Do m=number_model_parameters+1,number_of_parameters

            write(string,'(i2)') m-number_model_parameters

            write(15,*) 'alpha_'//trim(string)//' = ', bestfit(m)

        End Do
            
    End If

    If (doing_R11_analysis) then

       If ( ( use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. use_MW_as_anchor) then
    
          print *,'USE OF THREE ANCHORS SIMULTANEOUSLY NOT IMPLEMENTED YET'

          stop

       Else If ( ( use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then

          print *,'NGC4258+LMC NOT IMPLEMENTED YET'

          stop

       Else If ( ( use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. use_MW_as_anchor ) then

          print *,'NGC4258+MW NOT IMPLEMENTED YET'

          stop
          
       Else If ( ( .not.use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. use_MW_as_anchor ) then

          print *,'MW+LMC NOT IMPLEMENTED YET'

          stop

       Else If ( ( .not.use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. use_MW_as_anchor ) then

          print *,'MW NOT IMPLEMENTED YET'

          stop

       Else If ( ( .not.use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then

          print *,'LMC NOT IMPLEMENTED YET'

          stop

       Else If ( ( use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then

          If (use_metallicity) then 

             If (use_H_band) then

                print *,'H BAND NOT IMPLEMENTED INCLUDING METALLICITY DEPENDENCE'
                     
                stop
                
             Else

                write(15,*) 'MEANS FOR THE SAMPLES ARE : '

                write(15,*) 'mu01 = ', means(1)

                write(15,*) 'mu02 = ', means(2)

                write(15,*) 'mu03 = ', means(3)

                write(15,*) 'mu04 = ', means(4)

                write(15,*) 'mu05 = ', means(5)

                write(15,*) 'mu06 = ', means(6)

                write(15,*) 'mu07 = ', means(7)

                write(15,*) 'mu08 = ', means(8)

                write(15,*) 'mu04258 = ', means(9)

                write(15,*) 'zpw4258 = ', means(10)

                write(15,*) 'bw = ', means(11)

                write(15,*) 'm0v4258 = ', means(12)

                write(15,*) 'Zw = ', means(13)

             End If

          Else

             If (use_H_band) then

                write(15,*) 'MEANS FOR THE SAMPLES ARE : '

                write(15,*) 'bH = ', means(number_of_parameters)

             Else

                print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE'
                     
                stop
                        
             End If

          End If

       Else

          print *, 'USER MUST SET TRUE AT LEAST ONE ANCHOR DISTANCE IN FIDUCIAL MODULE'

          stop

       End If
    
    Else

       write(15,*) 'MEANS FOR THE SAMPLES ARE : '

       write(15,*) 'A = ', means(1)

       write(15,*) 'bw = ', means(2)

       !write(15,*) 'sigma_int = ', means(3)

    End If


    If (hyperparameters_as_mcmc .and. using_hyperparameters) then
    ! WRITING SAMPLE MEANS FOR HYPER-PARAMETERS
        Do m=number_model_parameters+1,number_of_parameters

            write(string,'(i2)') m-number_model_parameters

            write(15,*) 'alpha_'//trim(string)//' = ', means(m)

        End Do
            
    End If

    If (using_hyperparameters .and. .not.hyperparameters_as_mcmc) then

        write(15,*) 'Hyperparameters are :'

        If (doing_R11_analysis) then

           print *,'MUST IMPLEMENT EFFECTIVE HYPERPARAMETERS FOR R11 ANALYSIS'

        Else

           If (separate_dataA .and. include_dataA) then

              Do m=1,size(NameA)

                 If (using_jeffreys_prior) then

                    write(15,*) 'Point ', m,' in data set A = ', 1.d0/chi2A_i(bestfit(1),bestfit(2),prior_sigma_int,m)
                    
                 Else

                    If (chi2A_i(bestfit(1),bestfit(2),prior_sigma_int,m) .le. 1.d0 ) then

                       write(15,*) 'Point ', m,' in data set A = ', 1.d0

                    Else

                       write(15,*) 'Point ', m,' in data set A = ', 1.d0/chi2A_i(bestfit(1),bestfit(2),prior_sigma_int,m)

                    End If

                 End If

              End Do

           Else if (include_dataA .and. .not.separate_dataA) then

              !    write(15,*) 'For data set A = ', dble(size(NameA))/chi2A(bestfit(1),bestfit(2),bestfit(3))
              write(15,*) 'For data set A = ', dble(size(NameA))/chi2A(bestfit(1),bestfit(2),prior_sigma_int)

           End If

           If (separate_dataB .and. include_dataB) then

              Do m=1,size(NameB)

                 If (using_jeffreys_prior) then

                    write(15,*) 'Point ', m,' in data set B = ', 1.d0/chi2B_i(bestfit(1),bestfit(2),prior_sigma_int,m)

                 Else

                    If (chi2B_i(bestfit(1),bestfit(2),prior_sigma_int,m) .le. 1.d0 ) then

                       write(15,*) 'Point ', m,' in data set B = ', 1.d0

                    Else

                       write(15,*) 'Point ', m,' in data set B = ', 1.d0/chi2B_i(bestfit(1),bestfit(2),prior_sigma_int,m)

                    End If

                 End If

              End Do

           Else if (include_dataB .and. .not.separate_dataB) then

              !    write(15,*) 'For data set B = ', dble(size(NameB))/chi2B(bestfit(1),bestfit(2),bestfit(3))
              write(15,*) 'For data set B = ', dble(size(NameB))/chi2B(bestfit(1),bestfit(2),prior_sigma_int)

           End If

           If (separate_dataC .and. include_dataC) then

              Do m=1,size(NameC)

                 If (using_jeffreys_prior) then

                    write(15,*) 'Point ', m,' in data set C = ', 1.d0/chi2C_i(bestfit(1),bestfit(2),prior_sigma_int,m)

                 Else

                    If (chi2C_i(bestfit(1),bestfit(2),prior_sigma_int,m) .le. 1.d0) then

                       write(15,*) 'Point ', m,' in data set C = ', 1.d0

                    Else

                       write(15,*) 'Point ', m,' in data set C = ', 1.d0/chi2C_i(bestfit(1),bestfit(2),prior_sigma_int,m)

                    End If

                 End If

              End Do

           Else if (include_dataC .and. .not.separate_dataC) then

              !write(15,*) 'For data set C = ', dble(size(NameC))/chi2C(bestfit(1),bestfit(2),bestfit(3))
              write(15,*) 'For data set C = ', dble(size(NameC))/chi2C(bestfit(1),bestfit(2),prior_sigma_int)

           End If

           write(15,*) '\ln P(\vec{w},D) at the bestfit is ', log_Efstathiou_likelihood_hyperparameters(bestfit(1),bestfit(2),&
                prior_sigma_int)
           
        End If

     End If

     close(15)

    If ((.not. testing_Gaussian_likelihood) .and. (.not. using_hyperparameters) ) then

        deallocate (old_point,current_point,acceptance_probability,Name,Period,H,Sigma_m,V,II)

    End If

End Program mcmc




