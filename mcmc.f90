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

    Integer*4 :: m,n,i,j,q                                        ! INTERGER FOR SHORT LOOPS 
    Integer*4 :: seed1,seed2                                      ! SEEDS FOR RANDOM NUMBER GENERATOR 
    Integer*4 :: number_accepted_points,number_rejected_points    ! MCMC PARAMETERS
    Integer*4 :: weight                                           ! IT COUNTS THE NUMBER OF STEPS TAKEN BEFORE MOVING TO A NEW POINT IN MCMC 

    Real*4 :: average_acceptance_probability           ! SAVES ACCEPTANCE PROBABILITY 
    Real*4 :: genunf,gennor                            ! RANDOM UNIFOR DEVIATES 
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

        Covguess(1,1) = sigma_A**2 

        Covguess(2,2) = sigma_bw**2

        !Covguess(3,3) = sigma_sigma_int**2

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

                open(15,file='./output/execution_information_HP.txt')    ! OPEN FILE FOR EXECUTION INFORMATION
            
                write(15,*) 'WORKING WITH HYPER-PARAMETERS'     

                If (number_hyperparameters .ne. 0) then
                
                    write(15,*) 'NUMBER OF HYPER-PARAMETERS MUST BE ZERO WHEN SETTING FALSE "hyperparameters_as_mcmc" '
 
                    write(15,*) 'CHECK FIDUCIAL MODULE'
       
                    stop

                End If

            End If


        Else

            call read_data_Efstathiou(path_to_datafileAB)

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

            old_point(1) = prior_A         ! A 
            old_point(2) = prior_bw        ! bw 
            !old_point(3) = prior_sigma_int ! sigma_int 

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

            x_old(1) = genunf(real(prior_A - sigma_A),real(prior_A + sigma_A))         ! A
            x_old(2) = genunf(real(prior_bw - sigma_bw),real(prior_bw + sigma_bw)) ! bw
            !x_old(3) = genunf(real(-10.d0),real(0.d0)) ! log10(sigma_int)

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
            old_loglikelihood = log_Efstathiou_likelihood_hyperparameters(old_point(1),old_point(2),prior_sigma_int)
            
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

        write(16,*) 'A    A'

        write(16,*) 'bw    b_{w}'

        !write(16,*) 'sigma_int    \sigma_{int}'

        If (hyperparameters_as_mcmc) then
                ! WRITING PARAMNAMES FOR HYPER-PARAMETERS
                Do m=number_model_parameters+1,number_of_parameters

                    write(string,'(i2.2)') m-number_model_parameters

                    write(16,*) 'alpha_'//trim(string)//'    \alpha_{'//trim(string)//'}'

                    alpha_string(m-number_model_parameters) = 'alpha_'//trim(string)//'    '

                End Do
            
        End If

        close(16)

        write(17,*) 'A    0.    50. '

        write(17,*) 'bw    -20.    0. '

    !    write(17,*) 'sigma_int    1.e-10    1 '

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

    If (start_from_fiducial .and. (.not.testing_Gaussian_likelihood)) then

        write(15,*) '# FIDUCIAL MODEL IS (PARAMETERS ARE ORDERED AS IN CHAINS FILES) :', prior_A, prior_bw, prior_sigma_int

        write(15,'(a37,es18.10)') '# ln(L/L_max) AT THE FIDUCIAL MODEL :', old_loglikelihood

    End If

    If (hyperparameters_as_mcmc) then

        write(13,*) '# WEIGHT   -ln(L/L_{max})    A    bw   ', alpha_string(1:number_hyperparameters)
 
    Else

        !write(13,*) '# Weight   -ln(L/L_{max})    A    bw    sigma_int ' 
        write(13,*) '# WEIGHT   -ln(L/L_{max})    A    bw ' 

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

            plausibility(1) = (x_new(1) .le. real(0.d0)) .or. (x_new(1) .ge. real(5.d1))
            plausibility(2) = (x_new(2) .le. real(-2.d1)) .or. (x_new(2) .ge. real(0.d0))
            !plausibility(3) =  (x_new(3) .gt. real(0.d0)) .or. (x_new(3) .lt. real(-10.d0))    ! limit log10(sigma_int)

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

                    non_plausible_parameters = .true.

                    exit

                Else if (n .eq. number_of_parameters) then

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

                current_loglikelihood = log_Efstathiou_likelihood_hyperparameters(current_point(1),current_point(2),prior_sigma_int)

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
 
                write(14,*) weight,-old_loglikelihood,old_point(1:number_of_parameters)

            Else ! WRITE OUT INFORMATION IN DEFINITE FILE

                write(13,*) weight,-old_loglikelihood,old_point(1:number_of_parameters)

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

                    write(14,*) weight,-old_loglikelihood,old_point(1:number_of_parameters)

                Else ! WRITE OUT INFORMATION TO DEFINITE FILE

                    write(13,*) weight,-old_loglikelihood,old_point(1:number_of_parameters)

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

                call system('cd analyzer; python analyze_HP.py')

            End If

        Else

            call system('cd analyzer; python analyze.py')

        End If    

    End If

    call read_bestfit_mcmc(bestfit)

    call read_means_mcmc(means)

    write(15,*) 'BESTFIT IS : '

    write(15,*) 'A = ', bestfit(1)

    write(15,*) 'bw = ', bestfit(2)

    !write(15,*) 'sigma_int = ', bestfit(3)

    If (hyperparameters_as_mcmc .and. using_hyperparameters) then
    ! WRITING BESTFIT FOR HYPER-PARAMETERS
        Do m=number_model_parameters+1,number_of_parameters

            write(string,'(i2)') m-number_model_parameters

            write(15,*) 'alpha_'//trim(string)//' = ', bestfit(m)

        End Do
            
    End If

    write(15,*) 'MEANS FOR THE SAMPLES ARE : '

    write(15,*) 'A = ', means(1)

    write(15,*) 'bw = ', means(2)

!write(15,*) 'sigma_int = ', means(3)

    If (hyperparameters_as_mcmc .and. using_hyperparameters) then
    ! WRITING SAMPLE MEANS FOR HYPER-PARAMETERS
        Do m=number_model_parameters+1,number_of_parameters

            write(string,'(i2)') m-number_model_parameters

            write(15,*) 'alpha_'//trim(string)//' = ', means(m)

        End Do
            
    End If

    If (using_hyperparameters .and. .not.hyperparameters_as_mcmc) then

        write(15,*) 'Hyperparameters are :'

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

    close(15)

    If ((.not. testing_Gaussian_likelihood) .and. (.not. using_hyperparameters) ) then

        deallocate (old_point,current_point,acceptance_probability,Name,Period,H,Sigma_m,V,II)

    End If

End Program mcmc




