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

    Logical :: c1,c2,c4,c5,c6,c7,non_plausible_parameters   ! IT CONTROLS PLAUSIBLE VALUES OF COSMOLOGICAL PARAMETERS 

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

        open(15,file='./output/execution_information.txt')    ! OPEN FILE FOR EXECUTION INFORMATION 

    Else
        ! SETTING COVARIANCE MATRIX
        Covguess = 0.d0  

        Covguess(1,1) = sigma_A**2 

        Covguess(2,2) = sigma_bw**2

        !Covguess(3,3) = sigma_sigma_int**2
        ! COVARIANCE MATRIX SET

        jumping_factor = 2.38d0/sqrt(dble(number_of_parameters))!*1.d-3    !    MODIFY ACCORDING TO WANTED ACCEPTANCE PROBABILITY

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

            open(15,file='./output/execution_information_HP.txt')    ! OPEN FILE FOR EXECUTION INFORMATION

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
    
            Do m=1,number_of_parameters

                If (m .eq. 3) then

                    x_old(m) = real(log10(old_point(m)))    ! CONVERT TO log10(sigma_int)

                else

                    x_old(m) = real(old_point(m))

                End If

            End Do

        Else

            x_old(1) = genunf(real(prior_A - sigma_A),real(prior_A + sigma_A))         ! A
            x_old(2) = genunf(real(prior_bw - sigma_bw),real(prior_bw + sigma_bw)) ! bw
            !x_old(3) = genunf(real(-10.d0),real(0.d0)) ! log10(sigma_int)

            Do m=1,number_of_parameters

                If (m .eq. 3) then

                    old_point(m) = 10**(dble(x_old(m)))    !    CONVERT TO sigma_int

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

            !old_loglikelihood = log_Efstathiou_likelihood_hyperparameters(old_point(1),old_point(2),old_point(3))    !    Compute initial likelihood value
            old_loglikelihood = log_Efstathiou_likelihood_hyperparameters(old_point(1),old_point(2),prior_sigma_int)    !    Compute initial likelihood value

            open(16,file='./output/mcmc_final_output_HP.paramnames')    !    OPEN FILE WITH PARAMETER NAMES NEEDED BY GETDIST

            open(17,file='./output/mcmc_final_output_HP.ranges')    !    OPEN FILE WITH HARD BOUNDS NEEDED BY GETDIST

        Else
            ! OPEN FILE TO STORE MCMC COMPUTATION
            open(13,file='./output/mcmc_final_output.txt')

            write(15,*) 'NOT WORKING WITH HYPER-PARAMETERS' 

            write(15,*) 'COMPUTING LOG_LIKELIHOOD FOR INITIAL POINT'

            !old_loglikelihood = log_Efstathiou_likelihood(old_point(1),old_point(2),old_point(3))    !    Compute initial likelihood value
            old_loglikelihood = log_Efstathiou_likelihood(old_point(1),old_point(2),prior_sigma_int)    !    Compute initial likelihood value

            open(16,file='./output/mcmc_final_output.paramnames')    !    OPEN FILE WITH PARAMETER NAMES NEEDED BY GETDIST

            open(17,file='./output/mcmc_final_output.ranges')    !    OPEN FILE WITH HARD BOUNDS NEEDED BY GETDIST 

        End If

        write(16,*) 'A    A'

        write(16,*) 'bw    b_{w}'

        !write(16,*) 'sigma_int    \sigma_{int}'

        close(16)

        write(17,*) 'A    0.    50. '

        write(17,*) 'bw    -20.    0. '

    !    write(17,*) 'sigma_int    1.e-10    1 '

        close(17)

    End If

    ! OPEN TEMPORARY FILE TO SAVE CHAIN
    open(14,file='./output/mcmc_output.txt')     

    write(13,*) '# NUMBER OF ITERATIONS IN MCMC : ', number_iterations - steps_taken_before_definite_run

    If (start_from_fiducial .and. (.not.testing_Gaussian_likelihood)) then

        write(15,*) '# Fiducial model is (parameters ordered as below) :', prior_A, prior_bw, prior_sigma_int

        write(15,'(a37,es18.10)') '# ln(L/L_max) at the fiducial model :', old_loglikelihood

    End If

    !write(13,*) '# Weight   -ln(L/L_{max})    A    bw    sigma_int ' 
    write(13,*) '# WEIGHT   -ln(L/L_{max})    A    bw ' 

    ! LOOP TO EXPLORE PARAMETER SPACE STARTS HERE
Do m=1,number_iterations

    !######################################################################################################
    ! Generate new point in parameter space from a multivariate normal distribution. We use RANLIB library.  
    ! Be careful with x_old and old_point definitions 
    !######################################################################################################

    Do q=1,number_iterations 

        If (testing_Gaussian_likelihood) then

            call setgmn(x_old,real(jumping_factor*Covgauss),number_of_parameters,parm) ! used if testing Gaussian likelihood
 
            call genmn(parm,x_new,work)

            exit

        else

            call setgmn(x_old,real(jumping_factor*Covguess),number_of_parameters,parm) 

            call genmn(parm,x_new,work)

        End If

        c1 = (x_new(1) .le. real(0.d0)) .or. (x_new(1) .ge. real(5.d1))
        c2 = (x_new(2) .le. real(-2.d1)) .or. (x_new(2) .ge. real(0.d0))
!        c4 =  (x_new(3) .gt. real(0.d0)) .or. (x_new(3) .lt. real(-10.d0))    ! limit log10(sigma_int)
        !c5 = .false. !(x_new(5) .lt. real(0.d0)).or.(x_new(5).gt.real(85.d0))
        !c6 = .false. !x_new(6) .lt. real(0.d0)
        !c7 = .false. !x_new(7) .le. real(0.d0)
!        non_plausible_parameters = ((c1 .or. c2) .or. (c4 .or. c5)) .or. (c6 .or. c7) 
        non_plausible_parameters = c1 .or. c2

        If (non_plausible_parameters .and. (q .ne. number_iterations)) then

            call genmn(parm,x_new,work)

        else if (q .eq. number_iterations) then

            write(15,*) 'Loop to generate multivariate Gaussian deviate hit maximum number of iterations '

            stop

        else 

            exit

        End If

    End Do
    
    Do n=1,number_of_parameters

        If (n .eq. 3) then

            If (testing_Gaussian_likelihood) then

                current_point(n) = dble(x_new(n)) ! used if testing Gaussian likelihood

            else

                current_point(n) = 10**(dble(x_new(n))) ! Converting log10(sigma_int) to sigma_int 

            End If

        else

            current_point(n) = dble(x_new(n))

        End If

    End Do

    If (testing_Gaussian_likelihood) then

        Go to 3        ! used if testing Gaussian likelihood

    End If

    !#############################################################
    ! Evaluate log_likelihood for current point in parameter space
    !#############################################################

    If (using_hyperparameters) then    

!        current_loglikelihood = log_Efstathiou_likelihood_hyperparameters(current_point(1),current_point(2),current_point(3))
        current_loglikelihood = log_Efstathiou_likelihood_hyperparameters(current_point(1),current_point(2),prior_sigma_int)

    Else

!        current_loglikelihood = log_Efstathiou_likelihood(current_point(1),current_point(2),current_point(3))
        current_loglikelihood = log_Efstathiou_likelihood(current_point(1),current_point(2),prior_sigma_int)

    End If

    3 If (testing_Gaussian_likelihood) then 
    
          current_loglikelihood = log_Gaussian_likelihood(current_point) ! used if testing Gaussian likelihood

      End If

    !######################################################################################################
    ! Decide whether or not the current_point in parameter space becomes old_point in parameter space 
    !######################################################################################################

    If (current_loglikelihood .ge. old_loglikelihood) then ! It accepts the current point

        number_accepted_points = number_accepted_points + 1    ! Used to compute acceptance rate

        acceptance_probability(m) = min(1.d0,dexp(current_loglikelihood - old_loglikelihood))    

        If (m .le. steps_taken_before_definite_run) then
 
            write(14,*) weight,-old_loglikelihood,old_point(1),old_point(2)!,old_point(3) 

        else

            write(13,*) weight,-old_loglikelihood,old_point(1),old_point(2)!,old_point(3)

        End If
       
        weight = 1    

        old_loglikelihood = current_loglikelihood
        
        Do i=1,number_of_parameters 

            old_point(i) = current_point(i)

            If (i .eq. 3) then

                If (testing_Gaussian_likelihood) then

                    x_old(i) = real(old_point(i)) ! used if testing Gaussian likelihood

                else

                    x_old(i) = log10(real(old_point(i))) ! converting sigma_int to log10(sigma_int)

                End If 

            else

                x_old(i) = real(old_point(i))

            End If

        End Do
   
    else 

        random_uniform = dble(genunf(real(0.),real(1.)))
 
        If ( random_uniform .le. dexp(current_loglikelihood-old_loglikelihood)) then 
            ! It accetps the current point 

            number_accepted_points = number_accepted_points + 1    ! Used to compute acceptance rate

            acceptance_probability(m) = min(1.d0,dexp(current_loglikelihood - old_loglikelihood))    

            If (m .le. steps_taken_before_definite_run) then

                write(14,*) weight,-old_loglikelihood,old_point(1),old_point(2)!,old_point(3)

            else

                write(13,*) weight,-old_loglikelihood,old_point(1),old_point(2)!,old_point(3)

            End If

            weight = 1

            old_loglikelihood = current_loglikelihood

            Do i=1,number_of_parameters 

                old_point(i) = current_point(i)

                If (i .eq. 3) then

                    If (testing_Gaussian_likelihood) then

                        x_old(i) = real(old_point(i)) ! used when testing Gaussian likelihood

                    Else

                        x_old(i) = real(log10(old_point(i))) ! converting sigma_int to log10(sigma_int)

                    End If

                else

                    x_old(i) = real(old_point(i))

                End If

            End Do

        else   ! The code rejects the current point 

            If (m .gt. steps_taken_before_definite_run) then

                number_rejected_points = number_rejected_points + 1            

            End If

            acceptance_probability(m) = min(1.d0,dexp(current_loglikelihood - old_loglikelihood))    

            weight = weight + 1

            Do i=1,number_of_parameters 

                If (i .eq. 3) then

                    If (testing_Gaussian_likelihood) then

                        x_old(i) = real(old_point(i)) ! Used if testing Gaussian likelihood

                    Else

                        x_old(i) = real(log10(old_point(i))) ! convert sigma_int to log10(sigma_int)

                    End If

                else

                    x_old(i) = real(old_point(i))

                End If

            End Do

        End If

    End If

    !###################################################################################################
    ! Compute average acceptance probability and update covariance matrix and jumping factor (if needed)
    !###################################################################################################

    If ((mod(m,jumping_factor_update) .eq. 0) .and. (m .le. steps_taken_before_definite_run) ) then

        average_acceptance_probability = sum(acceptance_probability(m-jumping_factor_update+1:m))&
        /real(jumping_factor_update)

!        write(15,*) 'Current average acceptance probability ',average_acceptance_probability
        
        If (average_acceptance_probability .lt. 0.1) then 

            jumping_factor = jumping_factor*(1.d0 - step_size_changes)    !    Decreasing step size

        Else if (average_acceptance_probability .gt. 0.4) then

            jumping_factor = jumping_factor*(1.d0 + step_size_changes)    !    Increasing step size 

        End If

        If (testing_Gaussian_likelihood) then

            If ( mod(m,covariance_matrix_update) .eq. 0 ) then

                call system('cd output; python compute_covariance_matrix_Gaussian.py')
     
                close(14)

                call system('rm ./output/mcmc_output.txt')

                call read_covariance_matrix_mcmc(Covgauss)

!                write(15,*) 'Iteration ', m

!                write(15,*) 'Current covariance matrix ', Covgauss

                open(14,file='./output/mcmc_output.txt')

                jumping_factor = 2.38d0/sqrt(dble(number_of_parameters))

            End If

        Else   

            If ( mod(m,covariance_matrix_update) .eq. 0 ) then

                call system('cd output; python compute_covariance_matrix.py')
     
                close(14)

                call system('rm ./output/mcmc_output.txt')

                call read_covariance_matrix_mcmc(Covguess)

!                write(15,*) 'Iteration ', m

!                write(15,*) 'Current covariance matrix ', Covguess

                open(14,file='./output/mcmc_output.txt')

                jumping_factor = 2.38d0/sqrt(dble(number_of_parameters))

            End If

        End If

    End If

    !#########################################
    ! Loop to sample parameter space ends here
    !#########################################

End Do

!#######################################################
! Write last informations and close not needed data file 
!#######################################################

write(15,*) 'Number of rejected points: ', number_rejected_points

write(15,*) 'Acceptance ratio ', dble(number_iterations - steps_taken_before_definite_run - number_rejected_points)/&
dble(number_iterations - steps_taken_before_definite_run)

close(13)

close(14)

!###############################################################################
! Analyze samples, make figures, compute bestfit and hyperparameters (if needed)
!###############################################################################

If (using_hyperparameters) then

    call system('cd analyzer; python analyze_HP.py')

Else

    call system('cd analyzer; python analyze.py')

End If

call read_bestfit_mcmc(bestfit)

call read_means_mcmc(means)

write(15,*) 'Bestfit is :'

write(15,*) 'A = ', bestfit(1)

write(15,*) 'bw = ', bestfit(2)

!write(15,*) 'sigma_int = ', bestfit(3)

write(15,*) 'Means for the samples are :'

write(15,*) 'A = ', means(1)

write(15,*) 'bw = ', means(2)

!write(15,*) 'sigma_int = ', means(3)

If (using_hyperparameters) then

    write(15,*) 'Hyperparameters are :'

    If (separate_dataA .and. include_dataA) then

        Do m=1,size(NameA)

            If (using_jeffreys_prior) then

                write(15,*) 'Point ', m,' in data set A = ', 1.d0/chi2A_i(bestfit(1),bestfit(2),prior_sigma_int,m)

            Else

                If (chi2A_i(bestfit(1),bestfit(2),prior_sigma_int,m) .lt. 1.d0 ) then

                    write(15,*) 'Point ', m,' in data set A = ', 0.d0

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

                If (chi2B_i(bestfit(1),bestfit(2),prior_sigma_int,m) .lt. 1.d0 ) then

                    write(15,*) 'Point ', m,' in data set B = ', 0.d0

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

                If (chi2C_i(bestfit(1),bestfit(2),prior_sigma_int,m) .lt. 1.d0) then

                    write(15,*) 'Point ', m,' in data set C = ', 0.d0

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




