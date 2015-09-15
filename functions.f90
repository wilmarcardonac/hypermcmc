module functions

    Implicit none
     
contains

subroutine read_data_Efstathiou(path_to_datafile)
    use arrays
    Implicit none
    Integer*4 :: arrays_dimension,p
    Integer :: stat
    character(len=*) :: path_to_datafile

    open(11,file=path_to_datafile)

    arrays_dimension = 0

    Do 

        read(11,*,iostat=stat)

        If (stat .ne. 0) then

            exit

        Else

            arrays_dimension = arrays_dimension + 1 

        End If

    End Do

    close(11)

    allocate (Name(1:arrays_dimension),Period(1:arrays_dimension),H(1:arrays_dimension),&
    Sigma_m(1:arrays_dimension),V(1:arrays_dimension),II(1:arrays_dimension),stat=status1)

    open(11,file=path_to_datafile)

    Do p=1,arrays_dimension

        read(11,*) Name(p),Period(p),H(p),Sigma_m(p),V(p),II(p)

    End Do

    close(11)

end subroutine read_data_Efstathiou

subroutine read_data_EfstathiouA(path_to_datafile)
    use arrays
    Implicit none
    Integer*4 :: arrays_dimension,p
    Integer :: stat
    character(len=*) :: path_to_datafile

    open(11,file=path_to_datafile)

    arrays_dimension = 0

    Do 

        read(11,*,iostat=stat)

        If (stat .ne. 0) then

            exit

        Else

            arrays_dimension = arrays_dimension + 1 

        End If

    End Do

    close(11)

    allocate (NameA(1:arrays_dimension),PeriodA(1:arrays_dimension),HA(1:arrays_dimension),&
    Sigma_mA(1:arrays_dimension),VA(1:arrays_dimension),IIA(1:arrays_dimension),stat=status1)

    open(11,file=path_to_datafile)

    Do p=1,arrays_dimension

        read(11,*) NameA(p),PeriodA(p),HA(p),Sigma_mA(p),VA(p),IIA(p)

    End Do

    close(11)

end subroutine read_data_EfstathiouA

subroutine read_data_EfstathiouB(path_to_datafile)
    use arrays
    Implicit none
    Integer*4 :: arrays_dimension,p
    Integer :: stat
    character(len=*) :: path_to_datafile

    open(11,file=path_to_datafile)

    arrays_dimension = 0

    Do 

        read(11,*,iostat=stat)

        If (stat .ne. 0) then

            exit

        Else

            arrays_dimension = arrays_dimension + 1 

        End If

    End Do

    close(11)

    allocate (NameB(1:arrays_dimension),PeriodB(1:arrays_dimension),HB(1:arrays_dimension),&
    Sigma_mB(1:arrays_dimension),VB(1:arrays_dimension),IIB(1:arrays_dimension),stat=status1)

    open(11,file=path_to_datafile)

    Do p=1,arrays_dimension

        read(11,*) NameB(p),PeriodB(p),HB(p),Sigma_mB(p),VB(p),IIB(p)

    End Do

    close(11)

end subroutine read_data_EfstathiouB

subroutine read_data_EfstathiouC(path_to_datafile)
    use arrays
    Implicit none
    Integer*4 :: arrays_dimension,p
    Integer :: stat
    character(len=*) :: path_to_datafile

    open(11,file=path_to_datafile)

    arrays_dimension = 0

    Do 

        read(11,*,iostat=stat)

        If (stat .ne. 0) then

            exit

        Else

            arrays_dimension = arrays_dimension + 1 

        End If

    End Do

    close(11)

    allocate (NameC(1:arrays_dimension),PeriodC(1:arrays_dimension),HC(1:arrays_dimension),&
    Sigma_mC(1:arrays_dimension),VC(1:arrays_dimension),IIC(1:arrays_dimension),stat=status1)

    open(11,file=path_to_datafile)

    Do p=1,arrays_dimension

        read(11,*) NameC(p),PeriodC(p),HC(p),Sigma_mC(p),VC(p),IIC(p)

    End Do

    close(11)

end subroutine read_data_EfstathiouC

function observed_wesenheit_magnitude(H_band_magnitude,V_magnitude,I_magnitude)    !    It computes Wesenheit magnitudes: equation (1) in 
    Implicit none    !    published version of 1311.3461.
    Real*8 :: observed_wesenheit_magnitude,H_band_magnitude,V_magnitude,I_magnitude

    observed_wesenheit_magnitude =  H_band_magnitude - 0.41d0*(V_magnitude - I_magnitude)

end function observed_wesenheit_magnitude

function wesenheit_magnitude(A,bw,P)    !    It computes equation (2) in published version of 1311.3461
    Implicit none
    Real*8 :: wesenheit_magnitude,A,bw,P

    wesenheit_magnitude = A + bw*(log10(P) - 1.d0)

end function wesenheit_magnitude

function log_Efstathiou_likelihood(A,bw,sigma_int)    !    It computes equation (3) in published version of 1311.3461
    use arrays
    use fiducial
    Implicit none
    Real*8 :: log_Efstathiou_likelihood,A,bw,sigma_int,chi2,normalization
    Integer*4 :: m

    chi2 = 0.d0

    normalization = 0.d0

    Do m=1, size(Name)

        chi2 = ( observed_wesenheit_magnitude(H(m),V(m),II(m)) - wesenheit_magnitude(A,bw,Period(m)) )**2/&
        ( Sigma_m(m)**2 + sigma_int**2 ) + chi2

        normalization = log(Sigma_m(m)**2 + sigma_int**2) + normalization 
   
    End Do

    If ((abs(chi2) .ge. 0.d0) .and. (normalization**2 .ge. 0.d0 )) then

        log_Efstathiou_likelihood = - chi2/2.d0 - normalization/2.d0

    Else

        log_Efstathiou_likelihood = -1.d10
   
    End If

end function log_Efstathiou_likelihood

function new_chi2(chi2)
    use fiducial
    Implicit none
    Real*8 :: chi2,new_chi2

    If (chi2 .eq. 0.d0) then
  
        new_chi2 = sqrt(2.d0/Pi/9.d0)
  
    Else 

        new_chi2 = erf(sqrt(chi2/2.d0))/chi2**(3.d0/2.d0) - sqrt(2.d0/Pi)/chi2*exp(-chi2/2.d0)

    End If

end function new_chi2

function log_Efstathiou_likelihoodA(A,bw,sigma_int)    !    It computes equation (3) in published version of 1311.3461
    use arrays
    use fiducial
    Implicit none
    Real*8 :: log_Efstathiou_likelihoodA,A,bw,sigma_int,normalizationA
    Integer*4 :: m

    If (separate_dataA) then
    
        log_Efstathiou_likelihoodA = 0.d0

        Do m=1,size(NameA)
 
            If (using_jeffreys_prior) then

                log_Efstathiou_likelihoodA = -log(chi2A_i(A,bw,sigma_int,m))/2.d0 + log(N_tilde_A_i(sigma_int,m))&
                + log_Efstathiou_likelihoodA

            Else

                log_Efstathiou_likelihoodA = log(new_chi2(chi2A_i(A,bw,sigma_int,m))) + log(N_tilde_A_i(sigma_int,m))&
                + log_Efstathiou_likelihoodA

            End If

        End Do

        If ( abs(log_Efstathiou_likelihoodA) .ge. 0.d0 ) then

            continue

        Else 

            log_Efstathiou_likelihoodA = -1.d10

        End If

    Else

        normalizationA = 0.d0

        Do m=1, size(NameA)

            normalizationA = log(Sigma_mA(m)**2 + sigma_int**2) + normalizationA
   
        End Do

        If ((abs(chi2A(A,bw,sigma_int)) .ge. 0.d0) .and. (normalizationA**2 .ge. 0.d0)) then

            log_Efstathiou_likelihoodA = - dble(size(NameA))*log(chi2A(A,bw,sigma_int))/2.d0 - normalizationA/2.d0

        Else

            log_Efstathiou_likelihoodA = -1.d10
   
        End If

    End If

end function log_Efstathiou_likelihoodA

function chi2A(A,bw,sigma_int)    !    It computes equation (3) in published version of 1311.3461
    use arrays
    use fiducial
    Implicit none
    Real*8 :: A,bw,sigma_int,chi2A
    Integer*4 :: m

    chi2A = 0.d0

    Do m=1, size(NameA)

        chi2A = ( observed_wesenheit_magnitude(HA(m),VA(m),IIA(m)) - wesenheit_magnitude(A,bw,PeriodA(m)) )**2/&
        ( Sigma_mA(m)**2 + sigma_int**2 ) + chi2A

    End Do

end function chi2A

function chi2A_i(A,bw,sigma_int,m)    !    It computes equation (3) in published version of 1311.3461
    use arrays
    use fiducial
    Implicit none
    Real*8 :: A,bw,sigma_int,chi2A_i
    Integer*4 :: m

    chi2A_i = ( observed_wesenheit_magnitude(HA(m),VA(m),IIA(m)) - wesenheit_magnitude(A,bw,PeriodA(m)) )**2/&
    ( Sigma_mA(m)**2 + sigma_int**2 ) 

end function chi2A_i

function N_tilde_A_i(sigma_int,m)    !    It computes equation (3) in published version of 1311.3461
    use arrays
    use fiducial
    Implicit none
    Real*8 :: sigma_int,N_tilde_A_i
    Integer*4 :: m

    N_tilde_A_i = 1.d0/sqrt( Sigma_mA(m)**2 + sigma_int**2 ) 

end function N_tilde_A_i

function chi2B_i(A,bw,sigma_int,m)    !    It computes equation (3) in published version of 1311.3461
    use arrays
    use fiducial
    Implicit none
    Real*8 :: A,bw,sigma_int,chi2B_i
    Integer*4 :: m

    chi2B_i = ( observed_wesenheit_magnitude(HB(m),VB(m),IIB(m)) - wesenheit_magnitude(A,bw,PeriodB(m)) )**2/&
    ( Sigma_mB(m)**2 + sigma_int**2 ) 

end function chi2B_i

function N_tilde_B_i(sigma_int,m)    !    It computes equation (3) in published version of 1311.3461
    use arrays
    use fiducial
    Implicit none
    Real*8 :: sigma_int,N_tilde_B_i
    Integer*4 :: m

    N_tilde_B_i = 1.d0/sqrt( Sigma_mB(m)**2 + sigma_int**2 ) 

end function N_tilde_B_i

function chi2C_i(A,bw,sigma_int,m)    !    It computes equation (3) in published version of 1311.3461
    use arrays
    use fiducial
    Implicit none
    Real*8 :: A,bw,sigma_int,chi2C_i
    Integer*4 :: m

    chi2C_i = ( observed_wesenheit_magnitude(HC(m),VC(m),IIC(m)) - wesenheit_magnitude(A,bw,PeriodC(m)) )**2/&
    ( Sigma_mC(m)**2 + sigma_int**2 ) 

end function chi2C_i

function N_tilde_C_i(sigma_int,m)    !    It computes equation (3) in published version of 1311.3461
    use arrays
    use fiducial
    Implicit none
    Real*8 :: sigma_int,N_tilde_C_i
    Integer*4 :: m

    N_tilde_C_i = 1.d0/sqrt( Sigma_mC(m)**2 + sigma_int**2 ) 

end function N_tilde_C_i

function log_Efstathiou_likelihoodB(A,bw,sigma_int)    !    It computes equation (3) in published version of 1311.3461
    use arrays
    use fiducial
    Implicit none
    Real*8 :: log_Efstathiou_likelihoodB,A,bw,sigma_int,normalizationB
    Integer*4 :: m

    If (separate_dataB) then

        log_Efstathiou_likelihoodB = 0.d0
    
        Do m=1,size(NameB)

            If (using_jeffreys_prior) then

                log_Efstathiou_likelihoodB = -log(chi2B_i(A,bw,sigma_int,m))/2.d0 + log(N_tilde_B_i(sigma_int,m))&
                + log_Efstathiou_likelihoodB

            Else

                log_Efstathiou_likelihoodB = log(new_chi2(chi2B_i(A,bw,sigma_int,m))) + log(N_tilde_B_i(sigma_int,m))&
                + log_Efstathiou_likelihoodB

            End If

        End Do

        If ( abs(log_Efstathiou_likelihoodB) .ge. 0.d0 ) then

            continue

        Else 

            log_Efstathiou_likelihoodB = -1.d10

        End If

    Else

        normalizationB = 0.d0

        Do m=1, size(NameB)

            normalizationB = log(Sigma_mB(m)**2 + sigma_int**2) + normalizationB 
   
        End Do

        If ((abs(chi2B(A,bw,sigma_int)) .ge. 0.d0) .and. (normalizationB**2 .ge. 0.d0)) then

            log_Efstathiou_likelihoodB = -dble(size(NameB))*log(chi2B(A,bw,sigma_int))/2.d0 - normalizationB/2.d0

        Else

            log_Efstathiou_likelihoodB = -1.d10
   
        End If

    End If

end function log_Efstathiou_likelihoodB

function chi2B(A,bw,sigma_int)    !    It computes equation (3) in published version of 1311.3461
    use arrays
    use fiducial
    Implicit none
    Real*8 :: A,bw,sigma_int,chi2B
    Integer*4 :: m

    chi2B = 0.d0

    Do m=1, size(NameB)

        chi2B = ( observed_wesenheit_magnitude(HB(m),VB(m),IIB(m)) - wesenheit_magnitude(A,bw,PeriodB(m)) )**2/&
        ( Sigma_mB(m)**2 + sigma_int**2 ) + chi2B

    End Do

end function chi2B


function log_Efstathiou_likelihoodC(A,bw,sigma_int)    !    It computes equation (3) in published version of 1311.3461
    use arrays
    use fiducial
    Implicit none
    Real*8 :: log_Efstathiou_likelihoodC,A,bw,sigma_int,normalizationC
    Integer*4 :: m

    If (separate_dataC) then
    
        log_Efstathiou_likelihoodC = 0.d0

        Do m=1,size(NameC)

            If (using_jeffreys_prior) then

                log_Efstathiou_likelihoodC = -log(chi2C_i(A,bw,sigma_int,m))/2.d0 + log(N_tilde_C_i(sigma_int,m))&
                + log_Efstathiou_likelihoodC

            Else

                log_Efstathiou_likelihoodC = log(new_chi2(chi2C_i(A,bw,sigma_int,m))) + log(N_tilde_C_i(sigma_int,m))&
                + log_Efstathiou_likelihoodC

            End If

        End Do

        If ( abs(log_Efstathiou_likelihoodC) .ge. 0.d0 ) then

            continue

        Else 

            log_Efstathiou_likelihoodC = -1.d10

        End If

    Else

        normalizationC = 0.d0

        Do m=1, size(NameC)

            normalizationC = log(Sigma_mC(m)**2 + sigma_int**2) + normalizationC 
   
        End Do

        If ((abs(chi2C(A,bw,sigma_int)) .ge. 0.d0) .and. (normalizationC**2 .ge. 0.d0)) then

            log_Efstathiou_likelihoodC = -dble(size(NameC))*log(chi2C(A,bw,sigma_int))/2.d0 - normalizationC/2.d0

        Else

            log_Efstathiou_likelihoodC = -1.d10
   
        End If

    End If

end function log_Efstathiou_likelihoodC

function chi2C(A,bw,sigma_int)    !    It computes equation (3) in published version of 1311.3461
    use arrays
    use fiducial
    Implicit none
    Real*8 :: A,bw,sigma_int,chi2C
    Integer*4 :: m

    chi2C = 0.d0

    Do m=1, size(NameC)

        chi2C = ( observed_wesenheit_magnitude(HC(m),VC(m),IIC(m)) - wesenheit_magnitude(A,bw,PeriodC(m)) )**2/&
        ( Sigma_mC(m)**2 + sigma_int**2 ) + chi2C

    End Do

end function chi2C

function log_Efstathiou_likelihood_hyperparameters(A,bw,sigma_int)
    use arrays
    use fiducial 
    Implicit none
    Real*8 :: log_Efstathiou_likelihood_hyperparameters,A,bw,sigma_int
 
    If ((include_dataA .and. include_dataB) .and. include_dataC) then

        log_Efstathiou_likelihood_hyperparameters = log_Efstathiou_likelihoodA(A,bw,sigma_int) + &
        log_Efstathiou_likelihoodB(A,bw,sigma_int) + log_Efstathiou_likelihoodC(A,bw,sigma_int)

    Else if ( (include_dataA .and. include_dataB) .and. .not.include_dataC ) then

        log_Efstathiou_likelihood_hyperparameters = log_Efstathiou_likelihoodA(A,bw,sigma_int) + &
        log_Efstathiou_likelihoodB(A,bw,sigma_int) 

    Else if ( (include_dataA .and. .not.include_dataB) .and. include_dataC ) then

        log_Efstathiou_likelihood_hyperparameters = log_Efstathiou_likelihoodA(A,bw,sigma_int) + &
        log_Efstathiou_likelihoodC(A,bw,sigma_int)

    Else if ( (include_dataA .and. .not.include_dataB) .and. .not.include_dataC ) then

        log_Efstathiou_likelihood_hyperparameters = log_Efstathiou_likelihoodA(A,bw,sigma_int) 

    Else if ( (.not.include_dataA .and. include_dataB) .and. .not.include_dataC ) then

        log_Efstathiou_likelihood_hyperparameters = log_Efstathiou_likelihoodB(A,bw,sigma_int) 

    Else if ( (.not.include_dataA .and. .not.include_dataB) .and. include_dataC ) then

        log_Efstathiou_likelihood_hyperparameters = log_Efstathiou_likelihoodC(A,bw,sigma_int)

    Else if ( (.not.include_dataA .and. include_dataB) .and. include_dataC ) then

        log_Efstathiou_likelihood_hyperparameters = log_Efstathiou_likelihoodB(A,bw,sigma_int) + &
        log_Efstathiou_likelihoodC(A,bw,sigma_int)

    Else 

        print *,'Data must be included from fiducial module in order to compute likelihood '

        stop

    End If
    
end function log_Efstathiou_likelihood_hyperparameters

function log_Gaussian_likelihood(array)
    use fiducial
    Implicit none
    Integer*4 :: index
    Real*8 :: log_Gaussian_likelihood
    Real*8,dimension(number_of_parameters) :: array
    log_Gaussian_likelihood = 0.d0
    Do index=1,number_of_parameters
        log_Gaussian_likelihood = array(index)**2 + log_Gaussian_likelihood
    End Do
    log_Gaussian_likelihood = -log_Gaussian_likelihood/2.d0
end function log_Gaussian_likelihood

!function compute_determinant(A)
!    use fiducial
!    Implicit none
!    Integer*4 :: INFO,index
!    Real*8,dimension(max(1,nbins),nbins) :: A
!    Integer*4,dimension(min(nbins,nbins)) :: IPIV
!    Real*8,dimension(max(1,max(1,nbins))) :: WORK

!    Real*8 :: det,sgn,compute_determinant

!    call dgetrf(nbins,nbins,A,nbins,IPIV,INFO)

!    det = 1.d0
!    Do index=1,nbins

!        det = det*A(index,index)
!    End Do 

!    sgn = 1.d0
!    Do index=1,nbins
!        If (IPIV(index) .ne. index) then
!            sgn = -sgn
!        End If
!    End Do
 
!    compute_determinant = sgn*det
!end function compute_determinant

subroutine read_covariance_matrix_mcmc(matrix)
    use fiducial
    Implicit none
    Real*8,dimension(number_of_parameters,number_of_parameters) :: matrix
    Integer*4 :: index1
    open(12,file='./output/covariance_matrix.txt')
    Do index1=1,number_of_parameters
        read(12,*) matrix(index1,1),matrix(index1,2)!,matrix(index1,3)
    End Do
    close(12)
end subroutine read_covariance_matrix_mcmc

subroutine read_bestfit_mcmc(vector)
    use fiducial
    Implicit none
    Real*8,dimension(number_of_parameters) :: vector
    Integer*4 :: index1
    open(12,file='./output/bestfit.txt')
    Do index1=1,number_of_parameters
        read(12,*) vector(index1)
    End Do
    close(12)
end subroutine read_bestfit_mcmc

subroutine read_means_mcmc(vector)
    use fiducial
    Implicit none
    Real*8,dimension(number_of_parameters) :: vector
    Integer*4 :: index1
    open(12,file='./output/means.txt')
    Do index1=1,number_of_parameters
        read(12,*) vector(index1)
    End Do
    close(12)
end subroutine read_means_mcmc


!subroutine inverting_matrix()
!    use fiducial
!    use arrays
!    Implicit none
!    Integer*4 :: l,M,N,LDA,INFO,LWORK,i,j
!    Real*8,dimension(max(1,nbins*(nbins+1)/2),nbins*(nbins+1)/2) :: A
!    Integer*4,dimension(min(nbins*(nbins+1)/2,nbins*(nbins+1)/2)) :: IPIV
!    Real*8,dimension(max(1,max(1,nbins*(nbins+1)/2))) :: WORK
!    M = nbins*(nbins+1)/2
!    N = M
!    LDA = max(1,M)
!    LWORK = max(1,N)

!    Do l=lmin,lmax
!        Do i=1,M
!            Do j=1,M
!                A(i,j) = cov_l_IP(l,i,j)
!            End Do
!        End Do

!        call dgetrf(M,N,A,LDA,IPIV,INFO)

!        call dgetri(N,A,LDA,IPIV,WORK,LWORK,INFO)

!        Do i=1,M
!            Do j=1,M
!                inv_cov_l_IP(l,i,j) = A(i,j)
!            End Do
!        End Do
!    End Do
    
!end subroutine inverting_matrix

end module functions
