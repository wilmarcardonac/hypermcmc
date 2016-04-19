module functions

    Implicit none
     
contains

subroutine read_table2_R11(path_to_datafile)
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

    allocate (Field(1:arrays_dimension),ID(1:arrays_dimension),PeriodR11(1:arrays_dimension),&
    VIR11(1:arrays_dimension),F160WR11(1:arrays_dimension),eF160WR11(1:arrays_dimension),&
    OHR11(1:arrays_dimension),stat=status1)

    open(11,file=path_to_datafile)

    Do p=1,arrays_dimension

        read(11,'(a5,i8,f9.3,4f6.2)') Field(p),ID(p),PeriodR11(p),VIR11(p),F160WR11(p),eF160WR11(p),OHR11(p)

    End Do

    close(11)

end subroutine read_table2_R11

subroutine read_table3_R11(path_to_datafile)

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

    allocate (Fieldmvi(1:arrays_dimension),mvi5av(1:arrays_dimension),Sigma_mvi5av(1:arrays_dimension),stat=status1)

    open(11,file=path_to_datafile)

    Do p=1,arrays_dimension

        read(11,*) Fieldmvi(p),mvi5av(p),Sigma_mvi5av(p)

    End Do

    close(11)

end subroutine read_table3_R11

subroutine read_table2_LEEUWEN(path_to_datafile)

    use arrays
    Implicit none
    Integer*4 :: arrays_dimension,p
    Integer :: stat
    character(len=*) :: path_to_datafile
    Real*8,allocatable,dimension(:) :: LK

    open(11,file=path_to_datafile)

    read(11,*)

    read(11,*)

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

    allocate (FieldHipp(1:arrays_dimension),logP(1:arrays_dimension),Mw(1:arrays_dimension),&
         sigmaMw(1:arrays_dimension),LK(1:arrays_dimension),stat=status1)

    open(11,file=path_to_datafile)

    read(11,*)

    read(11,*)

    Do p=1,arrays_dimension

        read(11,*) FieldHipp(p),logP(p),Mw(p),sigmaMw(p),LK(p)

        MW(p) = Mw(p) + LK(p)

    End Do

    close(11)

end subroutine read_table2_LEEUWEN

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
    
    use fiducial
    Implicit none    !    published version of 1311.3461.
    Real*8 :: observed_wesenheit_magnitude,H_band_magnitude,V_magnitude,I_magnitude

    observed_wesenheit_magnitude =  H_band_magnitude - R*(V_magnitude - I_magnitude)

end function observed_wesenheit_magnitude

!FIX OBERVED_M_H FUNCTION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
function observed_m_H(W_band_magnitude,V_I_magnitude)    ! EQUATION (6) IN R09 
    
    use fiducial
    Implicit none
    Real*8 :: observed_m_H,W_band_magnitude,V_I_magnitude

    observed_m_H =  W_band_magnitude + R*V_I_magnitude

end function observed_m_H

function observed_m_W(H_band_magnitude,V_I_magnitude)    ! EQUATION (6) IN R09 
    
    use fiducial
    Implicit none
    Real*8 :: observed_m_W,H_band_magnitude,V_I_magnitude

    observed_m_W =  H_band_magnitude - R*V_I_magnitude

end function observed_m_W

function wesenheit_magnitude(A,bw,P)    !    It computes equation (2) in published version of 1311.3461
    Implicit none
    Real*8 :: wesenheit_magnitude,A,bw,P

    wesenheit_magnitude = A + bw*(log10(P) - 1.d0)

end function wesenheit_magnitude

function P_L_relation_passband_H(zpH_j,zpH_ref,bH,P_ij) ! EQUATION (2) IN R09

    Implicit none
    Real*8 :: P_L_relation_passband_H,zpH_j,zpH_ref,bH,P_ij

    P_L_relation_passband_H = ( zpH_j - zpH_ref ) + zpH_ref + bH*log10(P_ij)

end function P_L_relation_passband_H

function P_L_relation_passband_H_2(mu0i,mu0_ref,zpW_ref,bW,Zw,OH_ij,P_ij) ! EQUATION (7) IN R09

    Implicit none
    Real*8 :: P_L_relation_passband_H_2,mu0i,mu0_ref,zpW_ref,bW,P_ij,Zw,OH_ij

    P_L_relation_passband_H_2 = (mu0i - mu0_ref) + zpW_ref + bW*log10(P_ij) + Zw*OH_ij

end function P_L_relation_passband_H_2

function P_L_relation_passband_W(mu0i,mu0_ref,zpW_ref,bW,Zw,OH_ij,P_ij) ! EQUATION (7) IN R09

    Implicit none
    Real*8 :: P_L_relation_passband_W,mu0i,mu0_ref,zpW_ref,bW,P_ij,Zw,OH_ij

    P_L_relation_passband_W = (mu0i - mu0_ref) + zpW_ref + bW*log10(P_ij) + Zw*OH_ij

end function P_L_relation_passband_W

function P_L_relation_passband_W_2(zpW,bW,Zw,OH_ij,P_ij) ! EQUATION (7) IN R09

    Implicit none
    Real*8 :: P_L_relation_passband_W_2,zpW,bW,P_ij,Zw,OH_ij

    P_L_relation_passband_W_2 = zpW + bW*log10(P_ij) + Zw*OH_ij

end function P_L_relation_passband_W_2

function P_L_relation_passband_W_E14(mu0,M_w,bW,Zw,OH_ij,P_ij) ! EQUATION (7) IN R09

    Implicit none
    Real*8 :: P_L_relation_passband_W_E14,mu0,M_w,bW,P_ij,Zw,OH_ij

    P_L_relation_passband_W_E14 = mu0 + M_w + bW*(log10(P_ij) - 1.d0) + Zw*OH_ij

end function P_L_relation_passband_W_E14

function reddening_free_magnitude_SNIa(mu0i,H0,av) ! EQUATION (16) IN R09

    Implicit none
    Real*8 :: reddening_free_magnitude_SNIa,mu0i,H0,av

    reddening_free_magnitude_SNIa = mu0i + 5.d0*log10(H0) - 25.d0 - 5.d0*av

end function reddening_free_magnitude_SNIa

function log_Efstathiou_likelihood(A,bw,sigma_int)    !    It computes equation (3) in published version of 1311.3461
    use arrays
    use fiducial
    Implicit none
    Real*8 :: log_Efstathiou_likelihood,A,bw,sigma_int,chi2,normalization
    Integer*4 :: m

    If (using_hyperparameters .and. use_HP_per_cepheid) then

        log_Efstathiou_likelihood = 0.d0

        Do m=1,size(Name)
 
           If ( (Period(m) .gt. cepheid_lower_Period_limit) .and. (Period(m) .lt. cepheid_Period_limit)) then

              If (using_jeffreys_prior) then

                 log_Efstathiou_likelihood = -log(chi2_i(A,bw,sigma_int,m))/2.d0 + log(N_tilde_i(sigma_int,m))&
                      + log_Efstathiou_likelihood

              Else
               
                 log_Efstathiou_likelihood = log(new_chi2(chi2_i(A,bw,sigma_int,m))) + log(N_tilde_i(sigma_int,m))&
                      + log_Efstathiou_likelihood

              End If

            End If

         End Do

        If ( abs(log_Efstathiou_likelihood) .ge. 0.d0 ) then

           continue

        Else 

           log_Efstathiou_likelihood = -1.d10

        End If

    Else

       chi2 = 0.d0

       normalization = 0.d0

       Do m=1, size(Name)

          If ( (Period(m) .gt. cepheid_lower_Period_limit) .and. (Period(m) .lt. cepheid_Period_limit)) then

             chi2 = ( observed_wesenheit_magnitude(H(m),V(m),II(m)) - wesenheit_magnitude(A,bw,Period(m)) )**2/&
                  ( Sigma_m(m)**2 + sigma_int**2 ) + chi2

             normalization = log(Sigma_m(m)**2 + sigma_int**2) + normalization 

          End If

       End Do

       If ((abs(chi2) .ge. 0.d0) .and. (normalization**2 .ge. 0.d0 )) then

          log_Efstathiou_likelihood = - chi2/2.d0 - normalization/2.d0

       Else

          log_Efstathiou_likelihood = -1.d10
   
       End If

    End If

end function log_Efstathiou_likelihood

function log_R11_likelihood_H_NGC4258(mu0j,zpw_ref,bw,H0,Zw,av,sigma_int)    !    EQUATION (4) IN R09

    use arrays
    use fiducial
    Implicit none

    Real*8 :: log_R11_likelihood_H_NGC4258,zpw_ref,bw,H0,Zw,av,sigma_int!,normalizationA
    Real*8,dimension(number_of_hosts_galaxies) :: mu0j 
    Integer*4 :: m,index_host!,number_cepheid

    log_R11_likelihood_H_NGC4258 = 0.d0

    If (use_HP_per_host) then
 
       print *, 'NEED TO IMPLEMENT HP PER HOST IN LOG_R11_LIKELIHOOD_H_NGC4258'

       stop

!       If (using_jeffreys_prior) then

!          Do index_host=1,number_of_hosts_galaxies

!             normalizationA = 0.d0

!             number_cepheid = 0

!             Do m=1, size(Field)

!                If (host(index_host) .eq. Field(m)) then
                       
!                   number_cepheid = number_cepheid + 1

!                   normalizationA = log( eF160WR11(m)**2 + sigma_int**2 ) + normalizationA

!                End If

!             End Do

!             log_R11_likelihood_H_NGC4258 = - dble(number_cepheid)*log(chi2R11_W_host(mu0j(index_host),mu0j(9),&
!                  zpw_ref,bw,Zw,sigma_int,index_host))/2.d0 - normalizationA/2.d0 + log_R11_likelihood_H_NGC4258

!          End Do

!       Else

!          print *, 'UNIFORM PRIOR NOT IMPLEMENTED YET. NEED TO CODE EXPRESSION WITH GAMMA FUNCTIONS'

!          stop

!       End If

    Else

       If (use_HP_per_cepheid) then

          Do m=1,size(Field)

             Do index_host=1,number_of_hosts_galaxies
 
                If (host(index_host) .eq. Field(m)) then
    
                   If (using_jeffreys_prior) then

                      print *, 'IMPROPER JEFFREYS PRIOR LEADS TO SINGULARITIES AND THEREFORE IS NOT IMPLEMENTED'

                      stop

                   Else
                       
                      If (PeriodR11(m) .lt. cepheid_Period_limit) then

                         log_R11_likelihood_H_NGC4258 = log(new_chi2(chi2R11_H_2(mu0j(index_host),mu0j(9),zpw_ref,&
                              bw,Zw,sigma_int,m))) + log(N_tilde_R11_H(sigma_int,m)) + log_R11_likelihood_H_NGC4258

                      End If

                   End If

                End If

             End Do

          End Do

       Else

          Do m=1,size(Field)

             Do index_host=1,number_of_hosts_galaxies
 
                If (host(index_host) .eq. Field(m)) then
    
                   If (PeriodR11(m) .lt. cepheid_Period_limit) then

                      log_R11_likelihood_H_NGC4258 = -chi2R11_H_2(mu0j(index_host),mu0j(9),zpw_ref,bw,Zw,sigma_int,m)/2.d0 + &
                           log(N_tilde_R11_H(sigma_int,m)) + log_R11_likelihood_H_NGC4258

                   End If

                End If

             End Do

          End Do

       End If

    End If

    Do index_host=1,number_of_hosts_galaxies-1

       If (use_HP_in_SNIa) then

          log_R11_likelihood_H_NGC4258 = log(new_chi2(chi2R11_SNIa(mu0j(index_host),H0,av,index_host))) + &
               log(N_tilde_R11_SNIa(index_host)) + log_R11_likelihood_H_NGC4258  

       Else

          log_R11_likelihood_H_NGC4258 = -chi2R11_SNIa(mu0j(index_host),H0,av,index_host)/2.d0 + &
               log(N_tilde_R11_SNIa(index_host)) + log_R11_likelihood_H_NGC4258  

       End If
              
    End Do
        
    If (use_HP_in_av) then

       log_R11_likelihood_H_NGC4258 = log(new_chi2((a_v - av)**2/sigma_a_v**2)) - log(2.d0*Pi*sigma_a_v**2)/2.d0  +&
            log_R11_likelihood_H_NGC4258

    Else

       log_R11_likelihood_H_NGC4258 = -((a_v - av)**2/sigma_a_v**2 + log(2.d0*Pi*sigma_a_v**2) )/2.d0  + &
            log_R11_likelihood_H_NGC4258

    End If
        
    If (use_HP_in_anchor) then

       log_R11_likelihood_H_NGC4258 =  log(new_chi2(chi2R11_anchor_NGC4258(mu0j(9)))) - log(2.d0*Pi*sigma_mu_0_NGC4258**2)/2.d0 + &
            log_R11_likelihood_H_NGC4258

    Else

       log_R11_likelihood_H_NGC4258 =  -(chi2R11_anchor_NGC4258(mu0j(9)) + log(2.d0*Pi*sigma_mu_0_NGC4258**2))/2.d0 + &
            log_R11_likelihood_H_NGC4258

    End If

    If (use_prior_on_zpw4258) then

       log_R11_likelihood_H_NGC4258 = -((zpw_ref - prior_zpw4258)**2/sigma_zpw4258**2 + log(2.d0*Pi*sigma_zpw4258**2) )/2.d0  +&
            log_R11_likelihood_H_NGC4258

    Else 

       continue

    End If
        
    If ( abs(log_R11_likelihood_H_NGC4258) .ge. 0.d0 ) then

       continue

    Else 

       log_R11_likelihood_H_NGC4258 = -1.d10

    End If

end function log_R11_likelihood_H_NGC4258

function log_R11_likelihood_H(zpH,bH,sigma_int)    !    EQUATION (4) IN R09

    use arrays
    use fiducial
    Implicit none

    Real*8 :: log_R11_likelihood_H,bH,sigma_int
    Real*8,dimension(number_of_hosts_galaxies) :: zpH 
    Integer*4 :: m,index_host

    log_R11_likelihood_H = 0.d0

    Do m=1,size(Field)

       Do index_host=1,number_of_hosts_galaxies
 
          If (host(index_host) .eq. Field(m)) then
    
             If (using_jeffreys_prior) then

                print *, 'IMPROPER JEFFREYS PRIOR LEADS TO SINGULARITIES AND THEREFORE IS NOT IMPLEMENTED'

                stop

             Else

                If (PeriodR11(m) .lt. cepheid_Period_limit) then

                   log_R11_likelihood_H = log(new_chi2(chi2R11_H(zpH(index_host),zpH(9),bH,sigma_int,m))) + &
                        log(N_tilde_R11_H(sigma_int,m)) + log_R11_likelihood_H

                End If 

             End If

          End If

       End Do

    End Do

    If ( abs(log_R11_likelihood_H) .ge. 0.d0 ) then

       continue

    Else 

       log_R11_likelihood_H = -1.d10

    End If

end function log_R11_likelihood_H

function chi2R11_H(zpH_j,zpH_ref,bH,sigma_int,m)    !    It computes equation (3) in published version of 1311.3461

    use arrays
    use fiducial
    Implicit none

    Real*8 :: zpH_j,zpH_ref,bH,sigma_int,chi2R11_H
    Integer*4 :: m

    chi2R11_H = ( F160WR11(m) - P_L_relation_passband_H(zpH_j,zpH_ref,bH,PeriodR11(m)) )**2/&
    ( eF160WR11(m)**2 + sigma_int**2 ) 

end function chi2R11_H


function N_tilde_R11_H(sigma_int,m)    !    It computes equation (3) in published version of 1311.3461

    use arrays
    use fiducial
    Implicit none

    Real*8 :: sigma_int,N_tilde_R11_H
    Integer*4 :: m

    N_tilde_R11_H = 1.d0/sqrt( eF160WR11(m)**2 + sigma_int**2 ) 

end function N_tilde_R11_H

function log_R11_likelihood_W_LMC(mu0j,zpw_ref,bw,H0,Zw,av,acal,sigma_int,sigma_int_LMC)    !    EQUATION (4) IN R09

    use arrays
    use fiducial
    Implicit none

    Real*8 :: log_R11_likelihood_W_LMC,zpw_ref,bw,H0,Zw,av,acal,sigma_int,sigma_int_LMC,normalizationA
    Real*8,dimension(number_of_hosts_galaxies + 1) :: mu0j 
    Integer*4 :: m,index_host,number_cepheid

    log_R11_likelihood_W_LMC = 0.d0

    If (use_HP_per_host) then  ! R11 SAMPLE
 
       If (using_jeffreys_prior) then

          Do index_host=1,number_of_hosts_galaxies

             normalizationA = 0.d0

             number_cepheid = 0

             Do m=1, size(Field)

                If (host(index_host) .eq. Field(m)) then
                       
                   number_cepheid = number_cepheid + 1

                   normalizationA = log( eF160WR11(m)**2 + sigma_int**2 ) + normalizationA

                End If

             End Do

             log_R11_likelihood_W_LMC = - dble(number_cepheid)*log(chi2R11_W_host(mu0j(index_host),mu0j(10),&
                  zpw_ref,bw,Zw,sigma_int,index_host))/2.d0 - normalizationA/2.d0 + log_R11_likelihood_W_LMC

          End Do

       Else

          print *, 'UNIFORM PRIOR NOT IMPLEMENTED YET. NEED TO CODE EXPRESSION WITH GAMMA FUNCTIONS'

          stop

       End If

    Else

       If (use_HP_per_cepheid) then

          Do m=1,size(Field)

             Do index_host=1,number_of_hosts_galaxies
 
                If (host(index_host) .eq. Field(m)) then
    
                   If (using_jeffreys_prior) then

                      print *, 'IMPROPER JEFFREYS PRIOR LEADS TO SINGULARITIES AND THEREFORE IS NOT IMPLEMENTED'

                      stop

                   Else
                       
                      If (PeriodR11(m) .lt. cepheid_Period_limit) then

                         log_R11_likelihood_W_LMC = log(new_chi2(chi2R11_W(mu0j(index_host),mu0j(10),zpw_ref,bw,&
                              Zw,sigma_int,m))) + log(N_tilde_R11_W(sigma_int,m)) + log_R11_likelihood_W_LMC

                      End If

                   End If

                End If

             End Do

          End Do

       Else

          Do m=1,size(Field)

             Do index_host=1,number_of_hosts_galaxies
 
                If (host(index_host) .eq. Field(m)) then
    
                   If (PeriodR11(m) .lt. cepheid_Period_limit) then

                      log_R11_likelihood_W_LMC = -chi2R11_W(mu0j(index_host),mu0j(10),zpw_ref,bw,Zw,sigma_int,m)/2.d0 + &
                           log(N_tilde_R11_W(sigma_int,m)) + log_R11_likelihood_W_LMC

                   End If

                End If

             End Do

          End Do

       End If

    End If

    If (use_HP_per_cepheid) then ! LMC CEPHEID VARIABLES

       Do m=1,size(Name)

          If (using_jeffreys_prior) then

             print *, 'IMPROPER JEFFREYS PRIOR LEADS TO SINGULARITIES AND THEREFORE IS NOT IMPLEMENTED'

             stop

          Else
                       
             If (Period(m) .lt. cepheid_Period_limit) then

                log_R11_likelihood_W_LMC = log(new_chi2(chi2R11_W_LMC(mu0j(10),mu0j(10),zpw_ref,bw,Zw,sigma_int_LMC,m))) + &
                     log(N_tilde_R11_W_LMC(sigma_int_LMC,m)) + log_R11_likelihood_W_LMC
                      
             End If

          End If

       End Do

    Else

       Do m=1,size(Name)

          If (Period(m) .lt. cepheid_Period_limit) then

             log_R11_likelihood_W_LMC = -chi2R11_W_LMC(mu0j(10),mu0j(10),zpw_ref,bw,Zw,sigma_int_LMC,m)/2.d0 + &
                  log(N_tilde_R11_W_LMC(sigma_int_LMC,m)) + log_R11_likelihood_W_LMC

          End If

       End Do

    End If

    Do index_host=1,number_of_hosts_galaxies-1 ! SN Ia

       If (use_HP_in_SNIa) then

          log_R11_likelihood_W_LMC = log(new_chi2(chi2R11_SNIa(mu0j(index_host),H0,av,index_host))) + &
               log(N_tilde_R11_SNIa(index_host)) + log_R11_likelihood_W_LMC  

       Else

          log_R11_likelihood_W_LMC = -chi2R11_SNIa(mu0j(index_host),H0,av,index_host)/2.d0 + &
               log(N_tilde_R11_SNIa(index_host)) + log_R11_likelihood_W_LMC

       End If
              
    End Do
        
    If (use_HP_in_av) then ! a_v AND a_cal PARAMETERS

       log_R11_likelihood_W_LMC = log(new_chi2((a_v - av)**2/sigma_a_v**2)) - log(2.d0*Pi*sigma_a_v**2)/2.d0  + &
            log(new_chi2((a_cal - acal)**2/sigma_a_cal**2)) - log(2.d0*Pi*sigma_a_cal**2)/2.d0 + log_R11_likelihood_W_LMC

    Else

       log_R11_likelihood_W_LMC = -((a_v - av)**2/sigma_a_v**2 + log(2.d0*Pi*sigma_a_v**2) )/2.d0  - &
            ((a_cal - acal)**2/sigma_a_cal**2 + log(2.d0*Pi*sigma_a_cal**2) )/2.d0 + log_R11_likelihood_W_LMC

    End If
        
    If (use_HP_in_anchor) then ! ANCHOR

       log_R11_likelihood_W_LMC =  log(new_chi2(chi2R11_anchor_LMC(mu0j(10)))) - log(2.d0*Pi*sigma_mu_0_LMC**2)/2.d0 + &
            log_R11_likelihood_W_LMC

    Else

       log_R11_likelihood_W_LMC =  -(chi2R11_anchor_LMC(mu0j(10)) + log(2.d0*Pi*sigma_mu_0_LMC**2))/2.d0 + &
            log_R11_likelihood_W_LMC

    End If

    If (use_prior_on_Zw) then

       log_R11_likelihood_W_LMC = -((Zw - prior_Zw)**2/sigma_Zw_prior**2 + log(2.d0*Pi*sigma_Zw_prior**2) )/2.d0  +&
            log_R11_likelihood_W_LMC

    Else 

       continue

    End If

    If (use_prior_on_bw) then

       log_R11_likelihood_W_LMC = -((bw - prior_bw_from_LMC)**2/sigma_bw_prior**2 + log(2.d0*Pi*sigma_bw_prior**2) )/2.d0  +&
            log_R11_likelihood_W_LMC

    Else 

       continue

    End If
        
    If ( abs(log_R11_likelihood_W_LMC) .ge. 0.d0 ) then

       continue

    Else 

       log_R11_likelihood_W_LMC = -1.d10

    End If

end function log_R11_likelihood_W_LMC

function log_R11_likelihood_W_LMC_NGC4258(mu0j,M_w,bw,H0,Zw,av,acal,sigma_int,sigma_int_LMC)    !    EQUATION (4) IN R09

    use arrays
    use fiducial
    Implicit none

    Real*8 :: log_R11_likelihood_W_LMC_NGC4258,M_w,bw,H0,Zw,av,acal,sigma_int,sigma_int_LMC,normalizationA
    Real*8,dimension(number_of_hosts_galaxies + 1) :: mu0j 
    Integer*4 :: m,index_host,number_cepheid

    log_R11_likelihood_W_LMC_NGC4258 = 0.d0

    If (use_HP_per_host) then  ! R11 SAMPLE
 
       If (using_jeffreys_prior) then

          Do index_host=1,number_of_hosts_galaxies

             normalizationA = 0.d0

             number_cepheid = 0

             Do m=1, size(Field)

                If (host(index_host) .eq. Field(m)) then
                       
                   number_cepheid = number_cepheid + 1

                   normalizationA = log( eF160WR11(m)**2 + sigma_int**2 ) + normalizationA

                End If

             End Do

             log_R11_likelihood_W_LMC_NGC4258 = - dble(number_cepheid)*log(chi2R11_W_host_E14(mu0j(index_host),M_w,&
                  bw,Zw,sigma_int,index_host))/2.d0 - normalizationA/2.d0 + log_R11_likelihood_W_LMC_NGC4258

          End Do

       Else

          print *, 'UNIFORM PRIOR NOT IMPLEMENTED YET. NEED TO CODE EXPRESSION WITH GAMMA FUNCTIONS'

          stop

       End If

    Else

       If (use_HP_per_cepheid) then

          Do m=1,size(Field)

             Do index_host=1,number_of_hosts_galaxies
 
                If (host(index_host) .eq. Field(m)) then
    
                   If (using_jeffreys_prior) then

                      print *, 'IMPROPER JEFFREYS PRIOR LEADS TO SINGULARITIES AND THEREFORE IS NOT IMPLEMENTED'

                      stop

                   Else
                       
                      If (PeriodR11(m) .lt. cepheid_Period_limit) then

                         log_R11_likelihood_W_LMC_NGC4258 = log(new_chi2(chi2R11_W_E14(mu0j(index_host),M_w,bw,&
                              Zw,sigma_int,m))) + log(N_tilde_R11_W(sigma_int,m)) + log_R11_likelihood_W_LMC_NGC4258

                      End If

                   End If

                End If

             End Do

          End Do

       Else

          Do m=1,size(Field)

             Do index_host=1,number_of_hosts_galaxies
 
                If (host(index_host) .eq. Field(m)) then
    
                   If (PeriodR11(m) .lt. cepheid_Period_limit) then

                      log_R11_likelihood_W_LMC_NGC4258 = -chi2R11_W_E14(mu0j(index_host),M_w,bw,Zw,sigma_int,m)/2.d0 + &
                           log(N_tilde_R11_W(sigma_int,m)) + log_R11_likelihood_W_LMC_NGC4258

                   End If

                End If

             End Do

          End Do

       End If

    End If

    If (use_HP_per_cepheid) then ! LMC CEPHEID VARIABLES

       Do m=1,size(Name)

          If (using_jeffreys_prior) then

             print *, 'IMPROPER JEFFREYS PRIOR LEADS TO SINGULARITIES AND THEREFORE IS NOT IMPLEMENTED'

             stop

          Else
                       
             If (Period(m) .lt. cepheid_Period_limit) then

                log_R11_likelihood_W_LMC_NGC4258 = log(new_chi2(chi2R11_W_LMC_E14(mu0j(10),M_w,bw,Zw,sigma_int_LMC,m))) + &
                     log(N_tilde_R11_W_LMC(sigma_int_LMC,m)) + log_R11_likelihood_W_LMC_NGC4258
                      
             End If

          End If

       End Do

    Else

       Do m=1,size(Name)

          If (Period(m) .lt. cepheid_Period_limit) then

             log_R11_likelihood_W_LMC_NGC4258 = -chi2R11_W_LMC_E14(mu0j(10),M_w,bw,Zw,sigma_int_LMC,m)/2.d0 + &
                  log(N_tilde_R11_W_LMC(sigma_int_LMC,m)) + log_R11_likelihood_W_LMC_NGC4258

          End If

       End Do

    End If

    Do index_host=1,number_of_hosts_galaxies-1 ! SN Ia

       If (use_HP_in_SNIa) then

          log_R11_likelihood_W_LMC_NGC4258 = log(new_chi2(chi2R11_SNIa(mu0j(index_host),H0,av,index_host))) + &
               log(N_tilde_R11_SNIa(index_host)) + log_R11_likelihood_W_LMC_NGC4258  

       Else

          log_R11_likelihood_W_LMC_NGC4258 = -chi2R11_SNIa(mu0j(index_host),H0,av,index_host)/2.d0 + &
               log(N_tilde_R11_SNIa(index_host)) + log_R11_likelihood_W_LMC_NGC4258

       End If
              
    End Do
        
    If (use_HP_in_av) then ! a_v AND a_cal PARAMETERS

       log_R11_likelihood_W_LMC_NGC4258 = log(new_chi2((a_v - av)**2/sigma_a_v**2)) - log(2.d0*Pi*sigma_a_v**2)/2.d0  + &
            log(new_chi2((a_cal - acal)**2/sigma_a_cal**2)) - log(2.d0*Pi*sigma_a_cal**2)/2.d0 + log_R11_likelihood_W_LMC_NGC4258

    Else

       log_R11_likelihood_W_LMC_NGC4258 = -((a_v - av)**2/sigma_a_v**2 + log(2.d0*Pi*sigma_a_v**2) )/2.d0  - &
            ((a_cal - acal)**2/sigma_a_cal**2 + log(2.d0*Pi*sigma_a_cal**2) )/2.d0 + log_R11_likelihood_W_LMC_NGC4258

    End If
        
    If (use_HP_in_anchor) then ! ANCHOR

       log_R11_likelihood_W_LMC_NGC4258 =  log(new_chi2(chi2R11_anchor_LMC(mu0j(10)))) - log(2.d0*Pi*sigma_mu_0_LMC**2)/2.d0 + &
            log(new_chi2(chi2R11_anchor_NGC4258(mu0j(9)))) - log(2.d0*Pi*sigma_mu_0_NGC4258**2)/2.d0 + &
            log_R11_likelihood_W_LMC_NGC4258

    Else

       log_R11_likelihood_W_LMC_NGC4258 =  -(chi2R11_anchor_LMC(mu0j(10)) + log(2.d0*Pi*sigma_mu_0_LMC**2))/2.d0 - &
            (chi2R11_anchor_NGC4258(mu0j(9)) + log(2.d0*Pi*sigma_mu_0_NGC4258**2))/2.d0 + &
            log_R11_likelihood_W_LMC_NGC4258

    End If

    If (use_prior_on_Zw) then

       If (use_HP_in_Zw) then

          log_R11_likelihood_W_LMC_NGC4258 = log(new_chi2(chi2_Zw(Zw))) - log(2.d0*Pi*sigma_Zw_prior**2)/2.d0 + &
               log_R11_likelihood_W_LMC_NGC4258

       Else

          log_R11_likelihood_W_LMC_NGC4258 = -(chi2_Zw(Zw) + log(2.d0*Pi*sigma_Zw_prior**2) )/2.d0  +&
               log_R11_likelihood_W_LMC_NGC4258

       End If

    Else 

       continue

    End If
        
    If ( abs(log_R11_likelihood_W_LMC_NGC4258) .ge. 0.d0 ) then

       continue

    Else 

       log_R11_likelihood_W_LMC_NGC4258 = -1.d10

    End If

end function log_R11_likelihood_W_LMC_NGC4258

function log_R11_likelihood_W_LMC_MW_NGC4258(mu0j,M_w,bw,H0,Zw,av,acal,sigma_int,sigma_int_LMC)    !    EQUATION (4) IN R09

    use arrays
    use fiducial
    Implicit none

    Real*8 :: log_R11_likelihood_W_LMC_MW_NGC4258,M_w,bw,H0,Zw,av,acal,sigma_int,sigma_int_LMC,normalizationA
    Real*8 :: normalizationMW,chiMW
    Real*8,dimension(number_of_hosts_galaxies + 1) :: mu0j 
    Integer*4 :: m,index_host,number_cepheid

    log_R11_likelihood_W_LMC_MW_NGC4258 = 0.d0

    If (use_HP_per_host) then  ! R11 SAMPLE
 
       If (using_jeffreys_prior) then

          Do index_host=1,number_of_hosts_galaxies

             normalizationA = 0.d0

             number_cepheid = 0

             Do m=1, size(Field)

                If (host(index_host) .eq. Field(m)) then
                       
                   number_cepheid = number_cepheid + 1

                   normalizationA = log( eF160WR11(m)**2 + sigma_int**2 ) + normalizationA

                End If

             End Do

             log_R11_likelihood_W_LMC_MW_NGC4258 = - dble(number_cepheid)*log(chi2R11_W_host_E14(mu0j(index_host),M_w,&
                  bw,Zw,sigma_int,index_host))/2.d0 - normalizationA/2.d0 + log_R11_likelihood_W_LMC_MW_NGC4258

          End Do

       Else

          print *, 'UNIFORM PRIOR NOT IMPLEMENTED YET. NEED TO CODE EXPRESSION WITH GAMMA FUNCTIONS'

          stop

       End If

    Else

       If (use_HP_per_cepheid) then

          Do m=1,size(Field)

             Do index_host=1,number_of_hosts_galaxies
 
                If (host(index_host) .eq. Field(m)) then
    
                   If (using_jeffreys_prior) then

                      print *, 'IMPROPER JEFFREYS PRIOR LEADS TO SINGULARITIES AND THEREFORE IS NOT IMPLEMENTED'

                      stop

                   Else
                       
                      If (PeriodR11(m) .lt. cepheid_Period_limit) then

                         log_R11_likelihood_W_LMC_MW_NGC4258 = log(new_chi2(chi2R11_W_E14(mu0j(index_host),M_w,bw,&
                              Zw,sigma_int,m))) + log(N_tilde_R11_W(sigma_int,m)) + log_R11_likelihood_W_LMC_MW_NGC4258

                      End If

                   End If

                End If

             End Do

          End Do

       Else

          Do m=1,size(Field)

             Do index_host=1,number_of_hosts_galaxies
 
                If (host(index_host) .eq. Field(m)) then
    
                   If (PeriodR11(m) .lt. cepheid_Period_limit) then

                      log_R11_likelihood_W_LMC_MW_NGC4258 = -chi2R11_W_E14(mu0j(index_host),M_w,bw,Zw,sigma_int,m)/2.d0 + &
                           log(N_tilde_R11_W(sigma_int,m)) + log_R11_likelihood_W_LMC_MW_NGC4258

                   End If

                End If

             End Do

          End Do

       End If

    End If

    If (use_HP_per_cepheid) then ! LMC CEPHEID VARIABLES

       Do m=1,size(Name)

          If (using_jeffreys_prior) then

             print *, 'IMPROPER JEFFREYS PRIOR LEADS TO SINGULARITIES AND THEREFORE IS NOT IMPLEMENTED'

             stop

          Else
                       
             If (Period(m) .lt. cepheid_Period_limit) then

                log_R11_likelihood_W_LMC_MW_NGC4258 = log(new_chi2(chi2R11_W_LMC_E14(mu0j(10),M_w,bw,Zw,sigma_int_LMC,m))) + &
                     log(N_tilde_R11_W_LMC(sigma_int_LMC,m)) + log_R11_likelihood_W_LMC_MW_NGC4258
                      
             End If

          End If

       End Do

    Else

       Do m=1,size(Name)

          If (Period(m) .lt. cepheid_Period_limit) then

             log_R11_likelihood_W_LMC_MW_NGC4258 = -chi2R11_W_LMC_E14(mu0j(10),M_w,bw,Zw,sigma_int_LMC,m)/2.d0 + &
                  log(N_tilde_R11_W_LMC(sigma_int_LMC,m)) + log_R11_likelihood_W_LMC_MW_NGC4258

          End If

       End Do

    End If

    If (use_HP_per_MW_cepheid) then ! MW CEPHEID VARIABLES

       Do m=1,size(FieldHipp)

          If (using_jeffreys_prior) then

             print *, 'IMPROPER JEFFREYS PRIOR LEADS TO SINGULARITIES AND THEREFORE IS NOT IMPLEMENTED'

             stop

          Else
                       
             If (10**(logP(m)) .lt. cepheid_Period_limit) then

                log_R11_likelihood_W_LMC_MW_NGC4258 = log(new_chi2(chi2R11_W_MW(M_w,bw,Zw,prior_sigma_int_MW,m))) + &
                     log(N_tilde_R11_W_MW(prior_sigma_int_MW,m)) + log_R11_likelihood_W_LMC_MW_NGC4258
                      
             End If

          End If

       End Do

    Else

       If (use_HP_for_MW_dataset) then

             normalizationMW = 0.d0

             chiMW = 0.d0 

             Do m=1, size(FieldHipp)

                normalizationMW = log( sigmaMw(m)**2 + prior_sigma_int_MW**2 ) + normalizationMW

                chiMW = chi2R11_W_MW(M_w,bw,Zw,prior_sigma_int_MW,m) + chiMW

             End Do
             
             log_R11_likelihood_W_LMC_MW_NGC4258 = - dble(size(FieldHipp))*log(chiMW)/2.d0 - normalizationMW/2.d0&
                  + log_R11_likelihood_W_LMC_MW_NGC4258

       Else

          Do m=1,size(FieldHipp)

             If (10**(logP(m)) .lt. cepheid_Period_limit) then

                log_R11_likelihood_W_LMC_MW_NGC4258 = -chi2R11_W_MW(M_w,bw,Zw,prior_sigma_int_MW,m)/2.d0 + &
                     log(N_tilde_R11_W_MW(prior_sigma_int_MW,m)) + log_R11_likelihood_W_LMC_MW_NGC4258

             End If

          End Do

       End If

    End If

    Do index_host=1,number_of_hosts_galaxies-1 ! SN Ia

       If (use_HP_in_SNIa) then

          log_R11_likelihood_W_LMC_MW_NGC4258 = log(new_chi2(chi2R11_SNIa(mu0j(index_host),H0,av,index_host))) + &
               log(N_tilde_R11_SNIa(index_host)) + log_R11_likelihood_W_LMC_MW_NGC4258

       Else

          log_R11_likelihood_W_LMC_MW_NGC4258 = -chi2R11_SNIa(mu0j(index_host),H0,av,index_host)/2.d0 + &
               log(N_tilde_R11_SNIa(index_host)) + log_R11_likelihood_W_LMC_MW_NGC4258

       End If
              
    End Do
        
    If (use_HP_in_av) then ! a_v AND a_cal PARAMETERS

       log_R11_likelihood_W_LMC_MW_NGC4258 = log(new_chi2((a_v - av)**2/sigma_a_v**2)) - log(2.d0*Pi*sigma_a_v**2)/2.d0  + &
            log(new_chi2((a_cal - acal)**2/sigma_a_cal**2)) - log(2.d0*Pi*sigma_a_cal**2)/2.d0 + log_R11_likelihood_W_LMC_MW_NGC4258

    Else

       log_R11_likelihood_W_LMC_MW_NGC4258 = -((a_v - av)**2/sigma_a_v**2 + log(2.d0*Pi*sigma_a_v**2) )/2.d0  - &
            ((a_cal - acal)**2/sigma_a_cal**2 + log(2.d0*Pi*sigma_a_cal**2) )/2.d0 + log_R11_likelihood_W_LMC_MW_NGC4258

    End If
        
    If (use_HP_in_anchor) then ! ANCHOR

       log_R11_likelihood_W_LMC_MW_NGC4258 =  log(new_chi2(chi2R11_anchor_LMC(mu0j(10)))) - log(2.d0*Pi*sigma_mu_0_LMC**2)/2.d0 + &
            log(new_chi2(chi2R11_anchor_NGC4258(mu0j(9)))) - log(2.d0*Pi*sigma_mu_0_NGC4258**2)/2.d0 + &
            log_R11_likelihood_W_LMC_MW_NGC4258

    Else

       log_R11_likelihood_W_LMC_MW_NGC4258 =  -(chi2R11_anchor_LMC(mu0j(10)) + log(2.d0*Pi*sigma_mu_0_LMC**2))/2.d0 - &
            (chi2R11_anchor_NGC4258(mu0j(9)) + log(2.d0*Pi*sigma_mu_0_NGC4258**2))/2.d0 + &
            log_R11_likelihood_W_LMC_MW_NGC4258

    End If

    If (use_prior_on_Zw) then

       log_R11_likelihood_W_LMC_MW_NGC4258 = -((Zw - prior_Zw)**2/sigma_Zw_prior**2 + log(2.d0*Pi*sigma_Zw_prior**2) )/2.d0  +&
            log_R11_likelihood_W_LMC_MW_NGC4258

    Else 

       continue

    End If
        
    If ( abs(log_R11_likelihood_W_LMC_MW_NGC4258) .ge. 0.d0 ) then

       continue

    Else 

       log_R11_likelihood_W_LMC_MW_NGC4258 = -1.d10

    End If

end function log_R11_likelihood_W_LMC_MW_NGC4258

function log_likelihood_W_MW_alone(M_w,bw,sigma_int_MW)

    use arrays
    use fiducial
    Implicit none

    Real*8 :: log_likelihood_W_MW_alone,M_w,bw,sigma_int_MW,normalizationA
    Real*8 :: normalizationMW,chiMW
    Integer*4 :: m,index_host,number_cepheid

    log_likelihood_W_MW_alone = 0.d0

    If (use_HP_per_MW_cepheid) then ! MW CEPHEID VARIABLES

       Do m=1,size(FieldHipp)

          If (using_jeffreys_prior) then

             print *, 'IMPROPER JEFFREYS PRIOR LEADS TO SINGULARITIES AND THEREFORE IS NOT IMPLEMENTED'

             stop

          Else
                       
             If (10**(logP(m)) .lt. cepheid_Period_limit) then

                log_likelihood_W_MW_alone = log(new_chi2(chi2R11_W_MW(M_w,bw,0.d0,sigma_int_MW,m))) + &
                     log(N_tilde_R11_W_MW(sigma_int_MW,m)) + log_likelihood_W_MW_alone
                      
             End If

          End If

       End Do

    Else

       If (use_HP_for_MW_dataset) then

             normalizationMW = 0.d0

             chiMW = 0.d0 

             Do m=1, size(FieldHipp)

                normalizationMW = log( sigmaMw(m)**2 + sigma_int_MW**2 ) + normalizationMW

                chiMW = chi2R11_W_MW(M_w,bw,0.d0,sigma_int_MW,m) + chiMW

             End Do
             
             log_likelihood_W_MW_alone = - dble(size(FieldHipp))*log(chiMW)/2.d0 - normalizationMW/2.d0&
                  + log_likelihood_W_MW_alone

       Else

          Do m=1,size(FieldHipp)

             If (10**(logP(m)) .lt. cepheid_Period_limit) then

                log_likelihood_W_MW_alone = -chi2R11_W_MW(M_w,bw,0.d0,sigma_int_MW,m)/2.d0 + &
                     log(N_tilde_R11_W_MW(sigma_int_MW,m)) + log_likelihood_W_MW_alone

             End If

          End Do

       End If

    End If
        
    If ( abs(log_likelihood_W_MW_alone) .ge. 0.d0 ) then

       continue

    Else 

       log_likelihood_W_MW_alone = -1.d10

    End If

end function log_likelihood_W_MW_alone

function log_likelihood_W_MW_alone_2(M_w,bw,Zw,sigma_int_MW)

    use arrays
    use fiducial
    Implicit none

    Real*8 :: log_likelihood_W_MW_alone_2,M_w,bw,Zw,sigma_int_MW,normalizationA
    Real*8 :: normalizationMW,chiMW
    Integer*4 :: m,index_host,number_cepheid

    log_likelihood_W_MW_alone_2 = 0.d0

    If (use_HP_per_MW_cepheid) then ! MW CEPHEID VARIABLES

       Do m=1,size(FieldHipp)

          If (using_jeffreys_prior) then

             print *, 'IMPROPER JEFFREYS PRIOR LEADS TO SINGULARITIES AND THEREFORE IS NOT IMPLEMENTED'

             stop

          Else
                       
             If (10**(logP(m)) .lt. cepheid_Period_limit) then

                log_likelihood_W_MW_alone_2 = log(new_chi2(chi2R11_W_MW(M_w,bw,Zw,sigma_int_MW,m))) + &
                     log(N_tilde_R11_W_MW(sigma_int_MW,m)) + log_likelihood_W_MW_alone_2
                      
             End If

          End If

       End Do

    Else

       If (use_HP_for_MW_dataset) then

             normalizationMW = 0.d0

             chiMW = 0.d0 

             Do m=1, size(FieldHipp)

                normalizationMW = log( sigmaMw(m)**2 + sigma_int_MW**2 ) + normalizationMW

                chiMW = chi2R11_W_MW(M_w,bw,0.d0,sigma_int_MW,m) + chiMW

             End Do
             
             log_likelihood_W_MW_alone_2 = - dble(size(FieldHipp))*log(chiMW)/2.d0 - normalizationMW/2.d0&
                  + log_likelihood_W_MW_alone_2

       Else

          Do m=1,size(FieldHipp)

             If (10**(logP(m)) .lt. cepheid_Period_limit) then

                log_likelihood_W_MW_alone_2 = -chi2R11_W_MW(M_w,bw,Zw,sigma_int_MW,m)/2.d0 + &
                     log(N_tilde_R11_W_MW(sigma_int_MW,m)) + log_likelihood_W_MW_alone_2

             End If

          End Do

       End If

    End If

    If (use_prior_on_Zw) then

       log_likelihood_W_MW_alone_2 = -((Zw - prior_Zw)**2/sigma_Zw_prior**2 + log(2.d0*Pi*sigma_Zw_prior**2) )/2.d0  +&
            log_likelihood_W_MW_alone_2

    Else 

       continue

    End If

        
    If ( abs(log_likelihood_W_MW_alone_2) .ge. 0.d0 ) then

       continue

    Else 

       log_likelihood_W_MW_alone_2 = -1.d10

    End If

end function log_likelihood_W_MW_alone_2

function log_R11_likelihood_W_LMC_MW(mu0j,M_w,bw,H0,Zw,av,acal,sigma_int,sigma_int_LMC)    !    EQUATION (4) IN R09

    use arrays
    use fiducial
    Implicit none

    Real*8 :: log_R11_likelihood_W_LMC_MW,M_w,bw,H0,Zw,av,acal,sigma_int,sigma_int_LMC,normalizationA
    Real*8,dimension(number_of_hosts_galaxies + 1) :: mu0j 
    Integer*4 :: m,index_host,number_cepheid

    log_R11_likelihood_W_LMC_MW = 0.d0

    If (use_HP_per_host) then  ! R11 SAMPLE
 
       If (using_jeffreys_prior) then

          Do index_host=1,number_of_hosts_galaxies

             normalizationA = 0.d0

             number_cepheid = 0

             Do m=1, size(Field)

                If (host(index_host) .eq. Field(m)) then
                       
                   number_cepheid = number_cepheid + 1

                   normalizationA = log( eF160WR11(m)**2 + sigma_int**2 ) + normalizationA

                End If

             End Do

             log_R11_likelihood_W_LMC_MW = - dble(number_cepheid)*log(chi2R11_W_host_E14(mu0j(index_host),M_w,&
                  bw,Zw,sigma_int,index_host))/2.d0 - normalizationA/2.d0 + log_R11_likelihood_W_LMC_MW

          End Do

       Else

          print *, 'UNIFORM PRIOR NOT IMPLEMENTED YET. NEED TO CODE EXPRESSION WITH GAMMA FUNCTIONS'

          stop

       End If

    Else

       If (use_HP_per_cepheid) then

          Do m=1,size(Field)

             Do index_host=1,number_of_hosts_galaxies
 
                If (host(index_host) .eq. Field(m)) then
    
                   If (using_jeffreys_prior) then

                      print *, 'IMPROPER JEFFREYS PRIOR LEADS TO SINGULARITIES AND THEREFORE IS NOT IMPLEMENTED'

                      stop

                   Else
                       
                      If (PeriodR11(m) .lt. cepheid_Period_limit) then

                         log_R11_likelihood_W_LMC_MW = log(new_chi2(chi2R11_W_E14(mu0j(index_host),M_w,bw,&
                              Zw,sigma_int,m))) + log(N_tilde_R11_W(sigma_int,m)) + log_R11_likelihood_W_LMC_MW

                      End If

                   End If

                End If

             End Do

          End Do

       Else

          Do m=1,size(Field)

             Do index_host=1,number_of_hosts_galaxies
 
                If (host(index_host) .eq. Field(m)) then
    
                   If (PeriodR11(m) .lt. cepheid_Period_limit) then

                      log_R11_likelihood_W_LMC_MW = -chi2R11_W_E14(mu0j(index_host),M_w,bw,Zw,sigma_int,m)/2.d0 + &
                           log(N_tilde_R11_W(sigma_int,m)) + log_R11_likelihood_W_LMC_MW

                   End If

                End If

             End Do

          End Do

       End If

    End If

    If (use_HP_per_cepheid) then ! LMC CEPHEID VARIABLES

       Do m=1,size(Name)

          If (using_jeffreys_prior) then

             print *, 'IMPROPER JEFFREYS PRIOR LEADS TO SINGULARITIES AND THEREFORE IS NOT IMPLEMENTED'

             stop

          Else
                       
             If (Period(m) .lt. cepheid_Period_limit) then

                log_R11_likelihood_W_LMC_MW = log(new_chi2(chi2R11_W_LMC_E14(mu0j(10),M_w,bw,Zw,sigma_int_LMC,m))) + &
                     log(N_tilde_R11_W_LMC(sigma_int_LMC,m)) + log_R11_likelihood_W_LMC_MW
                      
             End If

          End If

       End Do

    Else

       Do m=1,size(Name)

          If (Period(m) .lt. cepheid_Period_limit) then

             log_R11_likelihood_W_LMC_MW = -chi2R11_W_LMC_E14(mu0j(10),M_w,bw,Zw,sigma_int_LMC,m)/2.d0 + &
                  log(N_tilde_R11_W_LMC(sigma_int_LMC,m)) + log_R11_likelihood_W_LMC_MW

          End If

       End Do

    End If

    If (use_HP_per_cepheid) then ! MW CEPHEID VARIABLES

       Do m=1,size(FieldHipp)

          If (using_jeffreys_prior) then

             print *, 'IMPROPER JEFFREYS PRIOR LEADS TO SINGULARITIES AND THEREFORE IS NOT IMPLEMENTED'

             stop

          Else
                       
             If (10**(logP(m)) .lt. cepheid_Period_limit) then

                log_R11_likelihood_W_LMC_MW = log(new_chi2(chi2R11_W_MW(M_w,bw,Zw,prior_sigma_int_MW,m))) + &
                     log(N_tilde_R11_W_MW(prior_sigma_int_MW,m)) + log_R11_likelihood_W_LMC_MW
                      
             End If

          End If

       End Do

    Else

       Do m=1,size(FieldHipp)

          If (10**(logP(m)) .lt. cepheid_Period_limit) then

             log_R11_likelihood_W_LMC_MW = -chi2R11_W_MW(M_w,bw,Zw,prior_sigma_int_MW,m)/2.d0 + &
                  log(N_tilde_R11_W_MW(prior_sigma_int_MW,m)) + log_R11_likelihood_W_LMC_MW

          End If

       End Do

    End If

    Do index_host=1,number_of_hosts_galaxies-1 ! SN Ia

       If (use_HP_in_SNIa) then

          log_R11_likelihood_W_LMC_MW = log(new_chi2(chi2R11_SNIa(mu0j(index_host),H0,av,index_host))) + &
               log(N_tilde_R11_SNIa(index_host)) + log_R11_likelihood_W_LMC_MW

       Else

          log_R11_likelihood_W_LMC_MW = -chi2R11_SNIa(mu0j(index_host),H0,av,index_host)/2.d0 + &
               log(N_tilde_R11_SNIa(index_host)) + log_R11_likelihood_W_LMC_MW

       End If
              
    End Do
        
    If (use_HP_in_av) then ! a_v AND a_cal PARAMETERS

       log_R11_likelihood_W_LMC_MW = log(new_chi2((a_v - av)**2/sigma_a_v**2)) - log(2.d0*Pi*sigma_a_v**2)/2.d0  + &
            log(new_chi2((a_cal - acal)**2/sigma_a_cal**2)) - log(2.d0*Pi*sigma_a_cal**2)/2.d0 + log_R11_likelihood_W_LMC_MW

    Else

       log_R11_likelihood_W_LMC_MW = -((a_v - av)**2/sigma_a_v**2 + log(2.d0*Pi*sigma_a_v**2) )/2.d0  - &
            ((a_cal - acal)**2/sigma_a_cal**2 + log(2.d0*Pi*sigma_a_cal**2) )/2.d0 + log_R11_likelihood_W_LMC_MW

    End If
        
    If (use_HP_in_anchor) then ! ANCHOR

       log_R11_likelihood_W_LMC_MW =  log(new_chi2(chi2R11_anchor_LMC(mu0j(10)))) - log(2.d0*Pi*sigma_mu_0_LMC**2)/2.d0 + &
            log_R11_likelihood_W_LMC_MW

    Else

       log_R11_likelihood_W_LMC_MW =  -(chi2R11_anchor_LMC(mu0j(10)) + log(2.d0*Pi*sigma_mu_0_LMC**2))/2.d0 + &
            log_R11_likelihood_W_LMC_MW

    End If

    If (use_prior_on_Zw) then

       log_R11_likelihood_W_LMC_MW = -((Zw - prior_Zw)**2/sigma_Zw_prior**2 + log(2.d0*Pi*sigma_Zw_prior**2) )/2.d0  +&
            log_R11_likelihood_W_LMC_MW

    Else 

       continue

    End If
        
    If ( abs(log_R11_likelihood_W_LMC_MW) .ge. 0.d0 ) then

       continue

    Else 

       log_R11_likelihood_W_LMC_MW = -1.d10

    End If

end function log_R11_likelihood_W_LMC_MW

function log_R11_likelihood_W_MW(mu0j,M_w,bw,H0,Zw,av,acal,sigma_int,sigma_int_MW)    !    EQUATION (4) IN R09

    use arrays
    use fiducial
    Implicit none

    Real*8 :: log_R11_likelihood_W_MW,M_w,bw,H0,Zw,av,acal,sigma_int,sigma_int_MW,normalizationA
    Real*8,dimension(number_of_hosts_galaxies) :: mu0j 
    Integer*4 :: m,index_host,number_cepheid

    log_R11_likelihood_W_MW = 0.d0

    If (use_HP_per_host) then  ! R11 SAMPLE
 
       If (using_jeffreys_prior) then

          Do index_host=1,number_of_hosts_galaxies

             normalizationA = 0.d0

             number_cepheid = 0

             Do m=1, size(Field)

                If (host(index_host) .eq. Field(m)) then
                       
                   number_cepheid = number_cepheid + 1

                   normalizationA = log( eF160WR11(m)**2 + sigma_int**2 ) + normalizationA

                End If

             End Do

             log_R11_likelihood_W_MW = - dble(number_cepheid)*log(chi2R11_W_host_E14(mu0j(index_host),&
                  M_w,bw,Zw,sigma_int,index_host))/2.d0 - normalizationA/2.d0 + log_R11_likelihood_W_MW

          End Do

       Else

          print *, 'UNIFORM PRIOR NOT IMPLEMENTED YET. NEED TO CODE EXPRESSION WITH GAMMA FUNCTIONS'

          stop

       End If

    Else

       If (use_HP_per_cepheid) then

          Do m=1,size(Field)

             Do index_host=1,number_of_hosts_galaxies
 
                If (host(index_host) .eq. Field(m)) then
    
                   If (using_jeffreys_prior) then

                      print *, 'IMPROPER JEFFREYS PRIOR LEADS TO SINGULARITIES AND THEREFORE IS NOT IMPLEMENTED'

                      stop

                   Else
                       
                      If (PeriodR11(m) .lt. cepheid_Period_limit) then

                         log_R11_likelihood_W_MW = log(new_chi2(chi2R11_W_E14(mu0j(index_host),M_w,bw,&
                              Zw,sigma_int,m))) + log(N_tilde_R11_W(sigma_int,m)) + log_R11_likelihood_W_MW

                      End If

                   End If

                End If

             End Do

          End Do

       Else

          Do m=1,size(Field)

             Do index_host=1,number_of_hosts_galaxies
 
                If (host(index_host) .eq. Field(m)) then
    
                   If (PeriodR11(m) .lt. cepheid_Period_limit) then

                      log_R11_likelihood_W_MW = -chi2R11_W_E14(mu0j(index_host),M_w,bw,Zw,sigma_int,m)/2.d0 + &
                           log(N_tilde_R11_W(sigma_int,m)) + log_R11_likelihood_W_MW

                   End If

                End If

             End Do

          End Do

       End If

    End If

    If (use_HP_per_MW_cepheid) then ! MW CEPHEID VARIABLES

       Do m=1,size(FieldHipp)

          If (using_jeffreys_prior) then

             print *, 'IMPROPER JEFFREYS PRIOR LEADS TO SINGULARITIES AND THEREFORE IS NOT IMPLEMENTED'

             stop

          Else
                       
             If (10**(logP(m)) .lt. cepheid_Period_limit) then

                log_R11_likelihood_W_MW = log(new_chi2(chi2R11_W_MW(M_w,bw,Zw,sigma_int_MW,m))) + &
                     log(N_tilde_R11_W_MW(sigma_int_MW,m)) + log_R11_likelihood_W_MW
                      
             End If

          End If

       End Do

    Else

       Do m=1,size(FieldHipp)

          If (10**(logP(m)) .lt. cepheid_Period_limit) then

             log_R11_likelihood_W_MW = -chi2R11_W_MW(M_w,bw,Zw,sigma_int_MW,m)/2.d0 + &
                  log(N_tilde_R11_W_MW(sigma_int_MW,m)) + log_R11_likelihood_W_MW

          End If

       End Do

    End If

    Do index_host=1,number_of_hosts_galaxies-1 ! SN Ia

       If (use_HP_in_SNIa) then

          log_R11_likelihood_W_MW = log(new_chi2(chi2R11_SNIa(mu0j(index_host),H0,av,index_host))) + &
               log(N_tilde_R11_SNIa(index_host)) + log_R11_likelihood_W_MW

       Else

          log_R11_likelihood_W_MW = -chi2R11_SNIa(mu0j(index_host),H0,av,index_host)/2.d0 + &
               log(N_tilde_R11_SNIa(index_host)) + log_R11_likelihood_W_MW

       End If
              
    End Do
        
    If (use_HP_in_av) then ! a_v AND a_cal PARAMETERS

       log_R11_likelihood_W_MW = log(new_chi2((a_v - av)**2/sigma_a_v**2)) - log(2.d0*Pi*sigma_a_v**2)/2.d0  + &
            log(new_chi2((a_cal - acal)**2/sigma_a_cal**2)) - log(2.d0*Pi*sigma_a_cal**2)/2.d0 + log_R11_likelihood_W_MW

    Else

       log_R11_likelihood_W_MW = -((a_v - av)**2/sigma_a_v**2 + log(2.d0*Pi*sigma_a_v**2) )/2.d0  - &
            ((a_cal - acal)**2/sigma_a_cal**2 + log(2.d0*Pi*sigma_a_cal**2) )/2.d0 + log_R11_likelihood_W_MW

    End If
        
    If (use_prior_on_Zw) then

       log_R11_likelihood_W_MW = -((Zw - prior_Zw)**2/sigma_Zw_prior**2 + log(2.d0*Pi*sigma_Zw_prior**2) )/2.d0  +&
            log_R11_likelihood_W_MW

    Else 

       continue

    End If
        
    If (use_prior_on_bw) then

       log_R11_likelihood_W_MW = -((bw - prior_bw_from_LMC)**2/sigma_bw_prior**2 + log(2.d0*Pi*sigma_bw_prior**2) )/2.d0  +&
            log_R11_likelihood_W_MW

    Else 

       continue

    End If

    If ( abs(log_R11_likelihood_W_MW) .ge. 0.d0 ) then

       continue

    Else 

       log_R11_likelihood_W_MW = -1.d10

    End If

end function log_R11_likelihood_W_MW

function log_R11_likelihood_W_MW_NGC4258(mu0j,M_w,bw,H0,Zw,av,acal,sigma_int,sigma_int_MW)    !    EQUATION (4) IN R09

    use arrays
    use fiducial
    Implicit none

    Real*8 :: log_R11_likelihood_W_MW_NGC4258,M_w,bw,H0,Zw,av,acal,sigma_int,sigma_int_MW,normalizationA
    Real*8,dimension(number_of_hosts_galaxies) :: mu0j 
    Integer*4 :: m,index_host,number_cepheid

    log_R11_likelihood_W_MW_NGC4258 = 0.d0

    If (use_HP_per_host) then  ! R11 SAMPLE
 
       If (using_jeffreys_prior) then

          Do index_host=1,number_of_hosts_galaxies

             normalizationA = 0.d0

             number_cepheid = 0

             Do m=1, size(Field)

                If (host(index_host) .eq. Field(m)) then
                       
                   number_cepheid = number_cepheid + 1

                   normalizationA = log( eF160WR11(m)**2 + sigma_int**2 ) + normalizationA

                End If

             End Do

             log_R11_likelihood_W_MW_NGC4258 = - dble(number_cepheid)*log(chi2R11_W_host_E14(mu0j(index_host),&
                  M_w,bw,Zw,sigma_int,index_host))/2.d0 - normalizationA/2.d0 + log_R11_likelihood_W_MW_NGC4258

          End Do

       Else

          print *, 'UNIFORM PRIOR NOT IMPLEMENTED YET. NEED TO CODE EXPRESSION WITH GAMMA FUNCTIONS'

          stop

       End If

    Else

       If (use_HP_per_cepheid) then

          Do m=1,size(Field)

             Do index_host=1,number_of_hosts_galaxies
 
                If (host(index_host) .eq. Field(m)) then
    
                   If (using_jeffreys_prior) then

                      print *, 'IMPROPER JEFFREYS PRIOR LEADS TO SINGULARITIES AND THEREFORE IS NOT IMPLEMENTED'

                      stop

                   Else
                       
                      If (PeriodR11(m) .lt. cepheid_Period_limit) then

                         log_R11_likelihood_W_MW_NGC4258 = log(new_chi2(chi2R11_W_E14(mu0j(index_host),M_w,bw,&
                              Zw,sigma_int,m))) + log(N_tilde_R11_W(sigma_int,m)) + log_R11_likelihood_W_MW_NGC4258

                      End If

                   End If

                End If

             End Do

          End Do

       Else

          Do m=1,size(Field)

             Do index_host=1,number_of_hosts_galaxies
 
                If (host(index_host) .eq. Field(m)) then
    
                   If (PeriodR11(m) .lt. cepheid_Period_limit) then

                      log_R11_likelihood_W_MW_NGC4258 = -chi2R11_W_E14(mu0j(index_host),M_w,bw,Zw,sigma_int,m)/2.d0 + &
                           log(N_tilde_R11_W(sigma_int,m)) + log_R11_likelihood_W_MW_NGC4258

                   End If

                End If

             End Do

          End Do

       End If

    End If

    If (use_HP_per_cepheid) then ! MW CEPHEID VARIABLES

       Do m=1,size(FieldHipp)

          If (using_jeffreys_prior) then

             print *, 'IMPROPER JEFFREYS PRIOR LEADS TO SINGULARITIES AND THEREFORE IS NOT IMPLEMENTED'

             stop

          Else
                       
             If (10**(logP(m)) .lt. cepheid_Period_limit) then

                log_R11_likelihood_W_MW_NGC4258 = log(new_chi2(chi2R11_W_MW(M_w,bw,Zw,sigma_int_MW,m))) + &
                     log(N_tilde_R11_W_MW(sigma_int_MW,m)) + log_R11_likelihood_W_MW_NGC4258
                      
             End If

          End If

       End Do

    Else

       Do m=1,size(FieldHipp)

          If (10**(logP(m)) .lt. cepheid_Period_limit) then

             log_R11_likelihood_W_MW_NGC4258 = -chi2R11_W_MW(M_w,bw,Zw,sigma_int_MW,m)/2.d0 + &
                  log(N_tilde_R11_W_MW(sigma_int_MW,m)) + log_R11_likelihood_W_MW_NGC4258

          End If

       End Do

    End If

    Do index_host=1,number_of_hosts_galaxies-1 ! SN Ia

       If (use_HP_in_SNIa) then

          log_R11_likelihood_W_MW_NGC4258 = log(new_chi2(chi2R11_SNIa(mu0j(index_host),H0,av,index_host))) + &
               log(N_tilde_R11_SNIa(index_host)) + log_R11_likelihood_W_MW_NGC4258

       Else

          log_R11_likelihood_W_MW_NGC4258 = -chi2R11_SNIa(mu0j(index_host),H0,av,index_host)/2.d0 + &
               log(N_tilde_R11_SNIa(index_host)) + log_R11_likelihood_W_MW_NGC4258

       End If
              
    End Do
        
    If (use_HP_in_av) then ! a_v AND a_cal PARAMETERS

       log_R11_likelihood_W_MW_NGC4258 = log(new_chi2((a_v - av)**2/sigma_a_v**2)) - log(2.d0*Pi*sigma_a_v**2)/2.d0  + &
            log(new_chi2((a_cal - acal)**2/sigma_a_cal**2)) - log(2.d0*Pi*sigma_a_cal**2)/2.d0 + log_R11_likelihood_W_MW_NGC4258

    Else

       log_R11_likelihood_W_MW_NGC4258 = -((a_v - av)**2/sigma_a_v**2 + log(2.d0*Pi*sigma_a_v**2) )/2.d0  - &
            ((a_cal - acal)**2/sigma_a_cal**2 + log(2.d0*Pi*sigma_a_cal**2) )/2.d0 + log_R11_likelihood_W_MW_NGC4258

    End If

    If (use_HP_in_anchor) then

       log_R11_likelihood_W_MW_NGC4258 =  log(new_chi2(chi2R11_anchor_NGC4258(mu0j(9)))) - &
            log(2.d0*Pi*sigma_mu_0_NGC4258**2)/2.d0 + log_R11_likelihood_W_MW_NGC4258

    Else

       log_R11_likelihood_W_MW_NGC4258 =  -(chi2R11_anchor_NGC4258(mu0j(9)) + log(2.d0*Pi*sigma_mu_0_NGC4258**2))/2.d0 + &
            log_R11_likelihood_W_MW_NGC4258

    End If
        
    If (use_prior_on_Zw) then

       log_R11_likelihood_W_MW_NGC4258 = -((Zw - prior_Zw)**2/sigma_Zw_prior**2 + log(2.d0*Pi*sigma_Zw_prior**2) )/2.d0  +&
            log_R11_likelihood_W_MW_NGC4258

    Else 

       continue

    End If
        
    If ( abs(log_R11_likelihood_W_MW_NGC4258) .ge. 0.d0 ) then

       continue

    Else 

       log_R11_likelihood_W_MW_NGC4258 = -1.d10

    End If

end function log_R11_likelihood_W_MW_NGC4258

function log_R11_likelihood_W(mu0j,zpw_ref,bw,H0,Zw,av,sigma_int)    !    EQUATION (4) IN R09

    use arrays
    use fiducial
    Implicit none

    Real*8 :: log_R11_likelihood_W,zpw_ref,bw,H0,Zw,av,sigma_int,normalizationA
    Real*8,dimension(number_of_hosts_galaxies) :: mu0j 
    Integer*4 :: m,index_host,number_cepheid

    log_R11_likelihood_W = 0.d0

    If (use_HP_per_host) then
 
       If (using_jeffreys_prior) then

          Do index_host=1,number_of_hosts_galaxies

             normalizationA = 0.d0

             number_cepheid = 0

             Do m=1, size(Field)

                If (host(index_host) .eq. Field(m)) then
                       
                   number_cepheid = number_cepheid + 1

                   normalizationA = log( eF160WR11(m)**2 + sigma_int**2 ) + normalizationA

                End If

             End Do

             log_R11_likelihood_W = - dble(number_cepheid)*log(chi2R11_W_host(mu0j(index_host),mu0j(9),&
                  zpw_ref,bw,Zw,sigma_int,index_host))/2.d0 - normalizationA/2.d0 + log_R11_likelihood_W

          End Do

       Else

          print *, 'UNIFORM PRIOR NOT IMPLEMENTED YET. NEED TO CODE EXPRESSION WITH GAMMA FUNCTIONS'

          stop

       End If

    Else

       If (use_HP_per_cepheid) then

          Do m=1,size(Field)

             Do index_host=1,number_of_hosts_galaxies
 
                If (host(index_host) .eq. Field(m)) then
    
                   If (using_jeffreys_prior) then

                      print *, 'IMPROPER JEFFREYS PRIOR LEADS TO SINGULARITIES AND THEREFORE IS NOT IMPLEMENTED'

                      stop

                   Else
                       
                      If (PeriodR11(m) .lt. cepheid_Period_limit) then

                         log_R11_likelihood_W = log(new_chi2(chi2R11_W(mu0j(index_host),mu0j(9),zpw_ref,bw,Zw,sigma_int,m))) + &
                              log(N_tilde_R11_W(sigma_int,m)) + log_R11_likelihood_W

                      End If

                   End If

                End If

             End Do

          End Do

       Else

          Do m=1,size(Field)

             Do index_host=1,number_of_hosts_galaxies
 
                If (host(index_host) .eq. Field(m)) then
    
                   If (PeriodR11(m) .lt. cepheid_Period_limit) then

                      log_R11_likelihood_W = -chi2R11_W(mu0j(index_host),mu0j(9),zpw_ref,bw,Zw,sigma_int,m)/2.d0 + &
                           log(N_tilde_R11_W(sigma_int,m)) + log_R11_likelihood_W

                   End If

                End If

             End Do

          End Do

       End If

    End If

    Do index_host=1,number_of_hosts_galaxies-1

       If (use_HP_in_SNIa) then

          log_R11_likelihood_W = log(new_chi2(chi2R11_SNIa(mu0j(index_host),H0,av,index_host))) + &
               log(N_tilde_R11_SNIa(index_host)) + log_R11_likelihood_W  

       Else

          log_R11_likelihood_W = -chi2R11_SNIa(mu0j(index_host),H0,av,index_host)/2.d0 + &
               log(N_tilde_R11_SNIa(index_host)) + log_R11_likelihood_W  

       End If
              
    End Do
        
    If (use_HP_in_av) then

       log_R11_likelihood_W = log(new_chi2((a_v - av)**2/sigma_a_v**2)) - log(2.d0*Pi*sigma_a_v**2)/2.d0  + log_R11_likelihood_W

    Else

       log_R11_likelihood_W = -((a_v - av)**2/sigma_a_v**2 + log(2.d0*Pi*sigma_a_v**2) )/2.d0  + log_R11_likelihood_W

    End If
        
    If (use_HP_in_anchor) then

       log_R11_likelihood_W =  log(new_chi2(chi2R11_anchor_NGC4258(mu0j(9)))) - log(2.d0*Pi*sigma_mu_0_NGC4258**2)/2.d0 + &
            log_R11_likelihood_W

    Else

       log_R11_likelihood_W =  -(chi2R11_anchor_NGC4258(mu0j(9)) + log(2.d0*Pi*sigma_mu_0_NGC4258**2))/2.d0 + &
            log_R11_likelihood_W

    End If

    If (use_prior_on_zpw4258) then

       log_R11_likelihood_W = -((zpw_ref - prior_zpw4258)**2/sigma_zpw4258**2 + log(2.d0*Pi*sigma_zpw4258**2) )/2.d0  +&
            log_R11_likelihood_W

    Else 

       continue

    End If

    If (use_prior_on_Zw) then

       log_R11_likelihood_W = -((Zw - prior_Zw)**2/sigma_Zw_prior**2 + log(2.d0*Pi*sigma_Zw_prior**2) )/2.d0  +&
            log_R11_likelihood_W

    Else 

       continue

    End If

    If (use_prior_on_bw) then

       log_R11_likelihood_W = -((bw - prior_bw_from_LMC)**2/sigma_bw_prior**2 + log(2.d0*Pi*sigma_bw_prior**2) )/2.d0  +&
            log_R11_likelihood_W

    Else 

       continue

    End If
        
    If ( abs(log_R11_likelihood_W) .ge. 0.d0 ) then

       continue

    Else 

       log_R11_likelihood_W = -1.d10

    End If

end function log_R11_likelihood_W

function log_R11_likelihood_W_cepheids(mu0j,zpw_ref,bw,Zw,sigma_int)    !    EQUATION (4) IN R09

    use arrays
    use fiducial
    Implicit none

    Real*8 :: log_R11_likelihood_W_cepheids,zpw_ref,bw,Zw,sigma_int
    Real*8,dimension(number_of_hosts_galaxies) :: mu0j 
    Integer*4 :: m,index_host


    log_R11_likelihood_W_cepheids = 0.d0

    Do m=1,size(Field)

       Do index_host=1,number_of_hosts_galaxies
 
          If (host(index_host) .eq. Field(m)) then
    
             If (using_jeffreys_prior) then

                print *, 'IMPROPER JEFFREYS PRIOR LEADS TO SINGULARITIES AND THEREFORE IS NOT IMPLEMENTED'

                stop

             Else
                       
                log_R11_likelihood_W_cepheids = log(new_chi2(chi2R11_W(mu0j(index_host),mu0j(9),zpw_ref,bw,&
                     Zw,sigma_int,m))) + log(N_tilde_R11_W(sigma_int,m)) + log_R11_likelihood_W_cepheids

             End If

          End If

       End Do

    End Do
        
    If ( abs(log_R11_likelihood_W_cepheids) .ge. 0.d0 ) then

       continue

    Else 

       log_R11_likelihood_W_cepheids = -1.d10

    End If

end function log_R11_likelihood_W_cepheids

function log_likelihood_only_cepheids(galaxy,zpw,bw,Zw,sigma_int)    !    EQUATION (4) IN R09

    use arrays
    use fiducial
    Implicit none

    Character(len=5) :: galaxy
    Real*8 :: log_likelihood_only_cepheids,zpw,bw,Zw,sigma_int
    Integer*4 :: m,index_host


    log_likelihood_only_cepheids = 0.d0

    Do m=1,size(Field)

       If (galaxy .eq. Field(m)) then

          Do index_host=1,number_of_hosts_galaxies
 
             If (host(index_host) .eq. Field(m)) then
    
                If (using_jeffreys_prior) then

                   print *, 'IMPROPER JEFFREYS PRIOR LEADS TO SINGULARITIES AND THEREFORE IS NOT IMPLEMENTED'

                   stop

                Else
                       
                   log_likelihood_only_cepheids = log(new_chi2(chi2R11_W_2(zpw,bw,Zw,sigma_int,m))) + &
                     log(N_tilde_R11_W(sigma_int,m)) + log_likelihood_only_cepheids

                End If

             End If

          End Do
          
       End If

    End Do

    If ( abs(log_likelihood_only_cepheids) .ge. 0.d0 ) then

       continue

    Else 

       log_likelihood_only_cepheids = -1.d10

    End If

end function log_likelihood_only_cepheids

function chi2R11_W(mu0_j,mu0_ref,zpw_ref,bw,Zw,sigma_int,m)    !    It computes equation (3) in published version of 1311.3461

    use arrays
    use fiducial
    Implicit none

    Real*8 :: mu0_j,mu0_ref,zpw_ref,bw,Zw,sigma_int,chi2R11_W
    Integer*4 :: m

    chi2R11_W = ( observed_m_W(F160WR11(m),VIR11(m)) - &
         P_L_relation_passband_W(mu0_j,mu0_ref,zpw_ref,bw,Zw,OHR11(m),PeriodR11(m)) )**2/( eF160WR11(m)**2 + sigma_int**2 ) 

end function chi2R11_W

function chi2R11_W_E14(mu0_j,M_w,bw,Zw,sigma_int,m)    !    It computes equation (3) in published version of 1311.3461

    use arrays
    use fiducial
    Implicit none

    Real*8 :: mu0_j,M_w,bw,Zw,sigma_int,chi2R11_W_E14
    Integer*4 :: m

    chi2R11_W_E14 = ( observed_m_W(F160WR11(m),VIR11(m)) - &
         P_L_relation_passband_W_E14(mu0_j,M_w,bw,Zw,OHR11(m),PeriodR11(m)) )**2/( eF160WR11(m)**2 + sigma_int**2 ) 

end function chi2R11_W_E14

function chi2R11_H_2(mu0_j,mu0_ref,zpw_ref,bw,Zw,sigma_int,m)    !    It computes equation (3) in published version of 1311.3461

    use arrays
    use fiducial
    Implicit none

    Real*8 :: mu0_j,mu0_ref,zpw_ref,bw,Zw,sigma_int,chi2R11_H_2
    Integer*4 :: m

    chi2R11_H_2 = ( F160WR11(m) - &
         P_L_relation_passband_H_2(mu0_j,mu0_ref,zpw_ref,bw,Zw,OHR11(m),PeriodR11(m)) )**2/( eF160WR11(m)**2 + sigma_int**2 ) 

end function chi2R11_H_2

function chi2R11_W_LMC(mu0_j,mu0_ref,zpw_ref,bw,Zw,sigma_int,m)    !    It computes equation (3) in published version of 1311.3461

    use arrays
    use fiducial
    Implicit none

    Real*8 :: mu0_j,mu0_ref,zpw_ref,bw,Zw,sigma_int,chi2R11_W_LMC
    Integer*4 :: m

    chi2R11_W_LMC = ( observed_m_W(H(m),V(m)-II(m)) - &
         P_L_relation_passband_W(mu0_j,mu0_ref,zpw_ref,bw,Zw,meanOH_LMC,Period(m)) )**2/( Sigma_m(m)**2 + sigma_int**2 ) 

end function chi2R11_W_LMC

function chi2R11_W_LMC_E14(mu0_j,M_w,bw,Zw,sigma_int,m)    !    It computes equation (3) in published version of 1311.3461

    use arrays
    use fiducial
    Implicit none

    Real*8 :: mu0_j,M_w,bw,Zw,sigma_int,chi2R11_W_LMC_E14
    Integer*4 :: m

    chi2R11_W_LMC_E14 = ( observed_m_W(H(m),V(m)-II(m)) - &
         P_L_relation_passband_W_E14(mu0_j,M_w,bw,Zw,meanOH_LMC,Period(m)) )**2/( Sigma_m(m)**2 + sigma_int**2 ) 

end function chi2R11_W_LMC_E14

function chi2R11_W_MW(M_w,bw,Zw,sigma_int,m)    !    It computes equation (3) in published version of 1311.3461

    use arrays
    use fiducial
    Implicit none

    Real*8 :: M_w,bw,Zw,sigma_int,chi2R11_W_MW
    Integer*4 :: m

    chi2R11_W_MW = ( Mw(m) - P_L_relation_passband_W_E14(0.d0,M_w,bw,Zw,meanOH_MW,&
         10**(logP(m))) )**2/( sigmaMw(m)**2 + sigma_int**2 ) 

end function chi2R11_W_MW

function chi2R11_W_host(mu0_j,mu0_ref,zpw_ref,bw,Zw,sigma_int,m)    !    It computes equation (3) in published version of 1311.3461

    use arrays
    use fiducial
    Implicit none

    Real*8 :: mu0_j,mu0_ref,zpw_ref,bw,Zw,sigma_int,chi2R11_W_host
    Integer*4 :: m,n

    chi2R11_W_host = 0.d0

    Do n=1, size(Field)

       If (host(m) .eq. Field(n)) then

          chi2R11_W_host = ( observed_m_W(F160WR11(n),VIR11(n)) - &
            P_L_relation_passband_W(mu0_j,mu0_ref,zpw_ref,bw,Zw,OHR11(n),PeriodR11(n)) )**2/( eF160WR11(n)**2 + sigma_int**2 ) + &
            chi2R11_W_host

       End If

    End Do

end function chi2R11_W_host

function chi2R11_W_host_E14(mu0_j,M_w,bw,Zw,sigma_int,m)    !    It computes equation (3) in published version of 1311.3461

    use arrays
    use fiducial
    Implicit none

    Real*8 :: mu0_j,M_w,bw,Zw,sigma_int,chi2R11_W_host_E14
    Integer*4 :: m,n

    chi2R11_W_host_E14 = 0.d0

    Do n=1, size(Field)

       If (host(m) .eq. Field(n)) then

          chi2R11_W_host_E14 = ( observed_m_W(F160WR11(n),VIR11(n)) - &
            P_L_relation_passband_W_E14(mu0_j,M_w,bw,Zw,OHR11(n),PeriodR11(n)) )**2/( eF160WR11(n)**2 + sigma_int**2 ) + &
            chi2R11_W_host_E14

       End If

    End Do

end function chi2R11_W_host_E14

function chi2R11_W_2(zpw,bw,Zw,sigma_int,m)    !    It computes equation (3) in published version of 1311.3461

    use arrays
    use fiducial
    Implicit none

    Real*8 :: zpw,bw,Zw,sigma_int,chi2R11_W_2
    Integer*4 :: m

    chi2R11_W_2 = ( observed_m_W(F160WR11(m),VIR11(m)) - &
         P_L_relation_passband_W_2(zpw,bw,Zw,OHR11(m),PeriodR11(m)) )**2/( eF160WR11(m)**2 + sigma_int**2 ) 

end function chi2R11_W_2

function chi2R11_SNIa(mu0_j,H0,av,snia)    !    It computes equation (3) in published version of 1311.3461

    use arrays
    use fiducial
    Implicit none

    Real*8 :: mu0_j,H0,av,chi2R11_SNIa
    Integer*4 :: snia

    chi2R11_SNIa = ( mvi5av(snia) - 5.d0*a_v - reddening_free_magnitude_SNIa(mu0_j,H0,av) )**2/&
         (Sigma_mvi5av(snia)**2 + (5.d0*sigma_a_v)**2)

end function chi2R11_SNIa

function chi2R11_anchor_NGC4258(mu09)    !    It computes equation (3) in published version of 1311.3461

    use arrays
    use fiducial
    Implicit none

    Real*8 :: mu09,chi2R11_anchor_NGC4258

    chi2R11_anchor_NGC4258 = ( mu09 - mu_0_NGC4258)**2/sigma_mu_0_NGC4258**2

end function chi2R11_anchor_NGC4258

function chi2_Zw(Zw)    !    It computes equation (3) in published version of 1311.3461

    use fiducial

    Implicit none

    Real*8 :: Zw,chi2_Zw

    chi2_Zw = ( Zw - prior_Zw)**2/sigma_Zw_prior**2

end function chi2_Zw

function chi2R11_anchor_LMC(mu09)    !    It computes equation (3) in published version of 1311.3461

    use arrays
    use fiducial
    Implicit none

    Real*8 :: mu09,chi2R11_anchor_LMC

    chi2R11_anchor_LMC = ( mu09 - mu_0_LMC)**2/sigma_mu_0_LMC**2

end function chi2R11_anchor_LMC

function N_tilde_R11_W(sigma_int,m)    !    It computes equation (3) in published version of 1311.3461

    use arrays
    use fiducial
    Implicit none

    Real*8 :: sigma_int,N_tilde_R11_W
    Integer*4 :: m

    N_tilde_R11_W = 1.d0/sqrt( eF160WR11(m)**2 + sigma_int**2 ) 

end function N_tilde_R11_W

function N_tilde_R11_W_LMC(sigma_int,m)    !    It computes equation (3) in published version of 1311.3461

    use arrays
    use fiducial
    Implicit none

    Real*8 :: sigma_int,N_tilde_R11_W_LMC
    Integer*4 :: m

    N_tilde_R11_W_LMC = 1.d0/sqrt( Sigma_m(m)**2 + sigma_int**2 ) 

end function N_tilde_R11_W_LMC

function N_tilde_R11_W_MW(sigma_int,m)    !    It computes equation (3) in published version of 1311.3461

    use arrays
    use fiducial
    Implicit none

    Real*8 :: sigma_int,N_tilde_R11_W_MW
    Integer*4 :: m

    N_tilde_R11_W_MW = 1.d0/sqrt( sigmaMw(m)**2 + sigma_int**2 ) 

end function N_tilde_R11_W_MW

function N_tilde_R11_SNIa(snia)    !    It computes equation (3) in published version of 1311.3461

    use arrays
    use fiducial
    Implicit none

    Real*8 :: N_tilde_R11_SNIa
    Integer*4 :: snia

    N_tilde_R11_SNIa = 1.d0/sqrt(Sigma_mvi5av(snia)**2 + (5.d0*sigma_a_v)**2) 

end function N_tilde_R11_SNIa

function log_likelihood_hyperparameters_as_mcmc(A,bw,sigma_int,oldpoint)    !    It computes equation (3) in published version of 1311.3461
    use arrays
    use fiducial
    Implicit none
    Real*8 :: log_likelihood_hyperparameters_as_mcmc,A,bw,sigma_int
    Integer*4 :: m
    logical :: oldpoint

    If (hyperparameters_as_mcmc) then
    
        log_likelihood_hyperparameters_as_mcmc = 0.d0

        Do m=1,size(Name)
 
            If (using_jeffreys_prior) then

                If (oldpoint) then

                    log_likelihood_hyperparameters_as_mcmc = -log(old_point(number_model_parameters+m))/2.d0 - &
                    old_point(number_model_parameters+m)*chi2_i(A,bw,sigma_int,m)/2.d0 + log(N_tilde_i(sigma_int,m))&
                    + log_likelihood_hyperparameters_as_mcmc

                Else

                    log_likelihood_hyperparameters_as_mcmc = -log(current_point(number_model_parameters+m))/2.d0 - &
                    current_point(number_model_parameters+m)*chi2_i(A,bw,sigma_int,m)/2.d0 + log(N_tilde_i(sigma_int,m))&
                    + log_likelihood_hyperparameters_as_mcmc

                End If

            Else
                
                If (oldpoint) then

                    log_likelihood_hyperparameters_as_mcmc = log(old_point(number_model_parameters+m))/2.d0 - &
                    old_point(number_model_parameters+m)*chi2_i(A,bw,sigma_int,m)/2.d0 + log(N_tilde_i(sigma_int,m))&
                    + log_likelihood_hyperparameters_as_mcmc 

                Else

                    log_likelihood_hyperparameters_as_mcmc = log(current_point(number_model_parameters+m))/2.d0 - &
                    current_point(number_model_parameters+m)*chi2_i(A,bw,sigma_int,m)/2.d0 + log(N_tilde_i(sigma_int,m))&
                    + log_likelihood_hyperparameters_as_mcmc 

                End If

            End If

        End Do

        If ( abs(log_likelihood_hyperparameters_as_mcmc) .ge. 0.d0 ) then

            continue

        Else 

            log_likelihood_hyperparameters_as_mcmc = -1.d10

        End If

    Else

        print *,'CURRENTLY NOT USING HYPER-PARAMETERS AS MCMC PARAMETERS. CHECK FIDUCIAL MODULE '

        stop

    End If

end function log_likelihood_hyperparameters_as_mcmc


function new_chi2(chi2)
    use fiducial
    Implicit none
    Real*8 :: chi2,new_chi2

    If ((chi2 .ge. 0.d0) .and. (chi2 .lt. 1.d-3)) then
  
        new_chi2 = sqrt(2.d0/Pi/9.d0) + 1/sqrt(2.d0*Pi)*(-chi2/5.d0 + chi2**2/28.d0 - chi2**3/216.d0)
  
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

function chi2_i(A,bw,sigma_int,m)    !    It computes equation (3) in published version of 1311.3461
    use arrays
    use fiducial
    Implicit none
    Real*8 :: A,bw,sigma_int,chi2_i
    Integer*4 :: m

    chi2_i = ( observed_wesenheit_magnitude(H(m),V(m),II(m)) - wesenheit_magnitude(A,bw,Period(m)) )**2/&
    ( Sigma_m(m)**2 + sigma_int**2 ) 

end function chi2_i

function N_tilde_A_i(sigma_int,m)    !    It computes equation (3) in published version of 1311.3461
    use arrays
    use fiducial
    Implicit none
    Real*8 :: sigma_int,N_tilde_A_i
    Integer*4 :: m

    N_tilde_A_i = 1.d0/sqrt( Sigma_mA(m)**2 + sigma_int**2 ) 

end function N_tilde_A_i

function N_tilde_i(sigma_int,m)    !    It computes equation (3) in published version of 1311.3461
    use arrays
    use fiducial
    Implicit none
    Real*8 :: sigma_int,N_tilde_i
    Integer*4 :: m

    N_tilde_i = 1.d0/sqrt( Sigma_m(m)**2 + sigma_int**2 ) 

end function N_tilde_i

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
    logical :: oldpoint
  
    If ((include_dataA .and. include_dataB) .and. include_dataC) then
 
        If (hyperparameters_as_mcmc) then

            If (A .eq. old_point(1)) then

                oldpoint = .true.

            Else if (A .eq. current_point(1)) then

                oldpoint = .false.

            Else 

                print *,'BESTFIT WILL BE IMPLEMENTED LATER'

                stop
  
            End If

            log_Efstathiou_likelihood_hyperparameters = log_likelihood_hyperparameters_as_mcmc(A,bw,&
            sigma_int,oldpoint)

        Else

            log_Efstathiou_likelihood_hyperparameters = log_Efstathiou_likelihoodA(A,bw,sigma_int) + &
            log_Efstathiou_likelihoodB(A,bw,sigma_int) + log_Efstathiou_likelihoodC(A,bw,sigma_int)

        End If

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

function compute_determinant(A)
    use fiducial
    Implicit none
    Integer*4 :: INFO,index
    Real*8,dimension(max(1,number_of_parameters),number_of_parameters) :: A
    Integer*4,dimension(min(number_of_parameters,number_of_parameters)) :: IPIV
!    Real*8,dimension(max(1,max(1,number_of_parameters))) :: WORK
    Real*8 :: det,sgn,compute_determinant

    call dgetrf(number_of_parameters,number_of_parameters,A,number_of_parameters,IPIV,INFO)

    det = 1.d0
    Do index=1,number_of_parameters

        det = det*A(index,index)
    End Do 

    sgn = 1.d0
    Do index=1,number_of_parameters
        If (IPIV(index) .ne. index) then
            sgn = -sgn
        End If
    End Do
 
    compute_determinant = sgn*det
end function compute_determinant

function compute_eigenvalues(A)
    use fiducial
    Implicit none
    Integer*4 :: INFO,index
    Integer*4,parameter :: LWORK = max(1,3*number_of_parameters-1)
    Real*8,dimension(max(1,number_of_parameters),number_of_parameters) :: A
    Real*8,dimension(max(1,LWORK)) :: WORK
    Real*8,dimension(number_of_parameters) :: W,compute_eigenvalues
    Character*1,parameter :: JOBZ = 'N'
    Character*1,parameter :: UPLO = 'U'


    call dsyev(JOBZ,UPLO,number_of_parameters,A,number_of_parameters,W,WORK,LWORK,INFO)


    If (INFO .eq. 0) then

        Do index=1,number_of_parameters
         
            compute_eigenvalues(index) = W(index)

        End Do

    Else

        print *,'ERROR COMPUTING EIGENVALUES'

    End If
 
end function compute_eigenvalues

subroutine read_covariance_matrix_mcmc(matrix1)
    use fiducial
    Implicit none
    Real*8,dimension(number_of_parameters,number_of_parameters) :: matrix,matrix1
    Integer*4 :: index1,INFO
    Integer*4,parameter :: LWORK = max(1,3*number_of_parameters-1)
    Real*8,dimension(max(1,LWORK)) :: WORK
    Real*8,dimension(number_of_parameters) :: W
    Character*1,parameter :: JOBZ = 'N'
    Character*1,parameter :: UPLO = 'U'
    Logical :: pos_def 
 
    open(12,file='./output/covariance_matrix.txt')

    Do index1=1,number_of_parameters

        read(12,*) matrix(index1,1:number_of_parameters)

    End Do

    close(12)

    call dsyev(JOBZ,UPLO,number_of_parameters,matrix,number_of_parameters,W,WORK,LWORK,INFO)

    If (INFO .eq. 0) then
 
        pos_def = .true.
        
        Do index1=1,number_of_parameters
         
            If (W(index1) .le. 0.d0) then

                pos_def = .false.

                exit

            End If

        End Do
      
        If (pos_def) then

            open(12,file='./output/covariance_matrix.txt')

            Do index1=1,number_of_parameters

                read(12,*) matrix1(index1,1:number_of_parameters)

            End Do

            close(12)

        Else

            print *,'COVARIANCE MATRIX IS NOT POSITIVE DEFINITE, KEEPING CURRENT COVARIANCE MATRIX'
            
        End If

    Else

        print *,'EIGENVALUES WERE NOT COMPUTED'

    End If

end subroutine read_covariance_matrix_mcmc

subroutine read_bestfit_mcmc(vector)
    use fiducial
    Implicit none
    Real*8,dimension(number_of_parameters) :: vector
    Integer*4 :: index1

    open(12,file='./output/chains/bestfit.txt')

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

    open(12,file='./output/chains/means.txt')

    Do index1=1,number_of_parameters

        read(12,*) vector(index1)

    End Do

    close(12)

end subroutine read_means_mcmc

subroutine set_covariance_matrix()

  use fiducial
  use arrays
  
  Implicit none

  Integer*4 :: m,n

  Logical :: cov_file, out_file

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

     open(UNIT_EXE_FILE,file= EXECUTION_INFORMATION)    ! OPEN FILE FOR EXECUTION INFORMATION 

  Else

     ! SETTING COVARIANCE MATRIX
     Covguess = 0.d0
  
     If (doing_R11_analysis) then  

        If (include_only_cepheids) then

           If (all_R11_hosts) then

              If (number_model_parameters .eq. 12) then
                          
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

                 Covguess(12,12) = sigma_Zw**2 

              Else
                        
                 print *,'WRONG NUMBER OF MODEL PARAMETERS. CHECK FIDUCIAL MODULE AND COMPARE WITH EQUATION (18) IN R09'

                 stop

              End If

           Else

              If (number_model_parameters .eq. 3) then
                          
                 Covguess(1,1) = sigma_zpw**2

                 Covguess(2,2) = sigma_bw**2

                 Covguess(3,3) = sigma_Zw**2

              Else
                        
                 print *,'WRONG NUMBER OF MODEL PARAMETERS. CHECK FIDUCIAL MODULE AND COMPARE WITH EQUATION (7) IN R09'

                 stop

              End If

           End If

        Else

           If ( ( use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. use_MW_as_anchor) then !!!HERE
    
              If (use_metallicity) then 

                 If (use_H_band) then
                     
                    print *,'H BAND NOT IMPLEMENTED USING LMC+MW+NGC4258 AS ANCHORS AND INCLUDING METALLICITY DEPENDENCE'
                     
                    stop
                     
                 Else

                    If (number_model_parameters .eq. 16) then
                          
                       Covguess(1,1) = sigma_mu1**2 

                       Covguess(2,2) = sigma_mu2**2 

                       Covguess(3,3) = sigma_mu3**2 

                       Covguess(4,4) = sigma_mu4**2 

                       Covguess(5,5) = sigma_mu5**2 

                       Covguess(6,6) = sigma_mu6**2 

                       Covguess(7,7) = sigma_mu7**2 

                       Covguess(8,8) = sigma_mu8**2 

                       Covguess(9,9) = sigma_mu9**2 

                       Covguess(10,10) = sigma_mu10**2 

                       Covguess(11,11) = sigma_Mw**2

                       Covguess(12,12) = sigma_bw**2 

                       Covguess(13,13) = sigma_H0**2 

                       Covguess(14,14) = sigma_Zw**2 
                          
                       Covguess(15,15) = sigma_a_v**2

                       Covguess(16,16) = sigma_a_cal**2

                    Else
                        
                       print *,'WRONG NUMBER OF MODEL PARAMETERS. CHECK FIDUCIAL MODULE'

                       stop

                    End If

                 End If

              Else

                 If (use_H_band) then

                    print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE FOR LMC+MW+NGC4258 AS ANCHORS'
                     
                    stop

                 Else

                    print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE FOR LMC+NGC4258+MW AS ANCHORS'
                     
                    stop

                 End If

              End If

           Else If ( ( use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then

              If (use_metallicity) then 

                 If (use_H_band) then
                     
                    print *,'H BAND NOT IMPLEMENTED USING LMC AND NGC4258 AS ANCHORS AND INCLUDING METALLICITY DEPENDENCE'
                     
                    stop
                     
                 Else

                    If (number_model_parameters .eq. 16) then
                          
                       Covguess(1,1) = sigma_mu1**2 

                       Covguess(2,2) = sigma_mu2**2 

                       Covguess(3,3) = sigma_mu3**2 

                       Covguess(4,4) = sigma_mu4**2 

                       Covguess(5,5) = sigma_mu5**2 

                       Covguess(6,6) = sigma_mu6**2 

                       Covguess(7,7) = sigma_mu7**2 

                       Covguess(8,8) = sigma_mu8**2 

                       Covguess(9,9) = sigma_mu9**2 

                       Covguess(10,10) = sigma_mu10**2 

                       Covguess(11,11) = sigma_Mw**2

                       Covguess(12,12) = sigma_bw**2 

                       Covguess(13,13) = sigma_H0**2 

                       Covguess(14,14) = sigma_Zw**2 

                       Covguess(15,15) = sigma_a_v**2

                       Covguess(16,16) = sigma_a_cal**2

                    Else
                          
                       print *,'WRONG NUMBER OF MODEL PARAMETERS. CHECK FIDUCIAL MODULE'

                       stop

                    End If

                 End If

              Else

                 If (use_H_band) then

                    print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE FOR LMC AND NGC4258 AS ANCHORS'

                    stop

                 Else

                    print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE FOR LMC AND NGC4258 AS ANCHORS'

                    stop

                 End If

              End If

           Else If ( ( use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. use_MW_as_anchor ) then

              If (use_metallicity) then 

                 If (use_H_band) then

                    print *,'H BAND NOT IMPLEMENTED USING MW+NGC4258 AS ANCHOR AND INCLUDING METALLICITY DEPENDENCE'

                    stop

                 Else

                    If (number_model_parameters .eq. 15) then

                       Covguess(1,1) = sigma_mu1**2 

                       Covguess(2,2) = sigma_mu2**2 

                       Covguess(3,3) = sigma_mu3**2 

                       Covguess(4,4) = sigma_mu4**2 

                       Covguess(5,5) = sigma_mu5**2 

                       Covguess(6,6) = sigma_mu6**2 

                       Covguess(7,7) = sigma_mu7**2 

                       Covguess(8,8) = sigma_mu8**2 

                       Covguess(9,9) = sigma_mu9**2 

                       Covguess(10,10) = sigma_Mw**2 

                       Covguess(11,11) = sigma_bw**2 

                       Covguess(12,12) = sigma_H0**2 

                       Covguess(13,13) = sigma_Zw**2 

                       Covguess(14,14) = sigma_a_v**2

                       Covguess(15,15) = sigma_a_cal**2

                    Else

                       print *,'WRONG NUMBER OF MODEL PARAMETERS. CHECK FIDUCIAL MODULE'

                       stop

                    End If

                 End If

              Else

                 If (use_H_band) then

                    print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE FOR MW+NGC4258 AS ANCHORS'

                    stop

                 Else

                    print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE FOR MW+NGC4258 AS ANCHORS'

                    stop

                 End If

              End If

           Else If ( ( .not.use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. use_MW_as_anchor ) then

              If (use_metallicity) then 

                 If (use_H_band) then

                    print *,'H BAND NOT IMPLEMENTED USING LMC AND MW AS ANCHORS AND INCLUDING METALLICITY DEPENDENCE'

                    stop

                 Else

                    If (number_model_parameters .eq. 16) then

                       Covguess(1,1) = sigma_mu1**2 

                       Covguess(2,2) = sigma_mu2**2 

                       Covguess(3,3) = sigma_mu3**2 

                       Covguess(4,4) = sigma_mu4**2 

                       Covguess(5,5) = sigma_mu5**2 

                       Covguess(6,6) = sigma_mu6**2 

                       Covguess(7,7) = sigma_mu7**2 

                       Covguess(8,8) = sigma_mu8**2 

                       Covguess(9,9) = sigma_mu9**2 

                       Covguess(10,10) = sigma_mu10**2 

                       Covguess(11,11) = sigma_Mw**2

                       Covguess(12,12) = sigma_bw**2 

                       Covguess(13,13) = sigma_H0**2 

                       Covguess(14,14) = sigma_Zw**2 

                       Covguess(15,15) = sigma_a_v**2

                       Covguess(16,16) = sigma_a_cal**2

                    Else

                       print *,'WRONG NUMBER OF MODEL PARAMETERS. CHECK FIDUCIAL MODULE'

                       stop

                    End If

                 End If

              Else

                 If (use_H_band) then

                    print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE FOR LMC AND MW AS ANCHORS'

                    stop

                 Else

                    print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE FOR LMC AND MW AS ANCHORS'

                    stop

                 End If

              End If

           Else If ( ( .not.use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. use_MW_as_anchor ) then

              If (use_metallicity) then 

                 If (use_H_band) then

                    print *,'H BAND NOT IMPLEMENTED USING MW AS ANCHOR AND INCLUDING METALLICITY DEPENDENCE'

                    stop

                 Else

                    If (number_model_parameters .eq. 17) then

                       Covguess(1,1) = sigma_mu1**2 

                       Covguess(2,2) = sigma_mu2**2 

                       Covguess(3,3) = sigma_mu3**2 

                       Covguess(4,4) = sigma_mu4**2 

                       Covguess(5,5) = sigma_mu5**2 

                       Covguess(6,6) = sigma_mu6**2 

                       Covguess(7,7) = sigma_mu7**2 

                       Covguess(8,8) = sigma_mu8**2 

                       Covguess(9,9) = sigma_mu9**2 

                       Covguess(10,10) = sigma_Mw**2 

                       Covguess(11,11) = sigma_bw**2 

                       Covguess(12,12) = sigma_H0**2 

                       Covguess(13,13) = sigma_Zw**2 

                       Covguess(14,14) = sigma_a_v**2

                       Covguess(15,15) = sigma_a_cal**2

                       Covguess(16,16) = sigma_sigma_int**2

                       Covguess(17,17) = sigma_sigma_int**2

                    Else

                       print *,'WRONG NUMBER OF MODEL PARAMETERS. CHECK FIDUCIAL MODULE'

                       stop

                    End If

                 End If

              Else

                 If (use_H_band) then

                    print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE FOR MW AS ANCHOR'

                    stop

                 Else

                    print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE FOR MW AS ANCHOR'

                    stop

                 End If

              End If

           Else If ( ( .not.use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then

              If (use_metallicity) then 

                 If (use_H_band) then

                    print *,'H BAND NOT IMPLEMENTED USING LMC AS ANCHOR AND INCLUDING METALLICITY DEPENDENCE'

                    stop

                 Else

                    If (number_model_parameters .eq. 16) then

                       Covguess(1,1) = sigma_mu1**2 

                       Covguess(2,2) = sigma_mu2**2 

                       Covguess(3,3) = sigma_mu3**2 

                       Covguess(4,4) = sigma_mu4**2 

                       Covguess(5,5) = sigma_mu5**2 

                       Covguess(6,6) = sigma_mu6**2 

                       Covguess(7,7) = sigma_mu7**2 

                       Covguess(8,8) = sigma_mu8**2 

                       Covguess(9,9) = sigma_mu9**2 

                       Covguess(10,10) = sigma_mu10**2 

                       Covguess(11,11) = sigma_zpwLMC**2

                       Covguess(12,12) = sigma_bw**2 

                       Covguess(13,13) = sigma_H0**2 

                       Covguess(14,14) = sigma_Zw**2 

                       Covguess(15,15) = sigma_a_v**2

                       Covguess(16,16) = sigma_a_cal**2

                    Else

                       print *,'WRONG NUMBER OF MODEL PARAMETERS. CHECK FIDUCIAL MODULE'

                       stop

                    End If

                 End If

              Else

                 If (use_H_band) then

                    print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE FOR LMC AS ANCHOR'

                    stop

                 Else

                    print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE FOR LMC AS ANCHOR'

                    stop

                 End If

              End If

           Else If ( ( use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then

              If (use_metallicity) then 

                 If (use_H_band) then

                    If (number_model_parameters .eq. 14) then ! ASSIGN VALUES TO COVARIANCE MATRIX. SINCE VALUES OF PARAMETERS DEPENDING 

                       Covguess(1,1) = sigma_mu1**2           ! ON H BAND ARE NOT EXPECTED TO BE TOO MUCH DIFFERENT FROM THOSE FOR W BAND 

                       Covguess(2,2) = sigma_mu2**2           ! I USE SAME VALUES AS FOR W BAND BELOW

                       Covguess(3,3) = sigma_mu3**2 

                       Covguess(4,4) = sigma_mu4**2 

                       Covguess(5,5) = sigma_mu5**2 

                       Covguess(6,6) = sigma_mu6**2 

                       Covguess(7,7) = sigma_mu7**2 

                       Covguess(8,8) = sigma_mu8**2 

                       Covguess(9,9) = sigma_mu9**2 

                       Covguess(10,10) = sigma_zpw**2

                       Covguess(11,11) = sigma_bw**2 

                       Covguess(12,12) = sigma_H0**2 

                       Covguess(13,13) = sigma_Zw**2 

                       Covguess(14,14) = sigma_a_v**2

                    Else

                       print *,'WRONG NUMBER OF MODEL PARAMETERS. CHECK FIDUCIAL MODULE AND COMPARE WITH EQUATION (18) IN R09'

                       stop

                    End If

                 Else

                    If (number_model_parameters .eq. 14) then

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

                       Covguess(12,12) = sigma_H0**2 

                       Covguess(13,13) = sigma_Zw**2 

                       Covguess(14,14) = sigma_a_v**2

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

           End If ! OF ANCHORS 

        End If ! OF INCLUDE ONLY CEPHEIDS
    
     Else

        If (fit_MW_cepheids_alone) then 

           If (use_metallicity) then

              Covguess(1,1) = sigma_Mw**2 

              Covguess(2,2) = sigma_bw**2

              Covguess(3,3) = sigma_Zw**2 

              Covguess(4,4) = sigma_sigma_int**2

           Else

              Covguess(1,1) = sigma_Mw**2 

              Covguess(2,2) = sigma_bw**2

              Covguess(3,3) = sigma_sigma_int**2

           End If

        Else

           Covguess(1,1) = sigma_A**2 

           Covguess(2,2) = sigma_bw**2

           !Covguess(3,3) = sigma_sigma_int**2

        End If
 
     End If ! OF R11 ANALYSIS
         
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

  End If

  inquire(file='./output/covariance_matrix.txt',exist=cov_file)  

  inquire(file='./output/mcmc_output.txt',exist=out_file)  
  
  If (cov_file) then

     call system('rm  ./output/covariance_matrix.txt')

  Else

     continue

  End If

  If (out_file) then

     call system('rm  ./output/mcmc_output.txt')

  Else

     continue

  End If

end subroutine set_covariance_matrix

subroutine read_data()

  use fiducial
  use arrays

  Implicit none

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

        call read_data_Efstathiou(path_to_datafileABC)

        call read_table2_LEEUWEN(path_to_table2_LEEUWEN)

     End If

     If (hyperparameters_as_mcmc) then

        open(UNIT_EXE_FILE,file=EXECUTION_INFORMATION)    ! OPEN FILE FOR EXECUTION INFORMATION

        write(UNIT_EXE_FILE,*) 'WORKING WITH HYPER-PARAMETERS'     

        write(UNIT_EXE_FILE,*) 'HYPER-PARAMETERS USED AS MCMC PARAMETERS'

        call read_data_Efstathiou(path_to_datafileABC)

        If (number_hyperparameters .ne. (size(NameA)+size(NameB)+size(NameC))) then

           write(UNIT_EXE_FILE,*) 'NUMBER OF HYPER-PARAMETERS MUST MATCH TOTAL NUMBER OF DATA POINTS ',&
                size(NameA)+size(NameB)+size(NameC)

           write(UNIT_EXE_FILE,*) 'CHECK FIDUCIAL MODULE'

           stop

        End If

     Else

        If (doing_R11_analysis) then

           open(UNIT_EXE_FILE,file=EXECUTION_INFORMATION)    ! OPEN FILE FOR EXECUTION INFORMATION

           write(UNIT_EXE_FILE,*) 'WORKING WITH HYPER-PARAMETERS AND DOING R11 ANALYSIS'

        Else

           open(UNIT_EXE_FILE,file=EXECUTION_INFORMATION)    ! OPEN FILE FOR EXECUTION INFORMATION

           write(UNIT_EXE_FILE,*) 'WORKING WITH HYPER-PARAMETERS BUT NOT DOING R11 ANALYSIS. DOING SECTION 2 OF EFSTATHIOU S PAPER'

        End If

        If (number_hyperparameters .ne. 0) then

           write(UNIT_EXE_FILE,*) 'NUMBER OF HYPER-PARAMETERS MUST BE ZERO WHEN SETTING FALSE "hyperparameters_as_mcmc" '

           write(UNIT_EXE_FILE,*) 'CHECK FIDUCIAL MODULE'

           stop

        End If

     End If

  Else

     call read_data_Efstathiou(path_to_datafileAB)

     If (include_table2_R11) then

        call read_table2_R11(path_to_table2_R11)

     End If

     open(UNIT_EXE_FILE,file=EXECUTION_INFORMATION)    ! OPEN FILE FOR EXECUTION INFORMATION 

  End If

end subroutine read_data

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
