Module arrays

    Integer :: status1,status2,status3,status4,status5,status6
    Integer*4,dimension(13) :: buff

    Real*8, allocatable, dimension(:) :: PeriodA,HA,Sigma_mA,VA,IIA,old_point
    Real*8, allocatable, dimension(:) :: PeriodB,HB,Sigma_mB,VB,IIB,current_point
    Real*8, allocatable, dimension(:) :: PeriodC,HC,Sigma_mC,VC,IIC
    Real*8, allocatable, dimension(:) :: Period,H,Sigma_m,V,II
    Real*8, allocatable, dimension(:,:,:,:,:) :: cov,inv_cov
    Real*4, allocatable, dimension(:) :: acceptance_probability

    character(len=10),allocatable,dimension(:) :: Name,NameA,NameB,NameC

End module arrays
