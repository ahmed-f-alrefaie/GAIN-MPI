!
module symmetry
  use fort_func
  implicit none
  public SymmetryInitialize,sym,max_irreps,Symm_initialize,Get_Nirreps

  integer, parameter :: sik         = selected_int_kind(4)       ! Small integers
  integer, parameter :: ik          = selected_int_kind(8)       ! "Normal" integers. This must map on
                                                                 ! C "int" type, or the DX interface won't
                                                                 ! work correctly.
  integer, parameter :: hik         = selected_int_kind(16)      ! "Pointer" integers - sufficient to store
                                                                 ! memory address
  integer, parameter :: drk         = selected_real_kind(12,25)  ! "Double" reals and complex (complexi? :-)
  integer, parameter :: rk          = selected_real_kind(12,25)  ! "Normal" reals and complex (complexi? :-)
  integer, parameter :: ark         = selected_real_kind(25,32)  ! "Accurate" reals and complex (complexi? :-)
  integer, parameter :: inp         = 5                          ! Output I/O channel
  integer, parameter :: out         = 6                          ! Output I/O channel
  integer, parameter :: nfilelegendre = 101                      ! Dump-outout channel for eigenfunction 

  ! universal constants
  real(drk), parameter :: planck     =  6.62606896e-27             ! Planck constant
  real(drk), parameter :: avogno     =  6.0221415E+23             ! Avogadro constant
  real(drk), parameter :: vellgt     =  2.99792458E+10            ! Speed of light constant
  real(drk), parameter :: boltz      =  1.380658E-16              ! Boltzmann constant
  real(drk), parameter :: bohr       =  0.529177249               ! a.u.
  
  real(rk)           :: safe_max                                 ! Largest number we want to work with
  real(rk)           :: safe_min                                 ! Smalles number we want to work with
  real(rk)           :: max_exp                                  ! Largest number OK for exponentiating
  real(ark)          :: small_a = epsilon(1.0_ark)                                  ! a positive model number that is almost 
  real(rk)           :: small_  = epsilon(1.0_rk)                                  ! a positive model number that is almost 
                                                                 ! negligible compared to unity in the current model 
  real(ark)          :: pi,twopi                                 ! Pi, 2*Pi
  real(ark)          :: epsil(3,3,3)                             ! epsil - antisymmetric tensor
  real(rk),parameter :: sqrt2 = 1.414213562373095048801689_rk    ! \sqrt{2}
  real(rk),parameter :: sqrt3 = 1.732050807568877293527446_rk    ! \sqrt{3}
  real(rk),parameter :: rad   = 57.295779513082320875_rk    ! radian = 180/pi
  integer, parameter :: cl          = 80                         ! Max character string length
  integer, parameter :: wl          = 500                       ! Very large max character string length 
  integer, parameter :: fititermax  = 200                        ! Max number of iteratons in different fittings 



   type MOrepres_arkT
      real(ark),pointer     :: repres(:,:,:)  
   end type MOrepres_arkT

  type  ScIIT
     integer(ik)          :: Noper     ! Number of operations in the CII operator
     integer(ik)          :: Nzeta     ! Number of symmetric elements taking into account the degeneracies 
     real(ark),pointer    :: ioper(:)  ! the operation number in the MS group
     integer(ik),pointer  :: coeff(:)  ! coefficients of the CII operator
     integer(ik),pointer  :: izeta(:)  ! symmetry indentification as a eigenvalues of the CII operator
  end type ScIIT

  type  SrepresT
     real(ark),pointer  :: repres(:,:)      ! matrix representation of the group 
  end type SrepresT

  type  SymmetryT
     character(len=cl)    :: group = 'C' ! The symmetry group 
     integer(ik)          :: Nrepresen   ! Number of irreduc. represent.
     integer(ik)          :: Noper       ! Number of operations
     integer(ik)          :: Nclasses   ! Number of classes
     integer,pointer :: Nelements(:)   ! Number of elements in a class
     real(rk),pointer, dimension(:,:)     :: characters! Character table
     type(SrepresT),pointer, dimension(:,:) :: irr     ! irreducible representaion 
     integer(ik),pointer, dimension(:)  :: degen       ! degeneracy
     character(len=4),pointer,dimension(:)  :: label  ! The symmetry label 
     integer(ik)          :: Maxdegen  = 1  ! Maximal degeneracy order
     integer(ik),pointer,dimension(:)  :: igenerator  ! address of the class generator in the sym%Ngroup list
     type(ScIIT)          :: CII            ! the elements of the CII operator 
     real(ark),pointer, dimension(:,:)   :: euler     ! rotational angles equivalent to the group operations
     integer(ik)          :: class_size_max = 8 ! current maximal class size 
     integer(ik)          :: N  = 1         ! The group order, currently desgined for Dnh where N is odd 
     integer(ik),pointer,dimension(:)  :: lquant      ! Store the value of the (vibrational) angular momentum 
     !
  end type SymmetryT




  type(SymmetryT) , save  :: sym
  integer(ik),parameter   :: max_irreps=100
  integer(ik),parameter   :: verbose_ = 3

contains 

  subroutine Get_Nirreps(J,Ntot,eigenvects) bind (c, name='c_GetNirreps')
  use,intrinsic :: iso_c_binding
  implicit none
	integer(4),intent(in)   :: J
	integer(4),intent(out)  :: Ntot(sym%Nrepresen)
	real(8),intent(out)	:: eigenvects(2*J+1,2*J+1)
	integer(ik)		::	Ncount,count_index(2*j+1,2*j+1)

	call MLrotsymmetry_generate(J,6,count_index,eigenvects,Ncount,Ntot)

  end subroutine Get_Nirreps

  subroutine Symm_initialize(input_string) bind (c, name='c_SymInit')
  use,intrinsic :: iso_c_binding
  implicit none
	!character (kind=c_char, len=1), dimension (10), intent (in) :: input_string
	character(len=cl) ::	sym_group
	 character (kind=c_char, len=1), dimension (80), intent (in) :: input_string
	integer(ik)	::	i	
	

	sym_group = " "
  	loop_string: do i=1, 80
 	     if ( input_string (i) == c_null_char ) then
        	 exit loop_string
      		else
         sym_group (i:i) = input_string (i)
      end if
   end do loop_string

	sym_group = 'TD(M)'
      safe_max = huge(1.0_rk)**(0.25_rk)
      safe_min = tiny(1.0_rk)**(0.25_rk)
      max_exp  = log(safe_max)
      pi       = 4.0_ark * atan2(1.0_ark,1.0_ark)
      twopi    = 2.0_ark * pi
      small_ = epsilon(1.0_rk )
      small_a= epsilon(1.0_ark)

      epsil = 0 
      epsil(1,2,3) = 1.0_ark
      epsil(1,3,2) =-1.0_ark
      epsil(2,1,3) =-1.0_ark
      epsil(2,3,1) = 1.0_ark
      epsil(3,1,2) = 1.0_ark
      epsil(3,2,1) =-1.0_ark
	call SymmetryInitialize(sym_group)
	!MLrotsymmetry_generate(J,5,count_index,eigenvects,Ncount,Ntot)

  end subroutine Symm_initialize

  subroutine GetNlevels(jI,jF,nLevelsI,nlevelsF) bind (c, name='c_GetNelevels')
    !
    use,intrinsic :: iso_c_binding
    implicit none 
    integer(ik), intent(in) :: jI,jF
    integer(ik), intent(out) :: nLevelsI,nlevelsF
    ! real (c_double), allocatable, dimension (:,:,:,:,:), target :: rot
    
    !
    integer(ik)         :: alloc
    integer(ik)         :: dimenI,nrootsI,idegI,ndeg,dimenf,nrootsF,idegF,irootI,irootF,dk
    integer(ik)         :: icountI,icountF,kI,kF,tauI,tauF,NcountI,NtotalI(sym%Nrepresen),NcountF,NtotalF(sym%Nrepresen)
    integer(ik)         :: isym,i,k0,tau0,sigmaI,sigmaF,n,indI,indF,dJ
    integer(ik),allocatable   :: count_indexI(:,:),rot_kI(:),rot_tauI(:),count_indexF(:,:),rot_kF(:),rot_tauF(:),NdegI(:),NdegF(:)
    real(rk),allocatable      :: eigenvectsI(:,:),eigenvectsF(:,:)
    real(rk),allocatable   :: mat(:,:,:),mat_(:,:)

    real(rk)           :: f3j,f_IF,f_I
    if(abs(jI-jF) > 1)  return

    ! In order to find the correlation of the vibrational quantum numbers at J=0 and J/=0 
    ! the J=0 quanta destibution has to be defined before. We assume here that this is the case, 
    ! bset_contr(1) is always reserved for J=0
    !
    !
    Ndeg = min(sym%Maxdegen,2)
    !
      !
    dimenI = 2*jI+1 
    nrootsI = dimenI
      !
    allocate(count_indexI(nrootsI,nrootsI),eigenvectsI(nrootsI,nrootsI),rot_kI(dimenI),rot_tauI(dimenI),NdegI(dimenI),stat=alloc)

      !
    rot_kI(1)   = 0
    rot_tauI(1) = mod(jI,2)
      !
    irootI = 1
    do k0 = 1,jI
     	do tau0 = 0,1
          !
          irootI = irootI + 1
          !
          rot_kI(irootI)  = k0
          rot_tauI(irootI)= tau0
          !
       enddo 
    enddo
      ! 
    call MLrotsymmetry_generate(jI,4,count_indexI,eigenvectsI,NcountI,NtotalI)
      icountI = 0
    do isym = 1,sym%Nrepresen
        !
      do i = 1,NtotalI(isym)
          !
        icountI = icountI + 1
          !
        NdegI(icountI) = sym%degen(isym)
          !
        enddo
    enddo
    nlevelsI  = icountI
         !
    
         !
    dimenF = 2*jF+1 
    nrootsF = dimenF
         !
    allocate(count_indexF(nrootsF,nrootsF),eigenvectsF(nrootsF,nrootsF),rot_kF(dimenF),rot_tauF(dimenF),NdegF(dimenF),stat=alloc)
    call MLrotsymmetry_generate(jF,4,count_indexF,eigenvectsF,NcountF,NtotalF)	

    icountF = 0
    do isym = 1,sym%Nrepresen
           !
      do i = 1,NtotalF(isym)
             !
          icountF = icountF + 1
          NdegF(icountF) = sym%degen(isym)
             !
             !do ideg = 1,Ndeg
             !  !
             !  iroot = count_index(icount,ideg)
             !  !
             !enddo
             !
      enddo
    enddo
    nlevelsF=icountF
    deallocate(count_indexF,eigenvectsF,rot_kF,rot_tauF,NdegF)
   
         !
 	!rot_ptr = c_loc (rot)
      !
   deallocate(count_indexI,eigenvectsI,rot_kI,rot_tauI,NdegI)

  end subroutine GetNlevels


  subroutine contraced_rotational_dipole(jI,jF,rot,nLevelsI,nlevelsF) bind (c, name='c_wignerGen')
    !
    use,intrinsic :: iso_c_binding
    implicit none 
    integer(ik), intent(in) :: jI,jF
    integer(ik), intent(inout) :: nLevelsI,nlevelsF
    real(rk), intent (out) :: rot(3,nLevelsI,nlevelsF,sym%Maxdegen,sym%Maxdegen)
    ! real (c_double), allocatable, dimension (:,:,:,:,:), target :: rot
    
    !
    integer(ik)         :: alloc,ijI,ijF,i1,i2
    integer(ik)         :: dimenI,nrootsI,idegI,ndeg,dimenf,nrootsF,idegF,irootI,irootF,dk
    integer(ik)         :: icountI,icountF,kI,kF,tauI,tauF,NcountI,NtotalI(sym%Nrepresen),NcountF,NtotalF(sym%Nrepresen)
    integer(ik)         :: isym,i,k0,tau0,sigmaI,sigmaF,n,indI,indF,dJ
    integer(ik),allocatable   :: count_indexI(:,:),rot_kI(:),rot_tauI(:),count_indexF(:,:),rot_kF(:),rot_tauF(:),NdegI(:),NdegF(:)
    real(rk),allocatable      :: eigenvectsI(:,:),eigenvectsF(:,:)
    real(rk),allocatable   :: mat(:,:,:),mat_(:,:),threej(:,:,:,:)

    real(rk)           :: f3j,f_IF,f_I
    if(abs(jI-jF) > 1)  return

    ! In order to find the correlation of the vibrational quantum numbers at J=0 and J/=0 
    ! the J=0 quanta destibution has to be defined before. We assume here that this is the case, 
    ! bset_contr(1) is always reserved for J=0
    !
    !
    Ndeg = min(sym%Maxdegen,2)
    !
      !
    dimenI = 2*jI+1 
    nrootsI = dimenI
      !
    allocate(count_indexI(nrootsI,nrootsI),eigenvectsI(nrootsI,nrootsI),rot_kI(dimenI),rot_tauI(dimenI),NdegI(dimenI),stat=alloc)

      !
    rot_kI(1)   = 0
    rot_tauI(1) = mod(jI,2)
      !
    irootI = 1
    do k0 = 1,jI
     	do tau0 = 0,1
          !
          irootI = irootI + 1
          !
          rot_kI(irootI)  = k0
          rot_tauI(irootI)= tau0
          !
       enddo 
    enddo
      ! 
    call MLrotsymmetry_generate(jI,4,count_indexI,eigenvectsI,NcountI,NtotalI)
      icountI = 0
    do isym = 1,sym%Nrepresen
        !
      do i = 1,NtotalI(isym)
          !
        icountI = icountI + 1
          !
        NdegI(icountI) = sym%degen(isym)
          !
        enddo
    enddo
    nlevelsI  = icountI

     allocate(threej(0:(jI+1),0:(jI+1),-1:1,-1:1), stat = alloc)

    do ijI = 0,jI+2
      do ijF = max(ijI-1,0),min(ijI+1,(jI+1))
        do kI = 0,ijI
          do kF = max(kI-1,0),min(kI+1,ijF)
            !
            !threej(jI, kI, jF - jI, kF - kI) = cg(jI, kI, jF - jI, kF - kI) *                                                 &
            !               (-1.0_rk) ** (jI - 1 + kF) / sqrt( real(2*jF + 1,rk) )
            !
            !dp = three_j(jI, 1, jF, kI, kF - kI, -kF)
            !
            !if (abs(dp-threej(jI, kI, jF - jI, kF - kI))>sqrt(small_)) then 
            !  write(out,"('a problem with threej for jI,jF,kI,kF = ',4i4,3es16.8)") jI,jF,kI,kF,dp,threej(jI, kI, jF - jI, kF - kI),dp-threej(jI, kI, jF - jI, kF - kI)
            !  stop 'a problem with threej'
            !endif
            !
            threej(jI, kI, jF - jI, kF - kI) = three_j(jI, 1, jF, kI, kF - kI, -kF)
	    !write(out,"('threeJ[',i4, ',' ,i4, ',' ,i4, ',' ,i4,']=',f12.6)") jI, kI, (jF - jI), (kF - kI),threej(jI, kI, jF - jI, kF - kI) 
	        !write(out, "( 'threej[' , i4 , ',' , i4 , ',' ,i4, ',' ,i4,'] = ',f12.6)") jI, kI, jF - jI, kF - kI,threej(jI, kI, jF - jI, kF - kI)
            !
          enddo
        enddo
      enddo
    enddo

      !
         !
         !
    dJ = jF-jI
         !
    
         !
    dimenF = 2*jF+1 
    nrootsF = dimenF
         !
    allocate(count_indexF(nrootsF,nrootsF),eigenvectsF(nrootsF,nrootsF),rot_kF(dimenF),rot_tauF(dimenF),NdegF(dimenF),stat=alloc)
     !     do i1=1,dimenI
!	do i2=1,dimenI
		
!		write(out,"('eigenvectsI(',i3,i3,')=',es14.3)") i1,i2,eigenvectsI(i1,i2)	
!	enddo
 !   enddo
         ! 
    call MLrotsymmetry_generate(jF,4,count_indexF,eigenvectsF,NcountF,NtotalF)
    !      do i1=1,dimenF
	!do i2=1,dimenF
		
	!	write(out,"('eigenvectsF(',i3,i3,')=',es14.3)") i1,i2,eigenvectsF(i1,i2)	
	!enddo
    !enddo
         !
    rot_kF(1)   = 0
    rot_tauF(1) = mod(jF,2)
         !
    irootF = 1
    do k0 = 1,jF
       do tau0 = 0,1
             !
         irootF = irootF + 1
             !
          rot_kF(irootF)  = k0
          rot_tauF(irootF)= tau0
             !
        enddo 
    enddo
         !
    icountF = 0
    do isym = 1,sym%Nrepresen
           !
      do i = 1,NtotalF(isym)
             !
          icountF = icountF + 1
          NdegF(icountF) = sym%degen(isym)
             !
             !do ideg = 1,Ndeg
             !  !
             !  iroot = count_index(icount,ideg)
             !  !
             !enddo
             !
      enddo
    enddo
    nlevelsF  = icountF
         !
    allocate(mat(3,dimenI,dimenF),mat_(dimenI,dimenF),stat=alloc)

         !
    !allocate(rot(3,nlevelsI,nlevelsF,sym%maxdegen,sym%maxdegen),stat=alloc)

         !
     rot = 0
         !
    mat = 0
         !
    do irootI  = 1,dimenI
            !
            kI   = rot_kI(irootI)
            tauI = rot_tauI(irootI)
            sigmaI = mod(kI, 3)*tauI
            !
            f_I = (-1.0_rk)**(sigmaI+kI)
            !
            do irootF  = 1,dimenF
               !
               kF   = rot_kF(irootF)
               tauF = rot_tauF(irootF)
               sigmaF = mod(kf, 3)*tauF
               !
               f_IF = f_I*(-1.0_rk)**(sigmaF)
               !
               dk = kI - kF
               !
               if (abs(dk)>1) cycle
               !
               f3j  =  threej(jI, kI, jF - jI, kF - kI)
               !We regenerate the threej each time because I'm lazy
	      ! f3j  =  three_j(jI, 1, jF, kI, kF - kI, -kF);
               !
               if (abs(f3j)<small_) cycle
               !
               if (kF == kI) then
                   !
                   mat(3,irootI,irootF) = real(tauF-tauI,rk)*f3j*f_IF
                   !
                   !mat(3,irootI,irootF) = f3j*f_IF
                   !
               elseif(tauF/=tauI) then 
                   !
                   mat(1,irootI,irootF) = real((kF-kI)*(tauF-tauI),rk)*f3j*f_IF
                   !
                   !mat(1,irootI,irootF) = f3j*f_IF
                   !
                   if (kI*kF /= 0) mat(1,irootI,irootF) = mat(1,irootI,irootF)/sqrt(2.0_rk)
                   !
               elseif(tauF==tauI) then 
                   !
                   mat(2,irootI,irootF) = -f3j*f_IF
                   !
                   !mat(2,irootI,irootF) = f3j
                   !
                   if (kI*kF /= 0) mat(2,irootI,irootF) = mat(2,irootI,irootF)/sqrt(2.0_rk)
                   !
               endif
               !
            enddo
            !
         enddo
         !
         do n = 1,3
           !
           mat_(:,:) = matmul( matmul( transpose(eigenvectsI),mat(n,:,:) ),eigenvectsF )
           !
           do icountI = 1,nlevelsI
             do idegI =1,NdegI(icountI)
               irootI = count_indexI(icountI,idegI)
               !
               do icountF = 1,nlevelsF
                 do idegF =1,NdegF(icountF)
                   irootF = count_indexF(icountF,idegF)
                   rot(n,icountI,icountF,idegI,idegF) = mat_(irootI,irootF)
		  ! write(out,"('wigner(',i4,',',i4,').rot(',i4,',',i4,',',i4,',',i4,',',i4,') = ',es14.3)") jI,dJ,n,icountI,icountF,idegI,idegF,rot(n,icountI,icountF,idegI,idegF)
                 enddo
               enddo
               !
             enddo
           enddo
           !
         enddo
         !

	

  deallocate(mat,mat_,count_indexF,eigenvectsF,rot_kF,rot_tauF,NdegF)
   deallocate(threej)
         !
 	!rot_ptr = c_loc (rot)
      !
  deallocate(count_indexI,eigenvectsI,rot_kI,rot_tauI,NdegI)
      
      !

    !
  end subroutine contraced_rotational_dipole

  subroutine MLrotsymmetry_generate(J,iverbose,count_index,eigenvects,Ncount,Ntot)
    !
    integer(ik),intent(in)   :: J,iverbose
    integer(ik),intent(out)  :: count_index(2*j+1,2*j+1),Ncount,Ntot(:)
    real(rk),   intent(out)  :: eigenvects(2*j+1,2*j+1)
    real(ark)                :: theta,phi,chi,coeff,vect_t(2*j+1)
    !
    integer(ik)              :: i1,k1,n1,i2,k2,n2,tau1,tau2,i,ioper,sigma,alloc,iroot,icount,ideg,nroots
    integer(ik)              :: izeta,irep,jroot,irep_t
    complex(ark)             :: temp,d
    complex(ark),allocatable :: R(:,:),T(:,:),Q(:,:)
    real(ark),allocatable    :: tmat(:,:,:)

    type(MOrepres_arkT),allocatable   :: irr(:)

    !
    integer(ik)              :: count_degen(2*j+1),Ntotal(1:sym%Nrepresen)
    integer(ik)              :: isym
    !
    if (iverbose>=5) write(out,"(/'MLrotsymmetry_generate/start')") 
    !
    ! trivial case of no symmetry
    !
  !  if (trim(molec%symmetry)=='C'.or.trim(molec%symmetry)=='C(M)') then
       !
       !gamma = 1
   !    !ideg = 1
    !   return
       !
    !endif
    !
    ! dimension of the problem
    !
    nroots  = 2*j+1
    !
    ! Obtain transformation-matrix 
    !
    allocate(tmat(sym%Noper,nroots,nroots),R(nroots,nroots),T(nroots,nroots),Q(nroots,nroots),stat=alloc)
    if (alloc/=0) stop 'MLrotsymmetry_generate - out of memory'

    !
    ! R is a transformation to the Wang functions
    !
    R = 0 
    !
    if (mod(j,2)==0) then 
      !
      r(1,1) = cmplx(1.0_ark,0.0_ark)
      !
    else
      !
      r(1,1) = cmplx(0.0_ark,-1.0_ark)
      !
    endif
    !
    i1 = 2
    do k1 = 1,j
       !
       sigma = mod(k1,3)
       !
       R(i1  ,i1  ) = cmplx(1.0_ark/sqrt(2.0_ark)           ,0.0_ark                               )
       R(i1  ,i1+1) = cmplx(0.0_ark                         ,-(-1.0_ark)**sigma/sqrt(2.0_ark)      )
       R(i1+1,i1  ) = cmplx((-1.0_ark)**(j+k1)/sqrt(2.0_ark),0.0_ark                               )
       R(i1+1,i1+1) = cmplx(0.0_ark                         ,(-1.0_ark)**(j+k1+sigma)/sqrt(2.0_ark))
       !
       i1 = i1 + 2
       !
    enddo
    !
    ! transformation to the wang rot. basis 
    !
    !
    do ioper = 1,sym%Noper
      !
      i1 = 0
      do k1 = 0,j
         do tau1 = 0,min(1,k1)
           !if (k1==0.and.mod(j+tau1,2)/=0) cycle
           i1 = i1 + 1
           !
           n1=k1*(-1)**tau1
           !
           i2 = 0
           do k2 = 0,j
              do tau2 = 0,min(1,k2)
                !if (k2==0.and.mod(j+tau2,2)/=0) cycle
                i2 = i2 + 1
                !
                n2=k2*(-1)**tau2
                !
                phi   = sym%euler(ioper,1)
                theta = sym%euler(ioper,2)
                chi   = sym%euler(ioper,3)
                !
                T(i2,i1) =  ddlmn_conj(j,n2,n1,phi,theta,chi)
                !
              enddo
           enddo
           !
         enddo
      enddo
      !
      Q = matmul( transpose( conjg(R) ),matmul( T,R ) )
      !
      tmat(ioper,:,:) = real(Q,ark)
      !
    enddo
    !
    allocate (irr(sym%Nrepresen),stat=alloc)
    !
    call MLsymmetrization_rot(nroots,tmat,iverbose,Ntotal)
    !
    do isym = 1,sym%Nrepresen
      !
      allocate (irr(isym)%repres(max(Ntotal(isym),1),sym%degen(isym),nroots),stat=alloc)

      !
      irr(isym)%repres = 0
      !
    enddo
    !
    call MLsymmetrization_rot(nroots,tmat,iverbose,Ntotal,irr)
    !
    iroot = 0
    icount = 0
    !
    do isym = 1,sym%Nrepresen
      !
      do i1 = 1,Ntotal(isym)
        !
        icount = icount + 1
        !
        do ideg = 1,sym%degen(isym)
          !
          iroot = iroot + 1
          ! 
          count_index(icount,ideg) = iroot
          count_degen(icount) = ideg
          !
          eigenvects(:,iroot) = irr(isym)%repres(i1,ideg,:)
          !
        enddo
        !
      enddo
      !
    enddo

    !do i1=1,2*j+1
	!do i2=1,2*j+1
	!	write(out,"('count_index(',i3,i3,') =',i3)") i1,i2,count_index(i1,i2)
	!	write(out,"('eigenvects(',i3,i3,')=',es14.3)") i1,i2,eigenvects(i1,i2)	
	!enddo
    !enddo
    !
    do ioper = 1,sym%Noper
      !
      Q(:,:) = matmul( ( transpose(eigenvects) ),matmul( tmat(ioper,:,:),(eigenvects) ) )
      !
      tmat(ioper,:,:) = real(Q)
      !
    enddo
    !
    Ncount = icount
    Ntot = Ntotal
    !
    deallocate(tmat,R,T,Q,irr)
    !
    if (iverbose>=5) write(out,"('MLrotsymmetry_generate/end')") 
    !
    !
  end subroutine MLrotsymmetry_generate

   !
   !                               Symmetrization.  
   !
   !
   ! Construction of irreducible representation - to project the rotational 
   ! basis set into the symmetrized representaion
   !

  subroutine MLsymmetrization_rot(dimen,transform,iverbose,Ntotal,irr)

   integer(ik),intent(in) :: dimen
   real(ark),intent(in)   :: transform(:,:,:)
   integer(ik),intent(in) :: iverbose
   integer(ik),intent(out):: Ntotal(1:sym%Nrepresen)
   type(MOrepres_arkT),optional,intent(out)  :: irr(1:sym%Nrepresen)
   !
   integer(ik)             :: alloc,isym,ioper,ideg,iclasses,Nirr(sym%Nrepresen)
   integer(ik)             :: i0,ielem_t
   real(ark)               :: chi(sym%Nclasses),Nirr_rk(sym%Nrepresen),f_t,g_t
   real(ark),allocatable   :: projector(:,:,:,:)
   real(ark),allocatable   :: vect(:,:),t_vect(:,:)
   !
   integer(ik)             :: Nelem,Ndeg,ielem,jelem,felem,try_elem,jdeg,kdeg,jsym
   !
   allocate (projector(sym%Nrepresen,sym%Noper,sym%maxdegen,sym%maxdegen),stat=alloc)

   !
   do isym =1,sym%Nrepresen
     !
     do ioper =1,sym%Noper
       !
       Ndeg = sym%degen(isym)
       !
       do ideg = 1,Ndeg
         !
         projector(isym,ioper,ideg,1:Ndeg) = sym%irr(isym,ioper)%repres(ideg,1:Ndeg)*real(sym%degen(isym),ark)/real(sym%Noper,ark)
         !
       enddo
     enddo 
   enddo 
   !
   ! We symmetrized the contracted basis set 
   ! by calculating irreducible represntations of the product of Nclasses the symmetrized components 
   ! Here the field 'represent' defines how the product of the contracted functions 
   ! transform with the group operations, i.e. the representation of basis set in the 
   ! contracted basis. 
   ! The field 'transform' will define the transformation from the contracted basis set 
   ! to the symmetrized basis set. 
   !
   ! first time only to estimate the total numbers (Ntotal) of irred. representations 
   !
   if (iverbose>=5.and..not.present(irr)) write(out,"('Counting the total number of irr. representations...')") 
   !
   Ntotal = 0 
   !
   if (iverbose>=5) then 
      write(out,"(//'Rotational irreducible representations')")
      write(out,"('  Find numbers of irr. representations...')")
   endif 
   !
   ! calculate the characters 
   !
   chi = 0
   ioper = 1
   !
   do iclasses =1,sym%Nclasses
     !
     do ideg = 1,dimen
       !
       chi(iclasses) = chi(iclasses) + transform(ioper,ideg,ideg)
       !
     enddo 
     !
     if (iverbose>=6) then 
       write(out,"('ioper = ',i5,'; character: ',f12.8)") ioper,chi(iclasses)
     endif 
     !
     ioper = ioper+sym%Nelements(iclasses)
     !
   enddo 
   !
   ! estimate the number of the irreducible representasions 
   !
   do isym =1,sym%Nrepresen
      !
      Nirr_rk(isym) = sum(real(sym%Nelements(:)*sym%characters(isym,:),ark)*chi(:))/real(sym%Noper,ark)
      Nirr(isym) = nint(Nirr_rk(isym),ik)
      !
   enddo 
   !
   if (all(Nirr(:)==0)) then 
     write(out,"('No symmetry defined, Nirr = 0')") 
     stop 'No symmetry defined, Nirr = 0'
   endif 
   !
   if (iverbose>=5) then 
     write(out,"('Number of irr. rep-s = ',30f12.4)") Nirr_rk(1:min(30,size(Nirr_rk)))
   endif 
   !
   Ntotal(:) = Nirr(:)
   !
   if (.not.present(irr)) return 
   !
   ! Now when we now the sizes we can allocate the matrices with all information 
   ! about the symmetry 
   !
   if (iverbose>=3) write(out,"('Generating the irr. representations...')") 
   !
   if (iverbose>=3) then 
      write(out,"(/'Total number of irr. representations: ',40i10)") Ntotal
   endif 
   !
   if (iverbose>=6) then 
      write(out,"(//'Construct IRREPS for each rotational level')")
   endif 
   !
   Nelem = dimen
   !
   allocate(vect(Nelem,Nelem),t_vect(Nelem,Nelem),stat=alloc)

   !
   do isym = 1,sym%Nrepresen
     !
     if (Nirr(isym)==0) cycle
     !
     Ndeg = sym%degen(isym)
     !
     ! Now we orthogonalize these vectors among themself, removing repeated ones
     !
     felem = 0
     try_elem = 0 
     !
     elem_loop2: do while (try_elem<Nelem.and.felem<Nirr(isym))
       !
       try_elem = try_elem + 1
       !
       do ideg = 1,Ndeg
         do jelem = 1,Nelem
           !
           t_vect(ideg,jelem) =  sum(projector(isym,1:sym%Noper,ideg,ideg)*transform(1:sym%Noper,try_elem,jelem))
           !
         enddo
       enddo
       !
       ielem_t = 0 
       !
       elem_loop: do while (ielem_t<Ndeg.and.felem<Nirr(isym))
         !
         ielem_t = ielem_t + 1
         !
         vect(1,1:Nelem) = t_vect(ielem_t,1:Nelem)
         !
         ! initial check for non-vanish 
         !
         f_t = sum(vect(1,:)**2)
         ! Continue for non-zero vectors
         !
         if (f_t>small_) then 
           !
           vect(1,:) = vect(1,:)/sqrt(f_t)
           !
           else
           !
           cycle elem_loop
           !
         endif
         !
         ! Schmidt orthogonalization
         !
         do ielem = 1,felem
           do  ideg= 1,sym%degen(isym)
             !
             f_t = sum(irr(isym)%repres(ielem,ideg,1:Nelem)*vect(1,1:Nelem))
             !
             vect(1,1:Nelem) = vect(1,1:Nelem)- irr(isym)%repres(ielem,ideg,1:Nelem)*f_t
             !
           enddo
         enddo
         !
         ! normalization 
         !
         f_t = sum(vect(1,:)**2)
         !
         ! Continue for non-zero vectors
         !
         if (f_t>small_) then 
           !
           vect(1,:) = vect(1,:)/sqrt(f_t)
           !
           else
           !
           cycle elem_loop
           !
         endif
         !
         ! Reconstructing other degenerate components
         !
         do jdeg= 2,Ndeg
           !
           vect(jdeg,1:Nelem) = 0
           !
           do ioper = 1,sym%Noper
            !
            do jelem = 1,Nelem
             !
             vect(jdeg,1:Nelem) =  vect(jdeg,1:Nelem) + &
                                   projector(isym,ioper,1,jdeg)*&
                                   transform(ioper,jelem,1:Nelem)*vect(1,jelem)
                                   !sym%irr(isym,ioper)%repres(ielem_t,jdeg)*&
                                   !real(sym%degen(isym),ark)/real(sym%Noper,ark)

            enddo

           enddo
           !
           !vect(jdeg,1:Nelem) = t_vect(jdeg,1:Nelem)
           !
         enddo
         !
         ! normalize 
         !
         do jdeg =2,Ndeg
           !
           f_t = sum(vect(jdeg,:)**2)
           !
           ! Continue for non-zero vectors
           !
           if (f_t>small_) then 
             !
             vect(jdeg,:) = vect(jdeg,:)/sqrt(f_t)
             !
           else
             !
             cycle elem_loop
             !
           endif
           !
           ! Orthogonalization with each other 
           !
           do kdeg = 1,jdeg-1
             !
             f_t = sum(vect(kdeg,:)*vect(jdeg,:))
             !
             vect(jdeg,:) = vect(jdeg,:)-vect(kdeg,:)*f_t
             !
           enddo
           !
         enddo
         !
         ! check if the found vector does transorm correctly with the group
         !         
         do ioper = 1,sym%Noper
           !
           do jdeg = 1,Ndeg
             !
             do jelem = 1,Nelem
               !
               f_t = sum( vect(jdeg,1:Nelem)*transform(ioper,1:Nelem,jelem) )
               !
               g_t = sum( sym%irr(isym,ioper)%repres(jdeg,1:Ndeg)*vect(1:Ndeg,jelem) )
               !
               ! Continue if this two quanttaties are the same 
               !
               if (abs(f_t-g_t)>100.0*sqrt(small_)) then 
                 !
                 continue
                 !
                 cycle elem_loop
                 !
               endif 
               !
             enddo
             !
           enddo
           !
         enddo
         !
         ! check if the found vector is really unique, by comparing also with other symmetries 
         !
         do jsym = 1,isym
           !
           do ielem = 1,min(Nirr(jsym),felem)
             !
             do kdeg =1,Ndeg
               !
               do jdeg =1,sym%degen(jsym)
                 !
                 if (all(abs(vect(kdeg,1:Nelem) - &
                             irr(jsym)%repres(ielem,jdeg,1:Nelem))<sqrt(small_))) cycle elem_loop
                 !
               enddo 
               !
             enddo 
             !
           enddo
           !
         enddo
         !
         ! if we are through all obstacles, we are finally there
         !
         felem = felem + 1
         !
         irr(isym)%repres(felem,1:Ndeg,1:Nelem) = vect(1:Ndeg,1:Nelem)
         !
        enddo elem_loop
        !
     enddo elem_loop2
     !
     ! Check if all representations have been found 
     !
     if (felem/=Nirr(isym)) then 
        write(out,"('degenerate_projectors: Not all irr. representations have been found for')")
        write(out,"(' isym = ',i4,' : only ',i7,' out of ',i7)") isym,felem,Nirr(isym)
        stop 'degenerate_projectors: Not all irr!. representations have been found'
     endif 
     !
  enddo 
  !
  deallocate (t_vect,vect,projector)

  !
  if (iverbose>=6) then 
    !
    do isym =1,sym%Nrepresen
      !
      do i0 = 1,Nirr(isym)
          !
          do ideg = 1,sym%degen(isym)
             write(out,"(' isym =',i4,' ideg =',i4,' irr. repres.: ',30f18.8)") &
                           isym,ideg,irr(isym)%repres(i0,ideg,1:dimen)
          enddo 
        !
      enddo
      !
    enddo
    !
  endif
  !
  end subroutine MLsymmetrization_rot


  function faclog(k)   result (v)
    integer(ik),intent(in) ::  k
    real(ark)              :: v 
    integer(ik) j

    v=0
    if(k>=2) then 
      do j=2,k
         v=v+log(real(j,ark))
      enddo 
    endif 
    
  end function faclog


  subroutine SymmetryInitialize(sym_group)
  character(len=cl),intent(in) :: sym_group
  integer(ik):: alloc,iclass,gamma,ioper,ielem,irepr,Nrot,irep,k,irot,N_Cn,ioper_,icn,NC2
  real(ark)  :: a,b,e,o,p2,p3,p4,p23,p43,phi,phi_n,factor,f_t,mat_t(2,2)
  character(len=4) :: Kchar
  !   
  sym%group=sym_group
  !
  select case(trim(sym_group))

  case("C(M)","C")

    sym%Nrepresen=1
    sym%Noper=1
    sym%Nclasses=1
    sym%CII%Noper = 0

    call simple_arrays_allocation

    sym%characters(1,1)=1.0_rk
    sym%degen=(/1/)
    sym%Nelements=(/1/)
    sym%label=(/'A'/)

    do gamma = 1,sym%Nrepresen
      ioper = 0
      do iclass = 1,sym%Nclasses
        do ielem = 1,sym%Nelements(iclass)
          !
          ioper = ioper + 1
          !
          allocate(sym%irr(gamma,ioper)%repres(sym%degen(gamma),sym%degen(gamma)),stat=alloc)
          !
          if (sym%degen(gamma)==1) then
             sym%irr(gamma,ioper)%repres(1,1)= sym%characters(gamma,iclass)
          endif 
          !
        enddo 
      enddo 
    enddo 


  case("CS(M)","CS")

    sym%Nrepresen=2
    sym%Noper=2
    sym%Nclasses=2
    sym%CII%Noper = 0

    call simple_arrays_allocation

    sym%characters= reshape( &
                                 (/ 1, 1,  & 
                                    1,-1  /),(/2,2/))
    sym%degen=(/1,1/)
    sym%Nelements=(/1,1/)
    sym%label=(/'A''','A"'/)

    call irr_allocation

  case("C2V(M)","C2V")

    sym%Nrepresen=4
    sym%Noper=4
    sym%Nclasses=4
    sym%CII%Noper = 0

    call simple_arrays_allocation

    sym%characters= reshape( &
                                  (/ 1, 1, 1, 1, &   ! A1
                                     1, 1,-1,-1, &   ! A2
                                     1,-1,-1, 1, &   ! B1
                                     1,-1, 1,-1 /),(/4,4/)) ! B2
    sym%degen=(/1,1,1,1/)
    sym%Nelements=(/1,1,1,1/)
    sym%label=(/'A1','A2','B1','B2'/)
    !
    o  = 0.0_ark
    p2 = 0.5_ark*pi
    p3 = 1.5_ark*pi
    !
    sym%euler( 1,:) = 0
    sym%euler( 2,:) = (/pi,o,o /)
    sym%euler( 3,:) = (/o,pi,o/)
    sym%euler( 4,:) = (/p2,pi,p3/)
    !
    call irr_allocation


  case("C2H(M)","C2H")

    sym%Nrepresen=4
    sym%Noper=4
    sym%Nclasses=4
    sym%CII%Noper = 0

    call simple_arrays_allocation

    sym%characters= reshape( &
                                  (/ 1, 1, 1, 1, &   ! Ag
                                     1, 1,-1,-1, &   ! Au
                                     1,-1,-1, 1, &   ! Bg
                                     1,-1, 1,-1 /),(/4,4/)) ! Bu
    sym%degen=(/1,1,1,1/)
    sym%Nelements=(/1,1,1,1/)
    sym%label=(/'Ag','Au','Bg','Bu'/)
    !
    o  = 0.0_ark
    p2 = 0.5_ark*pi
    p3 = 1.5_ark*pi
    !
    sym%euler( 1,:) = 0
    sym%euler( 3,:) = (/p2,pi,p3/)
    sym%euler( 2,:) = (/o,pi,o/)
    sym%euler( 4,:) = (/pi,o,o/)
    !
    ! 1324 and 1234 are working
    !
    call irr_allocation
    !
  case("CS(EM)")

    sym%Nrepresen=4
    sym%Noper=4
    sym%Nclasses=4
    sym%CII%Noper = 0

    call simple_arrays_allocation

    sym%characters= reshape( &     ! E  C2a  E* sab i
                                  (/ 1,  1,  1,  1, &   ! Ag 
                                     1,  1, -1, -1, &   ! Au 
                                     1, -1,  1, -1, &   ! Bg
                                     1, -1, -1,  1/),(/4 ,4/)) ! Bu
                                     

    sym%characters = transpose(sym%characters)
    sym%degen=(/1,1,1,1/)
    sym%Nelements=(/1,1,1,1/)
    sym%label=(/'Ag','Au','Bg','Bu'/)
    !
    o  = 0.0_ark
    p2 = 0.5_ark*pi
    p3 = 1.5_ark*pi
    !
    sym%euler( 1,:) = 0
    sym%euler( 2,:) = (/pi,o,o/) ! Rz
    sym%euler( 3,:) = (/p2,pi,p3/)
    sym%euler( 4,:) = (/o,pi,o/)
    !
    call irr_allocation
    !
  case("G4(M)","G4")

    sym%Nrepresen=4
    sym%Noper=4
    sym%Nclasses=4
    sym%CII%Noper = 0

    call simple_arrays_allocation

    sym%characters= reshape( &       ! E (12)(34)  E*  (12)(34)*
                                  (/ 1, 1, 1, 1, &   ! A+
                                     1, 1,-1,-1, &   ! A-
                                     1,-1,-1, 1, &   ! B+
                                     1,-1, 1,-1 /),(/4,4/)) ! B-
    sym%degen=(/1,1,1,1/)
    sym%Nelements=(/1,1,1,1/)
    sym%label=(/'A+','A-','B-','B+'/)
    o  = 0.0_ark
    p2 = 0.5_ark*pi
    p3 = 1.5_ark*pi
    !
    sym%euler( 1,:) = 0
    sym%euler( 2,:) = (/p2,pi,p3/)
    sym%euler( 3,:) = (/o,pi,o/)
    sym%euler( 4,:) = (/pi,o,o/)
    !
    call irr_allocation
    !
  case("G4(EM)")

    sym%Nrepresen=8
    sym%Noper=8
    sym%Nclasses=8
    sym%CII%Noper = 0

    call simple_arrays_allocation

    sym%characters= reshape( &     ! E   a   b  ab   E' E'a E'b E'ab 
                                  (/ 1,  1,  1,  1,  1,  1,  1,  1, &   ! Ags
                                     1,  1, -1, -1,  1,  1, -1, -1, &   ! Aus
                                     1, -1, -1,  1,  1, -1, -1,  1, &   ! Bgs
                                     1, -1,  1, -1,  1, -1,  1, -1, &   ! Bus
                                     1,  1, -1, -1, -1, -1,  1,  1, &   ! Agd 
                                     1,  1,  1,  1, -1, -1, -1, -1, &   ! Aud
                                     1, -1,  1, -1, -1,  1, -1,  1, &   ! Bgd
                                     1, -1, -1,  1, -1,  1,  1, -1 /),(/8 ,8/)) ! Bud
                                     

    sym%characters = transpose(sym%characters)
    sym%degen=(/1,1,1,1,1,1,1,1/)
    sym%Nelements=(/1,1,1,1,1,1,1,1/)
    sym%label=(/'Ags','Aus','Bgs','Bus','Agd','Aud','Bgd','Bud'/)
    !
    o  = 0.0_ark
    p2 = 0.5_ark*pi
    p3 = 1.5_ark*pi
    !
    sym%euler( 1,:) = 0
    sym%euler( 3,:) = (/p2,pi,p3/)
    sym%euler( 2,:) = (/o,pi,o/)
    sym%euler( 4,:) = (/pi,o,o/)
    sym%euler( 5,:) = (/pi,o,o/)
    sym%euler( 7,:) = (/o,pi,o/)
    sym%euler( 6,:) = (/p2,pi,p3/)
    sym%euler( 8,:) = 0


    !sym%euler( 1,:) = 0
    !sym%euler( 3,:) = (/p2,pi,p3/)
    !sym%euler( 4,:) = (/o,pi,o/)
    !sym%euler( 2,:) = (/pi,o,o/)
    !sym%euler( 7,:) = (/pi,o,o/)
    !sym%euler( 5,:) = (/o,pi,o/)
    !sym%euler( 6,:) = (/p2,pi,p3/)
    !sym%euler( 8,:) = 0


    !
    call irr_allocation

  case("D2H(M)")
    !
    sym%Nrepresen=8
    sym%Noper=8
    sym%Nclasses=8
    sym%CII%Noper = 0

    call simple_arrays_allocation

    sym%characters= reshape( &     ! E  C2a  C2b C2c E* sab sac sbc i
                                  (/ 1,  1,  1,  1,  1,  1,  1,  1, &   ! Ag 
                                     1,  1,  1,  1, -1, -1, -1, -1, &   ! Au 
                                     1,  1, -1, -1, -1, -1,  1,  1, &   ! B1g
                                     1,  1, -1, -1,  1,  1, -1, -1, &   ! B1u
                                     1, -1,  1, -1, -1,  1, -1,  1, &   ! B2g 
                                     1, -1,  1, -1,  1, -1,  1, -1, &   ! B2u
                                     1, -1, -1,  1,  1, -1, -1,  1, &   ! B3g
                                     1, -1, -1,  1, -1,  1,  1, -1 /),(/8 ,8/)) ! B3u
                                     

    sym%characters = transpose(sym%characters)
    sym%degen=(/1,1,1,1,1,1,1,1/)
    sym%Nelements=(/1,1,1,1,1,1,1,1/)
    sym%label=(/'Ag ','Au ','B1g','B1u','B2g','B2u','B3g','B3u'/)
    !
    o  = 0.0_ark
    p2 = 0.5_ark*pi
    p3 = 1.5_ark*pi
    !
    sym%euler( 1,:) = 0
    sym%euler( 2,:) = (/pi,o,o/) ! Rz
    sym%euler( 3,:) = (/o,pi,o/) ! Ry
    sym%euler( 4,:) = (/p2,pi,p3/) ! Rx
    !
    sym%euler( 5,:) = (/p2,pi,p3/)
    sym%euler( 6,:) = (/o,pi,o/)
    sym%euler( 7,:) = (/pi,o,o/)
    sym%euler( 8,:) = 0
    !
    call irr_allocation
    !
  case("D2H")

    sym%Nrepresen=8
    sym%Noper=8
    sym%Nclasses=8
    sym%CII%Noper = 0

    call simple_arrays_allocation

    sym%characters= reshape( &     ! E  C2z  C2y C2x i  sxy sxz syz  
                                  (/ 1,  1,  1,  1,  1,  1,  1,  1, &   ! Ag 
                                     1,  1,  1,  1, -1, -1, -1, -1, &   ! Au 
                                     1,  1, -1, -1,  1,  1, -1, -1, &   ! B1g
                                     1,  1, -1, -1, -1, -1,  1,  1, &   ! B1u
                                     1, -1,  1, -1,  1, -1,  1, -1, &   ! B2g 
                                     1, -1,  1, -1, -1,  1, -1,  1, &   ! B2u
                                     1, -1, -1,  1,  1, -1, -1,  1, &   ! B3g
                                     1, -1, -1,  1, -1,  1,  1, -1 /),(/8 ,8/)) ! B3u
                                     

    sym%characters = transpose(sym%characters)
    sym%degen=(/1,1,1,1,1,1,1,1/)
    sym%Nelements=(/1,1,1,1,1,1,1,1/)
    sym%label=(/'Ag ','Au ','B1g','B1u','B2g','B2u','B3g','B3u'/)
    !
    o  = 0.0_ark
    p2 = 0.5_ark*pi
    p3 = 1.5_ark*pi
    !
    sym%euler( 1,:) = 0
    sym%euler( 2,:) = (/pi,o,o/) ! Rz
    sym%euler( 3,:) = (/o,pi,o/) ! Ry
    sym%euler( 4,:) = (/p2,pi,p3/) ! Rx
    !
    sym%euler( 7,:) = (/p2,pi,p3/)
    sym%euler( 8,:) = (/o,pi,o/)
    sym%euler( 5,:) = (/pi,o,o/)
    sym%euler( 6,:) = 0
    !
    call irr_allocation    
    !
  case("C3V(M)","C3V")
    !
    sym%Nrepresen=3
    sym%Noper=6
    sym%Nclasses=3
    sym%CII%Noper = 0
    !
    call simple_arrays_allocation


    sym%characters= reshape( &      !A1 A2 E
                                  (/ 1, 1, 2, &  
                                     1, 1,-1, &  
                                     1,-1, 0 /),(/3,3/)) 
    sym%degen=(/1,1,2/)
    sym%Nelements=(/1,2,3/)
    sym%label=(/'A1','A2','E '/)
    !
    o  = 0.0_ark
    p2 = 0.5_ark*pi
    p3 = 2.0_ark/3.0_ark*pi
    p4 = 4.0_ark/3.0_ark*pi
    !
    sym%euler( 1,:) = 0
    sym%euler( 2,:) = (/p4,o ,o /)
    sym%euler( 3,:) = (/p3,o ,o /)
    sym%euler( 4,:) = (/o ,pi,o /)
    sym%euler( 5,:) = (/p3,pi,o /)
    sym%euler( 6,:) = (/p4,pi,o /)
    !
    sym%lquant(3) = 1
    !
    call irr_allocation
       !
       sym%irr(3,1)%repres = reshape((/1.0_ark,               0.0_ark, &
                                       0.0_ark,               1.0_ark/),(/2,2/))
       !
       sym%irr(3,2)%repres = reshape((/-0.5_ark,              -0.5_ark*sqrt(3.0_ark),&
                                        0.5_ark*sqrt(3.0_ark),-0.5_ark           /),(/2,2/))
       !
       sym%irr(3,3)%repres = reshape((/-0.5_ark,               0.5_ark*sqrt(3.0_ark),&
                                       -0.5_ark*sqrt(3.0_ark),-0.5_ark           /),(/2,2/))
       !
       sym%irr(3,4)%repres = reshape((/ 1.0_ark,               0.0_ark, &
                                        0.0_ark,              -1.0_ark/),(/2,2/))
       !
       sym%irr(3,5)%repres = reshape((/-0.5_ark,               0.5_ark*sqrt(3.0_ark),&
                                        0.5_ark*sqrt(3.0_ark), 0.5_ark            /),(/2,2/))
       !
       sym%irr(3,6)%repres = reshape((/-0.5_ark,              -0.5_ark*sqrt(3.0_ark),&
                                       -0.5_ark*sqrt(3.0_ark), 0.5_ark             /),(/2,2/))
       !
  case("D3H(M)","D3H")

    sym%Nrepresen=6
    sym%Noper=12
    sym%Nclasses=6
    sym%CII%Noper = 0

    call simple_arrays_allocation


    sym%characters= reshape( &      !A1 A2 E  A1 A2 E
                                  (/ 1, 1, 2, 1, 1, 2, &  
                                     1, 1,-1, 1, 1,-1, &  
                                     1,-1, 0, 1,-1, 0, &          
                                     1, 1, 2,-1,-1,-2, &  
                                     1, 1,-1,-1,-1, 1, &  
                                     1,-1, 0,-1, 1, 0 /),(/6,6/)) 
    sym%degen=(/1,1,2,1,1,2/)
    sym%Nelements=(/1,2,3,1,2,3/)
    sym%label=(/'A1''','A2''','E'' ','A1"','A2"','E" '/)
    !
    o   = 0.0_ark
    p2  = 0.5_ark*pi
    p3  = pi/3.0_ark
    p23 = 2.0_ark/3.0_ark*pi
    p43 = 4.0_ark/3.0_ark*pi
    !
    sym%euler( 1,:) = 0               ! (E)
    sym%euler( 2,:) = (/p23,o ,o /)   ! (132)
    sym%euler( 3,:) = (/p43,o ,o /)   ! (123)
    sym%euler( 4,:) = (/ pi,pi,o /)   ! (23)
    sym%euler( 5,:) = (/ p3,pi,o /)   ! (12)
    sym%euler( 6,:) = (/-p3,pi,o /)   ! (13)
    sym%euler( 7,:) = (/ pi, o,o /)   ! (E)*
    sym%euler( 8,:) = (/-p3, o,o /)   ! (132)*
    sym%euler( 9,:) = (/ p3, o,o /)   ! (123)*
    sym%euler(10,:) = (/  o,pi,o /)   ! (23)*
    sym%euler(11,:) = (/p43,pi,o /)   ! (12)*
    sym%euler(12,:) = (/p23,pi,o /)   ! (13)*
    !
    call irr_allocation
    !
    a = 0.5_ark ; b = 0.5_ark*sqrt(3.0_ark) ; e = 1.0_ark ; o = 0.0_ark
    !
    sym%irr(3,1)%repres = reshape((/ e, o,  &
                                     o, e/),(/2,2/))
    !
    sym%irr(3,2)%repres = reshape((/-a, b,  &
                                    -b,-a/),(/2,2/))
    !
    sym%irr(3,3)%repres = reshape((/-a,-b,  &
                                     b,-a/),(/2,2/))
    !
    sym%irr(3,4)%repres = reshape((/ e, o,  &
                                     o,-e/),(/2,2/))
    !
    sym%irr(3,5)%repres = reshape((/-a,-b,  &
                                    -b, a /),(/2,2/))
    !
    sym%irr(3,6)%repres = reshape((/-a, b,  &
                                     b, a  /),(/2,2/))
    !
    !
    sym%irr(3,7)%repres = reshape((/ e, o, &
                                     o, e/),(/2,2/))
    !
    sym%irr(3,8)%repres = reshape((/-a, b, &
                                    -b,-a/),(/2,2/))
    !
    sym%irr(3,9)%repres = reshape((/-a,-b, &
                                     b,-a/),(/2,2/))
    !
    sym%irr(3,10)%repres= reshape((/ e, o, &
                                     o,-e/),(/2,2/))
    !
    sym%irr(3,11)%repres= reshape((/-a,-b, &
                                    -b, a /),(/2,2/))
    !
    sym%irr(3,12)%repres= reshape((/-a, b, &
                                     b, a  /),(/2,2/))
    !
    sym%irr(6,1)%repres = reshape((/e, o,  &
                                    o, e/),(/2,2/))
    !
    sym%irr(6,2)%repres = reshape((/-a, b, &
                                    -b,-a/),(/2,2/))
    !
    sym%irr(6,3)%repres = reshape((/-a,-b, &
                                     b,-a/),(/2,2/))
    !
    sym%irr(6,4)%repres = reshape((/-e, o, &
                                     o, e/),(/2,2/))
    !
    sym%irr(6,5)%repres = reshape((/ a, b, &
                                     b,-a /),(/2,2/))
    !
    sym%irr(6,6)%repres = reshape((/ a,-b, &
                                    -b,-a  /),(/2,2/))
    !
    sym%irr(6,7)%repres = reshape((/-e, o, &
                                     o,-e/),(/2,2/))
    !
    sym%irr(6,8)%repres = reshape((/ a,-b, &
                                     b, a/),(/2,2/))
    !
    sym%irr(6,9)%repres = reshape((/ a, b, &
                                    -b, a/),(/2,2/))
    !
    sym%irr(6,10)%repres= reshape((/ e, o, &
                                     o,-e/),(/2,2/))
    !
    sym%irr(6,11)%repres= reshape((/-a,-b, &
                                    -b, a /),(/2,2/))
    !
    sym%irr(6,12)%repres= reshape((/-a, b, &
                                     b, a  /),(/2,2/))
    !
    sym%lquant(1:2) = 0
    sym%lquant(4:5) = 0
    sym%lquant(3) = 1
    sym%lquant(6) = 1
    !
  case("TD(M)","TD")

    sym%Nrepresen=5
    sym%Noper=24
    sym%Nclasses=5
    sym%CII%Noper = 6

    call simple_arrays_allocation

    sym%CII%coeff = (/1.0_ark,1.0_ark,4.0_ark,10.0_ark,1.0_ark,1.0_ark/)
    sym%CII%ioper = (/19,20,21,22,23,24/)



    sym%characters= reshape( &      !A1 A2 E  F1 F2
                                      (/ 1, 1, 2, 3, 3, &  
                                         1, 1,-1, 0, 0, &  
                                         1, 1, 2,-1,-1, &
                                         1,-1, 0, 1,-1, &
                                         1,-1, 0,-1, 1  /),(/5,5/)) 
    sym%degen=(/1,1,2,3,3/)
    sym%Nelements=(/1,8,3,6,6/)
    sym%label=(/'A1','A2','E ','F1','F2'/)

    sym%CII%Nzeta = sum(sym%degen(:))
    allocate(sym%CII%izeta(sym%CII%Nzeta),stat=alloc)
    !
    sym%CII%izeta = (/18,-18,12,-12,-14,-8,4,14,-4,8/)
    !
    o  = 0.0_ark
    p2 = 0.5_ark*pi
    p3 = 1.5_ark*pi
    !
    sym%euler( 1,:) = 0
    sym%euler( 2,:) = (/p3,p2,o /)
    sym%euler( 3,:) = (/pi,p2,p3/)
    sym%euler( 4,:) = (/ o,p2,p3/)
    sym%euler( 5,:) = (/p3,p2,pi/)
    sym%euler( 6,:) = (/p2,p2,o /)
    sym%euler( 7,:) = (/pi,p2,p2/)
    sym%euler( 8,:) = (/ o,p2,p2/)
    sym%euler( 9,:) = (/p2,p2,pi/)
    sym%euler(10,:) = (/p2,pi,p3/)
    sym%euler(11,:) = (/ o,pi,o /)
    sym%euler(12,:) = (/pi, o,o /)
    sym%euler(13,:) = (/ o,p3,o /)
    sym%euler(14,:) = (/ o,p2,o /)
    sym%euler(15,:) = (/p2, o,o /)
    sym%euler(16,:) = (/p3, o,o /)
    sym%euler(17,:) = (/p2,p3,p3/)
    sym%euler(18,:) = (/p2,p2,p3/)
    sym%euler(19,:) = (/p2,p2,p2/)
    sym%euler(20,:) = (/ o,p2,pi/)
    sym%euler(21,:) = (/ o,pi,p2/)
    sym%euler(22,:) = (/pi,pi,p2/)
    sym%euler(23,:) = (/pi,p2,o /)
    sym%euler(24,:) = (/p3,p2,p3/)
    !
    !sym%euler(:,1) = sym%euler(:,1)+pi*0.25_ark
    !sym%euler(:,2) = sym%euler(:,2)
    !sym%euler(:,3) = sym%euler(:,3)-pi*0.25_ark
    !
    call irr_allocation
    !
    sym%irr(3,1)%repres = reshape((/ 1.0_ark,               0.0_ark, &
                                     0.0_ark,               1.0_ark/),(/2,2/))
    !
    !sym%irr(3,2)%repres = reshape((/-0.5_ark,               0.5_ark*sqrt(3.0_ark),&
    !                                -0.5_ark*sqrt(3.0_ark),-0.5_ark           /),(/2,2/))
    ! working 
    sym%irr(3,2)%repres = reshape((/-0.5_ark,              -0.5_ark*sqrt(3.0_ark),&
                                     0.5_ark*sqrt(3.0_ark),-0.5_ark             /),(/2,2/))
    !
    !sym%irr(3,2)%repres = reshape((/-0.5_ark,              -0.5_ark*sqrt(3.0_ark),&
    !                                 0.5_ark*sqrt(3.0_ark),-0.5_ark           /),(/2,2/))
    !
    ! 
    !sym%irr(3,21)%repres = reshape((/-0.5_ark,              0.5_ark*sqrt(3.0_ark),&
    !                                  0.5_ark*sqrt(3.0_ark),0.5_ark             /),(/2,2/))
    !working    
    sym%irr(3,21)%repres = reshape((/-0.5_ark,              0.5_ark*sqrt(3.0_ark),&
                                      0.5_ark*sqrt(3.0_ark),0.5_ark             /),(/2,2/))

    !
    sym%irr(4,1)%repres = reshape((/ 1.0_ark, 0.0_ark, 0.0_ark, &
                                     0.0_ark, 1.0_ark, 0.0_ark, &
                                     0.0_ark, 0.0_ark, 1.0_ark/),(/3,3/))
    !
    !sym%irr(4,2)%repres = reshape((/ 0.0_ark,-1.0_ark, 0.0_ark, &
    !                                 0.0_ark, 0.0_ark,-1.0_ark, &
    !                                 1.0_ark, 0.0_ark, 0.0_ark/),(/3,3/))
    !
    !!sym%irr(4,2)%repres = reshape((/ 0.0_ark, 0.0_ark,-1.0_ark, &
    !!                                 1.0_ark, 0.0_ark, 0.0_ark, &
    !!                                 0.0_ark,-1.0_ark, 0.0_ark/),(/3,3/))
    ! working ->
    !
    sym%irr(4,2)%repres = reshape((/ 0.0_ark, 0.0_ark,-1.0_ark, &
                                     1.0_ark, 0.0_ark, 0.0_ark, &
                                     0.0_ark,-1.0_ark, 0.0_ark/),(/3,3/))

    !sym%irr(4,2)%repres = reshape((/ 0.0_ark, 0.0_ark, 1.0_ark, &
    !                                 1.0_ark, 0.0_ark, 0.0_ark, &
    !                                 0.0_ark, 1.0_ark, 0.0_ark/),(/3,3/))


    !sym%irr(4,2)%repres = transpose(reshape((/ 0.0_ark, 0.5_ark*sqrt(2.0_ark), 0.5_ark*sqrt(2.0_ark), &
    !                                 0.5_ark*sqrt(2.0_ark), 0.5_ark,-0.5_ark, &
    !                                -0.5_ark*sqrt(2.0_ark), 0.5_ark,-0.5_ark/),(/3,3/)))
    !
    !
    !sym%irr(4,21)%repres = reshape((/-1.0_ark, 0.0_ark, 0.0_ark, &
    !                                  0.0_ark, 0.0_ark,-1.0_ark, &
    !                                  0.0_ark,-1.0_ark, 0.0_ark/),(/3,3/))

    ! working->
    sym%irr(4,21)%repres = reshape((/  0.0_ark, 1.0_ark, 0.0_ark, &
                                       1.0_ark, 0.0_ark, 0.0_ark, &
                                       0.0_ark, 0.0_ark,-1.0_ark/),(/3,3/))


    !
    ! roman
    !

     !  sym%irr(4,2)%repres = reshape((/ 0.0_ark, 0.0_ark,-1.0_ark, &
     !                                   1.0_ark, 0.0_ark, 0.0_ark, &
     !                                   0.0_ark,-1.0_ark, 0.0_ark/),(/3,3/))
     !  !
     !  sym%irr(4,21)%repres = reshape((/ 0.0_ark, 0.0_ark,-1.0_ark, &
     !                                    0.0_ark,-1.0_ark, 0.0_ark, &
     !                                   -1.0_ark, 0.0_ark, 0.0_ark/),(/3,3/))




    !
    !sym%irr(4,21)%repres = reshape((/-1.0_ark, 0.0_ark, 0.0_ark, &
    !                                  0.0_ark, 1.0_ark, 0.0_ark, &
    !                                  0.0_ark, 0.0_ark,-1.0_ark/),(/3,3/))
    
    !sym%irr(4,21)%repres = reshape((/ 0.0_ark,-1.0_ark, 0.0_ark, &
    !                                 -1.0_ark, 0.0_ark, 0.0_ark, &
    !!                                  0.0_ark, 0.0_ark, -1.0_ark/),(/3,3/))
    !


    !
    !!!sym%irr(4,21)%repres = reshape((/ 0.0_ark, 0.0_ark,-1.0_ark, &
    !!!                                  0.0_ark,-1.0_ark, 0.0_ark, &
    !!!                                 -1.0_ark, 0.0_ark, 0.0_ark/),(/3,3/))
    !
    sym%irr(5,1)%repres = reshape((/ 1.0_ark, 0.0_ark, 0.0_ark, &
                                     0.0_ark, 1.0_ark, 0.0_ark, &
                                     0.0_ark, 0.0_ark, 1.0_ark/),(/3,3/))
    ! working
    !sym%irr(5,2)%repres = reshape((/ 0.0_ark,-1.0_ark, 0.0_ark, &
    !                                 0.0_ark, 0.0_ark,-1.0_ark, &
    !                                 1.0_ark, 0.0_ark, 0.0_ark/),(/3,3/))
    !
    sym%irr(5,2)%repres = reshape((/ 0.0_ark, 0.0_ark, 1.0_ark, &
                                    -1.0_ark, 0.0_ark, 0.0_ark, &
                                     0.0_ark,-1.0_ark, 0.0_ark/),(/3,3/))
    ! working
    !sym%irr(5,21)%repres = reshape((/ 0.0_ark, 0.0_ark,-1.0_ark, &
    !                                  0.0_ark, 1.0_ark, 0.0_ark, &
    !                                 -1.0_ark, 0.0_ark, 0.0_ark/),(/3,3/))
    !
    sym%irr(5,21)%repres = reshape((/ 0.0_ark, 0.0_ark,-1.0_ark, &
                                      0.0_ark, 1.0_ark, 0.0_ark, &
                                     -1.0_ark, 0.0_ark, 0.0_ark/),(/3,3/))
    !
    !sym%irr(5,2)%repres = reshape((/ 0.0_ark, 0.0_ark,  1.0_ark, &
    !                                -1.0_ark, 0.0_ark, 0.0_ark, &
    !                                 0.0_ark,-1.0_ark, 0.0_ark/),(/3,3/))


    !
    ! (2)  = (123)
    ! (21) = (14)*
    do irepr=3,5
      sym%irr(irepr,3)%repres = matmul(sym%irr(irepr,2)%repres,sym%irr(irepr,2)%repres)   ! (132)
      sym%irr(irepr,13)%repres= matmul(sym%irr(irepr,21)%repres,sym%irr(irepr,2)%repres)  ! (1423)*
      sym%irr(irepr,8)%repres = matmul(sym%irr(irepr,13)%repres,sym%irr(irepr,21)%repres) ! (234)
      sym%irr(irepr,5)%repres = matmul(sym%irr(irepr,8)%repres,sym%irr(irepr,3)%repres)   ! (134)
      sym%irr(irepr,4)%repres = matmul(sym%irr(irepr,5)%repres,sym%irr(irepr,5)%repres)   ! (143)
      sym%irr(irepr,6)%repres = matmul(sym%irr(irepr,4)%repres,sym%irr(irepr,3)%repres)   ! (142)
      sym%irr(irepr,7)%repres = matmul(sym%irr(irepr,6)%repres,sym%irr(irepr,6)%repres)   ! (124)
      sym%irr(irepr,9)%repres = matmul(sym%irr(irepr,8)%repres,sym%irr(irepr,8)%repres)   ! (243)
      sym%irr(irepr,10)%repres= matmul(sym%irr(irepr,7)%repres,sym%irr(irepr,2)%repres)   ! (13)(24)
      sym%irr(irepr,11)%repres= matmul(sym%irr(irepr,8)%repres,sym%irr(irepr,2)%repres)   ! (12)(34)
      sym%irr(irepr,12)%repres= matmul(sym%irr(irepr,4)%repres,sym%irr(irepr,2)%repres)   ! (14)(23)
      sym%irr(irepr,14)%repres= matmul(sym%irr(irepr,3)%repres,sym%irr(irepr,21)%repres)  ! (1324)*
      sym%irr(irepr,15)%repres= matmul(sym%irr(irepr,11)%repres,sym%irr(irepr,21)%repres) ! (1243)*
      sym%irr(irepr,16)%repres= matmul(sym%irr(irepr,10)%repres,sym%irr(irepr,21)%repres) ! (1342)*
      sym%irr(irepr,17)%repres= matmul(sym%irr(irepr,9)%repres,sym%irr(irepr,21)%repres)  ! (1432)*
      sym%irr(irepr,18)%repres= matmul(sym%irr(irepr,2)%repres,sym%irr(irepr,21)%repres)  ! (1234)*
      sym%irr(irepr,19)%repres= matmul(sym%irr(irepr,5)%repres,sym%irr(irepr,21)%repres)  ! (13)* 
      sym%irr(irepr,20)%repres= matmul(sym%irr(irepr,7)%repres,sym%irr(irepr,21)%repres)  ! (12)*
      sym%irr(irepr,22)%repres= matmul(sym%irr(irepr,12)%repres,sym%irr(irepr,21)%repres) ! (23)*
      sym%irr(irepr,23)%repres= matmul(sym%irr(irepr,4)%repres,sym%irr(irepr,21)%repres)  ! (34)*
      sym%irr(irepr,24)%repres= matmul(sym%irr(irepr,6)%repres,sym%irr(irepr,21)%repres)  ! (24)*
    enddo
    sym%lquant(1:2) = 0
    sym%lquant(3) = 1
    sym%lquant(4:5) = 2
    !
  case("DNH(M)","DNH") ! D_infinity_H(M)
    !
    if (mod(sym%N,2)==1) then
       !  write(out,"('symmetry: currently Dnh is only for an even ggroup order N ',i8)") sym%N
       !  stop 'symmetry: illegal order of Dnh group '
       !endif
       !
       ! Number of rotations to test for < infinity 
       !
       Nrot = sym%N ! must be >=1
       !
       ! Number of Cn classes 
       N_Cn = sym%N/2
       !
       sym%Noper=2+4*N_Cn+2*Nrot
       sym%Nclasses=4+N_Cn*2
       sym%Nrepresen= 4+N_Cn*2
       sym%CII%Noper = 0
       !
       phi = 2.0_ark*pi/real(Nrot,ark)
       !
       call simple_arrays_allocation
       !
       ! E nrotxCinf nrotxsigmav i  nrotxSinf nrotxC'2
       !
       sym%characters(:,:) = 0
       !
       ! A1g,A1u,A2g,A2u:
       ! E
       sym%characters(1:4,1) = 1
       ! Cinf
       sym%characters(1:4,1+1:1+N_Cn) = 1
       ! sigmav
       sym%characters(1,1+N_Cn+1) = 1
       sym%characters(2,1+N_Cn+1) =-1
       sym%characters(3,1+N_Cn+1) = 1
       sym%characters(4,1+N_Cn+1) =-1
       ! i
       sym%characters(1,1+N_Cn+2) = 1
       sym%characters(2,1+N_Cn+2) = 1
       sym%characters(3,1+N_Cn+2) =-1
       sym%characters(4,1+N_Cn+2) =-1
       ! Sinf
       sym%characters(1,1+N_Cn+2+1:1+N_Cn+2+N_Cn) = 1
       sym%characters(2,1+N_Cn+2+1:1+N_Cn+2+N_Cn) = 1
       sym%characters(3,1+N_Cn+2+1:1+N_Cn+2+N_Cn) =-1
       sym%characters(4,1+N_Cn+2+1:1+N_Cn+2+N_Cn) =-1
       ! C'inf
       sym%characters(1,1+N_Cn+2+N_Cn+1) = 1
       sym%characters(2,1+N_Cn+2+N_Cn+1) =-1
       sym%characters(3,1+N_Cn+2+N_Cn+1) =-1
       sym%characters(4,1+N_Cn+2+N_Cn+1) = 1
       !
       !sym%label(1:4)=(/'A1g','A2g','A1u','A2u'/)
       !
       sym%label(1:4)=(/'A1''','A2''','A1"','A2"'/)
       !
       !sym%label=(/'A1''','A2''','E'' ','A1"','A2"','E" '/)
       !
       ! E1' E1" E2' E2" E3' E3" ....
       !
       sym%lquant(1:4) = 0 
       !
       irep = 4
       do k = 1,(sym%Nrepresen-4)/2
         !
         irep = irep + 1
         !
         write(Kchar, '(i4)') K
         !
         sym%label(irep  ) = 'E'//trim(adjustl(Kchar))//''''
         sym%label(irep+1) = 'E'//trim(adjustl(Kchar))//'"'
         !
         sym%characters(irep  ,1) = 2.0_ark
         sym%characters(irep+1,1) = 2.0_ark
         !
         sym%characters(irep  ,1+N_Cn+2) = 2.0_ark
         sym%characters(irep+1,1+N_Cn+2) =-2.0_ark
         !
         do irot = 1,N_Cn
           sym%characters(irep  ,1+irot)          = 2.0_ark*cos(phi*irot*k)
           sym%characters(irep+1,1+irot)          = 2.0_ark*cos(phi*irot*k)
           sym%characters(irep  ,1+N_Cn+2+irot)   = 2.0_ark*cos(phi*irot*k)
           sym%characters(irep+1,1+N_Cn+2+irot)   =-2.0_ark*cos(phi*irot*k)
           !
           sym%lquant(irep  ) = k
           sym%lquant(irep+1) = k
           !
         enddo
         !
         irep = irep + 1
         !
       enddo
       !
       sym%degen(:)   = 2
       sym%degen(1:4) = 1
       !
       sym%Nelements(1) = 1
       sym%Nelements(1+ 1:1+ N_Cn) = 2
       sym%Nelements(1+N_Cn+1) = Nrot
       sym%Nelements(1+N_Cn+2) = 1
       sym%Nelements(1+N_Cn+2+1:1+N_Cn+2+N_Cn) = 2
       sym%Nelements(1+N_Cn+2+N_Cn+1) = Nrot
       !
       o  = 0.0_ark
       p2 = 0.5_ark*pi
       p3 = 1.5_ark*pi
       !
       sym%euler(:,:) = 0
       !
       ioper = 1
       do irot = 1,N_Cn
         !
         sym%euler(1+ioper  ,:) = (/o, phi*irot,o/) ! Rz
         sym%euler(1+ioper+1,:) = (/o,-phi*irot,o/) ! Rz
         sym%euler(1+2*N_Cn+Nrot+1+ioper,:)   = (/o,pi+phi*irot,o/) ! Rz
         sym%euler(1+2*N_Cn+Nrot+1+ioper+1,:) = (/o,pi-phi*irot,o/) ! Rz
         !
         ioper = ioper + 2
       enddo
       !
       call irr_allocation
       !
       ! Generate irr-representations
       !
       irep = 4
       do k = 1,(sym%Nrepresen-4)/2
         !
         irep = irep + 1
         !
         ioper = 1
         do ioper = 1,sym%Noper
           !
           factor = 1.0_ark
           !
           if (ioper==1) then ! E 
             !
             !
             sym%irr(irep,ioper)%repres(1,1) = 1.0_ark
             sym%irr(irep,ioper)%repres(1,2) = 0.0_ark
             !
             sym%irr(irep,ioper)%repres(2,1) = 0.0_ark
             sym%irr(irep,ioper)%repres(2,2) = 1.0_ark
             !
             sym%irr(irep+1,ioper)%repres(1,1) = 1.0_ark
             sym%irr(irep+1,ioper)%repres(1,2) = 0.0_ark
             !
             sym%irr(irep+1,ioper)%repres(2,1) = 0.0_ark
             sym%irr(irep+1,ioper)%repres(2,2) = 1.0_ark
             !
           elseif (ioper<=1+2*N_Cn) then !  Cinf
             !
             ioper_ =ioper-1 
             irot = (ioper_+1)/2
             !
             phi_n = phi*irot*k
             !
             ! Second oper in a class is with negative phi
             if ( mod(ioper_,2)==0 ) phi_n = -phi_n
             !
             sym%irr(irep,ioper)%repres(1,1) = cos(phi_n)
             sym%irr(irep,ioper)%repres(1,2) =-sin(phi_n)
             !
             sym%irr(irep,ioper)%repres(2,1) = sin(phi_n)
             sym%irr(irep,ioper)%repres(2,2) = cos(phi_n)
             !
             sym%irr(irep+1,ioper)%repres(1,1) = cos(phi_n)
             sym%irr(irep+1,ioper)%repres(1,2) =-sin(phi_n)
             !
             sym%irr(irep+1,ioper)%repres(2,1) = sin(phi_n)
             sym%irr(irep+1,ioper)%repres(2,2) = cos(phi_n)
             !
           elseif (ioper<=1+2*N_Cn+Nrot) then !  sigmav
             !
             irot = ioper-(1+2*N_Cn)
             !
             phi_n = phi*irot*k
             !
             sym%irr(irep,ioper)%repres(1,1) = cos(phi_n)
             sym%irr(irep,ioper)%repres(1,2) = sin(phi_n)
             !
             sym%irr(irep,ioper)%repres(2,1) = sin(phi_n)
             sym%irr(irep,ioper)%repres(2,2) =-cos(phi_n)
             !
             sym%irr(irep+1,ioper)%repres(1,1) = cos(phi_n)
             sym%irr(irep+1,ioper)%repres(1,2) =-sin(phi_n)
             !
             sym%irr(irep+1,ioper)%repres(2,1) =-sin(phi_n)
             sym%irr(irep+1,ioper)%repres(2,2) =-cos(phi_n)
             !
           elseif (ioper==1+2*N_Cn+Nrot+1) then ! i
             !
             phi_n = 0
             factor = -1.0_ark
             !
             sym%irr(irep,ioper)%repres(1,1) = 1.0_ark
             sym%irr(irep,ioper)%repres(1,2) = 0.0_ark
             !
             sym%irr(irep,ioper)%repres(2,1) = 0.0_ark
             sym%irr(irep,ioper)%repres(2,2) = 1.0_ark
             !
             sym%irr(irep+1,ioper)%repres(1,1) =-1.0_ark
             sym%irr(irep+1,ioper)%repres(1,2) = 0.0_ark
             !
             sym%irr(irep+1,ioper)%repres(2,1) = 0.0_ark
             sym%irr(irep+1,ioper)%repres(2,2) =-1.0_ark
             !
           elseif (ioper<=1+2*N_Cn+Nrot+1+2*N_Cn) then !  Sinf
             !
             ioper_ = ioper-(1+2*N_Cn+Nrot+1)
             !
             irot = (ioper_+1)/2
             !
             phi_n = phi*irot*k
             !
             ! Second oper in a class is with negative phi
             if (mod(ioper_,2)==0)  phi_n = -phi_n
             !
             sym%irr(irep,ioper)%repres(1,1) = cos(phi_n)
             sym%irr(irep,ioper)%repres(1,2) =-sin(phi_n)
             !
             sym%irr(irep,ioper)%repres(2,1) = sin(phi_n)
             sym%irr(irep,ioper)%repres(2,2) = cos(phi_n)
             !
             sym%irr(irep+1,ioper)%repres(1,1) =-cos(phi_n)
             sym%irr(irep+1,ioper)%repres(1,2) = sin(phi_n)
             !
             sym%irr(irep+1,ioper)%repres(2,1) =-sin(phi_n)
             sym%irr(irep+1,ioper)%repres(2,2) =-cos(phi_n)
             !
           elseif (ioper<=1+2*N_Cn+Nrot+1+2*N_Cn+Nrot) then !  C'2
             !
             irot = ioper-(1+2*N_Cn+Nrot+1+2*N_Cn)
             !
             phi_n  = phi*irot*k
             !
             sym%irr(irep,ioper)%repres(1,1) = cos(phi_n)
             sym%irr(irep,ioper)%repres(1,2) = sin(phi_n)
             !
             sym%irr(irep,ioper)%repres(2,1) = sin(phi_n)
             sym%irr(irep,ioper)%repres(2,2) =-cos(phi_n)
             !
             sym%irr(irep+1,ioper)%repres(1,1) =-cos(phi_n)
             sym%irr(irep+1,ioper)%repres(1,2) =-sin(phi_n)
             !
             sym%irr(irep+1,ioper)%repres(2,1) =-sin(phi_n)
             sym%irr(irep+1,ioper)%repres(2,2) = cos(phi_n)
             !
           else
             !
             stop  'symmetry: illegal ioper'
             !
           endif
           !
         enddo
         !
         irep = irep + 1
         !
       enddo
       !
    else  
       ! sym%N is even
       ! 
       Nrot = sym%N  ! Number of equivalent rotations
       NC2 = sym%N/2 ! Number of orthog. C2' axes
       !
       ! Number of Cn classes without C2
       N_Cn = sym%N/2-1
       !
       ! 1xE, 2xN_CnxCn, C2, C2', C2" ... 
       !
       sym%Noper=2*(1+2*N_Cn+1+2*NC2)
       sym%Nclasses=8+N_Cn*2
       sym%Nrepresen= 8+N_Cn*2
       sym%CII%Noper = 0
       !
       f_t = (-1.0_ark)**(Nrot/2)
       !
       phi = 2.0_ark*pi/real(Nrot,ark)
       !
       call simple_arrays_allocation
       !
       ! E nrotxCinf nrotxsigmav i  nrotxSinf nrotxC'2
       !
       sym%characters(:,:) = 0.0_ark
       !
       ! A1g,A2g,B1g,B2g,A1u,A2u,B1u,B2u:
       ! E
       sym%characters(1:8,1) = 1.0_ark
       !
       ! Cn - rotations
       !
       do icn=1,N_Cn
         !
         sym%characters(1,1+icn) = 1.0_ark
         sym%characters(2,1+icn) = 1.0_ark
         sym%characters(3,1+icn) =(-1.0_ark)**icn
         sym%characters(4,1+icn) =(-1.0_ark)**icn
         sym%characters(5,1+icn) = 1.0_ark
         sym%characters(6,1+icn) = 1.0_ark
         sym%characters(7,1+icn) =(-1.0_ark)**icn
         sym%characters(8,1+icn) =(-1.0_ark)**icn
         !
       enddo
       !
       ! C2
       sym%characters(1,2+N_Cn) = 1.0_ark
       sym%characters(2,2+N_Cn) = 1.0_ark
       sym%characters(3,2+N_Cn) = f_t
       sym%characters(4,2+N_Cn) = f_t
       sym%characters(5,2+N_Cn) = 1.0_ark
       sym%characters(6,2+N_Cn) = 1.0_ark
       sym%characters(7,2+N_Cn) = f_t
       sym%characters(8,2+N_Cn) = f_t
       ! nxC2'
       sym%characters(1,3+N_Cn) = 1.0_ark
       sym%characters(2,3+N_Cn) =-1.0_ark
       sym%characters(3,3+N_Cn) = 1.0_ark
       sym%characters(4,3+N_Cn) =-1.0_ark
       sym%characters(5,3+N_Cn) = 1.0_ark
       sym%characters(6,3+N_Cn) =-1.0_ark
       sym%characters(7,3+N_Cn) = 1.0_ark
       sym%characters(8,3+N_Cn) =-1.0_ark
       ! nxC2"
       sym%characters(1,4+N_Cn) = 1.0_ark
       sym%characters(2,4+N_Cn) =-1.0_ark
       sym%characters(3,4+N_Cn) =-1.0_ark
       sym%characters(4,4+N_Cn) = 1.0_ark
       sym%characters(5,4+N_Cn) = 1.0_ark
       sym%characters(6,4+N_Cn) =-1.0_ark
       sym%characters(7,4+N_Cn) =-1.0_ark
       sym%characters(8,4+N_Cn) = 1.0_ark
       ! i
       sym%characters(1,5+N_Cn) = 1.0_ark
       sym%characters(2,5+N_Cn) = 1.0_ark
       sym%characters(3,5+N_Cn) = 1.0_ark
       sym%characters(4,5+N_Cn) = 1.0_ark
       sym%characters(5,5+N_Cn) =-1.0_ark
       sym%characters(6,5+N_Cn) =-1.0_ark
       sym%characters(7,5+N_Cn) =-1.0_ark
       sym%characters(8,5+N_Cn) =-1.0_ark
       !
       ! Sn
       do icn=1,N_Cn
         sym%characters(1,5+N_Cn+icn) = 1.0_ark
         sym%characters(2,5+N_Cn+icn) = 1.0_ark
         sym%characters(3,5+N_Cn+icn) =(-1.0_ark)**icn*f_t
         sym%characters(4,5+N_Cn+icn) =(-1.0_ark)**icn*f_t
         sym%characters(5,5+N_Cn+icn) =-1.0_ark
         sym%characters(6,5+N_Cn+icn) =-1.0_ark
         sym%characters(7,5+N_Cn+icn) =-(-1.0_ark)**icn*f_t
         sym%characters(8,5+N_Cn+icn) =-(-1.0_ark)**icn*f_t
       enddo
       ! sigmah
       sym%characters(1,6+2*N_Cn) = 1.0_ark
       sym%characters(2,6+2*N_Cn) = 1.0_ark
       sym%characters(3,6+2*N_Cn) = f_t
       sym%characters(4,6+2*N_Cn) = f_t
       sym%characters(5,6+2*N_Cn) =-1.0_ark
       sym%characters(6,6+2*N_Cn) =-1.0_ark
       sym%characters(7,6+2*N_Cn) =-f_t
       sym%characters(8,6+2*N_Cn) =-f_t
       ! sigmav
       sym%characters(1,7+2*N_Cn) = 1.0_ark
       sym%characters(2,7+2*N_Cn) =-1.0_ark
       sym%characters(3,7+2*N_Cn) = f_t
       sym%characters(4,7+2*N_Cn) =-f_t
       sym%characters(5,7+2*N_Cn) =-1.0_ark
       sym%characters(6,7+2*N_Cn) = 1.0_ark
       sym%characters(7,7+2*N_Cn) =-f_t
       sym%characters(8,7+2*N_Cn) = f_t
       ! sigmad
       sym%characters(1,8+2*N_Cn) = 1.0_ark
       sym%characters(2,8+2*N_Cn) =-1.0_ark
       sym%characters(3,8+2*N_Cn) =-f_t
       sym%characters(4,8+2*N_Cn) = f_t
       sym%characters(5,8+2*N_Cn) =-1.0_ark
       sym%characters(6,8+2*N_Cn) = 1.0_ark
       sym%characters(7,8+2*N_Cn) = f_t
       sym%characters(8,8+2*N_Cn) =-f_t
       !
       !sym%label(1:4)=(/'A1g','A2g','A1u','A2u'/)
       !
       sym%label(1:8)=(/'A1g','A2g','B1g','B2g','A1u','A2u','B1u','B2u'/)
       !
       ! E1g E1u E2g E2u E3g E3u ....
       !
       sym%lquant(1:8) = 0 
       !
       irep = 8
       do k = 1,(sym%Nrepresen-8)/2
         !
         irep = irep + 1
         !
         write(Kchar, '(i4)') K
         !
         sym%label(irep  ) = 'E'//trim(adjustl(Kchar))//'g'
         sym%label(irep+1) = 'E'//trim(adjustl(Kchar))//'u'
         !
         sym%characters(irep  ,1) = 2.0_ark
         sym%characters(irep+1,1) = 2.0_ark
         !
         sym%characters(irep  ,5+N_Cn) = 2.0_ark
         sym%characters(irep+1,5+N_Cn) =-2.0_ark
         !
         do icn = 1,N_Cn
           !
           sym%characters(irep  ,1+icn)          = 2.0_ark*cos(phi*icn*k)
           sym%characters(irep+1,1+icn)          = 2.0_ark*cos(phi*icn*k)
           !
           sym%characters(irep  ,5+N_Cn+icn)   =-2.0_ark*cos(phi*icn*k)
           sym%characters(irep+1,5+N_Cn+icn)   = 2.0_ark*cos(phi*icn*k)
           !
           sym%lquant(irep  ) = k
           sym%lquant(irep+1) = k
           !
         enddo
         !
         sym%characters(irep  ,2+N_Cn) = 2.0_ark*(-1)**irep
         sym%characters(irep+1,2+N_Cn) = 2.0_ark*(-1)**irep
         !
         sym%characters(irep  ,6+2*N_Cn) = 2.0_ark*(-1)**irep
         sym%characters(irep+1,6+2*N_Cn) =-2.0_ark*(-1)**irep
         !
         irep = irep + 1
         !
       enddo
       !
       sym%degen(:)   = 2
       sym%degen(1:8) = 1
       !
       sym%Nelements(1) = 1
       sym%Nelements(1+ 1:1+ N_Cn) = 2
       sym%Nelements(2+N_Cn) = 1
       sym%Nelements(3+N_Cn) = NC2
       sym%Nelements(4+N_Cn) = NC2
       sym%Nelements(5+N_Cn) = 1
       sym%Nelements(5+N_Cn+1:5+2*N_Cn) = 2
       sym%Nelements(6+2*N_Cn) = 1
       sym%Nelements(7+2*N_Cn) = NC2
       sym%Nelements(8+2*N_Cn) = NC2
       !
       o  = 0.0_ark
       p2 = 0.5_ark*pi
       p3 = 1.5_ark*pi
       !
       sym%euler(:,:) = 0
       !
       !if (job%rotsym_do) then 
       !  write(out,"('Rotsym=euler is not implemented for this symmetry yet')")
       !  stop 'Rotsym=euler is not implemented for this symmetry yet'
       !endif
       !
       ioper = 1
       do irot = 1,N_Cn
         !
         ! not tested yet! 
         !
         sym%euler(1+ioper  ,:) = (/o, phi*irot,o/) ! Rz
         sym%euler(1+ioper+1,:) = (/o,-phi*irot,o/) ! Rz
         sym%euler(5+N_Cn+ioper,:)   = (/o,pi+phi*irot,o/) ! Rz
         sym%euler(5+N_Cn+ioper+1,:) = (/o,pi-phi*irot,o/) ! Rz
         !
         ioper = ioper + 2
         !
       enddo
       !
       call irr_allocation
       !
       ! Generate irr-representations
       !
       irep = 8
       do k = 1,(sym%Nrepresen-8)/2
         !
         irep = irep + 1
         !
         do ioper = 1,sym%Noper
           !
           if (ioper==1) then ! E 
             !
             sym%irr(irep,ioper)%repres(1,1) = 1.0_ark
             sym%irr(irep,ioper)%repres(1,2) = 0.0_ark
             !
             sym%irr(irep,ioper)%repres(2,1) = 0.0_ark
             sym%irr(irep,ioper)%repres(2,2) = 1.0_ark
             !
             sym%irr(irep+1,ioper)%repres(1,1) = 1.0_ark
             sym%irr(irep+1,ioper)%repres(1,2) = 0.0_ark
             !
             sym%irr(irep+1,ioper)%repres(2,1) = 0.0_ark
             sym%irr(irep+1,ioper)%repres(2,2) = 1.0_ark
             !
           elseif (ioper<=1+2*N_Cn) then !  Cn x 2 x(n/2-1)
             !
             ioper_ =ioper-1 
             irot = (ioper_+1)/2
             !
             phi_n = phi*irot*k
             !
             ! Second oper in a class is with negative phi
             if ( mod(ioper_,2)==0 ) phi_n = -phi_n
             !
             sym%irr(irep,ioper)%repres(1,1) = cos(phi_n)
             sym%irr(irep,ioper)%repres(1,2) =-sin(phi_n)
             !
             sym%irr(irep,ioper)%repres(2,1) = sin(phi_n)
             sym%irr(irep,ioper)%repres(2,2) = cos(phi_n)
             !
             sym%irr(irep+1,ioper)%repres(1,1) = cos(phi_n)
             sym%irr(irep+1,ioper)%repres(1,2) = sin(phi_n)
             !
             sym%irr(irep+1,ioper)%repres(2,1) =-sin(phi_n)
             sym%irr(irep+1,ioper)%repres(2,2) = cos(phi_n)
             !
           elseif (ioper<=1+2*N_Cn+1) then !  C2 only once
             !
             ioper_ =ioper-1 
             irot = (ioper_+1)/2
             !
             phi_n = phi*irot*k
             !
             ! Second oper in a class is with negative phi
             if ( mod(ioper_,2)==0 ) phi_n = -phi_n
             !
             sym%irr(irep,ioper)%repres(1,1) = cos(phi_n)
             sym%irr(irep,ioper)%repres(1,2) =-sin(phi_n)
             !
             sym%irr(irep,ioper)%repres(2,1) = sin(phi_n)
             sym%irr(irep,ioper)%repres(2,2) = cos(phi_n)
             !
             sym%irr(irep+1,ioper)%repres(1,1) = cos(phi_n)
             sym%irr(irep+1,ioper)%repres(1,2) = sin(phi_n)
             !
             sym%irr(irep+1,ioper)%repres(2,1) =-sin(phi_n)
             sym%irr(irep+1,ioper)%repres(2,2) = cos(phi_n)
             !
           elseif (ioper<=2+2*N_Cn+NC2) then !  C2'
             !
             irot =ioper-(2+2*N_Cn)
             !
             phi_n = phi*irot*2.0_ark
             !
             sym%irr(irep,ioper)%repres(1,1) =-cos(phi_n)
             sym%irr(irep,ioper)%repres(1,2) =-sin(phi_n)
             !
             sym%irr(irep,ioper)%repres(2,1) =-sin(phi_n)
             sym%irr(irep,ioper)%repres(2,2) = cos(phi_n)
             !
             sym%irr(irep+1,ioper)%repres(1,1) = cos(phi_n)
             sym%irr(irep+1,ioper)%repres(1,2) = sin(phi_n)
             !
             sym%irr(irep+1,ioper)%repres(2,1) = sin(phi_n)
             sym%irr(irep+1,ioper)%repres(2,2) =-cos(phi_n)
             !
           elseif (ioper<=2+2*N_Cn+2*NC2) then !  C2"
             !
             irot =ioper-(2+2*N_Cn+NC2)
             !
             phi_n =(-phi*0.5+phi*irot)*2.0_ark
             !
             sym%irr(irep,ioper)%repres(1,1) = cos(phi_n)
             sym%irr(irep,ioper)%repres(1,2) = sin(phi_n)
             !
             sym%irr(irep,ioper)%repres(2,1) = sin(phi_n)
             sym%irr(irep,ioper)%repres(2,2) =-cos(phi_n)
             !
             sym%irr(irep+1,ioper)%repres(1,1) = cos(phi_n)
             sym%irr(irep+1,ioper)%repres(1,2) = sin(phi_n)
             !
             sym%irr(irep+1,ioper)%repres(2,1) = sin(phi_n)
             sym%irr(irep+1,ioper)%repres(2,2) =-cos(phi_n)
             !
           elseif (ioper==3+2*N_Cn+2*NC2) then ! i
             !
             sym%irr(irep,ioper)%repres(1,1) = 1.0_ark
             sym%irr(irep,ioper)%repres(1,2) = 0.0_ark
             !
             sym%irr(irep,ioper)%repres(2,1) = 0.0_ark
             sym%irr(irep,ioper)%repres(2,2) = 1.0_ark
             !
             sym%irr(irep+1,ioper)%repres(1,1) =-1.0_ark
             sym%irr(irep+1,ioper)%repres(1,2) = 0.0_ark
             !
             sym%irr(irep+1,ioper)%repres(2,1) = 0.0_ark
             sym%irr(irep+1,ioper)%repres(2,2) =-1.0_ark
             !
           elseif (ioper<=3+4*N_Cn+2*NC2) then !  2xnxSn
             !
             ioper_ =ioper-(3+2*N_Cn+2*NC2)
             irot = (ioper_+1)/2
             !
             phi_n = phi*irot*k
             !
             ! Second oper in a class is with negative phi
             if ( mod(ioper_,2)==0 ) phi_n = -phi_n
             !
             sym%irr(irep,ioper)%repres(1,1) =-cos(phi_n)
             sym%irr(irep,ioper)%repres(1,2) =-sin(phi_n)
             !
             sym%irr(irep,ioper)%repres(2,1) = sin(phi_n)
             sym%irr(irep,ioper)%repres(2,2) =-cos(phi_n)
             !
             sym%irr(irep+1,ioper)%repres(1,1) = cos(phi_n)
             sym%irr(irep+1,ioper)%repres(1,2) =-sin(phi_n)
             !
             sym%irr(irep+1,ioper)%repres(2,1) = sin(phi_n)
             sym%irr(irep+1,ioper)%repres(2,2) = cos(phi_n)
             !
           elseif (ioper<=4+4*N_Cn+2*NC2) then !  sigma_h
             !
             sym%irr(irep,ioper)%repres(1,1) =-1.0_ark
             sym%irr(irep,ioper)%repres(1,2) = 0.0_ark
             !
             sym%irr(irep,ioper)%repres(2,1) = 0.0_ark
             sym%irr(irep,ioper)%repres(2,2) =-1.0_ark
             !
             sym%irr(irep+1,ioper)%repres(1,1) = 1.0_ark
             sym%irr(irep+1,ioper)%repres(1,2) = 0.0_ark
             !
             sym%irr(irep+1,ioper)%repres(2,1) = 0.0_ark
             sym%irr(irep+1,ioper)%repres(2,2) = 1.0_ark
             !
           elseif (ioper<=4+4*N_Cn+3*NC2) then !  sigmav
             !
             irot = ioper-(4+4*N_Cn+2*NC2)
             !
             phi_n = phi*irot*2.0_ark
             !
             sym%irr(irep,ioper)%repres(1,1) = cos(phi_n)
             sym%irr(irep,ioper)%repres(1,2) = sin(phi_n)
             !
             sym%irr(irep,ioper)%repres(2,1) = sin(phi_n)
             sym%irr(irep,ioper)%repres(2,2) =-cos(phi_n)
             !
             sym%irr(irep+1,ioper)%repres(1,1) = cos(phi_n)
             sym%irr(irep+1,ioper)%repres(1,2) = sin(phi_n)
             !
             sym%irr(irep+1,ioper)%repres(2,1) = sin(phi_n)
             sym%irr(irep+1,ioper)%repres(2,2) =-cos(phi_n)
             !
           elseif (ioper<=4+4*N_Cn+4*NC2) then !  sigmad
             !
             irot = ioper-(4+4*N_Cn+3*NC2)
             !
             phi_n = (-phi*0.5+phi*irot)*2.0_ark
             !
             sym%irr(irep,ioper)%repres(1,1) =-cos(phi_n)
             sym%irr(irep,ioper)%repres(1,2) =-sin(phi_n)
             !
             sym%irr(irep,ioper)%repres(2,1) =-sin(phi_n)
             sym%irr(irep,ioper)%repres(2,2) = cos(phi_n)
             !
             sym%irr(irep+1,ioper)%repres(1,1) = cos(phi_n)
             sym%irr(irep+1,ioper)%repres(1,2) = sin(phi_n)
             !
             sym%irr(irep+1,ioper)%repres(2,1) = sin(phi_n)
             sym%irr(irep+1,ioper)%repres(2,2) =-cos(phi_n)
             !
           else
             !
             stop  'symmetry: illegal ioper'
             !
           endif
           !
         enddo
         !
         irep = irep + 1
         !
       enddo


       do irep = 9,sym%Nrepresen
         !
         mat_t  = matmul(sym%irr(irep,2)%repres,sym%irr(irep,2)%repres) ! C4xC4 = C2
         mat_t  = matmul(sym%irr(irep,2)%repres,mat_t)                  ! C4xC2 = C4-2
         mat_t  = matmul(sym%irr(irep,2)%repres,sym%irr(irep,5)%repres) ! C2" = C4(2)xC2'(1)
         mat_t  = matmul(sym%irr(irep,3)%repres,sym%irr(irep,5)%repres) !
         mat_t  = matmul(sym%irr(irep,4)%repres,sym%irr(irep,12)%repres)-sym%irr(irep, 9)%repres ! i
         mat_t  = matmul(sym%irr(irep,2)%repres,sym%irr(irep,12)%repres)-sym%irr(irep,10)%repres ! Sn(2)
         mat_t  = matmul(sym%irr(irep,3)%repres,sym%irr(irep,12)%repres)-sym%irr(irep,11)%repres ! Sn(2)

         mat_t  = matmul(sym%irr(irep,3)%repres,sym%irr(irep,12)%repres)-sym%irr(irep,10)%repres ! Sn(2)
         mat_t  = matmul(sym%irr(irep,2)%repres,sym%irr(irep,12)%repres)-sym%irr(irep,11)%repres ! Sn(2)


         mat_t  = matmul(sym%irr(irep,5)%repres,sym%irr(irep,12)%repres)-sym%irr(irep,13)%repres ! sigmav1
         mat_t  = matmul(sym%irr(irep,6)%repres,sym%irr(irep,12)%repres)-sym%irr(irep,14)%repres ! sigmav2 
         mat_t  = matmul(sym%irr(irep,7)%repres,sym%irr(irep,12)%repres)-sym%irr(irep,15)%repres ! sigmad1
         mat_t  = matmul(sym%irr(irep,8)%repres,sym%irr(irep,12)%repres)-sym%irr(irep,16)%repres ! sigmad2
         !
         continue
         !
       enddo



    endif
    !
  case("DINFTYH(M)") ! D_infinity_H(M)

    ! Number of rotations to test for < infinity 
    !
    Nrot = 37 ! must be >=1
    !
    sym%Noper=6+2*Nrot
    sym%Nclasses=6
    sym%Nrepresen= max_irreps  ! must be even and >=4
    sym%CII%Noper = 0
    !
    phi = 2.0_ark*pi/real(Nrot,ark)
    !
    call simple_arrays_allocation
    !
    ! E nrotxCinf nrotxsigmav i  nrotxSinf nrotxC'2
    !
    sym%characters(:,:) = 0
    !
    ! A1g,A1u,A2g,A2u:
    ! E
    sym%characters(1:4,1) = 1
    ! Cinf
    sym%characters(1:4,2) = 1
    ! sigmav
    sym%characters(1,3) = 1
    sym%characters(2,3) = 1
    sym%characters(3,3) =-1
    sym%characters(4,3) =-1
    ! i
    sym%characters(1,4) = 1
    sym%characters(2,4) =-1
    sym%characters(3,4) = 1
    sym%characters(4,4) =-1
    ! Sinf
    sym%characters(1,5) = 1
    sym%characters(2,5) =-1
    sym%characters(3,5) = 1
    sym%characters(4,5) =-1
    ! C'inf
    sym%characters(1,6) = 1
    sym%characters(2,6) =-1
    sym%characters(3,6) =-1
    sym%characters(4,6) = 1
    !
    sym%label(1:4)=(/'A1g','A2g','A1u','A2u'/)
    !
    ! E_1g E_1u E_2g E_2u E_3g E_3u .... 
    !
    irep = 4
    do k = 1,(sym%Nrepresen-4)/2
      !
      irep = irep + 1
      !
      write(Kchar, '(i4)') K
      !
      sym%label(irep  ) = 'E'//trim(adjustl(Kchar))//'g'
      sym%label(irep+1) = 'E'//trim(adjustl(Kchar))//'u'
      !
      sym%characters(irep  ,1) = 2.0_ark
      sym%characters(irep+1,1) = 2.0_ark
      !
      sym%characters(irep  ,4) = 2.0_ark
      sym%characters(irep+1,4) =-2.0_ark
      !
      !do irot = 1,Nrot
      !  sym%characters(irep  ,2)          = 2.0_ark*cos(phi*irot*k)
      !  sym%characters(irep+1,2)          = 2.0_ark*cos(phi*irot*k)
      !  sym%characters(irep  ,5) = 2.0_ark*cos(phi*irot*k)*(-1)**k
      !  sym%characters(irep+1,5) = 2.0_ark*cos(phi*irot*k)*(-1)**(k+1)
      !enddo
      !
      sym%characters(irep  ,2)          = 2.0_ark*cos(phi*k)
      sym%characters(irep+1,2)          = 2.0_ark*cos(phi*k)
      sym%characters(irep  ,5) = 2.0_ark*cos(phi*k)*(-1)**k
      sym%characters(irep+1,5) = 2.0_ark*cos(phi*k)*(-1)**(k+1)
      !
      irep = irep + 1
      !
    enddo
    !
    sym%degen(:)   = 2
    sym%degen(1:4) = 1
    !
    sym%Nelements = (/1,2,Nrot,1,2,Nrot/)
    !
    o  = 0.0_ark
    p2 = 0.5_ark*pi
    p3 = 1.5_ark*pi
    !
    sym%euler(:,:) = 0
    !
    !do irot = 1,Nrot
      sym%euler(2,:) = (/o, phi,o/) ! Rz
      sym%euler(3,:) = (/o,-phi,o/) ! Rz
      sym%euler(3+Nrot+2,:) = (/o,pi+phi,o/) ! Rz
      sym%euler(3+Nrot+2,:) = (/o,pi-phi,o/) ! Rz
    !enddo
    !
    call irr_allocation        !
  case default

    write(out,"('symmetry: undefined symmetry group ',a)") trim(sym_group)
    stop 'symmetry: undefined symmetry group '

  end select
  !
  if (max_irreps<sym%Nrepresen) then 
    !
    write(out,"('symmetry: number of elements in _select_gamma_ is too small: ',i8)") 100 ! size(job%select_gamma)
    stop 'symmetry: size of _select_gamma_ is too small'
    !
  endif 
  !
  sym%maxdegen = maxval(sym%degen(:),dim=1)

  !
  ! store the address of the group generator from ioper = 1..Noper list 
  !
  ioper = 1
  !
  do iclass = 1,sym%Nclasses
    !
    sym%igenerator(iclass) = ioper
    ioper = ioper + sym%Nelements(iclass)
    !
  enddo
  
  !
  ! check 
  call check_characters_and_representation
  !

  contains

   subroutine simple_arrays_allocation

    integer(ik) :: alloc,nCII
    integer,pointer	:: Nelements(:)
    nCII = max(1,sym%CII%Noper)
    allocate(Nelements(sym%Nrepresen))
    allocate(sym%Nelements(sym%Nrepresen)  )
	!allocate(sym%Nelements(sym%Nrepresen))
    sym%Nelements = 0
    !
    allocate (sym%characters(sym%Nrepresen,sym%Nclasses))
    allocate (sym%irr(sym%Nrepresen,sym%Noper))
    allocate (sym%degen(sym%Nrepresen),sym%label(sym%Nrepresen))

    allocate (sym%igenerator(sym%Nclasses))

    allocate (sym%CII%ioper(nCII))
    allocate (sym%CII%coeff(nCII))
    allocate (sym%euler(sym%Noper,3))
    allocate (sym%lquant(sym%Nrepresen))
	alloc = 0

    if (alloc/=0) stop 'simple_arrays_allocation - out of memory'
    !
    sym%CII%coeff = 0
    sym%CII%ioper = 1
    sym%euler = 0
    sym%lquant = 0 
    !
   end subroutine simple_arrays_allocation



   subroutine irr_allocation

    integer(ik) :: gamma,ioper,iclass,ielem,alloc

    do gamma = 1,sym%Nrepresen
      ioper = 0
      do iclass = 1,sym%Nclasses
        do ielem = 1,sym%Nelements(iclass)
          !
          ioper = ioper + 1
          !
          allocate(sym%irr(gamma,ioper)%repres(sym%degen(gamma),sym%degen(gamma)),stat=alloc)
          !
          if (sym%degen(gamma)==1) then
             sym%irr(gamma,ioper)%repres(1,1)= sym%characters(gamma,iclass)
          endif 
          !
        enddo 
      enddo 
    enddo 


   if (alloc/=0) then
       write (out,"(' symmetryInitialize ',i9,' error trying to allocate symmetry')") alloc
      stop 'symmetryInitialize, symmetries - out of memory'
   end if


   end subroutine irr_allocation


 end subroutine symmetryInitialize



   subroutine check_characters_and_representation

    integer(ik) :: igamma,jgamma,ioper,ideg,iclass,ielem
    real(ark)   :: temp

    do igamma = 1,sym%Nrepresen
      do jgamma = igamma,sym%Nrepresen
        !
        temp = sum(real(sym%Nelements(:),ark)*sym%characters(igamma,:)*sym%characters(jgamma,:))
        !
        if (igamma/=jgamma.and.abs(temp)>sqrt(small_)) then 
          write (out,"(' check_characters_and_representation: not orhogonal for igamma,jgamma = ',2i4,' -> ',g18.6)") igamma,jgamma,temp
          stop 'check_characters_and_representation: not orhogonal'
        endif
        !
        if (igamma==jgamma.and.abs(temp-sym%Noper)>sqrt(small_)) then 
          write (out,"(' check_characters_and_representation: dot procut ',f16.2,' for isym = ',i4,' /= size of the group ',f16.0)") igamma,temp,sym%Noper
          stop 'check_characters_and_representation: not orhogonal'
        endif
        !
      enddo 
    enddo 

    !
    ! Check characters
    !
    if (verbose_>=5) write(out,"('Irrep matrices:')")
    !
    ioper = 0
    do iclass = 1,sym%Nclasses
      !
      do ielem = 1,sym%Nelements(iclass)
        !
        ioper = ioper + 1
        !
        do igamma = 1,sym%Nrepresen
          !
          temp = 0
          !
          do ideg = 1,sym%degen(igamma)
            !
            temp = temp + sym%irr(igamma,ioper)%repres(ideg,ideg)
            !
          enddo
          !
          if (verbose_>=5) then
            write(out,"('igamma,iclass,ioper = ',3i6)") igamma,iclass,ioper
            do ideg = 1,sym%degen(igamma)
          !    write(out,"(<sym%degen(igamma)>f18.8)") sym%irr(igamma,ioper)%repres(ideg,:)
            enddo
          endif
          !
          if (abs(temp-sym%characters(igamma,iclass))>sqrt(small_)) then 
            write (out,"(' symmetry: character and representation do not agree ',2f18.8,', igamma,iclass,ioper = ',3i5)") sym%characters(igamma,iclass),temp,igamma,iclass,ioper
            !stop 'symmetry: character and representation do not agree'
          endif
          !
        enddo
        !
      enddo
      !
    enddo


   end subroutine check_characters_and_representation

   function dlmn(j,n,m,t,error) result (d)
     integer(ik),intent(in) :: j,n,m
     integer(ik),intent(out) :: error
     real(ark),intent(in)   :: t
     real(ark)              :: ss,prefactor,d,f,g,g1,g2,x,sinx,cosx,gh,h,hs,cos_sin,hlog
     integer(ik)            :: s
     
         error = 0
         s  =0
         ss =0
         !
         x = 0.5_ark*t ; cosx = cos(x) ; sinx = sin(x)
         !
         g2 = 0.5_ark*(faclog(j-m)+faclog(j+m)+faclog(j-n)+faclog(j+n))
         !
         do s = 0,min(j-m,j-n)
            !
            if (j-m-s>=0.and.j-n-s>=0.and.m+n+s>=0) then
              !
              !f=1/((j-m-s)!*(j-n-s)!*s!*(m+n+s)!);
              !
              g1 = faclog(j-m-s)+faclog(j-n-s)+faclog(s)+faclog(m+n+s)
              !
              g = g2-g1
              !
              cos_sin = cosx**(2*s+m+n)*sinx**(2*j-2*s-m-n)*(-1.0_ark)**s
              hs = sign(1.0_ark,cos_sin)
              h = abs(cos_sin)
              !
              if (h>0) then
                !
                hlog = log(h)
                !
                gh = g + hlog
                !
                if (gh>max_exp) then 
                  !write(out,"('dlmn error:  the value is too large for exp:',g18.9)")  g
                  error = 1
                  d = 0
                  return 
                  !stop 'dlmn error:  the value is too large for exp'
                endif
                !
                f = exp(gh)*hs
                !
                ss = ss + f
                !
              endif
              !
              !f = exp(-g)
              !
            endif
            !
         enddo
         !
         !g = faclog(j-m)+faclog(j+m)+faclog(j-n)+faclog(j+n)
         !
         !prefactor=exp(g)
         !
         d=(-1.0_ark)**(j-n)*ss  !*sqrt(prefactor)
         !
   end function dlmn



   !
   ! D function conjugated
   !
   function ddlmn_conj(j,n,m,phi,theta,chi,info) result (d)
     integer(ik),intent(in) :: j,m,n
     integer(ik),optional,intent(out) :: info
     real(ark),intent(in)   :: phi,chi,theta
     complex(ark)           :: f,d
     real(ark)              :: prefactor,x,y
     integer(ik)            :: error
     
         !
         prefactor=(-1.0_ark)**(n-m)*dlmn(j,-n,-m,theta,error)
         !
         x = real(n,ark)*phi ; y = real(m,ark)*chi
         !
         f=cmplx(cos(x+y),sin(x+y))
         !
         d = f*prefactor
         !
         if (present(info)) info = error
         if (.not.present(info).and.error/=0) then
            info = error
            write(out,"('dlmn error:  the value is too large for exp ')")
            stop 'dlmn error:  the value is too large for exp'
         endif
         !
   end function ddlmn_conj



end module symmetry



