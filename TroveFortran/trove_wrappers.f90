!
module trove_wrappers
  use fort_func
  use accuracy
  use symmetry
  implicit none
  public SymmetryInitialize,sym,max_irreps,Symm_initialize,Get_Nirreps



   type MOrepres_arkT
      real(ark),pointer     :: repres(:,:,:)  
   end type MOrepres_arkT


contains 

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
	call accuracyInitialize
	call SymmetryInitialize(sym_group)
  end subroutine Symm_initialize

  subroutine Symm_GetNrepres(Nrepres) bind (c, name='c_SymGetNrepres')
    use,intrinsic :: iso_c_binding
    implicit none
	integer(ik),intent(out)	::	Nrepres	
	Nrepres = sym%Nrepresen


  end subroutine Symm_GetNrepres

  subroutine Symm_GetSymDegen(degen,symn) bind (c, name='c_GetSymDegen')
    use,intrinsic :: iso_c_binding
    implicit none
	integer(ik),intent(in) ::      symn
	integer(ik),intent(out)	::	degen
	degen = sym%degen(symn+1)


  end subroutine Symm_GetSymDegen


  subroutine Get_Nirreps(J,Ntot,eigenvects) bind (c, name='c_GetNirreps')
  use,intrinsic :: iso_c_binding
  implicit none
	integer(4),intent(in)   :: J
	integer(4),intent(out)  :: Ntot(sym%Nrepresen)
	real(8),intent(out)	:: eigenvects(2*J+1,2*J+1)
	integer(ik)		::	Ncount,count_index(2*j+1,2*j+1)

	call MLrotsymmetry_generate(J,6,count_index,eigenvects,Ncount,Ntot)

  end subroutine Get_Nirreps


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




end module trove_wrappers



