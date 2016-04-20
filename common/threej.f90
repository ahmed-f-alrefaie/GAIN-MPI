module FORT_FUNC


	implicit none
		public three_j
  integer, parameter :: silk         = selected_int_kind(4)       ! Small integers
  integer, parameter :: ilk          = selected_int_kind(8)       ! "Normal" integers. This must map on
                                                                 ! C "int" type, or the DX interface won't
                                                                 ! work correctly.
  integer, parameter :: hilk         = selected_int_kind(16)      ! "Pointer" integers - sufficient to store
                                                                 ! memory address
  integer, parameter :: drlk         = selected_real_kind(12,25)  ! "Double" reals and complex (complexi? :-)
  integer, parameter :: rlk          = selected_real_kind(12,25)  ! "Normal" reals and complex (complexi? :-)
  integer, parameter :: arlk         = selected_real_kind(25,32)  ! "Accurate" reals and complex (complexi? :-)
  contains
   function three_j(j1,j2,j3,k1,k2,k3) bind (c, name='c_three_j')
  use,intrinsic :: iso_c_binding
  implicit none

      real(rlk) :: three_j
      integer(ilk),intent(in) :: j1,j2,j3,k1,k2,k3
      real(rlk) :: a,b,c,al,be,ga
      !
      integer(ilk):: newmin,newmax,new,iphase
      real(rlk)   :: delta,clebsh,minus
      real(rlk)   :: term,term1,term2,term3,summ,dnew,term4,term5,term6,delta_log,term16,termlog

      a = j1
      b = j2
      c = j3
      al= k1
      be= k2
      ga= k3

      three_j=0
!
!     (j1+j2).ge.j and j.ge.abs(a-b)    -m=m1+m2    j1,j2,j.ge.0
!     abs(m1).le.j1    abs(m2).le.j2   abs(m).le.j
!
      if(c.gt.a+b) return
      if(c.lt.abs(a-b)) return
      if(a.lt.0.or.b.lt.0.or.c.lt.0) return
      if(a.lt.abs(al).or.b.lt.abs(be).or.c.lt.abs(ga)) return
      if(-1.0_rlk*ga.ne.al+be) return
!
!
!     compute delta(abc)
!
!     delta=sqrt(fakt(a+b-c)*fakt(a+c-b)*fakt(b+c-a)/fakt(a+b+c+1.0_rlk))
      delta_log = faclogf(a+b-c)+faclogf(a+c-b)+faclogf(b+c-a)-faclogf(a+b+c+1.0_rlk)
      !
      delta=sqrt(exp(delta_log)) 
!
!
      !term1=fakt(a+al)*fakt(a-al)
      !term2=fakt(b-be)*fakt(b+be)
      !term3=fakt(c+ga)*fakt(c-ga)
      !
      !term=sqrt( (2.0_rlk*c+1.0_rlk)*term1*term2*term3 )
      !
      !
      term1=faclogf(a+al)+faclogf(a-al)
      term2=faclogf(b-be)+faclogf(b+be)
      term3=faclogf(c+ga)+faclogf(c-ga)
      !
      termlog = ( term1+term2+term3+delta_log )*0.5_rlk
 
      term=sqrt( (2.0_rlk*c+1.0_rlk) )
!
!
!     now compute summation term
!
!     sum to get summation in eq(2.34) of brink and satchler.  sum until
!     a term inside factorial goes negative.  new is index for summation
!     .  now find what the range of new is.
!
!
      newmin=idnint(max((a+be-c),(b-c-al),0.0_rlk))
      newmax=idnint(min((a-al),(b+be),(a+b-c)))
!
!
      summ=0
!
!
      do new=newmin,newmax
        !
        dnew=real(new,rlk)
        !
        term4=faclogf(a-al-dnew)+faclogf(c-b+al+dnew)
        term5=faclogf(b+be-dnew)+faclogf(c-a-be+dnew)
        term6=faclogf(dnew)+faclogf(a+b-c-dnew)
        !
        term16=termlog-(term4+term5+term6)
        !
        summ=summ+(-1.0_rlk)**new*exp(term16)
        !
      enddo
!
!     so clebsch-gordon <j1j2m1m2ljm> is clebsh
!
      clebsh=term*summ ! /sqrt(10.0_rlk)
!
!     convert clebsch-gordon to three_j
!
      iphase=idnint(a-b-ga)
      minus = -1.0_rlk
      if (mod(iphase,2).eq.0) minus = 1.0_rlk
      three_j=minus*clebsh/term

!     threej=(-1.d0)**(iphase)*clebsh/sqrt(2.0_rlk*c+1.d0)
!
!
   end function three_j



   function fakt(a) result (f)

      real(rlk),intent(in) :: a
      real(rlk)            :: ax,f
      integer(ilk)         :: i,ic
!

!
      ax=a
      f=1.0_rlk
      if(abs(ax)<1.d-24) return
      f=.1d0
      if(ax.lt.0.d0) then 
         write (*,"(1h0,' fkt.err  negative argument for functi on fakt. argument = ',e12.5)") ax
         stop 'fkt.err  negative argument'
      endif 
      !
      ic=idnint(ax)
      ax=ax/10.0_rlk
      f=ax
      do  i=1,ic-1
        f=f*(ax-real(i,rlk)*0.1_rlk)
      enddo

    end function fakt

  function faclogf(a)   result (v)
    real(rlk),intent(in) ::  a
    real(arlk)              :: v 
    integer(ilk) j,k

    v=0
    k=nint(a)
    if(k>=2) then 
      do j=2,k
         v=v+log(real(j,arlk))
      enddo 
    endif 
    
  end function faclogf
end module FORT_FUNC
