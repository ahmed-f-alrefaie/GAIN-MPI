module FORT_FUNC
        use accuracy

	implicit none
		public three_j
  contains
   function three_j(j1,j2,j3,k1,k2,k3) bind (c, name='c_three_j')
  use,intrinsic :: iso_c_binding
  implicit none

      real(rk) :: three_j
      integer(ik),intent(in) :: j1,j2,j3,k1,k2,k3
      real(rk) :: a,b,c,al,be,ga
      !
      integer(ik):: newmin,newmax,new,iphase
      real(rk)   :: delta,clebsh,minus
      real(rk)   :: term,term1,term2,term3,summ,dnew,term4,term5,term6,delta_log,term16,termlog

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
      if(-1.0_rk*ga.ne.al+be) return
!
!
!     compute delta(abc)
!
!     delta=sqrt(fakt(a+b-c)*fakt(a+c-b)*fakt(b+c-a)/fakt(a+b+c+1.0_rk))
      delta_log = faclogf(a+b-c)+faclogf(a+c-b)+faclogf(b+c-a)-faclogf(a+b+c+1.0_rk)
      !
      delta=sqrt(exp(delta_log)) 
!
!
      !term1=fakt(a+al)*fakt(a-al)
      !term2=fakt(b-be)*fakt(b+be)
      !term3=fakt(c+ga)*fakt(c-ga)
      !
      !term=sqrt( (2.0_rk*c+1.0_rk)*term1*term2*term3 )
      !
      !
      term1=faclogf(a+al)+faclogf(a-al)
      term2=faclogf(b-be)+faclogf(b+be)
      term3=faclogf(c+ga)+faclogf(c-ga)
      !
      termlog = ( term1+term2+term3+delta_log )*0.5_rk
 
      term=sqrt( (2.0_rk*c+1.0_rk) )
!
!
!     now compute summation term
!
!     sum to get summation in eq(2.34) of brink and satchler.  sum until
!     a term inside factorial goes negative.  new is index for summation
!     .  now find what the range of new is.
!
!
      newmin=idnint(max((a+be-c),(b-c-al),0.0_rk))
      newmax=idnint(min((a-al),(b+be),(a+b-c)))
!
!
      summ=0
!
!
      do new=newmin,newmax
        !
        dnew=real(new,rk)
        !
        term4=faclogf(a-al-dnew)+faclogf(c-b+al+dnew)
        term5=faclogf(b+be-dnew)+faclogf(c-a-be+dnew)
        term6=faclogf(dnew)+faclogf(a+b-c-dnew)
        !
        term16=termlog-(term4+term5+term6)
        !
        summ=summ+(-1.0_rk)**new*exp(term16)
        !
      enddo
!
!     so clebsch-gordon <j1j2m1m2ljm> is clebsh
!
      clebsh=term*summ ! /sqrt(10.0_rk)
!
!     convert clebsch-gordon to three_j
!
      iphase=idnint(a-b-ga)
      minus = -1.0_rk
      if (mod(iphase,2).eq.0) minus = 1.0_rk
      three_j=minus*clebsh/term

!     threej=(-1.d0)**(iphase)*clebsh/sqrt(2.0_rk*c+1.d0)
!
!
   end function three_j



   function fakt(a) result (f)

      real(rk),intent(in) :: a
      real(rk)            :: ax,f
      integer(ik)         :: i,ic
!

!
      ax=a
      f=1.0_rk
      if(abs(ax)<1.d-24) return
      f=.1d0
      if(ax.lt.0.d0) then 
         write (*,"(1h0,' fkt.err  negative argument for functi on fakt. argument = ',e12.5)") ax
         stop 'fkt.err  negative argument'
      endif 
      !
      ic=idnint(ax)
      ax=ax/10.0_rk
      f=ax
      do  i=1,ic-1
        f=f*(ax-real(i,rk)*0.1_rk)
      enddo

    end function fakt

  function faclogf(a)   result (v)
    real(rk),intent(in) ::  a
    real(ark)              :: v 
    integer(ik) j,k

    v=0
    k=nint(a)
    if(k>=2) then 
      do j=2,k
         v=v+log(real(j,ark))
      enddo 
    endif 
    
  end function faclogf

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

end module FORT_FUNC
