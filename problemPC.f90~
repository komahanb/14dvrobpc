      program problemKriging

        use dimKrig,only:probtype,id_proc
!
      implicit none
!
!     include the Ipopt return codes
!
      include 'IpReturnCodes.inc'
      include 'mpif.h'
!
!     Size of the problem (number of variables and equality constraints)
!
      integer     N,     M,     NELE_JAC,     NELE_HESS,      IDX_STY
      parameter  (N = 16, M = 1, NELE_JAC = 16, NELE_HESS = 50)
      parameter  (IDX_STY = 1 )
!
!     Space for multipliers and constraints
!
      double precision LAM(M)
      double precision G(M)
!
!     Vector of variables
!
      double precision X(N)
!
!     Vector of lower and upper bounds
!
      double precision X_L(N), X_U(N), Z_L(N), Z_U(N)
      double precision G_L(M), G_U(M)
!
!     Private data for evaluation routines
!     This could be used to pass double precision and integer arrays untouched
!     to the evaluation subroutines EVAL_*
!
      double precision DAT(2000)
      integer IDAT(2000)
!
!     Place for storing the Ipopt Problem Handle
!
      integer*8 IPROBLEM
      integer*8 IPCREATE
!
      integer IERR
      integer IPSOLVE, IPADDSTROPTION
      integer IPADDNUMOPTION, IPADDINTOPTION
      integer IPOPENOUTPUTFILE
!
      double precision F,sigmax(N),Lifttarget,pi
      integer i,kprob

      double precision  infbound
      parameter        (infbound = 1.d+20)
!
!     The following are the Fortran routines for computing the model
!     functions and their derivatives - their code can be found further
!     down in this file.
!
      external EV_F, EV_G, EV_GRAD_F, EV_JAC_G, EV_HESS, ITER_CB

      call MPI_START

  !**********************************************************************

  pi=4.0*atan(1.0) ! constant for later use (visible globally)

  !======================================
  !(1)    Set initial point and bounds:
  !=====================================
  
  !
  ! (1)  Set initial point, bounds and sigmax for geometric variables (epistemic):
  !
  
  do i=1,14
     sigmax(i)=0.0025
     X(i) =0.0
     X_L(i) = -0.01
     X_U(i) = 0.01
  end do

  X_L(1) = -0.00125
  X_U(1) = 0.00125
  sigmax(1)=0.00125

  X_L(2) = -0.00125
  X_U(2) = 0.00125
  sigmax(2)=0.00125

  X_L(13) = -0.00125
  X_U(13) = 0.00125
  sigmax(13)=0.00125

  X_L(14) = -0.00125
  X_U(14) = 0.00125
  sigmax(14)=0.00125

  !X(1:7)=X_L(1:7)-sigmax(1:7)
  !X(8:14)=X_L(8:14)-sigmax(8:14)
  !X_L(1:7)=X(1:7)
  !X_L(8:14)=X(8:14)

  !
  !     Set initial point, bounds and sigmax for flow variables (aleatory):
  !

  sigmax(N-1)=0.1
  X(N-1) = 2.0    !alpha in degrees
  X_L(N-1) = 0.0
  X_U(N-1) = 4.0

  pi=4.0*atan(1.0)
  X(N-1) = X(N-1)*pi/180.0
  sigmax(N-1)=sigmax(N-1)*pi/180.0
  X_L(N-1) = X_L(N-1)*pi/180.0
  X_U(N-1) = X_U(N-1)*pi/180.0


  sigmax(N)=0.01
  X(N) = 0.56  !Minf
  X_L(N) = 0.1
  X_U(N) = 0.78
  
  !===================================================================
  !(2)     Integer Settings and store into IDAT (check for size above)
  !===================================================================
  
  kprob=1  
  probtype(:)=1

  IDAT(1)=kprob
  IDAT(2)=0
  IDAT(3:N+2)=probtype(1:N)

  !===============================================
  !(3)     Setup std dev and store in to dat(1:N)
  !===============================================

  Lifttarget=0.9  !0.268482186143556
  DAT(1)=Lifttarget
  do i=2,N+1
     DAT(i)=sigmax(i-1)
  end do

  !====================
  !(4)     Constraints
  !====================

  do i=1,M
     G_L(i)=-infbound
     G_U(i)=0.0
  end do

  !===========================================================
  !(5)    Other constants to be passed to the surrogate call
  !===========================================================
  
  pi=4.0*atan(1.0) ! constant for later use

!!$
!!$  !Problem data and other constants
!!$  dat(1000+1)=10.0 !height ref
!!$  dat(1000+2)=1.0e7 !E
!!$  dat(1000+3)=0.1 !gamma
!!$  dat(1000+4)=45.0*pi/180.0
!!$  dat(1000+5)=20000.0
!!$  ! Max constraint values
!!$  dat(1000+6)=5000.0    ! psi tensile_sigma1_max=dat(6)      
!!$  dat(1000+7)=20000.0    ! psi tensile_sigma2_max=dat(7)
!!$  dat(1000+8)=5000.0    ! psi tensile_sigma3_max=dat(8)
!!$  !Compressive
!!$  dat(1000+9)=5000.0    ! psi comp_sigma1_max=dat(9)
!!$  dat(1000+10)=20000.0   ! psi comp_sigma2_max=dat(10)
!!$  dat(1000+11)=5000.0   ! psi comp_sigma3_max=dat(11)
!!$  !Displacement
!!$  dat(1000+12)=0.005    ! in  max_u_disp=dat(12)
!!$  dat(1000+13)=0.005    ! in  max_v_disp=dat(12)
!!$  dat(1000+14)=1.0      ! Factor of safety

  dat(1000+20)=6       ! filenum for PC

!===========================================================================
!
!     First create a handle for the Ipopt problem (and read the options
!     file)
!
      IPROBLEM = IPCREATE(N, X_L, X_U, M, G_L, G_U, NELE_JAC, NELE_HESS,IDX_STY, EV_F, EV_G, EV_GRAD_F, EV_JAC_G, EV_HESS)
      if (IPROBLEM.eq.0) then
         write(*,*) 'Error creating an Ipopt Problem handle.'
         call stop_all
      endif
!
!     Open output files
!

      if (id_proc.eq.0) open(unit=76,file='Opt.his',form='formatted',status='replace')

      IERR = IPOPENOUTPUTFILE(IPROBLEM, 'IPOPT.OUT', 5)
      if (IERR.ne.0 ) then
         write(*,*) 'Error opening the Ipopt output file.'
         goto 9000
      endif
!!
!!     Set a callback function to give you control once per iteration.
!!     You can use it if you want to generate some output, or to stop
!!     the optimization early.
!!
      call IPSETCALLBACK(IPROBLEM, ITER_CB)

!
!     Call optimization routine
!
      if (id_proc.eq.0) then
          IERR = IPADDINTOPTION(IPROBLEM, 'print_level', 0)
          if (IERR.ne.0 ) goto 9990
      else
         IERR = IPADDINTOPTION(IPROBLEM, 'print_level', 0)
         if (IERR.ne.0 ) goto 9990
      end if

      IERR = IPSOLVE(IPROBLEM, X, G, F, LAM, Z_L, Z_U, IDAT, DAT)


!===================
! (6)  Output:
!==================


      if (id_proc.eq.0) then

         if( IERR.eq.IP_SOLVE_SUCCEEDED .or. IERR.eq.5) then
            write(*,*)
            write(*,*) 'The solution was found.'
            write(*,*)
         else
            write(*,*)
            write(*,*) 'An error occoured.'
            write(*,*) 'The error code is ',IERR
            write(*,*)
         endif
         
         write(*,*) 'The final value of the objective function is ',F
         write(*,*)
         write(*,*) 'The optimal values of X are:'
         write(*,*)
         X(N-1)=X(N-1)*180.0/pi
         do i = 1, N
            write(*,*) 'X  (',i,') = ',X(i)
         enddo
         write(*,*)
         write(*,*) 'The multipliers for the equality constraints are:'
         write(*,*)
         do i = 1, M
            write(*,*) 'LAM(',i,') = ',LAM(i)
         enddo
         write(*,*)
         write(*,*) 'Mean drag and its variance:',DAT(N+2),DAT(N+3)

      end if
         
!
 9000 continue
!
!     Clean up
!
      call IPFREE(IPROBLEM)
      
      if (id_proc.eq.0) close(76)

      call stop_all
!
 9990 continue
      write(*,*) 'Error setting an option'
      goto 9000

    end program problemKriging
!
! =============================================================================
!
!                    Computation of objective function
!
! =============================================================================
!
      subroutine EV_F(N, X, NEW_X, F, IDAT, DAT, IERR)
        use dimKrig,only:probtype,id_proc
      implicit none
      integer N, NEW_X,I
      double precision F, X(N),sigmax(N),fmeantmp,fvartmp,fmeanprimetmp(n),fvarprimetmp(n)
      double precision DAT(*),dfdD(n),dfdDD(n,n),v(n)
      double precision fmin,fmax,gradmin(N),gradmax(N),gtol,low(N-2),up(N-2),Xsave(N),lifttarget
      integer IDAT(*),kprob,NMC
      integer IERR
      integer::myflag(10) 

      if (id_proc.eq.0) print *,'Calc obj',X

      kprob=IDAT(1)
      probtype(1:N)=IDAT(3:N+2)

      Lifttarget=DAT(1)
      do i=1,N
         sigmax(i)=DAT(i+1)
         Xsave(i)=X(i)
      end do 

!!$      gtol=1e-4
!!$
!!$      low(1:N-2)=X(1:N-2)-sigmax(1:N-2)
!!$      up(1:N-2)=X(1:N-2)+sigmax(1:N-2)
!!$
!!$      call optimize(N-2,X,N,fmax,gradmax,low,up,gtol,.true.,.false.,10)

      !---- MEAN OF worst OBJECTIVE FUNCTION

      NMC=100000
      
      call Krigingestimate(2,N,x,sigmax,22,0,DAT(1001:1020),5,2,9,0,probtype,myflag,fmeantmp,fvartmp,fmeanprimetmp,fvarprimetmp)


      !!      call Krigingestimate(2,X,N,sigmax,fmeantmp,fvartmp,fmeanprimetmp,fvarprimetmp,NMC,0)

      !fvartmp=0.0
      !fvarprimetmp(:)=0.0
      !call Eulersolve(X,N,0,fmeantmp,dfdD,dfdDD,1,v,0)
      !fmeanprimetmp(:)=dfdD(:)

  
!---- COMBINED OBJECTIVE FUNCTION
      F=fmeantmp+fvartmp
      
      DAT(N+2)=fmeantmp
      DAT(N+3)=fvartmp

      if (id_proc.eq.0) then
      
         write(*,*) 'Design variables:',x(1:N)
         write(*,*) 'Mean drag value and variance:',fmeantmp,fvartmp
     
      end if

!---- OBJECTIVE FUNCTION gradient and x value

      do i=1,N-2
         DAT(N+4+i)=0.0
      end do
      DAT(N+4+N-1)=fmeanprimetmp(N-1)+fvarprimetmp(N-1)
      DAT(N+4+N)=fmeanprimetmp(N)+fvarprimetmp(N)

      do i=1,N
         DAT(2*N+4+i)=Xsave(i)
         X(i)=Xsave(i)
      end do

      IERR = 0
      return
      end

!
! =============================================================================
!
!                     Computation of constraints
!
! =============================================================================
!
      subroutine EV_G(N, X, NEW_X, M, G, IDAT, DAT, IERR)
        use dimKrig,only:probtype,id_proc
        implicit none
      integer N, NEW_X, M
      double precision G(M), X(N), sigmax(N), cmean(M), cstd(M), fmeantmp, fvartmp
      double precision DAT(*),Lifttarget,dfdD(n),dfdDtmp(n),dfdDD(n,n),v(n),fmeanprimetmp(n),fvarprimetmp(n)
      double precision fmin,fmax,gradmin(N),gradmax(N),gtol,low(N-2),up(N-2),Xsave(N),dc(M,N)
      integer IDAT(*),kprob,NMC
      integer IERR, i, j, cnt,myflag(10)

      if (id_proc.eq.0) print *,'Calc con',X

      kprob=IDAT(1)
      probtype(1:N)=IDAT(3:N+2)

      Lifttarget=DAT(1)
      do i=1,N
         sigmax(i)=DAT(i+1)
         Xsave(i)=X(i)
      end do 

!!$      gtol=1e-4
!!$
!!$      low(1:N-2)=X(1:N-2)-sigmax(1:N-2)
!!$      up(1:N-2)=X(1:N-2)+sigmax(1:N-2)

!---- Worst MEAN and GRADIENT OF INEQUALITY CONSTRAINT

!!$      call optimize(N-2,X,N,fmin,gradmin,low,up,gtol,.false.,.false.,14)
      
      dc(:,:)=0.0

      NMC=100000

      call Krigingestimate(2,N,x,sigmax,22,4,DAT(1001:1020),5,2,9,0,probtype,myflag,fmeantmp,fvartmp,fmeanprimetmp,fvarprimetmp)

!!      call Krigingestimate(2,X,N,sigmax,fmeantmp,fvartmp,fmeanprimetmp,fvarprimetmp,NMC,4)

      !fvartmp=0.0
      !fvarprimetmp(:)=0.0
      !call Eulersolve(x,n,0,fmeantmp,dfdD,dfdDD,1,v,4)
      !fmeanprimetmp(:)=dfdD(:) 


      cmean(1)=Lifttarget-fmeantmp
      cstd(1)=sqrt(fvartmp)

!---- COMBINED INEQUALITY CONSTRAINTS

      G(1)=cmean(1)+kprob*cstd(1)

     
      if (id_proc.eq.0) then
      
         write(*,*) 'Design variables:',x(1:N)
         write(*,*) 'Constraint value, mean lift and its standard deviation:',G(1),fmeantmp,cstd(1)
   
      end if

!---- INEQUALITY CONSTRAINTS gradient

      do j=N-1,N !can be 1,N ?
         dc(1,j)=-fmeanprimetmp(j)
         if (fvartmp.ne.0.0) then
            dc(1,j)=dc(1,j)+kprob*fvarprimetmp(j)/(2.0*sqrt(fvartmp))
         endif
      end do

      cnt=0
      do i=1,M
         do j=1,N
            cnt=cnt+1
            DAT(4*N+4+cnt)=dc(i,j)
         end do
      end do

      do i=1,N
         DAT(3*N+4+i)=Xsave(i)
         X(i)=Xsave(i)
      end do

      IERR = 0
      return
      end

!
! =============================================================================
!
!                Computation of gradient of objective function
!
! =============================================================================
!
      subroutine EV_GRAD_F(N, X, NEW_X, GRAD, IDAT, DAT, IERR)
        use dimKrig,only:probtype,id_proc
        implicit none
      integer N, NEW_X,i,j,k
      double precision GRAD(N),GRADtmp(N), X(N), sigmax(N), fmeantmp, fvartmp
      double precision DAT(*),objtmp,gradvar(n),dfdDD(n,n),v(n),fmeanprimetmp(n),fvarprimetmp(n)
      double precision fmin,fmax,gradmin(N),gradmax(N),gtol,low(N-2),up(N-2),Xsave(N),lifttarget
      integer IDAT(*),kprob,NMC
      integer IERR,myflag(10)
      logical samex
           
      !if (id_proc.eq.0) print *,'Calc obj grad',X

      samex=.true.
      do i=1,N
         if (x(i).ne.DAT(2*N+4+i)) samex=.false. 
      end do

      kprob=IDAT(1)
      probtype(1:N)=IDAT(3:N+2)

      Lifttarget=DAT(1)
      do i=1,N
         sigmax(i)=DAT(i+1)
         Xsave(i)=X(i)
      end do
      
      if (samex) then

         !if (id_proc.eq.0) print *,'Samex in obj'

         !---- TOTAL GRADIENT OF OBJECTIVE FUNCTION
         do i=1,n
            GRAD(i)=DAT(N+4+i)
         end do
         
      else

         !if (id_proc.eq.0) print *,'Not Samex in obj'

         NMC=100000

         call Krigingestimate(2,N,x,sigmax,22,0,DAT(1001:1020),5,2,9,0,probtype,myflag,fmeantmp,fvartmp,fmeanprimetmp,fvarprimetmp)

!!         call Krigingestimate(2,X,N,sigmax,fmeantmp,fvartmp,fmeanprimetmp,fvarprimetmp,NMC,0)

         !fvartmp=0.0
         !fvarprimetmp(:)=0.0
         !call Eulersolve(x,n,0,fmeantmp,GRADtmp,dfdDD,1,v,0)
         !fmeanprimetmp(:)=GRADtmp(:) 

         GRAD(N-1:N)=fmeanprimetmp(N-1:N)+fvarprimetmp(N-1:N)
   
      end if

      !---- GRADIENT OF worst OBJECTIVE FUNCTION

      call Eulersolve(Xsave,N,0,objtmp,GRADtmp,dfdDD,1,v,0)
      
      GRAD(1:N-2)=GRADtmp(1:N-2)

      do i=1,N
         X(i)=Xsave(i)
      end do

      !if (id_proc.eq.0) print *,'Obj Gradient',GRAD(1:3)

      IERR = 0
      return
      end

!
! =============================================================================
!
!                Computation of Jacobian of constraints
!
! =============================================================================
!
      subroutine EV_JAC_G(TASK, N, X, NEW_X, M, NZ, ACON, AVAR, JAC,IDAT, DAT, IERR)
        use dimKrig,only:probtype,id_proc
        implicit none
      integer TASK, N, NEW_X, M, NZ
      integer ACON(NZ), AVAR(NZ), I, J, K
      double precision X(N), JAC(NZ),sigmax(N), objtmp, fmeantmp, fvartmp  
      double precision DAT(*),dfdD(n),dfdDD(n,n),v(n),Varc(M),dc(M,N),fmeanprimetmp(n),fvarprimetmp(n)
      double precision fmin,fmax,gradmin(N),gradmax(N),gtol,low(N-2),up(N-2),Xsave(N),lifttarget
      integer IDAT(*)
      integer IERR, kprob, NMC, cnt
      logical samex
      integer::myflag(10) 

      if( TASK.eq.0 ) then 

         !if (id_proc.eq.0) print *,'Initialize Jacobian'
         !
         !     structure of Jacobian:
         !
         do i=1,N
            ACON(i)=1
            AVAR(i)=i
         end do

      else

         !if (id_proc.eq.0) print *,'Calc con grad',X

         samex=.true.
         do i=1,N
            if (x(i).ne.DAT(3*N+4+i)) samex=.false. 
         end do

         kprob=IDAT(1)
         probtype(1:N)=IDAT(3:N+2)

         Lifttarget=DAT(1)
         do i=1,N
            sigmax(i)=DAT(i+1)
            Xsave(i)=X(i)
         end do

         if (samex) then

            !if (id_proc.eq.0) print *,'Samex in con'

            cnt=0
            do i=1,M
               do j=1,N
                  cnt=cnt+1
                  dc(i,j)=DAT(4*N+4+cnt)
               end do
            end do
            

         else

            !if (id_proc.eq.0) print *,'Not samex in con'

            NMC=100000
            
            dc(:,:)=0.0

            call Krigingestimate(2,N,x,sigmax,22,4,DAT(1001:1020),5,2,9,0,probtype,myflag,fmeantmp,fvartmp,fmeanprimetmp,fvarprimetmp)

!!            call Krigingestimate(2,X,N,sigmax,fmeantmp,fvartmp,fmeanprimetmp,fvarprimetmp,NMC,4)

            !fvartmp=0.0
            !fvarprimetmp(:)=0.0
            !call Eulersolve(x,n,0,fmeantmp,dfdD,dfdDD,1,v,4)
            !fmeanprimetmp(:)=dfdD(:) 
               
            do j=N-1,N
               dc(1,j)=-fmeanprimetmp(j)
               if (fvartmp.ne.0.0) then
                  dc(1,j)=dc(1,j)+kprob*fvarprimetmp(j)/(2.0*sqrt(fvartmp))
               endif
            end do


         end if

         call Eulersolve(Xsave,N,0,objtmp,dfdD,dfdDD,1,v,4)

         jac(1:N-2)=-dfdD(1:N-2)

         jac(N-1:N)=dc(1,N-1:N)


         do i=1,N
            X(i)=Xsave(i)
         end do
           
         !if (id_proc.eq.0) print *,'Cons Gradients',jac(1:6)

      end if


      IERR = 0
      return
      end
!
! =============================================================================
!
!                Computation of Hessian of Lagrangian
!
! =============================================================================
!
      subroutine EV_HESS(TASK, N, X, NEW_X, OBJFACT, M, LAM, NEW_LAM,NNZH, IRNH, ICNH, HESS, IDAT, DAT, IERR)
      implicit none
      integer TASK, N, NEW_X, M, NEW_LAM, NNZH, i, ir
      double precision X(N), OBJFACT, LAM(M), HESS(NNZH), sigmax(N)
      integer IRNH(NNZH), ICNH(NNZH)
      double precision DAT(*)
      integer IDAT(*), kprob
      integer IERR
      double precision  rho, L, sigmay, pi, p, E, Lifttarget

      Lifttarget=DAT(1)

      kprob=IDAT(1)
      do i=1,N
         sigmax(i)=DAT(i+1)
      end do


      if( TASK.eq.0 ) then
!
!     structure of sparse Hessian (lower triangle):
!
         IERR = 1

      else

         IERR = 1

      endif

      return
      end











!
! =============================================================================
!
!                   Callback method called once per iteration
!
! =============================================================================
!
      subroutine ITER_CB(ALG_MODE, ITER_COUNT,OBJVAL, INF_PR, INF_DU,MU, DNORM, REGU_SIZE, ALPHA_DU, ALPHA_PR, LS_TRIAL, IDAT,DAT, ISTOP)
        use dimKrig,only:probtype,id_proc
        implicit none
      integer ALG_MODE, ITER_COUNT, LS_TRIAL
      double precision OBJVAL, INF_PR, INF_DU, MU, DNORM, REGU_SIZE
      double precision ALPHA_DU, ALPHA_PR
      double precision DAT(*)
      integer IDAT(*)
      integer ISTOP
!
!     You can put some output here
!
      
      if (id_proc.eq.0) then

         if (ITER_COUNT .eq.0) then
            write(*,*) 
            write(*,*) 'iter    objective      ||grad||        inf_pr          inf_du         lg(mu)'
            write(76,*) 'iter    objective      ||grad||        inf_pr          inf_du         lg(mu)'
         end if

         write(*,'(i5,5e15.7)') ITER_COUNT,OBJVAL,DNORM,INF_PR,INF_DU,MU
         write(76,'(i5,5e15.7)') ITER_COUNT,OBJVAL,DNORM,INF_PR,INF_DU,MU

      end if
!
!     And set ISTOP to 1 if you want Ipopt to stop now.  Below is just a
!     simple example.
!
      if (ITER_COUNT .gt. 1 .and. DNORM.le.1D-04) ISTOP = 1

      return
      end
!====================



    subroutine epigrads(fct,fctindx,dim,ndimt,xtmp,xstdt,ftmp,dftmp)
      use omp_lib
      !    use dimKrig, only: DS,fctindx,reusesamples,ndimt,xavgt,xstdt
      implicit none
      integer :: DIM,ndimt,fct,fctindx
      !   real*8 :: x(DIM),f,df(DIM),d2f(DIM,DIM),scal,scal2,prd,time,time2

      real*8,intent(in)  :: xtmp(ndimt),xstdt(ndimt)
      real*8,intent(out) :: dftmp(ndimt)
      real*8::ftmp

      real*8 :: gtol,low(ndimt-DIM),up(ndimt-DIM)

      gtol=1e-6
      low(1:ndimt-DIM)=xtmp(1:ndimt-DIM)
      up(1:ndimt-DIM)=xtmp(1:ndimt-DIM)

      if(fctindx.eq.0) then

         call optimize(ndimt-DIM,xtmp,ndimt,ftmp,dftmp,low,up,gtol,.true.,.false.,fctindx)

      else if (fctindx.eq.4) then
         
         call optimize(ndimt-DIM,xtmp,ndimt,ftmp,dftmp,low,up,gtol,.false.,.false.,fctindx)

      else 

         stop'wrong fct indx'

      end if

      return
    end subroutine epigrads
 
