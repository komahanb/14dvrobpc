subroutine CalcstuffBFGS(X,ndvart,fobj,dfdD,fct,fctindx)
  implicit none

  integer  :: ndvart,fct,fctindx
  double precision :: X(ndvart),fobj,dfdD(ndvart),x3
  double precision ::  rho, L, sigmay, pi, p, E, Fs  
!  print*,X,fct,fctindx,ndvart
 ! print*,'x:',X 
  call get_f(ndvart,12,x,fobj)
  call get_df(ndvart,12,x,dfDD)
!stop
  return
end subroutine CalcstuffBFGS
