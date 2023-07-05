program main

	use globals, only: problem,seed
	use myproblem, only: inip

	implicit none
	
	! LOCAL SCALARS
	integer :: n,m,outiter,nfev,ngev,inform,solver
	logical :: scaleF,checkder
	real(kind=4) :: time
	real(kind=8) :: theta
	
	! LOCAL ARRAYS
	real(kind=8),allocatable :: x(:),l(:),u(:)
	logical,allocatable :: strconvex(:)
	
!	Choose the solver
	
!	solver = 1: Algorithm 1 (Global BFGS)
!	solver = 2: BFGS-Wolfe
!	solver = 3: Cautious BFGS-Armijo
	
	solver = 1
	
!	Choose the problem to be solved
	
	problem = 'AP1'

	seed = 123456.0d0

	call inip(n,m,x,l,u,strconvex,scaleF,checkder)
	
!	Call the solver
	
	if ( solver == 1 ) call BFGS(n,m,x,l,u,strconvex,scaleF,checkder,outiter,time,nfev,ngev,theta,inform,1)
	if ( solver == 2 ) call BFGS(n,m,x,l,u,strconvex,scaleF,checkder,outiter,time,nfev,ngev,theta,inform,2)
	if ( solver == 3 ) call cautiousBFGS(n,m,x,l,u,scaleF,checkder,outiter,time,nfev,ngev,theta,inform)
	
end program main
