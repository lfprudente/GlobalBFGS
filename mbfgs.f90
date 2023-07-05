subroutine mbfgsupdate(n,m,x,xprev,JF,JFprev,lambda,B)

implicit none
	
	! SCALAR ARGUMENTS 
	integer,intent(in) :: n,m
	
	! ARRAY ARGUMENTS
	real(kind=8), intent(in) :: x(n),xprev(n),lambda(m),JF(m,n),JFprev(m,n)
	real(kind=8), intent(inout) :: B(n,n,m)
	
	! LOCAL SCALARS
	integer :: j,ind
	real(kind=8) :: eta,r,sTy,sTBs,eps,norm2dlambda
	
	! LOCAL ARRAYS
	real(kind=8) :: s(n),y(n),dlambda(n),Bs(n),BssTB(n,n),yyT(n,n),ysTB(n,n),BsyT(n,n)
	
	eps = 1.0d-1
	
	! Define s
	
	s = x - xprev
	
	! Compute dlambda
	
	dlambda(:) = 0.0d0
	do j = 1,m
		dlambda = dlambda + lambda(j) * JFprev(j,:)
	end do
	
	norm2dlambda = norm2(dlambda)
	
	do ind = 1,m
		
		! Define y 
		
		y = JF(ind,:) - JFprev(ind,:)
		
		! Compute eta
		
		eta = dot_product( s, y ) / norm2(s)**2
		
		! Compute r
		
		r = max( - eta, 0.0d0 ) + eps * norm2dlambda
		
		! Compute gamma (actually, update y)
		
		y = y + r * s
		
		! Compute <s,y>
		
		sTy = dot_product( s, y )
		
		! Compute Bs
		
		Bs = matmul( B(:,:,ind), s )
		
		! Compute Bss^TB
		
		do j = 1,n
			BssTB(:,j) = Bs(j) * Bs
		end do
		
		! Compute yy^T
		
		do j = 1,n
			yyT(:,j) = y(j) * y
		end do
		
		! Compute sTBs
		
		sTBs = dot_product( s, Bs )
		
		! Update B

		B(:,:,ind) = B(:,:,ind) - BssTB / sTBs + yyT / sTy

	end do	

end subroutine mbfgsupdate
