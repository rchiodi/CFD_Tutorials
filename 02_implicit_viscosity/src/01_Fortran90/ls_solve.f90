module ls_solve
use precision
implicit none

! Dummy module

end module ls_solve

subroutine ls_solve_thomas(n,a,b,c,d,sol)
	use ls_solve
	implicit none
	integer, intent(in) :: n
	real(DP), dimension(n), intent(inout) :: a,b,c,d
	real(DP), dimension(n), intent(out) :: sol
	
	integer :: k
	real(DP) :: m
	
	do k = 2,n
		m = a(k)/b(k-1)
		b(k) = b(k)-m*c(k-1)
		d(k) = d(k)-m*d(k-1)
	end do
	
	sol(n) = d(n)/b(n)
	do k = n-1,1,-1
		sol(k) = (d(k)-c(k)*sol(k+1))/b(k)
	end do
	
	return
end subroutine
