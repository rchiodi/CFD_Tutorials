! ==================================================================== !
! The pressure required to accelerate the fluid to having
! a solenoidal velocity field is found through solving Ax=b,
! where A is the operator matrix representing the discrete gradient
! operator, and b is the divergence of the original velocity field at
! time n^*. In reality, we will not form the matrices in order to
! reduce the memory use. For now, the Gauss-Seidel method will be used
! to solve the system.
! ==================================================================== !

module pressure
	use precision
   implicit none

	real(DP), dimension(:,:), allocatable :: RHS
	real(DP) :: Linf_fin, L2_fin

end module pressure

! Initialize variables needed for the pressure solver
subroutine pressure_init
	use pressure
	use indices
	use scratch
	implicit none

	! Allocate the RHS for the pressure solution
	allocate(RHS(imin:imax,jmin:jmax))

	return
end subroutine pressure_init

! Solve pressure Poisson equation to get pressure
! that leads to solenoidal velocity field
subroutine pressure_solve
	use data
	use operators
	use pressure
	implicit none
	
	real(DP) :: Linf_err, L2_err, Linf_init, L2_init
	real(DP), parameter :: ptol = 1.0e-4_DP
	integer :: piter
	integer, parameter :: piter_max = 10000
	integer, parameter :: check_interval = 20
	
	! Calculate the RHS (1/dt* div(\vec{U}))
	call pressure_RHS
	
	! Solve Ax = b with Gauss-Seidel
	Linf_err = huge(1.0_DP)
	L2_err = huge(1.0_DP)
	call pressure_solve_divfree_check(Linf_init, L2_init)
	if (Linf_init .lt. ptol) then
		Linf_fin = Linf_init
		L2_fin = L2_init
		return
	end if

   ! Perform pressure iterations
	do piter = 1, piter_max

		call pressure_solve_gauss_seidel
		! Only calculate the error every check_interval iterations
		if(mod(piter,check_interval) .eq. 0) then
			call pressure_solve_divfree_check(Linf_err, L2_err)
		end if
	
		! If the Linf(lap(P)<ptol), we've converged
		if(Linf_err .lt. ptol) exit
	end do
  ! Print to screen if max iterations reached
	if(piter.eq.piter_max+1) print*,'max iterations reached'
	! Store final Linf_err to export to screen
	Linf_fin = Linf_err
	L2_fin = L2_err
		
	return
end subroutine pressure_solve


! Calculate the right hand side of the pressure Poisson equation
! This is what we'll want to force to zero to satisfy the divergence
! free condition
subroutine pressure_RHS
	use pressure
	use operators
	use data
	use indices
	use time_info
	implicit none
	
	integer :: i, j
	real(DP) :: idt
	
	idt = 1.0_DP/dt

	do j = jmin, jmax
		do i = imin, imax
         ! 1/dt * div{ \vec{U} }
			RHS(i,j) = idt*(sum(cnt_divx(:,i,j)*U(i+div_m:i+div_p,j))+ &
									sum(cnt_divy(:,i,j)*V(i,j+div_m:j+div_p)) )
		end do
	end do
	
	return
end subroutine pressure_RHS

! Update P at each location according to Gauss-Seidel method
subroutine pressure_solve_gauss_seidel
	use data
	use operators
	use pressure
	use indices
	implicit none
	
	integer :: i,j, st
	real(DP) :: Rx

	do j = jmin, jmax
		do i = imin, imax
			Rx = 0.0_DP
			do st = lap_m,lap_p,2
				Rx = Rx + lapx(st,i,j)*P(i+st,j)+lapy(st,i,j)*P(i,j+st)
			end do
			P(i,j) = (RHS(i,j)-Rx)/(lapx(0,i,j)+lapy(0,i,j))
		end do
	end do
	
	return
end subroutine pressure_solve_gauss_seidel

! Calculate the L infinity and L2 norm of the residual
! to check for adequate convergence
subroutine pressure_solve_divfree_check(Linf_err, L2_err)
	use pressure
	use data
	use operators
	use indices
	implicit none
	
	real(DP), intent(out) :: Linf_err, L2_err
	real(DP) :: res
	integer :: i, j

	! Check that lap(P) = 0
	Linf_err = -huge(1.0_DP)
	L2_err = 0.0_DP
	do j = jmin, jmax
		do i = imin, imax
			res = abs(sum(lapx(:,i,j)*P(i+lap_m:i+lap_p,j)) + &
								sum(lapy(:,i,j)*P(i,j+lap_m:j+lap_p))-RHS(i,j))
			Linf_err = max(res,Linf_err)
			L2_err = L2_err + res*res
		end do
	end do
	L2_err = sqrt(L2_err)

	return
end subroutine pressure_solve_divfree_check
