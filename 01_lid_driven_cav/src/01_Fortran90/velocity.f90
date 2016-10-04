module velocity
! Dummy module
end module velocity

subroutine velocity_init

! No purpose as of now

	return
end subroutine velocity_init

subroutine velocity_solver
	use data
	implicit none
	
	! Store last times velocity
	Uold = U
	Vold = V
	
	! Advance velocity in time
	call velocity_solver_conv_u
	call velocity_solver_visc_u
	call velocity_solver_conv_v
	call velocity_solver_visc_v
	
	return
end subroutine velocity_solver

subroutine velocity_solver_conv_u
	use data
	use operators
	use scratch
	implicit none
	
	integer :: i,j
	
	! Calculating div(U \vec{U})
	
	! Loop and compute convective fluxes at left and bottom of U-cell
	! Flux is U \vec{U}
	do j = jmin,jmax+1
		do i = imin, imax+1
			Fx(i,j) = sum(int_x_c(:,i,j)*Uold(i-1+intxc_m:i-1+intxc_p,j))*&
								sum(int_x_c(:,i,j)*Uold(i-1+intxc_m:i-1+intxc_p,j))
								
			Fy(i,j) = sum(int_c_x(:,i,j)*Vold(i+intcx_m:i+intcx_p,j))*&
								sum(int_c_y(:,i,j)*Uold(i,j+intcy_m:j+intcy_p))
		end do
	end do
	
	! Take divergence with Fx and Fy for each cell, update U from n&
	! to n+1/2. U(imin,:) is on the boundary
	do j = jmin,jmax
		do i = imin+1,imax
			U(i,j) = Uold(i,j)+dt*(&
												 sum(cnt_divx(:,i,j)*Fx(i+div_m:i+div_p,j))&
												 sum(cnt_divy(:,i,j)*Fy(i,j+div_m:j+div_p))&
																																	 )
		end do
	end do
		
	return
end velocity_solver_conv_u

subroutine velocity_solver_visc_u
	use data
	use operators
	use scratch
	implicit none
	
	integer :: i,j

	! Calculating mu/rho * lap(U)
	
	! Loop and compute viscous fluxes at left and bottom of U-cell
	! Flux is mu*grad(U)
	do j = jmin, jmax+1
		do i = imin, imax+1
			Fx(i,j) = sum(int_x_c(:,i,j)*Uold(i+intxc_m:i+intxc_p,j))
		end do
	end do
	
	
	
	return
end subroutine velocity_solver_visc_u

