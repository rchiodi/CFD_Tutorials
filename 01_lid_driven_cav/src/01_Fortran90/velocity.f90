module velocity
! Dummy module
end module velocity

subroutine velocity_init

! No purpose as of now

	return
end subroutine velocity_init

subroutine velocity_solve
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
end subroutine velocity_solve

subroutine velocity_solver_conv_u
	use data
	use operators
	use scratch
	use time_info
	use indices
	implicit none
	
	integer :: i,j
	
	! Calculating div(U \vec{U})
	
	! Loop and compute convective fluxes at left and bottom of U-cell
	! Flux is U \vec{U}
	Fx = 0.0_DP
	Fy = 0.0_DP
	do j = jmin,jmax+1
		do i = imin, imax+1
			Fx(i,j) = -sum(int_x_c(:,i,j)*Uold(i-1+intxc_m:i-1+intxc_p,j))*&
								 sum(int_x_c(:,i,j)*Uold(i-1+intxc_m:i-1+intxc_p,j))
								
			Fy(i,j) = -sum(int_c_x(:,i,j)*Vold(i+intcx_m:i+intcx_p,j))*&
								 sum(int_c_y(:,i,j)*Uold(i,j+intcy_m:j+intcy_p))
		end do
	end do
	
	! Take divergence with Fx and Fy for each cell, update U from n&
	! to n+1/2. U(imin,:) is on the boundary and will be set with a BC
	do j = jmin,jmax
		do i = imin+1,imax
			U(i,j) = Uold(i,j)+dt*(&
												 sum(cnt_divx(:,i,j)*Fx(i+div_m:i+div_p,j))+&
												 sum(cnt_divy(:,i,j)*Fy(i,j+div_m:j+div_p)) &
																																	 )
		end do
	end do
		
	return
end subroutine velocity_solver_conv_u

subroutine velocity_solver_visc_u
	use data
	use parameters
	use operators
	use scratch
	use time_info
	use indices
	implicit none
	
	integer :: i,j

	! Calculating mu/rho * lap(U)
			
	! Loop and compute viscous fluxes at left and bottom of U-cell
	! Flux is mu/rho*grad(U)
	Fx = 0.0_DP
	Fy = 0.0_DP
	do j = jmin, jmax+1
		do i = imin, imax+1
			Fx(i,j) = sum(cnt_gradx(:,i-1,j)*Uold(i-1+grad_m:i-1+grad_p,j))
			Fy(i,j) = sum(cnt_grady(:,i,j-1)*Uold(i,j-1+grad_m:j-1+grad_p))
		end do
	end do
	! Multiply fluxes by mu/rho
	Fx = mu/rho*Fx
	Fy = mu/rho*Fy
		
	! Take divergence with Fx and Fy for each cell, update U to include
	! the viscous part. U(imin,:) is on the boundary, and will be set 
	! with a BC
	do j = jmin,jmax
		do i = imin+1,imax
			U(i,j) = U(i,j)+dt*(&
							  	sum(cnt_divx(:,i,j)*Fx(i+div_m:i+div_p,j))+&
					 			  sum(cnt_divy(:,i,j)*Fy(i,j+div_m:j+div_p)) &
																														)
		end do
	end do
	
	return
end subroutine velocity_solver_visc_u

subroutine velocity_solver_conv_v
	use data
	use operators
	use scratch
	use time_info
	use indices
	implicit none
	
	integer :: i,j
	
	! Calculating div(V \vec{U})
	
	! Loop and compute convective fluxes at left and bottom of V-cell
	! Flux is V \vec{U}
	Fx = 0.0_DP
	Fy = 0.0_DP
	do j = jmin,jmax+1
		do i = imin, imax+1
			Fx(i,j) = -sum(int_c_y(:,i,j)*Uold(i,j+intcy_m:j+intcy_p))*&
								 sum(int_c_x(:,i,j)*Vold(i+intcx_m:i+intcx_p,j))
								
			Fy(i,j) = -sum(int_y_c(:,i,j-1)*Vold(i,j-1+intyc_m:j-1+intyc_p))*&
								 sum(int_y_c(:,i,j-1)*Vold(i,j-1+intyc_m:j-1+intyc_p))
		end do
	end do
	
	! Take divergence with Fx and Fy for each cell, update U from n&
	! to n+1/2. U(imin,:) is on the boundary and will be set with a BC
	do j = jmin+1,jmax
		do i = imin,imax
			V(i,j) = Vold(i,j)+dt*(&
												 sum(cnt_divx(:,i,j)*Fx(i+div_m:i+div_p,j))+&
												 sum(cnt_divy(:,i,j)*Fy(i,j+div_m:j+div_p)) &
																																	 )
		end do
	end do
		
	return
end subroutine velocity_solver_conv_v

subroutine velocity_solver_visc_v
	use data
	use parameters
	use operators
	use scratch
	use time_info
	use indices
	implicit none
	
	integer :: i,j

	! Calculating mu/rho * lap(V)
			
	! Loop and compute viscous fluxes at left and bottom of U-cell
	! Flux is mu*grad(V)
	Fx = 0.0_DP
	Fy = 0.0_DP
	do j = jmin, jmax+1
		do i = imin, imax+1
			Fx(i,j) = sum(cnt_gradx(:,i-1,j)*Vold(i-1+grad_m:i-1+grad_p,j))
			Fy(i,j) = sum(cnt_grady(:,i,j-1)*Vold(i,j-1+grad_m:j-1+grad_p))
		end do
	end do
	! Multiply fluxes by mu
	Fx = mu/rho*Fx
	Fy = mu/rho*Fy
	
	! Take divergence with Fx and Fy for each cell, update U to include
	! the viscous part. U(imin,:) is on the boundary, and will be set 
	! with a BC
	do j = jmin+1,jmax
		do i = imin,imax
			V(i,j) = V(i,j)+dt*(&
							  	sum(cnt_divx(:,i,j)*Fx(i+div_m:i+div_p,j))+&
					 			  sum(cnt_divy(:,i,j)*Fy(i,j+div_m:j+div_p)) &
																														)
		end do
	end do
	
	return
end subroutine velocity_solver_visc_v

! Update velocity to solenoidal velocity field,
! from time n^* to n+1
subroutine velocity_correct
	use data
	use operators
	use indices
	use time_info
	implicit none
	
	integer :: i,j
	

	! Update U velocity
	do j = jmin, jmax
		do i = imin+1, imax
			U(i,j) = U(i,j) - dt*sum(cnt_gradx(:,i-1,j)*&
															 P(i-1+grad_m:i-1+grad_p,j))
		end do
	end do

	! Update V velocity
	do j = jmin+1, jmax
		do i = imin, imax
			V(i,j) = V(i,j) - dt*sum(cnt_grady(:,i,j-1)*&
															 P(i,j-1+grad_m:j-1+grad_p))
		end do
	end do
	
	return
end subroutine velocity_correct


