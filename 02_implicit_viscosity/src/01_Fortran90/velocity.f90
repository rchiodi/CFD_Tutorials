! ================================================ !
! These routines update the velocities.
! The velocity is updated from Uold to U^*
! with the convective and viscous terms. After
! the pressure is solved for, the velocity is also
! then accelerated by the pressure gradient.
! ================================================ !

module velocity
! Dummy module
end module velocity

! Initialize velocity variables
subroutine velocity_init

! No purpose as of now

	return
end subroutine velocity_init

! Update velocity from U^n to U^*
! with the convective and viscous terms
subroutine velocity_solve
	use data
	use indices
	implicit none
	
	! Store last times velocity
	Uold = U
	Vold = V
	
	! Advance velocity in time
	!call velocity_solver_visc_u_explicit
	call velocity_solver_visc_u_implicit
	call velocity_solver_conv_u
	!call velocity_solver_visc_v_explicit
	call velocity_solver_visc_v_implicit
	call velocity_solver_conv_v
		
	return
end subroutine velocity_solve

! Update U velocity with viscous term explicitly
subroutine velocity_solver_visc_u_explicit
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
	! the viscous part. U(imin,:) is on the boundary, and is set to 
  ! a no-slip BC
	do j = jmin,jmax
		do i = imin+1,imax
			U(i,j) = Uold(i,j)+dt*(&
							  	sum(cnt_divx(:,i,j)*Fx(i+div_m:i+div_p,j))+&
					 			  sum(cnt_divy(:,i,j)*Fy(i,j+div_m:j+div_p)) &
																														)
		end do
	end do
	
	return
end subroutine velocity_solver_visc_u_explicit

! Update U velocity with viscous term implicitly
subroutine velocity_solver_visc_u_implicit
	use data
	use parameters
	use operators
	use scratch
	use time_info
	use indices
	implicit none
	
	! Lower, Main, and Upper Diagonals in Operator Matrix
	real(DP), dimension((nx+2)*(ny+2)) :: ld, md, ud
	! Solution vector	
	real(DP), dimension((nx+2)*(ny+2)) :: sol
	! Right Hand Side vector
	real(DP),dimension((nx+2)*(ny+2)) :: RHS
	integer :: i,j,ind,q
	real(DP) :: coeff
	
	! Alternating Direction Implicit so that
	! the A matrix is triadiagonal, instead of Pentadiagonal with
	! a 3 bands.
	
	! First step, implicit in x, explicit in y
	! Create diagonals of the A operator matrix and the RHS vector
	ind = 0
	coeff = 0.5_DP*dt*mu/rho
	sol = 0.0_DP
	ld = 0.0_DP
	md = 0.0_DP
	ud = 0.0_DP
	RHS = 0.0_DP
	do j = jmin-1,jmax+1
		do i = imin, imax+1
			ind = ind + 1
			if(i.eq.imin) then ! Dirichlet BC
				ld(ind) = 0.0_DP
				md(ind) = 1.0_DP
				ud(ind) = 0.0_DP
				RHS(ind)= Uold(i,j)
			else if(i.eq.imax+1) then ! Dirichlet BC
				ld(ind) = 0.0_DP
				md(ind) = 1.0_DP
				ud(ind) = 0.0_DP
				RHS(ind)= Uold(i,j)
			else if(j.eq.jmin-1) then ! Dirichlet BC
				ld(ind) = 0.0_DP
				md(ind) = 1.0_DP
				ud(ind) = 0.0_DP
				RHS(ind)= Uold(i,j)
			else if(j.eq.jmax+1) then ! Dirichlet BC
				ld(ind) = 0.0_DP
				md(ind) = 1.0_DP
				ud(ind) = 0.0_DP
				RHS(ind)= Uold(i,j)
			else
				ld(ind) = -coeff*cnt_divx(0,i,j)*cnt_gradx(0,i-1,j)
				md(ind) =  1.0_DP-coeff*(cnt_divx(0,i,j)*&
									 cnt_gradx(1,i-1,j)+cnt_divx(1,i,j)*cnt_gradx(0,i,j))
				ud(ind) = -coeff*cnt_divx(1,i,j)*cnt_gradx(1,i,j)
				RHS(ind)=  coeff*(cnt_divy(0,i,j)*&
									 sum(cnt_grady(:,i,j-1)*Uold(i,j-1+grad_m:j-1+grad_p)) &
									 +cnt_divy(1,i,j)*&
									 sum(cnt_grady(:,i,j)*Uold(i,j+grad_m:j+grad_p)))&
									 +Uold(i,j)
			end if
		end do
	end do
	
	! Solve the system using Thomas Algorithm
	call ls_solve_thomas(ind,ld,md,ud,RHS,sol)
		
	! Map sol to 2D plane again
	do q = 1,ind
		j = (q-1)/(nx+1)
		i = mod(q-1,nx+1)+1
		U(i,j) = sol(q)
	end do 
	
	! Second step, implicit in y, explicit in x
	! Create diagonals of the A operator matrix and the RHS vector
	! Indexing order switched to maintain tight tridiagonal system
	ind = 0
	coeff = 0.5_DP*dt*mu/rho
	sol = 0.0_DP
	ld = 0.0_DP
	md = 0.0_DP
	ud = 0.0_DP
	RHS = 0.0_DP
	do i = imin, imax+1
		do j = jmin-1,jmax+1
			ind = ind + 1
			if(i.eq.imin) then ! Dirichlet BC
				ld(ind) = 0.0_DP
				md(ind) = 1.0_DP
				ud(ind) = 0.0_DP
				RHS(ind)= U(i,j)
			else if(i.eq.imax+1) then ! Dirichlet BC
				ld(ind) = 0.0_DP
				md(ind) = 1.0_DP
				ud(ind) = 0.0_DP
				RHS(ind)= U(i,j)
			else if(j.eq.jmin-1) then ! Dirichlet BC
				ld(ind) = 0.0_DP
				md(ind) = 1.0_DP
				ud(ind) = 0.0_DP
				RHS(ind)= U(i,j)
			else if(j.eq.jmax+1) then ! Dirichlet BC
				ld(ind) = 0.0_DP
				md(ind) = 1.0_DP
				ud(ind) = 0.0_DP
				RHS(ind)= U(i,j)
			else
				ld(ind) = -coeff*cnt_divy(0,i,j)*cnt_grady(0,i,j-1)
				md(ind) =  1.0_DP-coeff*(cnt_divy(0,i,j)*&
									 cnt_grady(1,i,j-1)+cnt_divy(1,i,j)*cnt_grady(0,i,j))
				ud(ind) = -coeff*cnt_divy(1,i,j)*cnt_grady(1,i,j)
				RHS(ind)=  coeff*(cnt_divx(0,i,j)*&
									 sum(cnt_gradx(:,i-1,j)*U(i-1+grad_m:i-1+grad_p,j)) &
									 +cnt_divx(1,i,j)*&
									 sum(cnt_gradx(:,i,j)*U(i+grad_m:i+grad_p,j)))&
									 +U(i,j)
			end if
		end do
	end do
	
	! Solve the system using Thomas Algorithm
	call ls_solve_thomas(ind,ld,md,ud,RHS,sol)

	! Map sol to 2D plane again
	do q = 1,ind
		j = mod(q-1,ny+2)
		i = (q-1)/(ny+2)+1
		U(i,j) = sol(q)
	end do 
		
	return
end subroutine velocity_solver_visc_u_implicit

! Update U velocity with convective term
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
	! to n^*. U(imin,:) is on the boundary and is set with no-slip BC
	do j = jmin,jmax
		do i = imin+1,imax
			U(i,j) = U(i,j)+dt*(&
												 sum(cnt_divx(:,i,j)*Fx(i+div_m:i+div_p,j))+&
												 sum(cnt_divy(:,i,j)*Fy(i,j+div_m:j+div_p)) &
																																	 )
		end do
	end do
		
	return
end subroutine velocity_solver_conv_u

! Update V velocity with viscous term explicitly
subroutine velocity_solver_visc_v_explicit
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
	! the viscous part. V(:,jmin) is on the boundary, and is set to
	! no-slip BC
	do j = jmin+1,jmax
		do i = imin,imax
			V(i,j) = Vold(i,j)+dt*(&
							  	sum(cnt_divx(:,i,j)*Fx(i+div_m:i+div_p,j))+&
					 			  sum(cnt_divy(:,i,j)*Fy(i,j+div_m:j+div_p)) &
																														)
		end do
	end do
	
	return
end subroutine velocity_solver_visc_v_explicit

! Update V velocity with viscous term implicitly
subroutine velocity_solver_visc_v_implicit
	use data
	use parameters
	use operators
	use scratch
	use time_info
	use indices
	implicit none
	
	! Lower, Main, and Upper Diagonals in Operator Matrix
	real(DP), dimension((nx+2)*(ny+2)) :: ld, md, ud
	! Solution vector	
	real(DP), dimension((nx+2)*(ny+2)) :: sol
	! Right Hand Side vector
	real(DP),dimension((nx+2)*(ny+2)) :: RHS
	integer :: i,j,ind,q
	real(DP) :: coeff
	
	! Alternating Direction Implicit so that
	! the A matrix is triadiagonal, instead of Pentadiagonal with
	! a large band.
	
	! First step, implicit in x, explicit in y
	! Create diagonals of the A operator matrix and the RHS vector
	ind = 0
	coeff = 0.5_DP*dt*mu/rho
	sol = 0.0_DP
	ld = 0.0_DP
	md = 0.0_DP
	ud = 0.0_DP
	RHS = 0.0_DP
	do j = jmin,jmax+1
		do i = imin-1, imax+1
			ind = ind + 1
			if(i.eq.imin-1) then ! Dirichlet BC
				ld(ind) = 0.0_DP
				md(ind) = 1.0_DP
				ud(ind) = 0.0_DP
				RHS(ind)= Vold(i,j)
			else if(i.eq.imax+1) then ! Dirichlet BC
				ld(ind) = 0.0_DP
				md(ind) = 1.0_DP
				ud(ind) = 0.0_DP
				RHS(ind)= Vold(i,j)
			else if(j.eq.jmin) then ! Dirichlet BC
				ld(ind) = 0.0_DP
				md(ind) = 1.0_DP
				ud(ind) = 0.0_DP
				RHS(ind)= Vold(i,j)
			else if(j.eq.jmax+1) then ! Dirichlet BC
				ld(ind) = 0.0_DP
				md(ind) = 1.0_DP
				ud(ind) = 0.0_DP
				RHS(ind)= Vold(i,j)
			else
				ld(ind) = -coeff*cnt_divx(0,i,j)*cnt_gradx(0,i-1,j)
				md(ind) =  1.0_DP-coeff*(cnt_divx(0,i,j)*&
									 cnt_gradx(1,i-1,j)+cnt_divx(1,i,j)*cnt_gradx(0,i,j))
				ud(ind) = -coeff*cnt_divx(1,i,j)*cnt_gradx(1,i,j)
				RHS(ind)=  coeff*(cnt_divy(0,i,j)*&
									 sum(cnt_grady(:,i,j-1)*Vold(i,j-1+grad_m:j-1+grad_p)) &
									 +cnt_divy(1,i,j)*&
									 sum(cnt_grady(:,i,j)*Vold(i,j+grad_m:j+grad_p)))&
									 +Vold(i,j)
			end if
		end do
	end do
	
	! Solve the system using Thomas Algorithm
	call ls_solve_thomas(ind,ld,md,ud,RHS,sol)
		
	! Map sol to 2D plane again
	do q = 1,ind
		j = (q-1)/(nx+2)+1
		i = mod(q-1,nx+2)
		V(i,j) = sol(q)
	end do 
	
	! Second step, implicit in y, explicit in x
	! Create diagonals of the A operator matrix and the RHS vector
	! Indexing order switched to maintain tight tridiagonal system
	ind = 0
	coeff = 0.5_DP*dt*mu/rho
	sol = 0.0_DP
	ld = 0.0_DP
	md = 0.0_DP
	ud = 0.0_DP
	RHS = 0.0_DP
	do i = imin-1, imax+1
		do j = jmin,jmax+1
			ind = ind + 1
			if(i.eq.imin-1) then ! Dirichlet BC
				ld(ind) = 0.0_DP
				md(ind) = 1.0_DP
				ud(ind) = 0.0_DP
				RHS(ind)= V(i,j)
			else if(i.eq.imax+1) then ! Dirichlet BC
				ld(ind) = 0.0_DP
				md(ind) = 1.0_DP
				ud(ind) = 0.0_DP
				RHS(ind)= V(i,j)
			else if(j.eq.jmin) then ! Dirichlet BC
				ld(ind) = 0.0_DP
				md(ind) = 1.0_DP
				ud(ind) = 0.0_DP
				RHS(ind)= V(i,j)
			else if(j.eq.jmax+1) then ! Dirichlet BC
				ld(ind) = 0.0_DP
				md(ind) = 1.0_DP
				ud(ind) = 0.0_DP
				RHS(ind)= V(i,j)
			else
				ld(ind) = -coeff*cnt_divy(0,i,j)*cnt_grady(0,i,j-1)
				md(ind) =  1.0_DP-coeff*(cnt_divy(0,i,j)*&
									 cnt_grady(1,i,j-1)+cnt_divy(1,i,j)*cnt_grady(0,i,j))
				ud(ind) = -coeff*cnt_divy(1,i,j)*cnt_grady(1,i,j)
				RHS(ind)=  coeff*(cnt_divx(0,i,j)*&
									 sum(cnt_gradx(:,i-1,j)*V(i-1+grad_m:i-1+grad_p,j)) &
									 +cnt_divx(1,i,j)*&
									 sum(cnt_gradx(:,i,j)*V(i+grad_m:i+grad_p,j)))&
									 +V(i,j)
			end if
		end do
	end do
	
	! Solve the system using Thomas Algorithm
	call ls_solve_thomas(ind,ld,md,ud,RHS,sol)

	! Map sol to 2D plane again
	do q = 1,ind
		j = mod(q-1,ny+1)+1
		i = (q-1)/(ny+1)
		V(i,j) = sol(q)
	end do 
		
	return
	
	
	return
end subroutine velocity_solver_visc_v_implicit

! Update V velocity with convective term
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
	! to n^*. V(:,jmin) is on the boundary and is set to no-slip BC
	do j = jmin+1,jmax
		do i = imin,imax
			V(i,j) = V(i,j)+dt*(&
												 sum(cnt_divx(:,i,j)*Fx(i+div_m:i+div_p,j))+&
												 sum(cnt_divy(:,i,j)*Fy(i,j+div_m:j+div_p)) &
																																	 )
		end do
	end do
		
	return
end subroutine velocity_solver_conv_v

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


