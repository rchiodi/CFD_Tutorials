! ================================================= !
! These subroutines print the important information
! to the screen while the simulation is running.
! Allows user to monitor the maximum velocity,
! convective and viscous CFLs, and the convergence
! of the pressure solver
! ================================================= !

module monitor

! Dummy module

end module monitor

subroutine monitor_init
	implicit none
	
	! Nothing done in this for now
	
	
	return
end subroutine monitor_init

! This subroutine writes out the CFLs, Velocity,
! and pressure convergence criteria to the screen
! to monitor
subroutine monitor_write
	use data
	use time_info
	use geometry
	use pressure
	use indices
	use parameters
	implicit none
	real(DP) :: Umax, Vmax
	real(DP) :: cCFL, vCFL
	
	! Calculate the max velocity in domain
	Umax = maxval(abs(U(imin:imax+1,jmin:jmax)))
	Vmax = maxval(abs(V(imin:imax,jmin:jmax+1)))
	
	! Calculate the convective and viscous CFLs
	cCFL = max(Umax,Vmax)*dt*idx
	vCFL = 2.0_DP*dt*mu*idx*idx/rho

	! Print to screen
	write(*,'((I10.1,1x),(ES10.4,1x),2(F8.4,1x),4(ES10.4,1x))') &
				titer, time, cCFL, vCFL, Umax, Vmax, Linf_fin, L2_fin

	return
end subroutine monitor_write
