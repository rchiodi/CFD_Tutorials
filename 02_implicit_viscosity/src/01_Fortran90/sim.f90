!> This is to solve the lid driven cavity case.
!! We are going to make a lot of assumptions and write a pretty
!! innacurate code, but it should be simpler to follow and still
!! produce reasonable results.
!! In the future, improvements will be made in accuracy, speed, and
!! capability.

module precision
!> Double (DP) and Single precision parameters
integer, parameter :: DP = kind(1.d0)
integer, parameter :: SP = kind(1.0)

end module precision

module time_info
	use precision
	integer :: titer
	real(DP) :: dt, time, tfin, write_freq

end module time_info

module parameters
	use precision
	real(DP) :: dw, dh, Ulid, rho, mu, init_cCFL, init_vCFL
	integer :: nx, ny ,ng
	character(len=80) :: fstart

end module parameters

module geometry
	use precision
	real(DP), dimension(:), allocatable :: x, y, xm, ym
	real(DP) :: dx, dy, idx, idy, area, iarea

end module geometry

module indices

	integer :: iming,imaxg,jming,jmaxg
	integer :: imin, imax, jmin, jmax

end module indices

module operators
	use precision
	! Derivatives for divergence of cell
  real(DP), dimension(:,:,:), allocatable :: cnt_divx, cnt_divy
	integer :: div_m, div_p
	
	! Derivatives for gradient of U or V cell
	real(DP), dimension(:,:,:), allocatable :: cnt_gradx, cnt_grady
	integer :: grad_m, grad_p
	
	! Derivatives for gradient of pressure (Neumann BC on boundaries used)
	real(DP), dimension(:,:,:), allocatable :: lapx, lapy
	integer :: lap_m, lap_p
	

	! Interpolation
	real(DP), dimension(:,:,:), allocatable :: int_x_c, int_y_c
	real(DP), dimension(:,:,:), allocatable :: int_c_x, int_c_y
	integer :: intxc_m, intxc_p, intyc_m, intyc_p, intcx_m, intcx_p
	integer :: intcy_m, intcy_p

end module operators

module data
	use precision
	real(DP), dimension(:,:), allocatable :: U, V, P
	real(DP), dimension(:,:), allocatable :: Uold, Vold

end module data

module scratch
	use precision
	real(DP), dimension(:,:), allocatable :: Fx, Fy
	real(DP), dimension(:,:), allocatable, target :: tmp1, tmp2

end module scratch

!> This just calls configuration, sets the initial data
!! and then starts the flow solver
program main

	call parameters_set
	call grid_make
	call operator_gen
	call data_init
	call flow_solve

	return
end program main

!> Parameters specified by user
!! Normally, this would be set through an input file or a GUI,
!! but for simplicity to start, these will be hardcoded here in this
!! function.
subroutine parameters_set
	use parameters
	use geometry
	use time_info
	implicit none
	! Domain width
	dw = 1.0_DP
	! Cells across the width
	nx = 50
	! Domain height
	dh = 1.0_DP
	! Cells across the height
	ny = 50
	! Lid velocity
	Ulid = 1.0_DP
	! Density of fluid
	rho = 100.0_DP
	! Dynamic viscosity of fluid
	mu = 1.0_DP
	
	! Number of ghost cells on each boundary
	ng = 2
	
	! Initial time step size
	init_cCFL = 0.85_DP !85% initial CFL
	init_vCFL = 0.5_DP!100.0_DP  ! Viscous CFL 
	! End time of simulation
	tfin = 50.0_DP
	! Time in between data dumps
	write_freq = 10.0_DP
	! File name used for data written to file
	fstart = 'data'

	return
end subroutine parameters_set
	
!> Sets up the configuration of the grid (x,y locations).
!! Rectangular grid is assumed, with a staggered cell arrangement
subroutine grid_make
	use parameters
	use geometry
	use indices
	implicit none
	
	integer :: i,j
		
	! Loop extents, g with ghost, otherwise in physical domain
	iming = 1-ng
	imaxg = nx+ng
	jming = 1-ng
	jmaxg = ny+ng
	imin = 1
	imax = nx
	jmin = 1
	jmax = ny
	
	! Allocate for physical domain, plus ng ghost cells on each side
	allocate(x(iming:imaxg))
	allocate(y(jming:jmaxg))
	allocate(xm(iming:imaxg))
	allocate(ym(jming:jmaxg))
	
	! Cell side lengths and area
	dx = dw/real(nx,DP)
	dy = dh/real(ny,DP)
	idx = 1.0_DP/dx
	idy = 1.0_DP/dy
	area = dx*dy
	iarea = 1.0_DP/area
	
	! Left cell faces
	do i = iming,imaxg
		x(i) = real(i-1,DP)*dx-0.5_DP*dw
	end do
	
	! Bottom cell faces
	do j = jming,jmaxg
		y(j) = real(j-1,DP)*dy-0.5_DP*dh
	end do
	
	! Cell centers
	do i = iming,imaxg
		xm(i) = x(i)+0.5_DP*dx
	end do
	do j = jming,jmaxg
		ym(j) = y(j)+0.5_DP*dy
	end do

	return
end subroutine grid_make

!> Creates the operators based on the mesh.
!! Operator is created for each cell, which will allow
!! imposing some boundary conditions through the operators themselves.
!! This is especially helpful for Neumann BCs (e.g. for pressure)
subroutine operator_gen
	use geometry
	use operators
	use indices
	implicit none
	
	integer :: i,j
	
	! Second order centered derivative, assume uniform mesh
	! Used to calculate divergence
	div_m = 0
	div_p = 1
	allocate(cnt_divx(div_m:div_p,iming-div_m:imaxg-div_p,jming-div_m:jmaxg-div_p))
	allocate(cnt_divy(div_m:div_p,iming-div_m:imaxg-div_p,jming-div_m:jmaxg-div_p))
	cnt_divx(0,:,:) = -1.0_DP
	cnt_divx(1,:,:) =  1.0_DP
	cnt_divy = cnt_divx*idy
	cnt_divx = cnt_divx*idx
	
	! Second order centered derivative, assume uniform mesh
	! Used to calculate gradient
	grad_m = 0
	grad_p = 1
	allocate(cnt_gradx(grad_m:grad_p,iming-grad_m:imaxg-grad_p,jming-grad_m:jmaxg-grad_p))
	allocate(cnt_grady(grad_m:grad_p,iming-grad_m:imaxg-grad_p,jming-grad_m:jmaxg-grad_p))
	cnt_gradx(0,:,:) = -1.0_DP
	cnt_gradx(1,:,:) =  1.0_DP
	cnt_grady = cnt_gradx*idy
	cnt_gradx = cnt_gradx*idx
	
	! Laplacian operator, assume uniform mesh
	! Used for evaluation of lap(P), assumes Neumann BC on domain boundary
	lap_m = -1
	lap_p =  1
	allocate(lapx(lap_m:lap_p,iming-lap_m:imaxg-lap_p,jming-lap_m:jmaxg-lap_p))
	allocate(lapy(lap_m:lap_p,iming-lap_m:imaxg-lap_p,jming-lap_m:jmaxg-lap_p))
	lapx(-1,:,:) =  1.0_DP
	lapx( 0,:,:) =  -2.0_DP
	lapx( 1,:,:) =  1.0_DP
	lapy = lapx
	! Add the Neumann BC on Boundaries
	lapx(-1,imin,:) = 0.0_DP; lapx(0,imin,:) = -1.0_DP ! Left wall
	lapx( 1,imax,:) = 0.0_DP; lapx(0,imax,:) = -1.0_DP ! Right wall
	lapy(-1,:,jmin) = 0.0_DP; lapy(0,:,jmin) = -1.0_DP ! Bottom wall
	lapy( 1,:,jmax) = 0.0_DP; lapy(0,:,jmax) = -1.0_DP ! Top wall
	! Divide by dx2
	lapy = lapy*idy*idy
	lapx = lapx*idx*idx
	
	! Interpolation from X-face to cell center (First order)
	intxc_m= 0
	intxc_p = 1
	allocate(int_x_c(intxc_m:intxc_p,iming-intxc_m:imaxg-intxc_p, &
																	 jming-intxc_m:jmaxg-intxc_p))
	int_x_c(0,:,:) = 0.5_DP
	int_x_c(1,:,:) = 0.5_DP

	! Interpolation from Y-face to cell center (First order)
	intyc_m = 0
	intyc_p = 1
	allocate(int_y_c(intyc_m:intyc_p,iming-intyc_m:imaxg-intyc_p, &
																	 jming-intyc_m:jmaxg-intyc_p))
	int_y_c(0,:,:) = 0.5_DP
	int_y_c(1,:,:) = 0.5_DP
		
	! Interpolation from cell center to X-face (First order)
	intcx_m = -1
	intcx_p =  0
	allocate(int_c_x(intcx_m:intcx_p,iming-intcx_m:imaxg-intcx_p, &
																	 jming-intcx_m:jmaxg-intcx_p))
	int_c_x(-1,:,:) = 0.5_DP
	int_c_x( 0,:,:) = 0.5_DP
	
	! Interpolation from cell center to Y-face (First order)
	intcy_m = -1
	intcy_p =  0
	allocate(int_c_y(intcy_m:intcy_p,imaxg-intcy_m:imaxg-intcy_p, &
																	 jming-intcy_m:jmaxg-intcy_p))
	int_c_y(-1,:,:) = 0.5_DP
	int_c_y( 0,:,:) = 0.5_DP
	
	return
end subroutine operator_gen

!> Initialize the data that will be used in the simulation
subroutine data_init
	use data
	use parameters
	use geometry
	use scratch
	use indices
	implicit none
	
	! Allocate the velocity in the x-direction, stored on X-faces
	allocate(U(iming:imaxg,jming:jmaxg))
	allocate(Uold(iming:imaxg,jming:jmaxg))
	! Allocate the velocity in the y-direction, stored on Y-faces
	allocate(V(iming:imaxg,jming:jmaxg))
	allocate(Vold(iming:imaxg,jming:jmaxg))
	! Allocate the pressure, stored in the cell center
	allocate(P(iming:imaxg,jming:jmaxg))
	
	! Scratch memory for fluxes
	allocate(Fx(iming:imaxg,jming:jmaxg))
	allocate(Fy(iming:imaxg,jming:jmaxg))
	! Other scratch memory arrays
	allocate(tmp1(iming:imaxg,jming:jmaxg))
	allocate(tmp2(iming:imaxg,jming:jmaxg))
	
	! Since this is lid driven cavity, initialize all to zero
	! except for the lid velocity.
	! Will use convention that cell (1,1) is in the bottom left.
	U = 0.0_DP
	V = 0.0_DP
	Uold = 0.0_DP
	Vold = 0.0_DP
	P = 1.0_DP
	
	U(:,jmax+1:jmaxg) = Ulid

	return
end subroutine data_init

!> Advance the simulation through time
subroutine flow_solve
	use parameters
	use time_info
	use geometry
	implicit none
	integer :: write_int

	! Initialize the different simulation models
	call velocity_init
	call pressure_init
	call monitor_init
	call dump_init
	
	! Start simulation at T = 0.0, with dt = init_dt
	time = 0.0_DP
	write_int = 1
	dt = init_cCFL*min(dx,dy)/Ulid
	dt = min(dt,init_vCFL*dx*dx*rho/(2.0_DP*mu))
	
	! Print initial monitor and write initial data to file
	call monitor_write
	call dump_data
	
	! Advance the simulation through time
	titer = 0
	do while(time.lt.tfin)
		! Update iteration counter
		titer = titer + 1
	  ! Advance U from n to n+1/2
		call velocity_solve
		! Solve for pressure to cause solenoidal velocity field
		call pressure_solve
		! Correct velocity to the solenoidal field
		call velocity_correct
		! Advance the simulation time
		time = time+dt
		
		! Monitor the simulation
		call monitor_write
		! Write out files for visualizing the simulation
		if(floor(time/write_freq) .ge. write_int) then
			write_int = write_int+1
			call dump_data
		end if
	end do

	return
end subroutine flow_solve
