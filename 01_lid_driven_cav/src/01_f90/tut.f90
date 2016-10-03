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

module parameter

	real(DP) :: dw, dh, Ulid, rho, mu
	integer :: nx, ny ,ng

end module parameter

module geometry

	real(DP), dimension(:), allocatable :: x, y, xm, ym
	real(DP) :: dx, dy, area, iarea
	integer :: iming,imaxg,jming,jmaxg
	integer :: imin, imax, jmin, jmax

end module geometry

module operators

	! Derivatives
	real(DP), dimension(:), allocatable :: sfdx, sbdx, scdx
	real(DP), dimension(:), allocatable :: sfdy, sbdy, scdy
	integer :: sfdm, sfdp, sbdm, smdp, scdm, scdp

	! Interpolation
	real(DP), dimension(:), allocatable :: sci_xc, sci_yc, sci_cx, sci_cy
	integer :: scim_xc, scip_xc, scim_yc, scip_yc, scim_cx, scip_cx
	integer :: scim_cy, scip_cy

end module operators

module data

	real(DP) :: U, V, P

end module data

!> This just calls configuration, sets the initial data
!! and then starts the flow solver
program main

	call parameter_set
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
subroutine parameter_set
		
	! Domain width
	dw = 1.0_DP
	! Cells across the width
	nx = 100
	! Domain height
	dh = 1.0_DP
	! Cells across the height
	ny = 100
	! Lid velocity
	Ulid = 1.0_DP
	! Density of fluid
	rho = 1.0_WP
	! Dynamic viscosity of fluid
	mu = 0.01_DP
	
	! Number of ghost cells on each boundary
	ng = 2

	return
end subroutine parameter_set
	
!> Sets up the configuration of the grid (x,y locations).
!! Rectangular grid is assumed, with a staggered cell arrangement
subroutine grid_make
	use parameter
	use geometry
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
		x(i) = real(i-1,DP)*dx-0.5_WP*dw
	end do
	
	! Bottom cell faces
	do j = jming,jmaxg
		y(j) = real(j-1,DP)*dy-0.5_WP*dh
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

!> Creates the operators based on the mesh
subroutine operator_gen
	use geometry
	use operators
	implicit none
	
	integer :: i
	
	! Second order forward derivative, assume uniform mesh
	sfdm = 0
	sfdp = 2
	allocate(sfdx(sfm:sfp))
	allocate(sfdy(sfm:sfp))
	sfdx(0) = -1.5_DP
	sfdx(1) =  2.0_DP
	sfdx(2) = -0.5_DP
	sfdy = sfdx*iarea*idy
	sfdx = sfdx*iarea*idx
	
	! Second order backward derivative, assume uniform mesh
	sbdm = -2
	sbdp =  0
	allocate(sbdx(sbm:sbp))
	allocate(sbdy(sbm:sbp))
	sbdx(0)  =  1.5_DP
	sbdx(-1) = -2.0_DP
	sbdx(-2) =  0.5_DP
	sbdy = sbdx*iarea*idy
	sbdx = sbdx*iarea*idx
	
	! Second order centered derivative, assume uniform mesh
	scdm = -1
	scdp = 1
	allocate(scdx(scm:scp))
	allocate(scdy(scm:scp))
	scdx(-1) = -0.5_DP
	scdx(0)  =  0.0_DP
	scdx(1)  =  0.5_DP
	scdy = scdx*iarea*idy
	scdx = scdx*iarea*idx
	
	! Interpolation from X-face to cell center (First order)
	scim_xc = 0
	scip_xc = 1
	allocate(sci_xc(scim:scip))
	scix(0) = 0.5_DP
	scix(1) = 0.5_DP

	! Interpolation from Y-face to cell center (First order)
	scim_yc = 0
	scip_yc = 1
	allocate(sci_yc(scim_yc:scip_yc))
	sci_yc(0) = 0.5_DP
	sci_yc(1) = 0.5_DP
	
	! Interpolation from cell center to X-face (First order)
	scim_cx = -1
	scip_cx =  0
	allocate(sci_cx(scim_cx:scip_cx))
	sci_cx(-1) = 0.5_DP
	sci_cx( 0) = 0.5_DP
	
	! Interpolation from cell center to X-face (First order)
	scim_cy = -1
	scip_cy =  0
	allocate(sci_cy(scim_cy:scip_cy))
	sci_cy(-1) = 0.5_DP
	sci_cy( 0) = 0.5_DP
	
	return
end subroutine operator_gen

!> Initialize the data that will be used in the simulation
subroutine data_init
	use data
	use geometry
	implicit none
	
	! Allocate the velocity in the x-direction, stored on X-faces
	allocate(U(iming:imaxg,jming:jmaxg))
	! Allocate the velocity in the y-direction, stored on Y-faces
	allocate(V(iming:imaxg,jming:jmaxg))
	! Allocate the pressure, stored in the cell center
	allocate(P(iming:imaxg,jming:jmaxg))
	
	! Since this is lid driven cavity, initialize all to zero
	! except for the lid velocity.
	! Will use convention that cell (1,1) is in the bottom left.
	U = 0.0_DP
	V = 0.0_DP
	P = 0.0_DP
	
	U(:,jmax+1:jmaxg) = Ulid

	return
end subroutine data_init

!> Advance the simulation through time
subroutine flow_solve

	return
end subroutine flow_solve
