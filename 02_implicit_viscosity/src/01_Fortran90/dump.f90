! ================================================= !
! These routines handle the dumping of the data
! for visualization. Since the velocity is
! interpolated to the cell centers, this cannot
! be used for restart. These current simulations
! are cheap enough that restarting should not
! be needed, so this will be added later.
! ================================================= !

module dump
	use precision

	! Main file list ID number
	integer :: filelist

	! Interpolating coefficients for mid-plane velocity
	integer :: imid
	real(DP) :: wx1, wx2

end module dump

! Initialize variables needed to dump data.
! Includes coefficients for midplane velocity
! interpolation.
subroutine dump_init
	use dump
	use indices
	use geometry
	implicit none
	
	integer :: i

	! Create the file that will hold the individual file names
	filelist=10
	open(unit=filelist,file="dlist.txt",status="REPLACE")
	close(filelist)
	
	! Create interpolating coefficient for mid-plane
	! Static mesh, so only need to do at beginning
	do i = imin, imax
		if(x(i).le.0.0_DP .and. x(i+1).gt.0.0_DP) then
			imid = i;
			wx1 = (0.0_DP-x(i))/(x(i+1)-x(i))
			wx2 = 1.0_DP-wx1
		end if
	end do

	return
end subroutine dump_init

! Routine to add data name to dlist.txt,
! Write the data file with formatting to allow
! use of pm3d in gnuplot, and the midplane velocity
subroutine dump_data
	use data
	use time_info
	use parameters
	use indices
	use dump
	use operators
	use geometry
	implicit none
	
	character(len=80) :: fname, buffer
	integer :: dwrite, mpwrite
	integer :: i,j
	
	! Write the file name tagged with the time
  write(buffer,'(ES12.3)') time
	fname = trim(fstart)//'_'//trim(adjustl(buffer))
	
	! Write this file name into a list of data files
	open(unit=filelist,file="dlist.txt",access="APPEND")
	write(filelist,"(A)") trim(adjustl(fname))
	close(filelist)
	
	! Open the time stamped data file and write the data to it
	dwrite = filelist+1
	open(unit=dwrite,file=trim(adjustl(fname)),status="REPLACE")
	! Write the xm, ym, cell centered U,V,and Pressure
	do j = jmin, jmax
			do i = imin, imax
				write(dwrite,"(5(ES14.6,1x))") xm(i), ym(j), &
				sum(int_x_c(:,i,j)*U(i+intxc_m:i+intxc_p,j)), &
				sum(int_y_c(:,i,j)*V(i,j+intyc_m:j+intyc_p)), &
				P(i,j)
			end do
		write(dwrite,"(A)") " "
	end do
	
	! Close the file
	close(dwrite)	
	
	! Dump the mid plane velocity non-dimensionalized by Ulid
	fname = 'U_midplane'
	mpwrite = filelist+2
	open(unit=mpwrite,file=trim(adjustl(fname)),status="REPLACE")
	do j = jmin,jmax
		write(mpwrite,"(2(ES14.6,1x))") ym(j),&
		(U(imid,j)*wx1+U(imid,j)*wx2)/Ulid 
	end do
	
	close(mpwrite)

	
	return
end subroutine dump_data
