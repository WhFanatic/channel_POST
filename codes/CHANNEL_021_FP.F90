
    

! caculate non-linear terms
! subroutine compute_convol not used, to avoid redundant operations of fft & ifft
subroutine fp
use global
implicit none

complex dh_x, dh_y, dh_z ! gradient of scalar field
common /scalar_gradient/ &
    dh_x (0:ndxhnpr-1, 0:ndy, 0:ndznpc-1, ns), &
    dh_y (0:ndxhnpr-1, 0:ndy, 0:ndznpc-1, ns), &
    dh_z (0:ndxhnpr-1, 0:ndy, 0:ndznpc-1, ns)

! local variables
real dh(ns) ! dt/dx--const

complex tempz1(0:ndxhnpr-1, 0:ndynpc, 0:ndz2-1)
complex tempz2(0:ndxhnpr-1, 0:ndynpc, 0:ndz2-1)
complex tempz3(0:ndxhnpr-1, 0:ndynpc, 0:ndz2-1)
complex tempz4(0:ndxhnpr-1, 0:ndynpc, 0:ndz2-1)
complex tempz5(0:ndxhnpr-1, 0:ndynpc, 0:ndz2-1)
complex tempz6(0:ndxhnpr-1, 0:ndynpc, 0:ndz2-1)
complex tempz7(0:ndxhnpr-1, 0:ndynpc, 0:ndz2-1)
complex tempz8(0:ndxhnpr-1, 0:ndynpc, 0:ndz2-1)
complex tempz9(0:ndxhnpr-1, 0:ndynpc, 0:ndz2-1)

complex tempx1(0:ndx2-1, 0:ndynpc, 0:ndz2npr-1)
complex tempx2(0:ndx2-1, 0:ndynpc, 0:ndz2npr-1)
complex tempx3(0:ndx2-1, 0:ndynpc, 0:ndz2npr-1)
complex tempx4(0:ndx2-1, 0:ndynpc, 0:ndz2npr-1)
complex tempx5(0:ndx2-1, 0:ndynpc, 0:ndz2npr-1)
complex tempx6(0:ndx2-1, 0:ndynpc, 0:ndz2npr-1)
complex tempx7(0:ndx2-1, 0:ndynpc, 0:ndz2npr-1)
complex tempx8(0:ndx2-1, 0:ndynpc, 0:ndz2npr-1)
complex tempx9(0:ndx2-1, 0:ndynpc, 0:ndz2npr-1)

complex tempz(0:ndz2-1), tempz0(0:ndz2-1)
complex tempx(0:ndx2-1), tempx0(0:ndx2-1)

integer i, j, k, is


tempx1 = 0.
tempx2 = 0.
tempx3 = 0.
tempx4 = 0.
tempx5 = 0.
tempx6 = 0.
tempx7 = 0.
tempx8 = 0.
tempx9 = 0.

tempz1 = 0.
tempz2 = 0.
tempz3 = 0.
tempz4 = 0.
tempz5 = 0.
tempz6 = 0.
tempz7 = 0.
tempz8 = 0.
tempz9 = 0.

!!! perform ifft in z direction for u, v, w, fx, fy, fz
! transpose division of field to ensure integrity in z direction
call transposeytoz( u , tempz1(0:ndxhnpr-1, 0:ndynpc, 0:ndz-1) )
call transposeytoz( v , tempz2(0:ndxhnpr-1, 0:ndynpc, 0:ndz-1) )
call transposeytoz( w , tempz3(0:ndxhnpr-1, 0:ndynpc, 0:ndz-1) )
call transposeytoz( fx, tempz4(0:ndxhnpr-1, 0:ndynpc, 0:ndz-1) )
call transposeytoz( fy, tempz5(0:ndxhnpr-1, 0:ndynpc, 0:ndz-1) )
call transposeytoz( fz, tempz6(0:ndxhnpr-1, 0:ndynpc, 0:ndz-1) )
! applying 3/2 rule for de-aliasing
tempz1(:, :, ndz+1:ndz2-1) = tempz1(:, :, ndzh+1:ndz-1)
tempz2(:, :, ndz+1:ndz2-1) = tempz2(:, :, ndzh+1:ndz-1)
tempz3(:, :, ndz+1:ndz2-1) = tempz3(:, :, ndzh+1:ndz-1)
tempz4(:, :, ndz+1:ndz2-1) = tempz4(:, :, ndzh+1:ndz-1)
tempz5(:, :, ndz+1:ndz2-1) = tempz5(:, :, ndzh+1:ndz-1)
tempz6(:, :, ndz+1:ndz2-1) = tempz6(:, :, ndzh+1:ndz-1)

tempz1(:, :, ndzh:ndz-1) = 0.
tempz2(:, :, ndzh:ndz-1) = 0.
tempz3(:, :, ndzh:ndz-1) = 0.
tempz4(:, :, ndzh:ndz-1) = 0.
tempz5(:, :, ndzh:ndz-1) = 0.
tempz6(:, :, ndzh:ndz-1) = 0.
! do ifft
do j = 0, ndynpc
do i = 0, ndxhnpr-1
    tempz = tempz1(i, j, :)
    call fftw_f77(my_plan_z_f, 1, tempz, 1, ndz2, tempz0, 0, 0)
    tempz1(i, j, :) = tempz
enddo
enddo

do j = 0, ndynpc
do i = 0, ndxhnpr-1
    tempz = tempz2(i, j, :)
    call fftw_f77(my_plan_z_f, 1, tempz, 1, ndz2, tempz0, 0, 0)
    tempz2(i, j, :) = tempz
enddo
enddo

do j = 0, ndynpc
do i = 0, ndxhnpr-1
    tempz = tempz3(i, j, :)
    call fftw_f77(my_plan_z_f, 1, tempz, 1, ndz2, tempz0, 0, 0)
    tempz3(i, j, :) = tempz
enddo
enddo

do j = 0, ndynpc
do i = 0, ndxhnpr-1
    tempz = tempz4(i, j, :)
    call fftw_f77(my_plan_z_f, 1, tempz, 1, ndz2, tempz0, 0, 0)
    tempz4(i, j, :) = tempz
enddo
enddo

do j = 0, ndynpc
do i = 0, ndxhnpr-1
    tempz = tempz5(i, j, :)
    call fftw_f77(my_plan_z_f, 1, tempz, 1, ndz2, tempz0, 0, 0)
    tempz5(i, j, :) = tempz
enddo
enddo

do j = 0, ndynpc
do i = 0, ndxhnpr-1
    tempz = tempz6(i, j, :)
    call fftw_f77(my_plan_z_f, 1, tempz, 1, ndz2, tempz0, 0, 0)
    tempz6(i, j, :) = tempz
enddo
enddo

!!! perform ifft in x direction for u, v, w, fx, fy, fz
! transpose division of field to ensure integrity in x direction
call transposeztox( tempz1, tempx1(0:ndxh-1, 0:ndynpc, 0:ndz2npr-1) )
call transposeztox( tempz2, tempx2(0:ndxh-1, 0:ndynpc, 0:ndz2npr-1) )
call transposeztox( tempz3, tempx3(0:ndxh-1, 0:ndynpc, 0:ndz2npr-1) )
call transposeztox( tempz4, tempx4(0:ndxh-1, 0:ndynpc, 0:ndz2npr-1) )
call transposeztox( tempz5, tempx5(0:ndxh-1, 0:ndynpc, 0:ndz2npr-1) )
call transposeztox( tempz6, tempx6(0:ndxh-1, 0:ndynpc, 0:ndz2npr-1) )
! applying 3/2 rule for de-aliasing
tempx1(ndx2-1:ndx+1:-1, :, :) = conjg( tempx1(1:ndxh-1, :, :) )
tempx2(ndx2-1:ndx+1:-1, :, :) = conjg( tempx2(1:ndxh-1, :, :) )
tempx3(ndx2-1:ndx+1:-1, :, :) = conjg( tempx3(1:ndxh-1, :, :) )
tempx4(ndx2-1:ndx+1:-1, :, :) = conjg( tempx4(1:ndxh-1, :, :) )
tempx5(ndx2-1:ndx+1:-1, :, :) = conjg( tempx5(1:ndxh-1, :, :) )
tempx6(ndx2-1:ndx+1:-1, :, :) = conjg( tempx6(1:ndxh-1, :, :) )
! do ifft
do k = 0, ndz2npr-1
do j = 0, ndynpc
    tempx = tempx1(:, j, k)
    call fftw_f77(my_plan_x_f, 1, tempx, 1,ndx2, tempx0, 0, 0)
    tempx1(:, j, k) = tempx
enddo
enddo

do k = 0, ndz2npr-1
do j = 0, ndynpc
    tempx = tempx2(:, j, k)
    call fftw_f77(my_plan_x_f, 1, tempx, 1, ndx2, tempx0, 0, 0)
    tempx2(:, j, k) = tempx
enddo
enddo

do k = 0, ndz2npr-1
do j = 0, ndynpc
    tempx = tempx3(:, j, k)
    call fftw_f77(my_plan_x_f, 1, tempx, 1, ndx2, tempx0, 0, 0)
    tempx3(:, j, k) = tempx
enddo
enddo

do k = 0, ndz2npr-1
do j = 0, ndynpc
    tempx = tempx4(:, j, k)
    call fftw_f77(my_plan_x_f, 1, tempx, 1, ndx2, tempx0, 0, 0)
    tempx4(:, j, k) = tempx
enddo
enddo

do k = 0, ndz2npr-1
do j = 0, ndynpc
    tempx = tempx5(:, j, k)
    call fftw_f77(my_plan_x_f, 1, tempx, 1, ndx2, tempx0, 0, 0)
    tempx5(:, j, k) = tempx
enddo
enddo

do k = 0, ndz2npr-1
do j = 0, ndynpc
    tempx = tempx6(:, j, k)
    call fftw_f77(my_plan_x_f, 1, tempx, 1, ndx2, tempx0, 0, 0)
    tempx6(:, j, k) = tempx
enddo
enddo

!!! get non-linear terms in physical space
tempx7 = real(tempx2) * real(tempx6) - real(tempx3) * real(tempx5)
tempx8 = real(tempx3) * real(tempx4) - real(tempx1) * real(tempx6)
tempx9 = real(tempx1) * real(tempx5) - real(tempx2) * real(tempx4)

!!! perform fft in x direction for fx, fy, fz
do k = 0, ndz2npr-1
do j = 0, ndynpc
    tempx = tempx7(:, j, k)
    call fftw_f77(my_plan_x_b, 1, tempx, 1, ndx2, tempx0, 0, 0)
    tempx7(0:ndxh-1, j, k) = tempx(0:ndxh-1) / ndx2
enddo
enddo

do k = 0, ndz2npr-1
do j = 0, ndynpc
    tempx = tempx8(:, j, k)
    call fftw_f77(my_plan_x_b, 1, tempx, 1, ndx2, tempx0, 0, 0)
    tempx8(0:ndxh-1, j, k) = tempx(0:ndxh-1) / ndx2
enddo
enddo

do k = 0, ndz2npr-1
do j = 0, ndynpc
    tempx = tempx9(:, j, k)
    call fftw_f77(my_plan_x_b, 1, tempx, 1, ndx2, tempx0, 0, 0)
    tempx9(0:ndxh-1, j, k) = tempx(0:ndxh-1) / ndx2
enddo
enddo

!!! perform fft in z direction for fx, fy, fz
call transposextoz( tempx7(0:ndxh-1, 0:ndynpc, 0:ndz2npr-1), tempz7(0:ndxhnpr-1, 0:ndynpc, 0:ndz2-1) )
call transposextoz( tempx8(0:ndxh-1, 0:ndynpc, 0:ndz2npr-1), tempz8(0:ndxhnpr-1, 0:ndynpc, 0:ndz2-1) )
call transposextoz( tempx9(0:ndxh-1, 0:ndynpc, 0:ndz2npr-1), tempz9(0:ndxhnpr-1, 0:ndynpc, 0:ndz2-1) )

do j = 0, ndynpc
do i = 0, ndxhnpr-1
    tempz = tempz7(i, j, :)
    call fftw_f77(my_plan_z_b, 1, tempz, 1, ndz2, tempz0, 0, 0)
    tempz7(i, j, 0:ndzh-1) = tempz(0:ndzh-1) / ndz2
    tempz7(i, j, ndzh+1:ndz-1) = tempz(ndz+1:ndz2-1) / ndz2
    tempz7(i, j, ndzh) = 0.
enddo
enddo

do j = 0, ndynpc
do i = 0, ndxhnpr-1
    tempz = tempz8(i, j, :)
    call fftw_f77(my_plan_z_b, 1, tempz, 1, ndz2, tempz0, 0, 0)
    tempz8(i, j, 0:ndzh-1) = tempz(0:ndzh-1) / ndz2
    tempz8(i, j, ndzh+1:ndz-1) = tempz(ndz+1:ndz2-1) / ndz2
    tempz8(i, j, ndzh) = 0.
enddo
enddo

do j = 0, ndynpc
do i = 0, ndxhnpr-1
    tempz = tempz9(i, j, :)
    call fftw_f77(my_plan_z_b, 1, tempz, 1, ndz2, tempz0, 0, 0)
    tempz9(i, j, 0:ndzh-1) = tempz(0:ndzh-1) / ndz2
    tempz9(i, j, ndzh+1:ndz-1) = tempz(ndz+1:ndz2-1) / ndz2
    tempz9(i, j, ndzh) = 0.
enddo
enddo

!!! transpose back to original state
call transposeztoy( tempz7(0:ndxhnpr-1, 0:ndynpc, 0:ndz-1), fx(0:ndxhnpr-1, 0:ndy, 0:ndznpc-1) )
call transposeztoy( tempz8(0:ndxhnpr-1, 0:ndynpc, 0:ndz-1), fy(0:ndxhnpr-1, 0:ndy, 0:ndznpc-1) )
call transposeztoy( tempz9(0:ndxhnpr-1, 0:ndynpc, 0:ndz-1), fz(0:ndxhnpr-1, 0:ndy, 0:ndznpc-1) )

!!! same procedures for scalars
if (ifscl .ge. 1) then
    do is = 1, ifscl
        tempz4 = 0.
	    tempz5 = 0.
	    tempz6 = 0.
	    
	    tempx4 = 0.
	    tempx5 = 0.
	    tempx6 = 0.

	    call chebp( dh_x(0, 0, 0, is), ndxhnpr, ndy, ndznpc )
	    call chebp( dh_y(0, 0, 0, is), ndxhnpr, ndy, ndznpc )
	    call chebp( dh_z(0, 0, 0, is), ndxhnpr, ndy, ndznpc )
	    call transposeytoz( dh_x(0:ndxhnpr-1, 0:ndy, 0:ndznpc-1, is), tempz4(0:ndxhnpr-1, 0:ndynpc, 0:ndz-1) )
	    call transposeytoz( dh_y(0:ndxhnpr-1, 0:ndy, 0:ndznpc-1, is), tempz5(0:ndxhnpr-1, 0:ndynpc, 0:ndz-1) )
	    call transposeytoz( dh_z(0:ndxhnpr-1, 0:ndy, 0:ndznpc-1, is), tempz6(0:ndxhnpr-1, 0:ndynpc, 0:ndz-1) )

        tempz4(0:ndxhnpr-1, :, ndzh) = 0.
        tempz5(0:ndxhnpr-1, :, ndzh) = 0.
        tempz6(0:ndxhnpr-1, :, ndzh) = 0.

        tempz4(0:ndxhnpr-1, :, ndz+1:ndz2-1) = tempz4(0:ndxhnpr-1, :, ndzh+1:ndz-1)
        tempz5(0:ndxhnpr-1, :, ndz+1:ndz2-1) = tempz5(0:ndxhnpr-1, :, ndzh+1:ndz-1)
        tempz6(0:ndxhnpr-1, :, ndz+1:ndz2-1) = tempz6(0:ndxhnpr-1, :, ndzh+1:ndz-1)

        tempz4(0:ndxhnpr-1, :, ndzh+1:ndz-1) = 0.
        tempz5(0:ndxhnpr-1, :, ndzh+1:ndz-1) = 0.
        tempz6(0:ndxhnpr-1, :, ndzh+1:ndz-1) = 0.
	    
	    do j = 0, ndynpc
        do i = 0, ndxhnpr-1
            tempz = tempz4(i, j, :)
            call fftw_f77(my_plan_z_f, 1, tempz, 1, ndz2, tempz0, 0, 0)
            tempz4(i, j, :) = tempz
        enddo
        enddo
        
	    do j = 0, ndynpc
        do i = 0, ndxhnpr-1
            tempz = tempz5(i, j, :)
            call fftw_f77(my_plan_z_f, 1, tempz, 1, ndz2, tempz0, 0, 0)
            tempz5(i, j, :) = tempz
        enddo
        enddo
        
        do j = 0, ndynpc
        do i = 0, ndxhnpr-1
            tempz = tempz6(i, j, :)
            call fftw_f77(my_plan_z_f, 1, tempz, 1, ndz2, tempz0, 0, 0)
            tempz6(i, j, :) = tempz
        enddo
        enddo
	    
	    call transposeztox( tempz4, tempx4(0:ndxh-1, 0:ndynpc, 0:ndz2npr-1) )
        call transposeztox( tempz5, tempx5(0:ndxh-1, 0:ndynpc, 0:ndz2npr-1) )
        call transposeztox( tempz6, tempx6(0:ndxh-1, 0:ndynpc, 0:ndz2npr-1) )
	    
	    tempx4(ndx2-1:ndx+1:-1, :, :) = conjg( tempx4(1:ndxh-1, :, :) )
        tempx5(ndx2-1:ndx+1:-1, :, :) = conjg( tempx5(1:ndxh-1, :, :) )
        tempx6(ndx2-1:ndx+1:-1, :, :) = conjg( tempx6(1:ndxh-1, :, :) )

        do k = 0, ndz2npr-1
        do j = 0, ndynpc
            tempx = tempx4(:, j, k)
            call fftw_f77(my_plan_x_f, 1, tempx, 1, ndx2, tempx0, 0, 0)
            tempx4(:, j, k) = tempx
        enddo
        enddo
        
        do k = 0, ndz2npr-1
        do j = 0, ndynpc
            tempx = tempx5(:, j, k)
            call fftw_f77(my_plan_x_f, 1, tempx, 1, ndx2, tempx0, 0, 0)
            tempx5(:, j, k) = tempx
        enddo
        enddo
        
        do k = 0, ndz2npr-1
        do j = 0, ndynpc
            tempx = tempx6(:, j, k)
            call fftw_f77(my_plan_x_f, 1, tempx, 1, ndx2, tempx0, 0, 0)
            tempx6(:, j, k) = tempx
        enddo
        enddo
        
        tempx7 = - real(tempx1) * real(tempx4) - real(tempx2) * real(tempx5) - real(tempx3) * real(tempx6)   
        tempx7 = tempx7 + real(tempx1) * dh(is)

        do k = 0, ndz2npr-1
        do j = 0, ndynpc
            tempx = tempx7(:, j, k)
            call fftw_f77(my_plan_x_b, 1, tempx, 1, ndx2, tempx0, 0, 0)
            tempx4(0:ndxh-1, j, k) = tempx(0:ndxh-1) / ndx2
        enddo
        enddo
        
        call transposextoz( tempx4(0:ndxh-1, 0:ndynpc, 0:ndz2npr-1), tempz4(0:ndxhnpr-1, 0:ndynpc, 0:ndz2-1) )
        
        do j = 0, ndynpc
        do i = 0, ndxhnpr-1
            tempz = tempz4(i, j, :)
            call fftw_f77(my_plan_z_b, 1, tempz, 1, ndz2, tempz0, 0, 0)
            tempz4(i, j, 0:ndzh-1) = tempz(0:ndzh-1) / ndz2
            tempz4(i, j, ndzh+1:ndz-1) = tempz(ndz+1:ndz2-1) / ndz2
            tempz4(i, j, ndzh) = 0.
        enddo
        enddo
	    
        call transposeztoy( tempz4(0:ndxhnpr-1, 0:ndynpc, 0:ndz-1), fh(0:ndxhnpr-1, 0:ndy, 0:ndznpc-1, is) )
    
    enddo
endif
	
end

    
! bulk force
subroutine bulkforce
use global
implicit none

if (impi_myid .eq. 0) then
    if (iff .eq. 1) fx(0, :, 0) = fx(0, :, 0) - dp
    if (iff .eq. 0) fx(0, :, 0) = fx(0, :, 0) - dp0
endif

end
    
    