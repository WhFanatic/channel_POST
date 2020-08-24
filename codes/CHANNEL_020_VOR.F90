
! calculate vorticity
subroutine vor
use global
implicit none

integer i, j, k, ii, kk
complex rhsu(0:ndy), rhsw(0:ndy)

do k = 0, ndznpc-1
    kk = - ( coord(1) * ndznpc + k )
    if (kk .lt. -ndzh) kk = kk + ndz
    do i = 0, ndxhnpr-1
        ii = - ( coord(0) * ndxhnpr + i )
        
        call diffy( u(i,:,k), rhsu )
        call diffy( w(i,:,k), rhsw )
        
        fx(i, :, k) = rhsw - (0., 1.) * bat * kk * v(i,:,k)
        fy(i, :, k) = (0., 1.) * ( bat * kk * u(i,:,k) - arf * ii * w(i,:,k) )
        fz(i, :, k) = (0., 1.) * arf * ii * v(i,:,k) - rhsu
    enddo	
enddo

!! add average rotation here if computing rotating channel
!if (impi_myid .eq. 0) then
!    fz(0, :, 0) = fz(0, :, 0) + Ro
!endif

end


! calculate gradient of scalars
subroutine cdh
use global
implicit none

! gradient of scalar field
complex dh_x, dh_y, dh_z
common /scalar_gradient/ &
    dh_x (0:ndxhnpr-1, 0:ndy, 0:ndznpc-1, ns), &
    dh_y (0:ndxhnpr-1, 0:ndy, 0:ndznpc-1, ns), &
    dh_z (0:ndxhnpr-1, 0:ndy, 0:ndznpc-1, ns)


!-------------------------------------------------------------------------------
!    local variable
!-------------------------------------------------------------------------------
complex th(0:ndxhnpr-1,0:ndy+1)
integer ii , k , kk , jj , i , is

do is = 1 , ns
    do k = 0 , ndznpc - 1		   
        kk = - ( coord(1) * ndznpc + k)
        if(kk .lt. -ndzh) kk = kk + ndz

	    th(: , ndy) = (0. , 0.)
	    th(: , ndy+1) = (0. , 0.)

	    do jj = ndy,2,-1
	        th(: , jj-1) = th(: , jj+1) + 2. * jj * h(: , jj , k , is)
	    enddo

        th(: , 0) = .5 * th(: , 2) + h(: , 1 , k , is)

	    dh_y(: , 0:ndy , k , is) = th(: , 0:ndy)
	    dh_z(: , : , k , is) = (0. , 1.) * bat  * kk  * h(: , : , k , is)

        do i = 0 , ndxhnpr-1
	        ii = -( coord(0) * ndxhnpr + i) 
		    dh_x(i , : , k , is) = (0. , 1.) * ii * arf * h(i , : , k , is)
	    enddo
    enddo
enddo	
   														 		
end
