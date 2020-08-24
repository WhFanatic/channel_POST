

subroutine mari2new
use global
implicit none

complex prhsf(0:ndxhnpr-1, 0:ndy, 0:ndznpc-1)
complex rhsp(0:ndy)
complex fpre(0:ndy)
complex u_temp(0:ndy), v_temp(0:ndy), w_temp(0:ndy)
complex &
    rhs2x (0:ndxhnpr-1, 0:ndy, 0:ndznpc-1), &
    rhs2y (0:ndxhnpr-1, 0:ndy, 0:ndznpc-1), &
    rhs2z (0:ndxhnpr-1, 0:ndy, 0:ndznpc-1)
integer i, j, k, ii, kk


u_temp = u(0, :, 0)
v_temp = v(0, :, 0)
w_temp = w(0, :, 0)

!calculate right hand side
do k = 0, ndznpc-1
do i = 0, ndxhnpr - 1
    if ( (impi_myid .eq. 0) .and. (i .eq. 0) .and. (k .eq. 0) ) then ! 00 mode
        v(0, :, 0) = 0.0
    else
        call ri2new(i, k, fpre)
        prhsf(i, :, k) = fpre
        call penta(ndy+1, ztp(:,i,k), gmp(:,i,k), afp(:,i,k), btp(:,i,k), qtp(:,i,k), fpre)
        p(i, :, k) = fpre
    endif
enddo
enddo

if ( (impi_myid .ne. 0) .or. (i .ne. 0) .or. (k .ne. 0) ) then ! not 00 mode
    ! updata velocity
    do k = 0, ndznpc-1
        kk = - ( coord(1) * ndznpc + k )
        if (kk .lt. -ndzh) kk = kk + ndz
        do i = 0, ndxhnpr-1
            ii = - ( coord(0) * ndxhnpr + i )
            
            call diffy( p(i,:,k), rhsp )
            
            do j = 0, ndy
                rhs2x(i, j, k) = u(i, j, k) - (0., 1.) * arf * ii * p(i, j, k) * dt
                rhs2y(i, j, k) = v(i, j, k) - rhsp(j) * dt
                rhs2z(i, j, k) = w(i, j, k) - (0., 1.) * bat * kk * p(i, j, k) * dt
            enddo
            
            do j = 0, ndy
                u(i, j, k) = rhs2x(i, j, k)
                v(i, j, k) = rhs2y(i, j, k)
                w(i, j, k) = rhs2z(i, j, k)
            enddo
        
        enddo
    enddo
endif

!if (ifmfu == 1) then
!    do k = 0, ndznpc-1
!        kk = - (coord(1) * ndznpc + k)
!        if (kk .lt. -ndzh) kk = kk + ndz
!        do i = 0, ndxhnpr-1
!            ii = - (coord(0) * ndxhnpr + i)
!            if ( (ii/=0) .and. (kk==0) ) then
!                u(i, :, k) = 0.0
!                v(i, :, k) = 0.0
!            endif
!        enddo
!    enddo
!endif

end subroutine




subroutine ri2new(i, k, fpre)
use global
implicit none

complex rhsv(0:ndy)
complex fpre(0:ndy), temp(0:ndy)

integer i, j, k, kk, ii

kk = - ( coord(1) * ndznpc + k )
if (kk .lt. -ndzh) kk = kk + ndz
ii = - ( coord(0) * ndxhnpr + i )

call diffy( v(i,:,k), rhsv )

do j = 1, ndy-1
    temp(j) = ( (0., 1.) * (bat * kk * w(i,j,k) + arf * ii * u(i,j,k)) + rhsv(j) ) / dt
enddo

do j = 1, ndy-1
    fpre(j) = secdera(j) * temp(j-1) + secderb(j) * temp(j) + secderc(j) * temp(j+1)
enddo

fpre(0) = 3. * pb0(i,k,1) - 3. * pb1(i,k,1) + pb2(i,k,1)
fpre(ndy)=3. * pb0(i,k,2) - 3. * pb1(i,k,2) + pb2(i,k,2)

fpre(0) = fpre(0) - cofm1(i,k) * fpre(1)
fpre(ndy)=fpre(ndy)-cofm2(i,k) * fpre(ndy-1)

end subroutine