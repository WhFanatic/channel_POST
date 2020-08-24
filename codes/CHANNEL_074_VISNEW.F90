

subroutine mari3new
use global
implicit none

integer i, j, k, ii, kk
complex fvisu(0:ndy), fvisv(0:ndy), fvisw(0:ndy)


do k = 0, ndznpc-1
do i = 0, ndxhnpr-1
    call ri3new(i, k, fvisu, fvisv, fvisw)
    call penta(ndy+1, ztv(:,i,k), gmv(:,i,k), afv(:,i,k), btv(:,i,k), qtv(:,i,k), fvisu)
    call penta(ndy+1, ztv(:,i,k), gmv(:,i,k), afv(:,i,k), btv(:,i,k), qtv(:,i,k), fvisv)
    call penta(ndy+1, ztv(:,i,k), gmv(:,i,k), afv(:,i,k), btv(:,i,k), qtv(:,i,k), fvisw)
    u(i, :, k) = fvisu
    v(i, :, k) = fvisv
    w(i, :, k) = fvisw
enddo
enddo

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



subroutine ri3new(i, k, fvisu, fvisv, fvisw)
use global
implicit none

complex fvisu(0:ndy), fvisv(0:ndy), fvisw(0:ndy)
integer i, j, k

do j = 1, ndy-1
    fvisu(j) = secdera(j) * u(i,j-1,k) + secderb(j) * u(i,j,k) + secderc(j) * u(i,j+1,k)
    fvisv(j) = secdera(j) * v(i,j-1,k) + secderb(j) * v(i,j,k) + secderc(j) * v(i,j+1,k)
    fvisw(j) = secdera(j) * w(i,j-1,k) + secderb(j) * w(i,j,k) + secderc(j) * w(i,j+1,k)
enddo

fvisu = fvisu * (-re/dt)
fvisv = fvisv * (-re/dt)
fvisw = fvisw * (-re/dt)

! no-slip boundary condition
fvisu(0) = 0.0
fvisv(0) = 0.0
fvisw(0) = 0.0
fvisu(ndy) = 0.0
fvisv(ndy) = 0.0
fvisw(ndy) = 0.0

end subroutine  	          