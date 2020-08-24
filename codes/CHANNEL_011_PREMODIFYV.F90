
subroutine premodifyv
use global
implicit none

!complex umm, vmm, wmm, ump, vmp, wmp
!common /velomodify_full/ &
!    umm(0:ndxhnpr-1, 0:ndy, 0:ndznpc-1), &
!    vmm(0:ndxhnpr-1, 0:ndy, 0:ndznpc-1), &
!    wmm(0:ndxhnpr-1, 0:ndy, 0:ndznpc-1), &
!    ump(0:ndxhnpr-1, 0:ndy, 0:ndznpc-1), &
!    vmp(0:ndxhnpr-1, 0:ndy, 0:ndznpc-1), &
!    wmp(0:ndxhnpr-1, 0:ndy, 0:ndznpc-1)

integer i, j, k
complex fvisup(0:ndy), fvisvp(0:ndy), fviswp(0:ndy)
complex fvisum(0:ndy), fvisvm(0:ndy), fviswm(0:ndy)

do k = 0, ndznpc-1
do i = 0, ndxhnpr-1

    call ri3pmmody(i, k, fvisup, fvisvp, fviswp, fvisum, fvisvm, fviswm)

    call penta(ndy+1, ztv(:,i,k), gmv(:,i,k), afv(:,i,k), btv(:,i,k), qtv(:,i,k), fvisup)
    call penta(ndy+1, ztv(:,i,k), gmv(:,i,k), afv(:,i,k), btv(:,i,k), qtv(:,i,k), fvisvp)
    call penta(ndy+1, ztv(:,i,k), gmv(:,i,k), afv(:,i,k), btv(:,i,k), qtv(:,i,k), fviswp)

    call penta(ndy+1, ztv(:,i,k), gmv(:,i,k), afv(:,i,k), btv(:,i,k), qtv(:,i,k), fvisum)
    call penta(ndy+1, ztv(:,i,k), gmv(:,i,k), afv(:,i,k), btv(:,i,k), qtv(:,i,k), fvisvm)
    call penta(ndy+1, ztv(:,i,k), gmv(:,i,k), afv(:,i,k), btv(:,i,k), qtv(:,i,k), fviswm)

    ump(i,:,k) = fvisup
    vmp(i,:,k) = fvisvp
    wmp(i,:,k) = fviswp

    umm(i,:,k) = fvisum
    vmm(i,:,k) = fvisvm
    wmm(i,:,k) = fviswm

enddo
enddo

end subroutine


subroutine ri3pmmody(i, k, fvisup, fvisvp, fviswp, fvisum, fvisvm, fviswm)
use global
implicit none

!complex pplus, pminus
!common /ppm/ pplus(0:ndxhnpr-1, 0:ndy, 0:ndznpc-1), pminus(0:ndxhnpr-1, 0:ndy, 0:ndznpc-1)

integer i, j, k, kk, ii
complex fvisum(0:ndy), fvisvm(0:ndy), fviswm(0:ndy)
complex fvisup(0:ndy), fvisvp(0:ndy), fviswp(0:ndy)

complex tempu(0:ndy), tempv(0:ndy), tempw(0:ndy)
complex dpplus(0:ndy), dpminus(0:ndy)


kk = -( coord(1) * ndznpc + k )
if ( kk .lt. -ndzh ) kk = kk + ndz
ii = -( coord(0) * ndxhnpr + i )


call diffy( pplus(i,:,k), dpplus )
do j = 1, ndy-1
    tempu(j) = (0., 1.) * (arf * ii * pplus(i,j,k))
    tempv(j) = dpplus(j)
    tempw(j) = (0., 1.) * (bat * kk * pplus(i,j,k))
enddo
do j = 1, ndy-1
    fvisup(j) = secdera(j) * tempu(j-1) + secderb(j) * tempu(j) + secderc(j) * tempu(j+1)
    fvisvp(j) = secdera(j) * tempv(j-1) + secderb(j) * tempv(j) + secderc(j) * tempv(j+1)
    fviswp(j) = secdera(j) * tempw(j-1) + secderb(j) * tempw(j) + secderc(j) * tempw(j+1)
enddo
fvisup = fvisup * re
fvisvp = fvisvp * re
fviswp = fviswp * re
fvisup(0) = 0.0
fvisvp(0) = 0.0
fviswp(0) = 0.0
fvisup(ndy) = 0.0
fvisvp(ndy) = 0.0
fviswp(ndy) = 0.0


call diffy( pminus(i,:,k), dpminus )
do j = 1, ndy-1
    tempu(j) = (0., 1.) * (arf * ii * pminus(i,j,k))
    tempv(j) = dpminus(j)
    tempw(j) = (0., 1.) * (bat * kk * pminus(i,j,k))
enddo
do j = 1, ndy-1
    fvisum(j) = secdera(j) * tempu(j-1) + secderb(j) * tempu(j) + secderc(j) * tempu(j+1)
    fvisvm(j) = secdera(j) * tempv(j-1) + secderb(j) * tempv(j) + secderc(j) * tempv(j+1)
    fviswm(j) = secdera(j) * tempw(j-1) + secderb(j) * tempw(j) + secderc(j) * tempw(j+1)
enddo
fvisum = fvisum * re
fvisvm = fvisvm * re
fviswm = fviswm * re
fvisum(0) = 0.0
fvisvm(0) = 0.0
fviswm(0) = 0.0
fvisum(ndy) = 0.0
fvisvm(ndy) = 0.0
fviswm(ndy) = 0.0  
  	          
endsubroutine