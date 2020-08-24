

subroutine premodifyp
use global
implicit none

integer i, j, k
complex rhsp(0:ndy), fprep(0:ndy), fprem(0:ndy)

do k = 0, ndznpc-1
do i = 0, ndxhnpr-1
    if ( (i.eq.0) .and. (k.eq.0) .and. (impi_myid .eq. 0) ) then
        pplus(i, :, k) = 0.5 * ( y + 1.0 )
        pminus(i, :, k) = -0.5 * ( y - 1.0 )
    else
        call ri2pmmody(i, k, fprep, fprem)
        call penta(ndy+1, ztp1(:,i,k), gmp1(:,i,k), afp1(:,i,k), btp1(:,i,k), qtp1(:,i,k), fprep)
        call penta(ndy+1, ztp1(:,i,k), gmp1(:,i,k), afp1(:,i,k), btp1(:,i,k), qtp1(:,i,k), fprem)
        pplus(i, :, k) = fprep
        pminus(i, :, k) = fprem
    endif
enddo
enddo

end subroutine
    

subroutine ri2pmmody(i, k, fprep, fprem)
use global
implicit none

integer i, j, k, kk, ii
complex rhsv(0:ndy)
complex fprep(0:ndy), fprem(0:ndy), temp(0:ndy)

kk = -( coord(1) * ndznpc + k )
if (kk .lt. -ndzh) kk = kk + ndz
ii = -( coord(0) * ndxhnpr + i )

do j = 1, ndy-1
    temp(j) = 0.0
enddo

do j = 1, ndy-1
    fprep(j) = secdera(j) * temp(j-1) + secderb(j) * temp(j) + secderc(j) * temp(j+1)
    fprem(j) = secdera(j) * temp(j-1) + secderb(j) * temp(j) + secderc(j) * temp(j+1)
enddo

fprep(0) = 0.0
fprem(0) = 1.0
fprep(ndy) = 1.0
fprem(ndy) = 0.0

end subroutine
