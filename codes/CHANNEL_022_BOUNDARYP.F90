
! caculate the 1st term of pb
subroutine pb1condition
use global
implicit none

integer i, k, ii, kk

do k = 0, ndznpc-1
    kk = - ( coord(1) * ndznpc + k )
    if (kk .lt. -ndzh) kk = kk + ndz
	do i = 0, ndxhnpr-1
        ii = - ( coord(0) * ndxhnpr + i )
        pb0(i, k, 1) = (0., 1.) / re * ( arf * ii * fz(i, 0, k) - bat * kk * fx(i, 0, k) )
        pb0(i, k, 2) = (0., 1.) / re * ( arf * ii * fz(i,ndy,k) - bat * kk * fx(i,ndy,k) )
    enddo
enddo
    
end


! caculate the 2nd term of pb   	
subroutine pb2condition
use global
implicit none

integer i, j, k, is

pb0(:, :, 1) = pb0(:, :, 1) + fy(:, 0 , :)
pb0(:, :, 2) = pb0(:, :, 2) + fy(:,ndy, :)

!if (impi_myid .eq. 0) then
!    if (iff .eq. 1) fx(0, :, 0) = fx(0, :, 0) - dp
!    if (iff .eq. 0) fx(0, :, 0) = fx(0, :, 0) - dp0
!endif

end
