
! correct flow rates
subroutine frcorr
use global
implicit none

if (impi_myid .eq. 0) u(0, :, 0) = u(0,:,0) + 1.5 * ( 1 - y**2 ) * ( 1. - fr/fr0 )

end


! correct flow rates	
subroutine hrcorr
use global
implicit none

real hr
common /heatrate/ hr(ns)

integer is

if (impi_myid .eq. 0) then
    do is = 1, ifscl
        h(0, 0, 0, is) = h(0, 0, 0, is) + 95./256.* ( 1. - hr(is) / fr0 )
        h(0, 2, 0, is) = h(0, 2, 0, is) - 25./64. * ( 1. - hr(is) / fr0 )
        h(0, 4, 0, is) = h(0, 4, 0, is) + 5. /256.* ( 1. - hr(is) / fr0 )
    enddo
endif

end