
! monitor flow rate on screen
subroutine wrfr
use global
implicit none

integer j
complex f(0:ndy)

if (impi_myid .eq. 0) then
    call integy( u(0,:,0), f )
    fr = real( f(ndy) )
endif

end
