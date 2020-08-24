
! compute tau_w (or equivalently, pressure gradient)
subroutine utp
use global
implicit none

integer j
complex rhsu(0:ndy)
real utu, utl
common /wallstress/ utu, utl

if (impi_myid .eq. 0) then
    call diffy( u(0,:,0), rhsu )
    utl = abs(rhsu(0)) / re
    utu = abs(rhsu(ndy)) / re
    dp = - 0.5  * (utu + utl) ! 2H * dp/dx = - ( tau_up + tau_low )
endif
 
end
