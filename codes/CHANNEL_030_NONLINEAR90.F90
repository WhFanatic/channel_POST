
! solve nonlinear term
subroutine mari1
use global
implicit none

! non-linear term
complex &
    rhs1x (0:ndxhnpr-1, 0:ndy, 0:ndznpc-1), &
    rhs1y (0:ndxhnpr-1, 0:ndy, 0:ndznpc-1), &
    rhs1z (0:ndxhnpr-1, 0:ndy, 0:ndznpc-1)

integer i, j, k, kk, ii

u = 3. * u0 - 1.5 * u1 + u2 / 3. + dt * (3. * fx - 3. * fx1 + fx2)
v = 3. * v0 - 1.5 * v1 + v2 / 3. + dt * (3. * fy - 3. * fy1 + fy2)
w = 3. * w0 - 1.5 * w1 + w2 / 3. + dt * (3. * fz - 3. * fz1 + fz2)

if (ifmfu == 1) then
    do k = 0, ndznpc-1
        kk = - (coord(1) * ndznpc + k)
        if (kk .lt. -ndzh) kk = kk + ndz
        do i = 0, ndxhnpr-1
            ii = - (coord(0) * ndxhnpr + i)
            if ( (ii/=0) .and. (kk==0) ) then
                u(i, :, k) = 0.0
                v(i, :, k) = 0.0
            endif
        enddo
    enddo
endif

end


! solve nonlinear term for scalars
subroutine mari1_scalar
use global
implicit none	
h = 3. * h0 - 1.5 * h1 + h2 / 3. + dt * (3. * fh - 3. * fh1 + fh2)
end
