
! displace old & new field after evolution
subroutine bdisplace
use global
implicit none

u2 = u1
v2 = v1
w2 = w1
p2 = p1
fx2 = fx1
fy2 = fy1
fz2 = fz1
pb2 = pb1

u1 = u0
v1 = v0
w1 = w0
p1 = p0
fx1 = fx
fy1 = fy
fz1 = fz
pb1 = pb0

u0 = u
v0 = v
w0 = w
p0 = p

if (ifscl .ge. 1) then
    h2  = h1
    h1  = h0
    h0  = h
    fh2 = fh1
    fh1 = fh
endif

end