
! initiate field from laminar
subroutine inifield0(time0)
use global
implicit none

real time0
integer i, j, k, is, ii, kk


time0 = 0.0
u = 0
v = 0
w = 0

if (impi_myid .eq. 0) u(0, :, 0) = 1.5 * ( 1 - y**2 )

do i = 0, ndxhnpr-1
do k = 0, ndznpc-1
    ii = coord(0) * ndxhnpr + i 
    kk = - ( coord(1) * ndznpc + k )
    if (kk < -ndzh) kk = kk + ndz 
    ! add initial disturbance
    if ( &
        ( (ii >= 1) .and. (ii <= 8) ) .and. &
        ( ((kk >= 1) .and. (kk <= 8)) .or. ((kk <= -1) .and. (kk >= -8)) ) &
    ) then
        u(i, :, k) = u(i,:,k) + 1.5 * ( 1 - y**2 ) / 1000.
        v(i, :, k) = v(i,:,k) + 1.5 * ( 1 - y**2 ) / 1000.
        w(i, :, k) = w(i,:,k) + 1.5 * ( 1 - y**2 ) / 1000.
    endif
enddo
enddo

! assign values to older time steps
u0 = u
u1 = u
u2 = u
v0 = v
v1 = v
v2 = v
w0 = w
w1 = w
w2 = w
 
call vor
call utp
call wrfr
call pb1condition
call fp
call bulkforce
call pb2condition

pb1 = pb0
pb2 = pb0

fx1 = fx
fy1 = fy
fz1 = fz
fx2 = fx
fy2 = fy
fz2 = fz

fh1 = fh
fh2 = fh

u = u0
v = v0
w = w0
h = h0

if (impi_myid .eq. 0) then
    write(*,*) 'Initiated field from laminar.'
    open(20, file = 'runtimedata/INITFIELD.DAT')
        do j = 0, ndy
            write(20,*) y(j), real( u(0,j,0) )
        enddo
    close(20)
endif

end


! initiate field from a velocity field
subroutine inifield1(time0)
use global
implicit none

real time0
integer i, j, k, ii, kk, is, recn
character*2 cis
integer iname
character*8 cname
integer ndx0, ndz0
complex &
    ur(0:ndxhnpr-1, 0:ndyr, 0:ndznpc-1), &
    vr(0:ndxhnpr-1, 0:ndyr, 0:ndznpc-1), &
    wr(0:ndxhnpr-1, 0:ndyr, 0:ndznpc-1) 
character :: filenames(3) = (/ 'U', 'V', 'W' /)
integer :: filenums(3) = (/ 15, 16, 17 /)

if (impi_myid .eq. 0) then
    open(255, file = 'XINDAT')
    do i = 1, 19
        read(255,*) ! skip useless lines
    enddo
    read(255,*) iname, ndx0, ndz0
    close(255)
    write(*,*) 'Reading initial field from one-step velocity field.'
    write(*,*) 'iname =', iname
    write(*,*) 'ndx0 =', ndx0
    write(*,*) 'ndy0 =', ndyr
    write(*,*) 'ndz0 =', ndz0
endif

call mpi_bcast(iname,1, mpi_integer, 0, mpi_comm_world, impi_errorinfo)
call mpi_bcast(ndx0, 1, mpi_integer, 0, mpi_comm_world, impi_errorinfo)
call mpi_bcast(ndz0, 1, mpi_integer, 0, mpi_comm_world, impi_errorinfo)

time0 = 0.0
ur = 0.0
vr = 0.0
wr = 0.0
u = 0.0
v = 0.0
w = 0.0
h = 0.0

write(cname, '(I8.8)') iname

do i = 1, 3
    open(filenums(i), file = 'initialfield/'//filenames(i)//cname//'.BIN', &
        access = 'DIRECT', &
        recl = 8, & ! read one complex number each time
        status = 'UNKNOWN', &
        form = 'BINARY' )
enddo

do j = 0, ndyr
do k = 0, ndznpc-1
do i = 0, ndxhnpr-1
    ii = coord(0) * ndxhnpr + i
    kk = coord(1) * ndznpc + k
    if (ii .ge. ndx0/2) cycle ! read in only the corresponding Fourier modes of direction x
    if (kk .lt. ndzh) then
        if (kk .ge. ndz0/2) cycle ! read in only the corresponding Fourier modes of direction z
        recn = (j+1) * (ndx0/2*ndz0) + kk * ndx0/2 + ii + 1 ! j+1 for skipping the info section
    else
        if ( (ndz-kk) .ge. ndz0/2 ) cycle
        recn = (j+1) * (ndx0/2*ndz0) + (kk-ndz+ndz0) * ndx0/2 + ii + 1
    endif
    read(filenums(1), rec = recn) ur(i,j,k)
    read(filenums(2), rec = recn) vr(i,j,k)
    read(filenums(3), rec = recn) wr(i,j,k)
enddo
enddo
enddo

close(filenums(1))
close(filenums(2))
close(filenums(3))

! for wall normal direction, read in corresponding chebyshev points
if (ndy .ne. ndyr) then
    call chebsr(ur, ndxhnpr, ndyr, ndznpc)
    call chebsr(vr, ndxhnpr, ndyr, ndznpc)
    call chebsr(wr, ndxhnpr, ndyr, ndznpc)

    u(:, 0:min(ndy,ndyr), :) = ur(:, 0:min(ndy,ndyr), :)
    v(:, 0:min(ndy,ndyr), :) = vr(:, 0:min(ndy,ndyr), :)
    w(:, 0:min(ndy,ndyr), :) = wr(:, 0:min(ndy,ndyr), :)

    call chebp(u, ndxhnpr, ndy, ndznpc)
    call chebp(v, ndxhnpr, ndy, ndznpc)
    call chebp(w, ndxhnpr, ndy, ndznpc)
else
    u(:, 0:min(ndy,ndyr), :) = ur(:, 0:min(ndy,ndyr), :)
    v(:, 0:min(ndy,ndyr), :) = vr(:, 0:min(ndy,ndyr), :)
    w(:, 0:min(ndy,ndyr), :) = wr(:, 0:min(ndy,ndyr), :)
endif

! assign values to older time steps
u0 = u
u1 = u
u2 = u
v0 = v
v1 = v
v2 = v
w0 = w
w1 = w
w2 = w

call vor

h0 = h
h1 = h
h2 = h

if (ifscl .ge. 1) call cdh
call utp
call wrfr
call pb1condition
call fp
call bulkforce
call pb2condition

pb1 = pb0
pb2 = pb0

fx1 = fx
fx2 = fx
fy1 = fy
fy2 = fy
fz1 = fz
fz2 = fz

fh1 = fh
fh2 = fh

u = u0
v = v0
w = w0
h = h0

if (impi_myid .eq. 0) then
    open(20, file = 'runtimedata/INITFIELD.DAT')
        do j = 0, ndy
            write(20,*) y(j), real( u(0,j,0) )
        enddo
    close(20)
endif

end


! initiate field from mid files for continuing computation                    
subroutine inifieldm(time0)
use global
implicit none

real time0
integer i, j, k, is, kk
character*2 cis

!-------------------------------------------------------------------------------
!    execute statement
!-------------------------------------------------------------------------------	    
open(20, file = 'runtimedata/TIME.MID')
open(21, file = 'runtimedata/MIDVELO.MID', &
    access = 'DIRECT', &
    recl = 9 * ilrecxnpy, & ! read a y column each time
    status = 'UNKNOWN', &
    form = 'BINARY')
open(22, file = 'runtimedata/MIDFX.MID', &
    access = 'DIRECT', &
    recl = 6 * ilrecxnpy, & ! read a y column each time
    status = 'UNKNOWN', &
    FORM = 'BINARY')
open(23, file = 'runtimedata/MIDPB.MID', &
    access = 'DIRECT', &
    recl = 2 * ilrecxnpz, & ! read a y column each time
    status = 'UNKNOWN', &
    form = 'BINARY')

if (impi_myid .eq. 0) read(20,*) time0
call mpi_bcast(time0, 1, mpi_real, 0, mpi_comm_world, impi_errorinfo)

do k = 0, ndznpc-1
do i = 0, ndxhnpr-1
    kk = (coord(1) * ndznpc + k) * ndxh + coord(0) * ndxhnpr + i + 1
    read(21, rec = kk) &
        ( u0(i,j,k), j = 0, ndy ), &
        ( v0(i,j,k), j = 0, ndy ), &
        ( w0(i,j,k), j = 0, ndy ), &
        ( u1(i,j,k), j = 0, ndy ), &
        ( v1(i,j,k), j = 0, ndy ), &
        ( w1(i,j,k), j = 0, ndy ), &
        ( u2(i,j,k), j = 0, ndy ), &
        ( v2(i,j,k), j = 0, ndy ), &
        ( w2(i,j,k), j = 0, ndy )
    read(22, rec = kk) &
        ( fx1(i,j,k), j = 0, ndy ), &
        ( fy1(i,j,k), j = 0, ndy ), &
        ( fz1(i,j,k), j = 0, ndy ), &
        ( fx2(i,j,k), j = 0, ndy ), &
        ( fy2(i,j,k), j = 0, ndy ), &
        ( fz2(i,j,k), j = 0, ndy )
    read(23, rec = kk) pb1(i,k,1), pb2(i,k,1)
    read(23, rec = kk + ndxh*ndz) pb1(i,k,2), pb2(i,k,2) 
enddo
enddo

close(20)
close(21)
close(22)
close(23)

u = u0
v = v0
w = w0

if (ifscl .ge. 1) then
    
    do is = 1, ifscl
        write(cis, '(I2.2)') is
    
        open (24, file = 'MIDSCALAR'//cis//'.MID', &
            access = 'DIRECT', &
            recl = 5 * ilrecxnpy, &
            status = 'UNKNOWN', &
            form = 'BINARY')

        do k = 0, ndznpc-1
        do i = 0, ndxhnpr-1
            kk = (coord(1) * ndznpc + k) * ndxh + coord(0) * ndxhnpr + i + 1
            read(24, rec = kk) &
                ( h0 (i,j,k,is), j = 0, ndy ), &
                ( h1 (i,j,k,is), j = 0, ndy ), &
                ( h2 (i,j,k,is), j = 0, ndy ), &
                ( fh1(i,j,k,is), j = 0, ndy ), &
                ( fh2(i,j,k,is), j = 0, ndy ) 
        enddo
        enddo

    	close(24)
    enddo
    
    h = h0
    
endif

end
