
    
subroutine write_files(iwt, ikt, time)
use global
implicit none

integer iwt, ikt
real time

if ( mod(iwt, nsutp) .eq. 0 ) call write_utp(iwt, time)
if ( mod(iwt, nsp) .eq. 0) call write_fields(iwt, time)
if ( (mod(iwt, nsprb) .eq. 0) .and. (ifprb .ge. 1) ) call write_probe(iwt, time)
if ( (mod(iwt, nsm) .eq. 0) .or. (ikt .eq. nt) ) call writemid(time)

if (ifscl .ge. 1) then
    if (mod(iwt, nsp) .eq. 0) call write_scalar(iwt, time)
    if ( (mod(iwt, nsm) .eq. 0) .or. (ikt .eq. nt) ) call writemid_scalar(time)
endif

end

    
! watch wall stress, flow rate, divergence and mean profile on screen
subroutine write_utp(iwt, time)
use global
implicit none

integer iwt
real time

real utu, utl
common /wallstress/ utu, utl
real maxdiv
integer maxdivi, maxdivj, maxdivk
common /maxdiv/ maxdiv
common /maxdivpos/ maxdivi, maxdivj, maxdivk

integer j
integer file_stat
real file_time

call div

if (impi_myid .eq. 0) then
    
    open(111, file = 'runtimedata/XUTP.DAT')
        do while(.true.)
            read(111, *, iostat = file_stat) file_time
            if ( (file_stat .lt. 0) .or. (file_time .ge. time) ) exit
        enddo
        backspace(111)
        write(111,*) time, -dp, utu, utl, re*sqrt(abs(dp)) ! time, average wall shear, upper wall shear, lower wall shear, friction Reynolds number
    close(111)
    
    open (41, file = 'runtimedata/XFR.DAT')
        do while(.true.)
            read(41, *, iostat = file_stat) file_time
            if ( (file_stat .lt. 0) .or. (file_time .ge. time) ) exit
        enddo
        backspace(41)
        write(41,*) time, fr
    close(41)

    open(42, file = 'runtimedata/DIV.DAT')
        do while(.true.)
            read(42, *, iostat = file_stat) file_time
            if ( (file_stat .lt. 0) .or. (file_time .ge. time) ) exit
        enddo
        backspace(42)
        write(42,*) time, maxdiv, maxdivi, maxdivj, maxdivk
    close(42)
    
    open(20, file = 'runtimedata/UMEAN.DAT')
        do j = 0, ndy
            write(20,*) y(j), real( u(0,j,0) )
        enddo
    close(20)
        
endif

end


! write velocity, pressure and their time derivative fields
subroutine write_fields(iwt, time)
use global
implicit none

integer iwt
real time

complex &
    uu (0:ndxhnpr-1, 0:ndy, 0:ndznpc-1), & ! convolutions of velocities
    vv (0:ndxhnpr-1, 0:ndy, 0:ndznpc-1), &
    ww (0:ndxhnpr-1, 0:ndy, 0:ndznpc-1), &
    uv (0:ndxhnpr-1, 0:ndy, 0:ndznpc-1), &
    vw (0:ndxhnpr-1, 0:ndy, 0:ndznpc-1), &
    wu (0:ndxhnpr-1, 0:ndy, 0:ndznpc-1), &
    ps0(0:ndxhnpr-1, 0:ndy, 0:ndznpc-1), & ! static pressure
    ps1(0:ndxhnpr-1, 0:ndy, 0:ndznpc-1), &
    ps2(0:ndxhnpr-1, 0:ndy, 0:ndznpc-1)

character*8 cname
character*5 :: filenames(14) = (/ &
    'U/U', 'V/V', 'W/W', 'P/P', &
    'UU/UU', 'VV/VV', 'WW/WW', 'UV/UV', 'VW/VW', 'WU/WU', &
    'UT/UT', 'VT/VT', 'WT/WT', 'PT/PT' /)
integer :: filenums(14) = (/ &
    15, 16, 17, 18, &
    19, 20, 21, 22, 23, 24, &
    25, 26, 27, 28 /)
integer i, j, k, kk, temp1, temp2, temp3, temp4

call compute_convol(u2, u2, uu)
call compute_convol(v2, v2, vv)
call compute_convol(w2, w2, ww)
ps2 = p2 - 0.5 * (uu + vv + ww)
if (impi_myid .eq. 0) ps2(0,:,0) = - vv(0,:,0) ! average pressure not computed in the program but assigned here

call compute_convol(u1, u1, uu)
call compute_convol(v1, v1, vv)
call compute_convol(w1, w1, ww)
ps1 = p1 - 0.5 * (uu + vv + ww)
if (impi_myid .eq. 0) ps1(0,:,0) = - vv(0,:,0)

call compute_convol(u, u, uu)
call compute_convol(v, v, vv)
call compute_convol(w, w, ww)
call compute_convol(u, v, uv)
call compute_convol(v, w, vw)
call compute_convol(w, u, wu)
ps0 = p0 - 0.5 * (uu + vv + ww)
if (impi_myid .eq. 0) ps0(0,:,0) = - vv(0,:,0)


write(cname, '(I8.8)') iwt

do i = 1, 14
    open(filenums(i), file = 'fielddata/'//trim(filenames(i))//cname//'.BIN', &
        access = 'DIRECT', &
        recl = 8 * ndxhnpr, & ! write an XZ plane each time
        status = 'UNKNOWN', &
        form = 'BINARY' )
    if (impi_myid .eq. 0) write(filenums(i), rec = 1) time ! information section
enddo

do j = 0, ndy
do k = 0, ndznpc-1
    kk = (j+1) * (npr*ndz) + (coord(1) * ndznpc + k) * npr + coord(0) + 1 ! j+1 for skipping the info section
    write(filenums(1), rec = kk) (u (i,j,k), i = 0, ndxhnpr-1)
    write(filenums(2), rec = kk) (v (i,j,k), i = 0, ndxhnpr-1)
    write(filenums(3), rec = kk) (w (i,j,k), i = 0, ndxhnpr-1)
    write(filenums(4), rec = kk) (ps0(i,j,k), i = 0, ndxhnpr-1)
    write(filenums(5), rec = kk) (uu(i,j,k), i = 0, ndxhnpr-1)
    write(filenums(6), rec = kk) (vv(i,j,k), i = 0, ndxhnpr-1)
    write(filenums(7), rec = kk) (ww(i,j,k), i = 0, ndxhnpr-1)
    write(filenums(8), rec = kk) (uv(i,j,k), i = 0, ndxhnpr-1)
    write(filenums(9), rec = kk) (vw(i,j,k), i = 0, ndxhnpr-1)
    write(filenums(10),rec = kk) (wu(i,j,k), i = 0, ndxhnpr-1)
    write(filenums(11),rec = kk) (( (1.5*u0(i,j,k)-2.0*u1(i,j,k)+0.5*u2(i,j,k)) / dt ), i = 0, ndxhnpr-1)
    write(filenums(12),rec = kk) (( (1.5*v0(i,j,k)-2.0*v1(i,j,k)+0.5*v2(i,j,k)) / dt ), i = 0, ndxhnpr-1)
    write(filenums(13),rec = kk) (( (1.5*w0(i,j,k)-2.0*w1(i,j,k)+0.5*w2(i,j,k)) / dt ), i = 0, ndxhnpr-1)
    write(filenums(14),rec = kk) (( (1.5*ps0(i,j,k)-2.0*ps1(i,j,k)+0.5*ps2(i,j,k)) / dt ), i = 0, ndxhnpr-1)
enddo
enddo

do i = 1, 14
    close(filenums(i))
enddo

call mpi_barrier( mpi_comm_world, impi_errorinfo )
if (impi_myid .eq. 0) write(*,*) 'Finish writing 14 whole fields at step ', iwt

! update probe status
if (impi_myid .eq. 0) then
    open(20, file = 'XINDAT')
        do i = 1, 3
            read(20, *) ! skip useless lines
        enddo
        read(20, *) temp1, temp2, temp3, temp4, nsprb
        do i = 1, 5
            read(20, *) ! skip useless lines
        enddo
        read(20, *) temp1, temp2, temp3, ifprb
        do i = 1, 7
            read(20, *) ! skip useless lines
        enddo
        read(20, *) j_probes(1:ifprb)
    close(20)
endif
call mpi_bcast( nsprb, 1, mpi_integer, 0, mpi_comm_world, impi_errorinfo )
call mpi_bcast( ifprb, 1, mpi_integer, 0, mpi_comm_world, impi_errorinfo )
call mpi_bcast( j_probes, ndy+1, mpi_integer, 0, mpi_comm_world, impi_errorinfo )

end


! write high-time-resolution data of slices specified in XINDAT
subroutine write_probe(iwt, time)
use global
implicit none

integer iwt
real time

complex &
    uu (0:ndxhnpr-1, 0:ndy, 0:ndznpc-1), & ! convolutions of velocities
    vv (0:ndxhnpr-1, 0:ndy, 0:ndznpc-1), &
    ww (0:ndxhnpr-1, 0:ndy, 0:ndznpc-1), &
    ps (0:ndxhnpr-1, 0:ndy, 0:ndznpc-1) ! static pressure
    
character*8 cname
character*3 :: filenames(4) = (/ 'U/U', 'V/V', 'W/W', 'P/P' /)
integer :: filenums(4) = (/ 15, 16, 17, 18 /)
integer i, j, k, kk, jprb
integer start_iwt, file_stat

! get start iwt
open(111, file = 'probedata/start_iwt')
read(111, *, iostat = file_stat) start_iwt
if (file_stat .lt. 0) then
    start_iwt = iwt
    write(111,*) start_iwt
endif
close(111)

! static pressure is solved for whole field here, which is not necessary. better to be solved for only j_probes
call compute_convol(u, u, uu)
call compute_convol(v, v, vv)
call compute_convol(w, w, ww)
ps = p - 0.5 * (uu + vv + ww) ! total pressure to static pressure
if (impi_myid .eq. 0) ps(0,:,0) = - vv(0,:,0) ! average pressure assigned here

            
do j = 1, ifprb
    jprb = j_probes(j)
    
    ! solve static pressure
    if (impi_myid .eq. 0) then
        ! gather data from other cores
        ! do ifft in z direction
        ! do ifft in x direction
        ! compute product
        ! do fft in x direction
        ! do fft in z direction
        ! broadcast to other cores
        ! save into uu, vv, ww
        ! get static pressure (including average pressure)
    endif
    
    
    write(cname, '(I8.8)') jprb
    
    do i = 1, 4
        open(filenums(i), file = 'probedata/'//trim(filenames(i))//cname//'.BIN', &
            access = 'DIRECT', &
            recl = 8 * ndxhnpr, &
            status = 'UNKNOWN', &
            form = 'BINARY' )
        if (impi_myid .eq. 0) write(filenums(i), rec = 1) start_iwt*dt, time ! information section
    enddo

    do k = 0, ndznpc-1
        kk = ((iwt-start_iwt)/nsprb + 1) * (npr*ndz) + (coord(1) * ndznpc + k) * npr + coord(0) + 1
        write(15, rec = kk) (u (i,jprb,k), i = 0,ndxhnpr-1)
        write(16, rec = kk) (v (i,jprb,k), i = 0,ndxhnpr-1)
        write(17, rec = kk) (w (i,jprb,k), i = 0,ndxhnpr-1)
        write(18, rec = kk) (ps(i,jprb,k), i = 0,ndxhnpr-1)
    enddo

    do i = 1, 4
        close(filenums(i))
    enddo
    
enddo

call mpi_barrier( mpi_comm_world, impi_errorinfo )
if (impi_myid .eq. 0) write(*,*) 'Finish probing ', ifprb, ' slices at step ', iwt

end


! write mid-file for continue computing
subroutine writemid(time)
use global
implicit none

real time

integer i, j, k, kk

if (impi_myid .eq. 0) write(*,*) 'Writing middle files...'

! write the major mid files
open(20, file = 'runtimedata/TIME.MID')
open(21, file = 'runtimedata/MIDVELO.MID', &
    access = 'DIRECT', &
    recl = 9 * ilrecxnpy, & ! write a y column each time
    status = 'UNKNOWN', &
    form = 'BINARY')
open(22, file = 'runtimedata/MIDFX.MID', &
    access = 'DIRECT', &
    recl = 6 * ilrecxnpy, & ! write a y column each time
    status = 'UNKNOWN', &
    FORM = 'BINARY')
open(23, file = 'runtimedata/MIDPB.MID', &
    access = 'DIRECT', &
    recl = 2 * ilrecxnpz, & ! write a y column each time
    status = 'UNKNOWN', &
    form = 'BINARY')

if (impi_myid .eq. 0) write(20,*) time

do k = 0, ndznpc-1
do i = 0, ndxhnpr-1
    kk = (coord(1) * ndznpc + k) * ndxh + coord(0) * ndxhnpr + i + 1
    write(21, rec = kk) &
        ( u0(i,j,k), j = 0, ndy ), &
        ( v0(i,j,k), j = 0, ndy ), &
        ( w0(i,j,k), j = 0, ndy ), &
        ( u1(i,j,k), j = 0, ndy ), &
        ( v1(i,j,k), j = 0, ndy ), &
        ( w1(i,j,k), j = 0, ndy ), &
        ( u2(i,j,k), j = 0, ndy ), &
        ( v2(i,j,k), j = 0, ndy ), &
        ( w2(i,j,k), j = 0, ndy )
    write(22, rec = kk) &
        ( fx1(i,j,k), j = 0, ndy ), &
        ( fy1(i,j,k), j = 0, ndy ), &
        ( fz1(i,j,k), j = 0, ndy ), &
        ( fx2(i,j,k), j = 0, ndy ), &
        ( fy2(i,j,k), j = 0, ndy ), &
        ( fz2(i,j,k), j = 0, ndy )
    write(23, rec = kk) pb1(i,k,1), pb2(i,k,1)
    write(23, rec = kk + ndxh*ndz) pb1(i,k,2), pb2(i,k,2) 
enddo
enddo

close(20)
close(21)
close(22)
close(23)

call mpi_barrier( mpi_comm_world, impi_errorinfo )
if (impi_myid .eq. 0) write(*,*) 'Mid files written.'      

! write backup mid files, in case system error occur during mid file writing
open(20, file = 'runtimedata/TIME1.MID')
open(21, file = 'runtimedata/MIDVELO1.MID', &
    access = 'DIRECT', &
    recl = 9 * ilrecxnpy, & ! write a y column each time
    status = 'UNKNOWN', &
    form = 'BINARY')
open(22, file = 'runtimedata/MIDFX1.MID', &
    access = 'DIRECT', &
    recl = 6 * ilrecxnpy, & ! write a y column each time
    status = 'UNKNOWN', &
    FORM = 'BINARY')
open(23, file = 'runtimedata/MIDPB1.MID', &
    access = 'DIRECT', &
    recl = 2 * ilrecxnpz, & ! write a y column each time
    status = 'UNKNOWN', &
    form = 'BINARY')

if (impi_myid .eq. 0) write(20,*) time

do k = 0, ndznpc-1
do i = 0, ndxhnpr-1
    kk = (coord(1) * ndznpc + k) * ndxh + coord(0) * ndxhnpr + i + 1
    write(21, rec = kk) &
        ( u0(i,j,k), j = 0, ndy ), &
        ( v0(i,j,k), j = 0, ndy ), &
        ( w0(i,j,k), j = 0, ndy ), &
        ( u1(i,j,k), j = 0, ndy ), &
        ( v1(i,j,k), j = 0, ndy ), &
        ( w1(i,j,k), j = 0, ndy ), &
        ( u2(i,j,k), j = 0, ndy ), &
        ( v2(i,j,k), j = 0, ndy ), &
        ( w2(i,j,k), j = 0, ndy )
    write(22, rec = kk) &
        ( fx1(i,j,k), j = 0, ndy ), &
        ( fy1(i,j,k), j = 0, ndy ), &
        ( fz1(i,j,k), j = 0, ndy ), &
        ( fx2(i,j,k), j = 0, ndy ), &
        ( fy2(i,j,k), j = 0, ndy ), &
        ( fz2(i,j,k), j = 0, ndy )
    write(23, rec = kk) pb1(i,k,1), pb2(i,k,1)
    write(23, rec = kk + ndxh*ndz) pb1(i,k,2), pb2(i,k,2) 
enddo
enddo

close(20)
close(21)
close(22)
close(23)

call mpi_barrier( mpi_comm_world, impi_errorinfo )
if (impi_myid .eq. 0) write(*,*) 'Backup mid files written.' 

end


    
    
    
    
    
! write scalars to file
subroutine write_scalar(iwt, time)
use global
implicit none

integer iwt
real time

character*5 cname
character*2 cis
integer i, j, k, is, kk

do is = 1, ifscl
    
    write(cis, '(I2.2)') is
    write(cname, '(I8.8)') iwt
    
    open(15, file = 'fielddata/'//cis//'S'//cname//'.BIN', &
        access = 'DIRECT', &
        recl = 8 * ndxhnpr, & ! write an XZ plane each time
        status = 'UNKNOWN', &
        form = 'BINARY' )

    if (impi_myid .eq. 0) write(*,*) 'Writing Scalar field ', cis//'S'//cname//'.BIN'
    if (impi_myid .eq. 0) write(15, rec = 1) time ! information section

    do j = 0, ndy
    do k = 0, ndznpc-1
        kk = (j+1) * (npr*ndz) + (coord(1) * ndznpc + k) * npr + coord(0) + 1 ! j+1 for skipping the info section
        write(15, rec = kk) (h(i,j,k,is), i = 0,ndxhnpr-1)
    enddo
    enddo

    close( 15 )       

enddo

end

! output mid-file of scalar for continue caculate 
subroutine  writemid_scalar(time)
use global
implicit none

real time
integer i, j, k, is, kk
character*2 cis

do is = 1, ifscl
    write(cis, '(I2.2)') is
    
    open (21, file = 'MIDSCALAR'//cis//'.MID', &
        access = 'DIRECT', &
        recl = 5 * ilrecxnpy, &
        status = 'UNKNOWN', &
        form = 'BINARY')
    
    if (impi_myid .eq. 0) write(*,*) 'Writing middle Scalar files...'

    do k = 0, ndznpc-1
    do i = 0, ndxhnpr-1
        kk = (coord(1) * ndznpc + k) * ndxh + coord(0) * ndxhnpr + i + 1
        write(21, rec = kk) &
            ( h0(i,j,k,is), j = 0, ndy ), &
            ( h1(i,j,k,is), j = 0, ndy ), &
            ( h2(i,j,k,is), j = 0, ndy ), &
            ( fh1(i,j,k,is), j = 0, ndy ), &
            ( fh2(i,j,k,is), j = 0, ndy ) 
    enddo
    enddo

    close(21)
enddo

end

    
    