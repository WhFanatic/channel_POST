
! read file from XINDAT
subroutine readinitialfile
use global
implicit none

integer mx, my, mz

if ( impi_myid .eq. 0 ) then
    open(20, file = 'XINDAT')
        read(20, *)
        read(20, *) re, arf, bat, dt, mx, my, mz
        read(20, *)
        read(20, *) nt, nsp, nsm, nsutp, nsprb
        read(20, *)
        read(20, *) ifc, iff
        read(20, *)
        read(20, *) fr0, dp0
        read(20, *)
        read(20, *) ifscl, ifsgs, ifmfu, ifprb
        read(20, *)
	    read(20, *) ifscalarbound 
	    read(20, *)
	    read(20, *) ifscalarsgs
	    read(20, *)
	    read(20, *) pr
        read(20, *)
        read(20, *) j_probes(1:ifprb)
    close(20)

    write(*,*) 'Re =', re
    write(*,*) 'Lx =', 2/arf, 'pi'
    write(*,*) 'Lz =', 2/bat, 'pi'
    write(*,*) 'dt =', dt
    write(*,*) 'Nx =', mx
    write(*,*) 'Ny =', my
    write(*,*) 'Nz =', mz
    write(*,*) 'nt =', nt
    write(*,*) 'nsp =', nsp
    write(*,*) 'nsm =', nsm
    write(*,*) 'nsutp =', nsutp
    write(*,*) 'nsprb =', nsprb
endif

call mpi_bcast( re, 1, mpi_real, 0, mpi_comm_world, impi_errorinfo )
call mpi_bcast( arf,1, mpi_real, 0, mpi_comm_world, impi_errorinfo )
call mpi_bcast( bat,1, mpi_real, 0, mpi_comm_world, impi_errorinfo )
call mpi_bcast( dt, 1, mpi_real, 0, mpi_comm_world, impi_errorinfo )
call mpi_bcast( mx, 1, mpi_integer, 0, mpi_comm_world, impi_errorinfo )
call mpi_bcast( my, 1, mpi_integer, 0, mpi_comm_world, impi_errorinfo )
call mpi_bcast( mz, 1, mpi_integer, 0, mpi_comm_world, impi_errorinfo )
call mpi_bcast( nt, 1, mpi_integer, 0, mpi_comm_world, impi_errorinfo )
call mpi_bcast( nsp,1, mpi_integer, 0, mpi_comm_world, impi_errorinfo )
call mpi_bcast( nsm,1, mpi_integer, 0, mpi_comm_world, impi_errorinfo )
call mpi_bcast( nsutp,1, mpi_integer, 0, mpi_comm_world, impi_errorinfo )
call mpi_bcast( nsprb,1, mpi_integer, 0, mpi_comm_world, impi_errorinfo )
call mpi_bcast( ifc,1, mpi_integer, 0, mpi_comm_world, impi_errorinfo )
call mpi_bcast( iff,1, mpi_integer, 0, mpi_comm_world, impi_errorinfo )
call mpi_bcast( fr0,1, mpi_real, 0, mpi_comm_world, impi_errorinfo )
call mpi_bcast( dp0,1, mpi_real, 0, mpi_comm_world, impi_errorinfo )
call mpi_bcast( ifscl, 1, mpi_integer, 0, mpi_comm_world, impi_errorinfo )
call mpi_bcast( ifsgs, 1, mpi_integer, 0, mpi_comm_world, impi_errorinfo )
call mpi_bcast( ifmfu, 1, mpi_integer, 0, mpi_comm_world, impi_errorinfo )
call mpi_bcast( ifprb, 1, mpi_integer, 0, mpi_comm_world, impi_errorinfo )
call mpi_bcast( ifscalarbound, ns, mpi_integer, 0, mpi_comm_world, impi_errorinfo )
call mpi_bcast( ifscalarsgs, ns, mpi_integer, 0, mpi_comm_world, impi_errorinfo )
call mpi_bcast( pr, ns, mpi_real, 0, mpi_comm_world, impi_errorinfo )
call mpi_bcast( j_probes, ndy+1, mpi_integer, 0, mpi_comm_world, impi_errorinfo )

! parameters check
if ( (mx .ne. ndx) .or. (my .ne. ndy) .or. (mz .ne. ndz) ) then
    write(*,*) 'Grid size dose not match !'
    stop
endif

do while (mod(mx, 2) .eq. 0)
    mx = mx / 2
enddo
do while (mod(mx, 3) .eq. 0)
    mx = mx / 3
enddo
do while (mod(mz, 2) .eq. 0)
    mz = mz / 2
enddo
do while (mod(mz, 3) .eq. 0)
    mz = mz / 3
enddo
if ( (mx .ne. 1) .or. (mz .ne. 1) .or. &
    (mod(ndxh,2*npr) .ne. 0) .or. (mod(ndxh,2*npc) .ne. 0) .or. &
    (mod(ndy, 2*npr) .ne. 0) .or. (mod(ndy, 2*npc) .ne. 0) .or. &
    (mod(ndz, 2*npr) .ne. 0) .or. (mod(ndz, 2*npc) .ne. 0) ) then
    write(*,*) 'Invalid grid size !'
    stop
endif

if ( ns .lt. ifscl ) then
    write(*,*) 'Error in number of scalar !'
    stop
endif

end


! initialzie constant variable
subroutine initializevariable
use global
implicit none

integer i, j, k ! position of the point in the sub-field on one processor
integer kk, ii ! position of the point in the whole field
real hm1(0:ndy), hm2(0:ndy), hp1(0:ndy), hp2(0:ndy), hp3, hm3 ! y intervals
real matfir(1:4,1:4), matsec(1:5,1:5), matin(1:7,1:7) ! matrixes for determining FD coefficients
real bfir(1:4), bsec(1:5), bin(1:7) ! vectors for determining FD coefficients (RHS of the linear equation)
real kpre(0:ndxhnpr-1, 0:ndznpc-1), kvis(0:ndxhnpr-1, 0:ndznpc-1)
real secrhsc(0:ndy, 1:5)
real &
    secdervisa(2:ndy, 0:ndxhnpr-1, 0:ndznpc-1), &
    secdervisb(1:ndy, 0:ndxhnpr-1, 0:ndznpc-1), &
    secdervisc(0:ndy, 0:ndxhnpr-1, 0:ndznpc-1), &
    secdervisd(0:ndy-1, 0:ndxhnpr-1, 0:ndznpc-1), &
    secdervise(0:ndy-2, 0:ndxhnpr-1, 0:ndznpc-1)
real &
    secderprea(2:ndy, 0:ndxhnpr-1, 0:ndznpc-1), &
    secderpreb(1:ndy, 0:ndxhnpr-1, 0:ndznpc-1), &
    secderprec(0:ndy, 0:ndxhnpr-1, 0:ndznpc-1), &
    secderpred(0:ndy-1, 0:ndxhnpr-1, 0:ndznpc-1), &
    secderpree(0:ndy-2, 0:ndxhnpr-1, 0:ndznpc-1)

!-------------------------------------------------------------------------------
!    execute statement
!-------------------------------------------------------------------------------

call prefft
!call iy_matrix
!call dy_matrix

!!!!!!!!!!!!!!!!!!!!!! set mesh grids !!!!!!!!!!!!!!!!!!!!!!!!!!!
do i = 0, ndx-1
    x(i) = 2.*pai / arf * i/ndx
enddo

do k = 0, ndz-1
    z(k) = 2.*pai / bat * k/ndz
enddo

do j = 0, ndy
    y(j) = -cos( j * pai / ndy )
enddo

!!!!!!!!!!!!!!!!!!!! calculate FD coefficients for compact 1st derivative !!!!!!!!!!!!!!!!!!!!!

firderb = 1.0

do j = 0, ndy
    select case(j)
        
    case(0)
        hp1(j) = y(j+1) - y(j)
        hp2(j) = y(j+2) - y(j)
        hp3 = y(j+3) - y(j)
        matfir(1, 1:4) = (/ 1.0, 1.0, 1.0, 0.0 /)
        matfir(2, 1:4) = (/ 0.0, hp1(j), hp2(j), -1.0 /)
        matfir(3, 1:4) = (/ 0.0, (hp1(j))**2.0, (hp2(j))**2.0, -2.0*hp1(j) /)
        matfir(4, 1:4) = (/ 0.0, (hp1(j))**3.0, (hp2(j))**3.0, -3.0*(hp1(j))**2.0 /)
        bfir(1:4) = (/ 0.0, 1.0, 0.0, 0.0 /)
        call gaussian_eli(4, matfir, bfir)
        firrhsc(j,3) = bfir(1)
        firrhsc(j,4) = bfir(2)
        firrhsc(j,5) = bfir(3)
        firderc(j) = bfir(4)

    case(1)
        hp1(j) = y(j+1) - y(j)
        hm1(j) = y(j-1) - y(j)
        matsec(1, 1:5) = (/ 1.0, 1.0, 1.0, 0.0, 0.0 /)
        matsec(2, 1:5) = (/ (hm1(j))**1, 0.0, (hp1(j))**1.0, -1.0, -1.0 /)
        matsec(3, 1:5) = (/ (hm1(j))**2, 0.0, (hp1(j))**2.0, -2.0*(hm1(j))**1.0, -2.0*(hp1(j))**1.0 /)
        matsec(4, 1:5) = (/ (hm1(j))**3, 0.0, (hp1(j))**3.0, -3.0*(hm1(j))**2.0, -3.0*(hp1(j))**2.0 /)
        matsec(5, 1:5) = (/ (hm1(j))**4, 0.0, (hp1(j))**4.0, -4.0*(hm1(j))**3.0, -4.0*(hp1(j))**3.0 /)
        bsec(1:5) = (/ 0.0, 1.0, 0.0, 0.0, 0.0 /)
        call gaussian_eli(5, matsec, bsec)
        firrhsc(j,2) = bsec(1)
        firrhsc(j,3) = bsec(2)
        firrhsc(j,4) = bsec(3)
        firdera(j) = bsec(4)
        firderc(j) = bsec(5)
   
    case(ndy-1)
        hp1(j) = y(j+1) - y(j)
        hm1(j) = y(j-1) - y(j)
        matsec(1, 1:5) = (/ 1.0, 1.0, 1.0, 0.0, 0.0 /)
        matsec(2, 1:5) = (/ (hm1(j))**1, 0.0, (hp1(j))**1.0, -1.0, -1.0 /)
        matsec(3, 1:5) = (/ (hm1(j))**2, 0.0, (hp1(j))**2.0, -2.0*(hm1(j))**1.0, -2.0*(hp1(j))**1.0 /)
        matsec(4, 1:5) = (/ (hm1(j))**3, 0.0, (hp1(j))**3.0, -3.0*(hm1(j))**2.0, -3.0*(hp1(j))**2.0 /)
        matsec(5, 1:5) = (/ (hm1(j))**4, 0.0, (hp1(j))**4.0, -4.0*(hm1(j))**3.0, -4.0*(hp1(j))**3.0 /)
        bsec(1:5) = (/ 0.0, 1.0, 0.0, 0.0, 0.0 /)
        call gaussian_eli(5, matsec, bsec)
        firrhsc(j,2) = bsec(1)
        firrhsc(j,3) = bsec(2)
        firrhsc(j,4) = bsec(3)
        firdera(j) = bsec(4)
        firderc(j) = bsec(5)
    
    case(ndy)
        hm1(j) = y(j-1) - y(j)
        hm2(j) = y(j-2) - y(j)
        hm3 = y(j-3) - y(j)
        matfir(1, 1:4) = (/ 1.0, 1.0, 1.0, 0.0 /)
        matfir(2, 1:4) = (/ 0.0, hm1(j), hm2(j), -1.0 /)
        matfir(3, 1:4) = (/ 0.0, (hm1(j))**2.0, (hm2(j))**2.0, -2.0*hm1(j) /)
        matfir(4, 1:4) = (/ 0.0, (hm1(j))**3.0, (hm2(j))**3.0, -3.0*(hm1(j))**2.0 /)
        bfir(1:4) = (/ 0.0, 1.0, 0.0, 0.0 /)
        call gaussian_eli(4, matfir, bfir)
        firrhsc(j,1) = bfir(3)
        firrhsc(j,2) = bfir(2)
        firrhsc(j,3) = bfir(1)
        firdera(j) = bfir(4)

    case default
        hp1(j) = y(j+1) - y(j)
        hp2(j) = y(j+2) - y(j)
        hm1(j) = y(j-1) - y(j)
        hm2(j) = y(j-2) - y(j)
        matin(1, 1:7) = (/ 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0 /)
        matin(2, 1:7) = (/ (hm2(j))**1.0, (hm1(j))**1.0, 0.0, (hp1(j))**1.0, (hp2(j))**1.0, -1.0, -1.0 /)
        matin(3, 1:7) = (/ (hm2(j))**2.0, (hm1(j))**2.0, 0.0, (hp1(j))**2.0, (hp2(j))**2.0, -2.0*(hm1(j))**1.0, -2.0*(hp1(j))**1.0 /)
        matin(4, 1:7) = (/ (hm2(j))**3.0, (hm1(j))**3.0, 0.0, (hp1(j))**3.0, (hp2(j))**3.0, -3.0*(hm1(j))**2.0, -3.0*(hp1(j))**2.0 /)
        matin(5, 1:7) = (/ (hm2(j))**4.0, (hm1(j))**4.0, 0.0, (hp1(j))**4.0, (hp2(j))**4.0, -4.0*(hm1(j))**3.0, -4.0*(hp1(j))**3.0 /)
        matin(6, 1:7) = (/ (hm2(j))**5.0, (hm1(j))**5.0, 0.0, (hp1(j))**5.0, (hp2(j))**5.0, -5.0*(hm1(j))**4.0, -5.0*(hp1(j))**4.0 /)
        matin(7, 1:7) = (/ (hm2(j))**6.0, (hm1(j))**6.0, 0.0, (hp1(j))**6.0, (hp2(j))**6.0, -6.0*(hm1(j))**5.0, -6.0*(hp1(j))**5.0 /)
        bin(1:7) = (/ 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
        call gaussian_eli(7, matin, bin)
        firrhsc(j,1) = bin(1)
        firrhsc(j,2) = bin(2)
        firrhsc(j,3) = bin(3)
        firrhsc(j,4) = bin(4)
        firrhsc(j,5) = bin(5)
        firdera(j) = bin(6)
        firderc(j) = bin(7)
        
    end select 
enddo

!!!!!!!!!!!!!!!!!!!!!!!!! calculate FD coefficients for compact 2nd derivative !!!!!!!!!!!!!!!!!!!!!!!!!!!!

! operator for d^2/dy^2
secderb = 1.0

do j = 1, ndy-1
    select case(j)

    case(1)
        hp1(j) = y(j+1) - y(j)
        hp2(j) = y(j+2) - y(j)
        hm1(j) = y(j-1) - y(j)
        matsec(1, 1:5) = (/ 1.0, 1.0, 1.0, 1.0, 0.0 /)
        matsec(2, 1:5) = (/ (hm1(j))**1, 0.0, (hp1(j))**1.0, (hp2(j))**1.0, 0.0 /)
        matsec(3, 1:5) = (/ (hm1(j))**2, 0.0, (hp1(j))**2.0, (hp2(j))**2.0, -2.0 /)
        matsec(4, 1:5) = (/ (hm1(j))**3, 0.0, (hp1(j))**3.0, (hp2(j))**3.0, -6.0*(hp1(j))**1.0 /)
        matsec(5, 1:5) = (/ (hm1(j))**4, 0.0, (hp1(j))**4.0, (hp2(j))**4.0, -12.0*(hp1(j))**2.0 /)
        bsec(1:5) = (/ 0.0, 0.0, 2.0, 0.0, 0.0 /)
        call gaussian_eli(5, matsec, bsec)
        secrhsc(j,2) = bsec(1)
        secrhsc(j,3) = bsec(2)
        secrhsc(j,4) = bsec(3)
        secrhsc(j,5) = bsec(4)
        secderc(j) = bsec(5)
        secdera(j) = 0.0
   
    case(ndy-1)
        hp1(j) = y(j+1) - y(j)
        hm1(j) = y(j-1) - y(j)
        hm2(j) = y(j-2) - y(j)
        matsec(1, 1:5) = (/ 1.0, 1.0, 1.0, 1.0, 0.0 /)
        matsec(2, 1:5) = (/ (hm2(j))**1, (hm1(j))**1, 0.0, (hp1(j))**1.0, 0.0 /)
        matsec(3, 1:5) = (/ (hm2(j))**2, (hm1(j))**2, 0.0, (hp1(j))**2.0, -2.0 /)
        matsec(4, 1:5) = (/ (hm2(j))**3, (hm1(j))**3, 0.0, (hp1(j))**3.0, -6.0*(hm1(j))**1.0 /)
        matsec(5, 1:5) = (/ (hm2(j))**4, (hm1(j))**4, 0.0, (hp1(j))**4.0, -12.0*(hm1(j))**2.0 /)
        bsec(1:5) = (/ 0.0, 0.0, 2.0, 0.0, 0.0 /)
        call gaussian_eli(5, matsec, bsec)
        secrhsc(j,1) = bsec(1)
        secrhsc(j,2) = bsec(2)
        secrhsc(j,3) = bsec(3)
        secrhsc(j,4) = bsec(4)
        secdera(j) = bsec(5)
        secderc(j) = 0.0

    case default
        hp1(j) = y(j+1) - y(j)
        hp2(j) = y(j+2) - y(j)
        hm1(j) = y(j-1) - y(j)
        hm2(j) = y(j-2) - y(j)
        matin(1, 1:7) = (/ 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0 /)
        matin(2, 1:7) = (/ (hm2(j))**1.0, (hm1(j))**1.0, 0.0, (hp1(j))**1.0, (hp2(j))**1.0, 0.0, 0.0 /)
        matin(3, 1:7) = (/ (hm2(j))**2.0, (hm1(j))**2.0, 0.0, (hp1(j))**2.0, (hp2(j))**2.0, -2.0, -2.0 /)
        matin(4, 1:7) = (/ (hm2(j))**3.0, (hm1(j))**3.0, 0.0, (hp1(j))**3.0, (hp2(j))**3.0, -6.0*(hm1(j))**1.0, -6.0*(hp1(j))**1.0 /)
        matin(5, 1:7) = (/ (hm2(j))**4.0, (hm1(j))**4.0, 0.0, (hp1(j))**4.0, (hp2(j))**4.0, -12.0*(hm1(j))**2.0, -12.0*(hp1(j))**2.0 /)
        matin(6, 1:7) = (/ (hm2(j))**5.0, (hm1(j))**5.0, 0.0, (hp1(j))**5.0, (hp2(j))**5.0, -20.0*(hm1(j))**3.0, -20.0*(hp1(j))**3.0 /)
        matin(7, 1:7) = (/ (hm2(j))**6.0, (hm1(j))**6.0, 0.0, (hp1(j))**6.0, (hp2(j))**6.0, -30.0*(hm1(j))**4.0, -30.0*(hp1(j))**4.0 /)
        bin(1:7) = (/ 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0 /)
        call gaussian_eli(7, matin, bin)
        secrhsc(j,1) = bin(1)
        secrhsc(j,2) = bin(2)
        secrhsc(j,3) = bin(3)
        secrhsc(j,4) = bin(4)
        secrhsc(j,5) = bin(5)
        secdera(j) = bin(6)
        secderc(j) = bin(7)

    end select 
enddo

! operator for ( \nabla^2 ) & ( \nabla^2 - \gamma_0 * Re / dt )
do k = 0, ndznpc-1
    kk = -( coord(1) * ndznpc + k )
    if (kk .lt. -ndzh) kk = kk + ndz
    do i = 0, ndxhnpr-1
        ii = -( coord(0) * ndxhnpr + i )
        
        kpre(i,k) = -( (arf * ii)**2 + (bat * kk)**2 )
        kvis(i,k) = -( (arf * ii)**2 + (bat * kk)**2 + 11./6. * re/dt )
        
        do j = 1, ndy-1
            secdervisb(j, i, k) = secrhsc(j,2) + kvis(i,k) * secdera(j)
            secdervisc(j, i, k) = secrhsc(j,3) + kvis(i,k) * secderb(j)
            secdervisd(j, i, k) = secrhsc(j,4) + kvis(i,k) * secderc(j)

            secderpreb(j, i, k) = secrhsc(j,2) + kpre(i,k) * secdera(j)
            secderprec(j, i, k) = secrhsc(j,3) + kpre(i,k) * secderb(j)
            secderpred(j, i, k) = secrhsc(j,4) + kpre(i,k) * secderc(j)

            select case(j)
            case (1)
                secdervise(j, i, k) = secrhsc(j,5)
                secderpree(j, i, k) = secrhsc(j,5)
            case (ndy-1)
                secdervisa(j, i, k) = secrhsc(j,1)
                secderprea(j, i, k) = secrhsc(j,1)
            case default
                secdervisa(j, i, k) = secrhsc(j,1)
                secderprea(j, i, k) = secrhsc(j,1)
                secdervise(j, i, k) = secrhsc(j,5)
                secderpree(j, i, k) = secrhsc(j,5)
            end select
        enddo
        
    enddo
enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!! manipulate coefficient matrix on boundaries to apply BC !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Neumann for pressure, Dirichlet for viscuous
do k = 0, ndznpc-1
do i = 0, ndxhnpr-1
    matfir(1, 1:4) = (/ 1.0, 1.0, 1.0, 1.0 /)
    matfir(2, 1:4) = (/ 0.0, hp1(0), hp2(0), hp3 /)
    matfir(3, 1:4) = (/ 0.0, (hp1(0))**2.0, (hp2(0))**2.0, hp3**2 /)
    matfir(4, 1:4) = (/ 0.0, (hp1(0))**3.0, (hp2(0))**3.0, hp3**3 /)
    bfir(1:4) = (/ 0.0, 1.0, 0.0, 0.0 /)
    call gaussian_eli(4, matfir, bfir)
    secderprec(0, i, k) = bfir(1) - bfir(4) / secderpree(1,i,k) * secderpreb(1,i,k)
    secderpred(0, i, k) = bfir(2) - bfir(4) / secderpree(1,i,k) * secderprec(1,i,k)
    secderpree(0, i, k) = bfir(3) - bfir(4) / secderpree(1,i,k) * secderpred(1,i,k)
    cofm1(i,k) = bfir(4) / secderpree(1,i,k)
    secdervisc(0, i, k) = 1.0
    secdervisd(0, i, k) = 0.0
    secdervise(0, i, k) = 0.0

    matfir(1, 1:4) = (/ 1.0, 1.0, 1.0, 1.0 /)
    matfir(2, 1:4) = (/ hm3, hm2(ndy), hm1(ndy), 0.0 /)
    matfir(3, 1:4) = (/ hm3**2, (hm2(ndy))**2.0, (hm1(ndy))**2.0, 0.0 /)
    matfir(4, 1:4) = (/ hm3**3, (hm2(ndy))**3.0, (hm1(ndy))**3.0, 0.0 /)
    bfir(1:4) = (/ 0.0, 1.0, 0.0, 0.0 /)
    call gaussian_eli(4, matfir, bfir)
    secderprea(ndy, i, k) = bfir(2) - bfir(1) / secderprea(ndy-1,i,k) * secderpreb(ndy-1,i,k)
    secderpreb(ndy, i, k) = bfir(3) - bfir(1) / secderprea(ndy-1,i,k) * secderprec(ndy-1,i,k)
    secderprec(ndy, i, k) = bfir(4) - bfir(1) / secderprea(ndy-1,i,k) * secderpred(ndy-1,i,k)
    cofm2(i,k) = bfir(1) / secderprea(ndy-1,i,k)
    secdervisa(ndy, i, k) = 0.0
    secdervisb(ndy, i, k) = 0.0
    secdervisc(ndy, i, k) = 1.0
enddo
enddo

! LU decomposition of penta-diagonal matrix into lower triangle (zt,gm,af) and upper triangle (1,bt,qt)
do k = 0, ndznpc-1
do i = 0, ndxhnpr-1
    
    afp(1, i, k) = secderprec(0,i,k)
    btp(1, i, k) = secderpred(0,i,k) / afp(1,i,k)
    gmp(2, i, k) = secderpreb(1,i,k)
    afp(2, i, k) = secderprec(1,i,k) - gmp(2,i,k) * btp(1,i,k)
    do j = 1, (ndy+1)-2
        qtp(j, i, k) = secderpree(j-1,i,k) / afp(j,i,k)
        btp(j+1, i, k) = ( secderpred(j,i,k) - gmp(j+1,i,k) * qtp(j,i,k) ) / afp(j+1,i,k)
        ztp(j+2, i, k) = secderprea(j+1,i,k)
        gmp(j+2, i, k) = secderpreb(j+1,i,k) - ztp(j+2,i,k) * btp(j,i,k)
        afp(j+2, i, k) = secderprec(j+1,i,k) - ztp(j+2,i,k) * qtp(j,i,k) - gmp(j+2,i,k) * btp(j+1,i,k)
    enddo

    afv(1, i, k) = secdervisc(0,i,k)
    btv(1, i, k) = secdervisd(0,i,k) / afv(1,i,k)
    gmv(2, i, k) = secdervisb(1,i,k)
    afv(2, i, k) = secdervisc(1,i,k) - gmv(2,i,k) * btv(1,i,k)
    do j = 1, (ndy+1)-2
        qtv(j, i, k) = secdervise(j-1,i,k) / afv(j,i,k)
        btv(j+1, i, k) = ( secdervisd(j,i,k) - gmv(j+1,i,k) * qtv(j,i,k) ) / afv(j+1,i,k)
        ztv(j+2, i, k) = secdervisa(j+1,i,k)
        gmv(j+2, i, k) = secdervisb(j+1,i,k) - ztv(j+2,i,k) * btv(j,i,k)
        afv(j+2, i, k) = secdervisc(j+1,i,k) - ztv(j+2,i,k) * qtv(j,i,k) - gmv(j+2,i,k) * btv(j+1,i,k)
    enddo
    
enddo
enddo


! Dirichlet for pressure
do k = 0, ndznpc-1
do i = 0, ndxhnpr-1
    secderprec(0, i, k) = 1.0
    secderpred(0, i, k) = 0.0
    secderpree(0, i, k) = 0.0

    secderprea(ndy, i, k) = 0.0
    secderpreb(ndy, i, k) = 0.0
    secderprec(ndy, i, k) = 1.0
enddo
enddo

! LU decomposition of penta-diagonal matrix into lower triangle (zt,gm,af) and upper triangle (1,bt,qt)
do k = 0, ndznpc-1
do i = 0, ndxhnpr-1
    afp1(1, i, k) = secderprec(0,i,k)
    btp1(1, i, k) = secderpred(0,i,k) / afp1(1,i,k)
    gmp1(2, i, k) = secderpreb(1,i,k)
    afp1(2, i, k) = secderprec(1,i,k) - gmp1(2,i,k) * btp1(1,i,k)
    do j = 1, (ndy+1)-2
        qtp1(j, i, k) = secderpree(j-1,i,k) / afp1(j,i,k)
        btp1(j+1, i, k) = ( secderpred(j,i,k) - gmp1(j+1,i,k) * qtp1(j,i,k) ) / afp1(j+1,i,k)
        ztp1(j+2, i, k) = secderprea(j+1,i,k)
        gmp1(j+2, i, k) = secderpreb(j+1,i,k) - ztp1(j+2,i,k) * btp1(j,i,k)
        afp1(j+2, i, k) = secderprec(j+1,i,k) - ztp1(j+2,i,k) * qtp1(j,i,k) - gmp1(j+2,i,k) * btp1(j+1,i,k)
    enddo
enddo
enddo

!!!!!!!!!!!!!!!!! prepare coefficients for divergence modification of pressure & velocity !!!!!!!!!!!!!!!!!!!!!
call premodifyp
call premodifyv
    
end

    
! create plans for fft
subroutine prefft
use global
implicit none

call fftw_f77_create_plan( my_plan_x_b, ndx2, fftw_backward, fftw_estimate + fftw_in_place )
call fftw_f77_create_plan( my_plan_z_b, ndz2, fftw_backward, fftw_estimate + fftw_in_place )
call fftw_f77_create_plan( my_plan_x_f, ndx2, fftw_forward , fftw_estimate + fftw_in_place )
call fftw_f77_create_plan( my_plan_z_f, ndz2, fftw_forward , fftw_estimate + fftw_in_place )
call fftw_f77_create_plan( my_plan_y_f, ndy , fftw_forward , fftw_estimate + fftw_in_place )
call fftw_f77_create_plan( my_plan_y_fr,ndyr, fftw_forward , fftw_estimate + fftw_in_place )

end

    
! generate the matrix for differentiation in y direction
subroutine dy_matrix
use global
implicit none

real ddy ! matrix for calculating d/dy
common /ddy_matrix/ ddy(0:ndy, 0:ndy)

integer j, jp, jj
real c(0:ndy), g(0:ndy, 0:ndy), dg(0:ndy+1)
complex ddy_temp(0:ndy, 0:ndy)

c = 1.
c(0) = 2.
c(ndy) = 2.

do j = 0, ndy
do jp = 0, ndy
    g(jp, j) = cos( pai * real(j * jp) / real(ndy) )       
enddo
enddo

do jp = 0, ndy
    do j = 0, ndy
        g(jp, j) =  2.0d0 / ( c(jp) * c(j) * real(ndy) ) * g(jp,j)
    enddo
    dg(ndy) = 0.
    dg(ndy+1) = 0.
    
    do jj = ndy, 2, -1
        dg(jj-1) = dg(jj+1) + 2.0d0 * jj * g(jp,jj)
    enddo
    dg(0) = 0.5d0 * dg(2) + g(jp,1)
    ddy_temp(jp, 0:ndy) = dg(0:ndy)
enddo

call chebp(ddy_temp, ndy+1, ndy, 1)

ddy = real(ddy_temp)

end

    
! generate the matrix for integration in y direction
subroutine iy_matrix
use global
implicit none

real inty ! matrix for calculating int(y)
common /inty_matrix/ inty(0:ndy)

integer j, jp, jj, n
real g(0:ndy, 0:ndy), t(0:ndy, 0:ndy) , c(0:ndy)

c = 1.
c(0) = 2.
c(ndy) = 2.

do j = 0, ndy
    do jp = 0, ndy
        t(jp, j) = cos( pai * real(j * jp) / real(ndy) )
    enddo
enddo

do jp = 0, ndy
    do j = 0, ndy
        g(jp, j) = 2. / ( c(jp) * c(j) * real(ndy) ) * t(j,jp)
    enddo
enddo
       
do jp = 0, ndy
    inty(jp) = 0.
    do j = 0, ndy, 2
        inty(jp) = inty(jp) + g(jp, j) * 2. / real(1 - j*j)
    enddo
enddo

end

