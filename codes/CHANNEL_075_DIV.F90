
subroutine div
use global
implicit none

integer i, j, k, ii, kk
complex rhsv(0:ndy)
integer maxdivi, maxdivj, maxdivk
real localdiv, maxdiv
integer recvbuf_maxdivi(np), recvbuf_maxdivj(np), recvbuf_maxdivk(np)
real recvbuf_maxdiv(np)
common /maxdiv/ maxdiv
common /maxdivpos/ maxdivi, maxdivj, maxdivk

maxdiv = 0.0

do k = 0, ndznpc-1
    kk = - ( coord(1) * ndznpc + k )
    if (kk .lt. -ndzh) kk = kk + ndz
    do i = 0, ndxhnpr-1
        ii = - ( coord(0) * ndxhnpr + i )
            
        call diffy( v(i,:,k), rhsv )
            
        do j = 0, ndy
            localdiv = abs( (0.,1.) * arf * ii * u(i,j,k) + rhsv(j) + (0.,1.) * bat * kk * w(i,j,k) )
            if ( localdiv .gt. maxdiv ) then
                maxdiv = localdiv
                maxdivi = coord(0) * ndxhnpr + i
                maxdivj = j
                maxdivk = coord(1) * ndznpc + k
            endif
        enddo
    enddo
enddo

! gather maxdiv from other cores and choose the largest one
call mpi_barrier( mpi_comm_world, impi_errorinfo )
call mpi_gather(maxdiv, 1, mpi_real, recvbuf_maxdiv, 1, mpi_real, 0, mpi_comm_world, impi_errorinfo)
call mpi_gather(maxdivi, 1, mpi_integer, recvbuf_maxdivi, 1, mpi_integer, 0, mpi_comm_world, impi_errorinfo)
call mpi_gather(maxdivj, 1, mpi_integer, recvbuf_maxdivj, 1, mpi_integer, 0, mpi_comm_world, impi_errorinfo)
call mpi_gather(maxdivk, 1, mpi_integer, recvbuf_maxdivk, 1, mpi_integer, 0, mpi_comm_world, impi_errorinfo)

if (impi_myid .eq. 0) then
    do i = 1, np
        if ( recvbuf_maxdiv(i) .gt. maxdiv ) then
            maxdiv = recvbuf_maxdiv(i)
            maxdivi = recvbuf_maxdivi(i)
            maxdivj = recvbuf_maxdivj(i)
            maxdivk = recvbuf_maxdivk(i)
        endif
    enddo
endif

end subroutine