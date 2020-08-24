program turbulentchannel
use global
implicit none

integer impi_numprocs, mpi_comm_cart

call mpi_init( impi_errorinfo )
call mpi_comm_rank( mpi_comm_world, impi_myid, impi_errorinfo )
call mpi_comm_size( mpi_comm_world, impi_numprocs, impi_errorinfo )
call mpi_cart_create( mpi_comm_world, 2, (/npr,npc/), (/.false.,.false./), .false., mpi_comm_cart, impi_errorinfo )
call mpi_cart_coords( mpi_comm_cart, impi_myid, 2, coord, impi_errorinfo )
call mpi_cart_sub( mpi_comm_cart, (/.true.,.false./), mpi_comm_col, impi_errorinfo ) ! create sub communication domains for cores in the same columns
call mpi_cart_sub( mpi_comm_cart, (/.false.,.true./), mpi_comm_row, impi_errorinfo ) ! create sub communication domains for cores in the same rows

! the number of process is right
if ( np .eq. impi_numprocs ) then
    call mpi_barrier( mpi_comm_world, impi_errorinfo )
    call readinitialfile ! read XINDAT
    call initializevariable ! set constant coefficients such as finite difference coefficients
    !call premodifyp ! prepare coefficients needed in divergence modification of pressure
    !call premodifyv ! prepare coefficients needed in divergence modification of velocity
	call mainprogram ! main loop
endif

write(*,*) 'Computing finished !'

call mpi_barrier( mpi_comm_world, impi_errorinfo )
call mpi_finalize( impi_errorinfo )

end


subroutine mainprogram
use global
implicit none

real time0, time
integer iwt, ikt

!-------------------------------------------------------------------------------
!    execute statement
!-------------------------------------------------------------------------------
if (ifc .eq. 0) call inifield0(time0) ! initiate field from laminar
if (ifc .eq. 1) call inifield1(time0) ! initiate field from a velocity field
if (ifc .eq. 2) call inifieldm(time0) ! initiate field from mid files for continuing computation

write(*,*) 'Core ', impi_myid, 'init finished.'

do ikt = 1, nt
    
    time = time0 + ikt * dt
    iwt  = int(anint( time / dt )) ! in continuing computing, iwt and time resume from last time, otherwise time0 = 0 and iwt = ikt
    
    if ( (impi_myid .eq. 0) .and. (mod(iwt, 100) .eq. 0) ) write(*,*) 'Computing step ', iwt
    
    call vor
    call utp
    call pb1condition
    call fp
    call bulkforce
    call pb2condition
    call mari1
    call mari2new
    call mari3new
    call wrfr
    call modifydiv
    
    if (iff .eq. 1) then
        call frcorr
	    call hrcorr
    endif

	call bdisplace
    
    call write_files(iwt, ikt, time)
    
    
enddo

end
