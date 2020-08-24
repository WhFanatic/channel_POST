
! change the directions of parallelization from (X & Z) to (X & Y)
subroutine transposeytoz(array_in, array_out)
use global
implicit none

complex array_in (0:ndxhnpr-1, 0:ndy, 0:ndznpc-1)
complex array_out(0:ndxhnpr-1, 0:ndynpc, 0:ndz-1)
complex temp1(0:ndxhnpr-1, 0:ndynpc, 0:ndznpc-1, 0:npc-1)
complex temp2(0:ndxhnpr-1, 0:ndynpc, 0:ndznpc-1, 0:npc-1)
integer i, st, en

st = 0
en = ndynpc
do i = 0, npc-1
    temp1(:, :, :, i) = array_in(:,st:en,:)
    st = en ! ndy must be exactly divided by npc, otherwise here's gonna be a problem
    en = st + ndynpc
enddo

call mpi_alltoall( &
    temp1, x1count, mpi_complex, &
    temp2, y1count, mpi_complex, &
    mpi_comm_row, impi_errorinfo )

st = 0
en = ndznpc-1
do i = 0, npc-1
    array_out(:, :, st:en) = temp2(:,:,:,i)
    st = en + 1
    en = st + ndznpc-1
enddo

end


! change the directions of parallelization from (X & Y) back to (X & Z)
subroutine transposeztoy(array_in, array_out)
use global
implicit none

complex array_in(0:ndxhnpr-1, 0:ndynpc, 0:ndz-1)
complex array_out(0:ndxhnpr-1, 0:ndy, 0:ndznpc-1)
complex temp1(0:ndxhnpr-1, 0:ndynpc, 0:ndznpc-1, 0:npc-1)
complex temp2(0:ndxhnpr-1, 0:ndynpc, 0:ndznpc-1, 0:npc-1)
integer i, st, en

st = 0
en = ndznpc-1
do i = 0, npc-1
    temp1(:, :, :, i) = array_in(:,:,st:en) 
    st = en + 1
    en = st + ndznpc-1
enddo

call mpi_alltoall( &
    temp1, y1count, mpi_complex, &
    temp2, x1count, mpi_complex, &
    mpi_comm_row, impi_errorinfo )

st = 0
en = ndynpc
do i= 0, npc-1
    array_out(:, st:en, :) = temp2(:,:,:,i) 
    st = en ! ndy must be exactly divided by npc, otherwise here's gonna be a problem
    en = st + ndynpc
enddo

end

    
! change the directions of parallelization from (X & Y) to (Z & Y)
subroutine transposeztox(array_in, array_out)
use global
implicit none

complex array_in (0:ndxhnpr-1, 0:ndynpc, 0:ndz2-1)
complex array_out(0:ndxh-1, 0:ndynpc, 0:ndz2npr-1)
complex temp1(0:ndxhnpr-1, 0:ndynpc, 0:ndz2npr-1, 0:npr-1)
complex temp2(0:ndxhnpr-1, 0:ndynpc, 0:ndz2npr-1, 0:npr-1)
integer i, st, en

st = 0
en = ndz2npr-1
do i = 0, npr-1
    temp1(:, :, :, i) = array_in(:,:,st:en) 
    st = en + 1
    en = st + ndz2npr-1
enddo

call mpi_alltoall( &
    temp1, x2count, mpi_complex, &
    temp2, z2count, mpi_complex, &
    mpi_comm_col, impi_errorinfo )

st = 0
en = ndxhnpr-1
do i = 0, npr-1
    array_out(st:en, :, :) = temp2(:,:,:,i) 
    st = en + 1
    en = st + ndxhnpr-1
enddo

end


! change the directions of parallelization from (Z & Y) back to (X & Y)
subroutine transposextoz(array_in, array_out)
use global
implicit none

complex array_in (0:ndxh-1, 0:ndynpc, 0:ndz2npr-1)
complex array_out(0:ndxhnpr-1, 0:ndynpc, 0:ndz2-1)
complex temp1(0:ndxhnpr-1, 0:ndynpc, 0:ndz2npr-1, 0:npr-1)
complex temp2(0:ndxhnpr-1, 0:ndynpc, 0:ndz2npr-1, 0:npr-1)
integer i, st, en

st = 0
en = ndxhnpr-1
do i= 0, npr-1
    temp1(:, :, :, i) = array_in(st:en,:,:)   
    st = en + 1
    en = st + ndxhnpr-1
enddo

call mpi_alltoall( &
    temp1, x2count, mpi_complex, &
    temp2, z2count, mpi_complex, &
    mpi_comm_col, impi_errorinfo )

st = 0
en = ndz2npr-1
do i = 0, npr-1
    array_out(:, :, st:en) = temp2(:,:,:,i) 
    st = en + 1
    en = st + ndz2npr-1
enddo

end