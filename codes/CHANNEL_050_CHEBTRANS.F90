
! transfer from physical space to Chebyshev spectral space    
subroutine chebs(a, nx, ny, nz)
use global
implicit none

integer nx, ny, nz
complex a(nx, 0:ny, nz)
complex b(nx, 0:ny-1)
real urr(nx, 0:ny), ucr(nx, 0:ny)
real uri(nx, 0:ny), uci(nx, 0:ny)
real c01(nx), c02(nx), cn1(nx), cn2(nx)
real xni, pin, ssi, alj, glj
integer nh, l2n, l2nm1, l2np1
complex temp(ny)
complex temp1(ny)
integer i, k, iz

xni= 1. / ny
nh = ny / 2
pin= pai / ny            

do iz = 1, nz
    b = 0.
    urr(:, :) = real( a(:,:,iz) )
    uri(:, :) = imag( a(:,:,iz) )
    b(:, 0) = cmplx( urr(:, 0), uri(:, 0) )
    b(:,nh) = cmplx( urr(:,ny), uri(:,ny) )
    c01(:) = (0.5 * ( urr(:,0 ) + urr(:,ny) ) ) + urr(:,ny-1 )
    cn1(:) = (0.5 * ( urr(:,0 ) + urr(:,ny) ) ) - urr(:,ny-1 )  
    c02(:) = (0.5 * ( uri(:,0 ) + uri(:,ny) ) ) + uri(:,ny-1 )
    cn2(:) = (0.5 * ( uri(:,0 ) + uri(:,ny) ) ) - uri(:,ny-1 )   

    l2n = 0
    do k = 1, nh-1
        l2n   = l2n + 2
        l2nm1 = l2n - 1
        l2np1 = l2n + 1   
        b(:, ny-k) = cmplx( urr(:,l2n) - uri(:,l2np1) + uri(:,l2nm1), uri(:,l2n) + urr(:,l2np1) - urr(:,l2nm1) )  
    enddo
   
    l2n = ny
    do k = nh+1 , ny-1
        l2n   = l2n - 2
        l2nm1 = l2n + 1
        l2np1 = l2n - 1  
        b(:, ny-k) = cmplx( urr(:,l2n) - uri(:,l2np1) + uri(:,l2nm1), uri(:,l2n) + urr(:,l2np1) - urr(:,l2nm1) ) 
    enddo

    do k = 1, ny-2, 2   
        c01(:) = c01(:) + urr(:,k) + urr(:,k+1)
        c02(:) = c02(:) + uri(:,k) + uri(:,k+1) 
        cn1(:) = cn1(:) - urr(:,k) + urr(:,k+1)
        cn2(:) = cn2(:) - uri(:,k) + uri(:,k+1)         
    enddo 

    do i = 1, nx
        temp(1:ny) = b(i, 0:ny-1)
        call fftw_f77(my_plan_y_f, 1, temp, 1, ny, temp1, 0, 0)
        b(i, 0:ny-1) = temp(1:ny)
    enddo
    do k = 1, nh - 1
        ssi = 0.25 / sin( k * pin )
        alj = xni * ( 0.5 + ssi )
        glj = xni * ( 0.5 - ssi ) 

        ucr(:,    k ) = alj * real(b(:,   k)) + glj * real(b(:,ny-k))
        ucr(:, ny-k ) = alj * real(b(:,ny-k)) + glj * real(b(:,   k))         
        uci(:,    k ) = alj * imag(b(:,   k)) + glj * imag(b(:,ny-k))      
        uci(:, ny-k ) = alj * imag(b(:,ny-k)) + glj * imag(b(:,   k))         
    enddo

    ucr(:,nh) = xni * real(b(:,nh))
    uci(:,nh) = xni * imag(b(:,nh))          
    ucr(:, 0) = xni * c01(:)
    uci(:, 0) = xni * c02(:)          
    ucr(:,ny) = xni * cn1(:)          
    uci(:,ny) = xni * cn2(:)  

    a(:, :, iz) = ucr(:,:) + (0., 1.) * uci(:,:)
enddo

end


! transfer from Chebyshev spectral space to physical space
subroutine chebp(a, nx, ny, nz)
use global
implicit none

integer nx, ny, nz
complex a(nx, 0:ny, nz)
complex b(nx, 0:ny-1)
real urr(nx, 0:ny), ucr(nx, 0:ny)
real uri(nx, 0:ny), uci(nx, 0:ny)
real u01(nx), u02(nx), un1(nx), un2(nx)
real xni, pin, ssi, alj, glj
integer nh, l2n, l2nm1, l2np1
complex temp(ny)
complex temp1(ny)
integer i, k, iz

nh = ny / 2 
pin= pai / ny

do iz = 1, nz 
    ucr(:, :) = real( a(:,:,iz) )
    uci(:, :) = imag( a(:,:,iz) )

    u01(:) = (ucr(:,0) + ucr(:,ny)) + ucr(:,ny-1)
    un1(:) = (ucr(:,0) + ucr(:,ny)) - ucr(:,ny-1)
    u02(:) = (uci(:,0) + uci(:,ny)) + uci(:,ny-1)
    un2(:) = (uci(:,0) + uci(:,ny)) - uci(:,ny-1)          

    b(:, 0)  = cmplx( 2. * ucr(:, 0) , 2. * uci(:, 0) )
    b(:, nh) = cmplx( 2. * ucr(:,ny) , 2. * uci(:,ny) ) 

    do k = 1, ny-2, 2
        u01(:) = u01(:) + ucr(:,k) + ucr(:,k+1)
        u02(:) = u02(:) + uci(:,k) + uci(:,k+1)      
        un1(:) = un1(:) - ucr(:,k) + ucr(:,k+1)
        un2(:) = un2(:) - uci(:,k) + uci(:,k+1)      
    enddo

    l2n = 0
    do k = 1, nh-1
        l2n   = l2n + 2
        l2nm1 = l2n - 1
        l2np1 = l2n + 1    
        b(:, ny-k) = cmplx( ucr(:,l2n) - uci(:,l2np1) + uci(:,l2nm1), uci(:,l2n) + ucr(:,l2np1) - ucr(:,l2nm1) )           
    enddo
    
    l2n = ny
    do k = nh+1, ny-1
        l2n   = l2n - 2
        l2nm1 = l2n + 1
        l2np1 = l2n - 1  
        b(:, ny-k) = cmplx( ucr(:,l2n) - uci(:,l2np1) + uci(:,l2nm1), uci(:,l2n) + ucr(:,l2np1) - ucr(:,l2nm1) )              
    enddo

    do i = 1 , nx
        temp(1:ny) = b(i,0:ny-1)
        call fftw_f77(my_plan_y_f, 1, temp, 1, ny, temp1, 0, 0)
        b(i, 0:ny-1) = temp(1:ny)
    enddo

    do k = 1, nh-1
        ssi = 0.125 / sin( k * pin )
        alj = 0.25 + ssi
        glj = 0.25 - ssi    
        urr( : ,      k ) = alj * real(b(:,   k)) + glj * real(b(:,ny-k))
        urr( : , ny - k ) = alj * real(b(:,ny-k)) + glj * real(b(:,   k))           
        uri( : ,      k ) = alj * imag(b(:,   k)) + glj * imag(b(:,ny-k))          
        uri( : , ny - k ) = alj * imag(b(:,ny-k)) + glj * imag(b(:,   k))     
    enddo
   
    urr(:, nh) = 0.5 * real(b(:,nh))
    uri(:, nh) = 0.5 * imag(b(:,nh))          

    urr(:, 0) = u01(:)
    urr(:,ny) = un1(:)                    

    uri(:, 0) = u02(:)          
    uri(:,ny) = un2(:)     

    a(:, :, iz) = cmplx( urr(:,:) , uri(:,:) )
enddo

end 


! transfer from physical space to Chebyshev spectral space for the length of ndyr (relevant in fft events)
subroutine chebsr(a, nx, ny, nz)
use global
implicit none

integer nx, ny, nz
complex a(nx, 0:ny, nz)
complex b(nx, 0:ny-1)
real urr(nx, 0:ny), ucr(nx, 0:ny)
real uri(nx, 0:ny), uci(nx, 0:ny)
real c01(nx), c02(nx), cn1(nx), cn2(nx)
real xni, pin, ssi, alj, glj
integer nh, l2n, l2nm1, l2np1
complex temp(ny)
complex temp1(ny)
integer i, k, iz

xni = 1. / ny
nh  = ny / 2
pin = pai / ny

do iz = 1, nz
    b = 0.
    urr(:, :) = real( a(:,:,iz) )
    uri(:, :) = imag( a(:,:,iz) )
    b(:, 0) = cmplx( urr(:, 0), uri(:, 0) )
    b(:,nh) = cmplx( urr(:,ny), uri(:,ny) )
    c01(:) = (0.5 * ( urr(:,0) + urr(:,ny) ) ) + urr(:,ny-1)
    cn1(:) = (0.5 * ( urr(:,0) + urr(:,ny) ) ) - urr(:,ny-1)  
    c02(:) = (0.5 * ( uri(:,0) + uri(:,ny) ) ) + uri(:,ny-1)
    cn2(:) = (0.5 * ( uri(:,0) + uri(:,ny) ) ) - uri(:,ny-1)   

    l2n = 0
    do k = 1, nh-1
        l2n   = l2n + 2
        l2nm1 = l2n - 1
        l2np1 = l2n + 1   
        b(:, ny-k) = cmplx( urr(:,l2n) - uri(:,l2np1) + uri(:,l2nm1), uri(:,l2n) + urr(:,l2np1) - urr(:,l2nm1) )  
    enddo
   
    l2n = ny
    do k = nh+1, ny-1
        l2n   = l2n - 2
        l2nm1 = l2n + 1
        l2np1 = l2n - 1  
        b(:, ny-k) = cmplx( urr(:,l2n) - uri(:,l2np1) + uri(:,l2nm1), uri(:,l2n) + urr(:,l2np1) - urr(:,l2nm1) ) 
    enddo

    do k = 1, ny-2, 2   
        c01(:) = c01(:) + urr(:,k) + urr(:,k+1)
        c02(:) = c02(:) + uri(:,k) + uri(:,k+1) 
        cn1(:) = cn1(:) - urr(:,k) + urr(:,k+1)
        cn2(:) = cn2(:) - uri(:,k) + uri(:,k+1)         
    enddo 

    do i = 1, nx
        temp(1:ny) = b(i,0:ny-1)
        call fftw_f77(my_plan_y_fr, 1, temp, 1, ny, temp1, 0, 0)
        b(i, 0:ny-1) = temp(1:ny)
    enddo
    do k = 1, nh-1
        ssi = 0.25 / sin( k * pin )
        alj = xni * ( 0.5 + ssi )
        glj = xni * ( 0.5 - ssi ) 

        ucr( : ,      k ) = alj * real(b( : ,      k)) + glj * real(b( : , ny - k))
        ucr( : , ny - k ) = alj * real(b( : , ny - k)) + glj * real(b( : ,      k))         
        uci( : ,      k ) = alj * imag(b( : ,      k)) + glj * imag(b( : , ny - k))      
        uci( : , ny - k ) = alj * imag(b( : , ny - k)) + glj * imag(b( : ,      k))         
    enddo

    ucr(:,nh) = xni * real(b(:,nh))
    uci(:,nh) = xni * imag(b(:,nh))          
    ucr(:, 0) = xni * c01(:)
    uci(:, 0) = xni * c02(:)          
    ucr(:,ny) = xni * cn1(:)          
    uci(:,ny) = xni * cn2(:)  

    a(:, :, iz) = ucr(:,:) + (0., 1.) * uci(:,:)
enddo

    end


! finite difference in y direction: f = da/dy, equation set [firder_abc][f] = [firrhsc][a]
subroutine diffy(a, f)
use global
implicit none

complex a(0:ndy), f(0:ndy)
integer j

do j = 0, ndy
    select case(j)
    case(0)
        f(j) = firrhsc(j,3) * a(j) + firrhsc(j,4) * a(j+1) + firrhsc(j,5) * a(j+2)
    case(1)
        f(j) = firrhsc(j,2) * a(j-1) + firrhsc(j,3) * a(j) + firrhsc(j,4) * a(j+1)
    case(ndy-1)
        f(j) = firrhsc(j,2) * a(j-1) + firrhsc(j,3) * a(j) + firrhsc(j,4) * a(j+1)
    case(ndy)
        f(j) = firrhsc(j,1) * a(j-2) + firrhsc(j,2) * a(j-1) + firrhsc(j,3) * a(j)
    case default
        f(j) = firrhsc(j,1) * a(j-2) + firrhsc(j,2) * a(j-1) + firrhsc(j,3) * a(j) + firrhsc(j,4) * a(j+1) + firrhsc(j,5) * a(j+2)
    end select
enddo

call triold(ndy+1, firdera, firderb, firderc, f)

end


! integral in y direction: f = \int_0^y{ a }, equation set [firrhsc][f] = [firder_abc][a]
subroutine integy(a, f)
use global
implicit none

complex a(0:ndy), f(0:ndy)
integer j

do j = 0, ndy
    select case(j)
    case(0)
        f(j) = a(j) * firderb(j) + a(j+1) * firderc(j)
    case(ndy)
        f(j) = a(j) * firderb(j) + a(j-1) * firdera(j)
    case default
        f(j) = a(j-1) * firdera(j) + a(j) * firderb(j) + a(j+1) * firderc(j)
    end select
enddo

call pentaold(ndy+1, firrhsc(2:ndy, 1), firrhsc(1:ndy, 2), firrhsc(0:ndy, 3), firrhsc(0:ndy-1, 4), firrhsc(0:ndy-2, 5), f)

f = f - f(0)

end
    
! compute convolution using pseudo spectral method, applying 3/2 rule for de-aliasing
subroutine compute_convol(a, b, ab)
use global
implicit none

complex &
    a (0:ndxhnpr-1, 0:ndy, 0:ndznpc-1), &
    b (0:ndxhnpr-1, 0:ndy, 0:ndznpc-1), &
    ab(0:ndxhnpr-1, 0:ndy, 0:ndznpc-1)

complex tempz1(0:ndxhnpr-1, 0:ndynpc, 0:ndz2-1)
complex tempz2(0:ndxhnpr-1, 0:ndynpc, 0:ndz2-1)
complex tempz3(0:ndxhnpr-1, 0:ndynpc, 0:ndz2-1)

complex tempx1(0:ndx2-1, 0:ndynpc, 0:ndz2npr-1)
complex tempx2(0:ndx2-1, 0:ndynpc, 0:ndz2npr-1)
complex tempx3(0:ndx2-1, 0:ndynpc, 0:ndz2npr-1)

complex tempz(0:ndz2-1), tempz0(0:ndz2-1)
complex tempx(0:ndx2-1), tempx0(0:ndx2-1)

integer i, j, k

tempx1 = 0.
tempx2 = 0.
tempx3 = 0.

tempz1 = 0.
tempz2 = 0.
tempz3 = 0.

!!! perform ifft from spectral space to physical space
! transpose division of field to ensure integrity in z direction
call transposeytoz( a , tempz1(0:ndxhnpr-1, 0:ndynpc, 0:ndz-1) )
call transposeytoz( b , tempz2(0:ndxhnpr-1, 0:ndynpc, 0:ndz-1) )
! applying 3/2 rule
tempz1(:, :, ndz+1:ndz2-1) = tempz1(:, :, ndzh+1:ndz-1)
tempz2(:, :, ndz+1:ndz2-1) = tempz2(:, :, ndzh+1:ndz-1)
tempz1(:, :, ndzh:ndz-1) = 0.
tempz2(:, :, ndzh:ndz-1) = 0.
! do ifft in z direction
do j = 0, ndynpc
do i = 0, ndxhnpr-1
    tempz = tempz1(i, j, :)
    call fftw_f77(my_plan_z_f, 1, tempz, 1, ndz2, tempz0, 0, 0)
    tempz1(i, j, :) = tempz
enddo
enddo
do j = 0, ndynpc
do i = 0, ndxhnpr-1
    tempz = tempz2(i, j, :)
    call fftw_f77(my_plan_z_f, 1, tempz, 1, ndz2, tempz0, 0, 0)
    tempz2(i, j, :) = tempz
enddo
enddo
! transpose division of field to ensure integrity in x direction
call transposeztox( tempz1, tempx1(0:ndxh-1, 0:ndynpc, 0:ndz2npr-1) )
call transposeztox( tempz2, tempx2(0:ndxh-1, 0:ndynpc, 0:ndz2npr-1) )
! complement the other half in x direction
tempx1(ndx2-1:ndx+1:-1, :, :) = conjg( tempx1(1:ndxh-1, :, :) )
tempx2(ndx2-1:ndx+1:-1, :, :) = conjg( tempx2(1:ndxh-1, :, :) )
! do ifft in x direction
do k = 0, ndz2npr-1
do j = 0, ndynpc
    tempx = tempx1(:, j, k)
    call fftw_f77(my_plan_x_f, 1, tempx, 1,ndx2, tempx0, 0, 0)
    tempx1(:, j, k) = tempx
enddo
enddo
do k = 0, ndz2npr-1
do j = 0, ndynpc
    tempx = tempx2(:, j, k)
    call fftw_f77(my_plan_x_f, 1, tempx, 1, ndx2, tempx0, 0, 0)
    tempx2(:, j, k) = tempx
enddo
enddo

!!! compute product in physical space
tempx3 = real(tempx1) * real(tempx2)

!!! perform fft from physical space to spectral space
! do fft in x direction
do k = 0, ndz2npr-1
do j = 0, ndynpc
    tempx = tempx3(:, j, k)
    call fftw_f77(my_plan_x_b, 1, tempx, 1, ndx2, tempx0, 0, 0)
    tempx3(0:ndxh-1, j, k) = tempx(0:ndxh-1) / ndx2
enddo
enddo
! transpose division of field to ensure integrity in z direction
call transposextoz( tempx3(0:ndxh-1, 0:ndynpc, 0:ndz2npr-1), tempz3(0:ndxhnpr-1, 0:ndynpc, 0:ndz2-1) )
! do fft in z direction
do j = 0, ndynpc
do i = 0, ndxhnpr-1
    tempz = tempz3(i, j, :)
    call fftw_f77(my_plan_z_b, 1, tempz, 1, ndz2, tempz0, 0, 0)
    tempz3(i, j, 0:ndzh-1) = tempz(0:ndzh-1) / ndz2
    tempz3(i, j, ndzh+1:ndz-1) = tempz(ndz+1:ndz2-1) / ndz2
    tempz3(i, j, ndzh) = 0.
enddo
enddo
! transpose division of field back to normal (integral in y direction)
call transposeztoy( tempz3(0:ndxhnpr-1, 0:ndynpc, 0:ndz-1), ab(0:ndxhnpr-1, 0:ndy, 0:ndznpc-1) )

end