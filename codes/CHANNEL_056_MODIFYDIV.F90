
subroutine modifydiv
use global
implicit none

!complex umm, vmm, wmm, ump, vmp, wmp
!complex pplus, pminus
!common /velomodify_full/ &
!    umm (0:ndxhnpr-1, 0:ndy, 0:ndznpc-1), &
!    vmm (0:ndxhnpr-1, 0:ndy, 0:ndznpc-1), &
!    wmm (0:ndxhnpr-1, 0:ndy, 0:ndznpc-1), &
!    ump (0:ndxhnpr-1, 0:ndy, 0:ndznpc-1), &
!    vmp (0:ndxhnpr-1, 0:ndy, 0:ndznpc-1), &
!    wmp (0:ndxhnpr-1, 0:ndy, 0:ndznpc-1)
!common /ppm/ &
!    pplus (0:ndxhnpr-1, 0:ndy, 0:ndznpc-1), &
!    pminus(0:ndxhnpr-1, 0:ndy, 0:ndznpc-1)

integer i, j, k
complex uc(0:ndy), vc(0:ndy), wc(0:ndy), pc(0:ndy)
complex localdivd, localdivu, localdivdm, localdivum, localdivdp, localdivup
complex aa, bb, cc, dd, ee, ff, deltap, deltam


do k = 0, ndznpc-1
do i = 0, ndxhnpr-1
    call modifydiv_coeff(i, k, localdivd, localdivu, localdivdm, localdivum, localdivdp, localdivup)
    
    aa = localdivup
    bb = localdivum
    cc = localdivdp
    dd = localdivdm
    ee =-localdivu
    ff =-localdivd

    deltap = (dd*ee - bb*ff) / (aa*dd - bb*cc)
    deltam = (aa*ff - cc*ee) / (aa*dd - bb*cc)

    pc = deltap * pplus(i,:,k) + deltam * pminus(i,:,k)
    uc = deltap * ump(i,:,k) + deltam * umm(i,:,k)
    vc = deltap * vmp(i,:,k) + deltam * vmm(i,:,k)
    wc = deltap * wmp(i,:,k) + deltam * wmm(i,:,k)

    u(i,:,k) = u(i,:,k) + uc
    v(i,:,k) = v(i,:,k) + vc
    w(i,:,k) = w(i,:,k) + wc
    p(i,:,k) = p(i,:,k) + pc
enddo
enddo

endsubroutine


subroutine modifydiv_coeff(i, k, localdivd, localdivu, localdivdm, localdivum, localdivdp, localdivup)
use global
implicit none

!complex umm, vmm, wmm, ump, vmp, wmp
!common /velomodify_full/ &
!    umm (0:ndxhnpr-1, 0:ndy, 0:ndznpc-1), &
!    vmm (0:ndxhnpr-1, 0:ndy, 0:ndznpc-1), &
!    wmm (0:ndxhnpr-1, 0:ndy, 0:ndznpc-1), &
!    ump (0:ndxhnpr-1, 0:ndy, 0:ndznpc-1), &
!    vmp (0:ndxhnpr-1, 0:ndy, 0:ndznpc-1), &
!    wmp (0:ndxhnpr-1, 0:ndy, 0:ndznpc-1)  

integer i, j, k, ii, kk
complex rhsv(0:ndy)
complex &
    localdivdp, localdivup, &
    localdivdm, localdivum, &
    localdivd , localdivu

kk = - ( coord(1) * ndznpc + k )
if (kk .lt. -ndzh) kk = kk + ndz
ii = - ( coord(0) * ndxhnpr + i )

call diffy( v(i,:,k), rhsv )
localdivd = (0.,1.) * arf * ii * u(i, 0 ,k) + (0.,1.) * bat * kk * w(i, 0 ,k) + rhsv(0)
localdivu = (0.,1.) * arf * ii * u(i,ndy,k) + (0.,1.) * bat * kk * w(i,ndy,k) + rhsv(ndy)

call diffy( vmm(i,:,k), rhsv )
localdivdm = (0.,1.) * arf * ii * umm(i, 0 ,k) + (0.,1.) * bat * kk * wmm(i, 0 ,k) + rhsv(0)
localdivum = (0.,1.) * arf * ii * umm(i,ndy,k) + (0.,1.) * bat * kk * wmm(i,ndy,k) + rhsv(ndy)

call diffy( vmp(i,:,k), rhsv )
localdivdp = (0.,1.) * arf * ii * ump(i, 0 ,k) + (0.,1.) * bat * kk * wmp(i, 0 ,k) + rhsv(0)
localdivup = (0.,1.) * arf * ii * ump(i,ndy,k) + (0.,1.) * bat * kk * wmp(i,ndy,k) + rhsv(ndy)

end subroutine