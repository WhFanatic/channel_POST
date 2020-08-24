! gaussian elimination: n matrix a and n vector b, solve b\a and return the result in b
subroutine gaussian_eli(n,a,b)

integer n,k,i,j
real a(1:n,1:n), b(1:n), l(1:n,1:n), solu(1:n), sum

do k = 1, n-1
    do i = k+1, n
        l(i,k) = a(i,k) / a(k,k)
    enddo

    do j = 1, k
        do i = k+1, n
            a(i,j) = 0.0
        enddo
    enddo 
  
    do i = k+1, n
        do j = k+1, n
            a(i,j) = a(i,j) - l(i,k) * a(k,j)
        enddo
        b(i) = b(i) - l(i,k) * b(k)
    enddo 
enddo

solu(n) = b(n) / a(n,n)

do i = n-1, 1, -1
    sum = 0.0
    do j = i+1, n
        sum = sum + a(i,j) * solu(j)
    enddo
    solu(i) = (b(i) - sum) / a(i,i)
enddo

b = solu

end subroutine


! solve penta equation with penta diagonal matrix coefficients directly: first do LU decomposition, then solve 2 equations
subroutine pentaold (n, a, b, c, d, e, f)
integer n
real a(3:n), b(2:n), c(1:n), d(1:n-1), e(1:n-2)
complex f(1:n)
real af(1:n), bt(1:n-1), gm(2:n), qt(1:n-2), zt(3:n)
complex y(1:n), solu(1:n)
integer i

af(1) = c(1)
bt(1) = d(1) / af(1)
gm(2) = b(2)
af(2) = c(2) - gm(2) * bt(1)
 
do i = 1, n-2
    qt(i) = e(i) / af(i)
    bt(i+1) = ( d(i+1) - gm(i+1) * qt(i) ) / af(i+1)
    zt(i+2) = a(i+2)
    gm(i+2) = b(i+2) - zt(i+2) * bt(i)
    af(i+2) = c(i+2) - zt(i+2) * qt(i) - gm(i+2) * bt(i+1)
enddo
 
y(1) = f(1) / af(1)
y(2) =(f(2) - gm(2) * y(1)) / af(2)
do i = 3, n
    y(i) = ( f(i) - zt(i) * y(i-2) - gm(i) * y(i-1) ) / af(i)
enddo
 
solu(n) = y(n)
solu(n-1) = y(n-1) - bt(n-1) * solu(n)
do i = n-2, 1, -1
    solu(i) = y(i) - qt(i) * solu(i+2) - bt(i) * solu(i+1)
enddo

f = solu

end subroutine


! solve penta equation that is pre-LU-decomposed
subroutine penta (n, zt, gm, af, bt, qt, f)
integer n
complex f(1:n)
real af(1:n), bt(1:n-1), gm(2:n), qt(1:n-2), zt(3:n)
complex y(1:n), solu(1:n)
integer i
 
y(1) = f(1) / af(1)
y(2) = ( f(2) - gm(2) * y(1) ) / af(2)
do i = 3, n
    y(i) = ( f(i) - zt(i) * y(i-2) - gm(i) * y(i-1) ) / af(i)
enddo

solu(n) = y(n)
solu(n-1) = y(n-1) - bt(n-1) * solu(n)
do i = n-2, 1, -1
    solu(i) = y(i) - qt(i) * solu(i+2) - bt(i) * solu(i+1)
enddo

f = solu

end subroutine


! solve tridiagonal matrix using chasing method
subroutine triold (n, a, b, c, f)
integer n
real a(2:n), b(1:n), c(1:n-1), bt(1:n-1)
complex f(1:n)
complex y(1:n), solu(1:n)
integer i
 
bt(1) = c(1) / b(1)
do i = 2, n-1
    bt(i) = c(i) / ( b(i) - a(i) * bt(i-1) )
enddo

y(1) = f(1) / b(1)
do i = 2, n
    y(i) = ( f(i) - a(i) * y(i-1) ) / ( b(i) - a(i) * bt(i-1) )
enddo

solu(n) = y(n)
do i = n-1, 1, -1
    solu(i) = y(i) - bt(i) * solu(i+1)
enddo
f = solu
end subroutine


!SUBROUTINE TRITRI(N,A,B,C,F)
! INTEGER N
! REAL FLT(2:N),UT(1:N)
! REAL A(2:N),B(1:N),C(1:N-1)
! COMPLEX F(1:N)
! COMPLEX Y(1:N)
! INTEGER I
! 
! UT(1)=B(1)
! DO I =2,N
!  FLT(I)=A(I)/UT(I-1)
!  UT(I)=B(I)-FLT(I)*C(I-1)
! ENDDO
! 
! Y(1)=F(1)
! DO I=2,N
! Y(I)=F(I)-FLT(I)*Y(I-1)
! ENDDO
! 
! F(N)=Y(N)/UT(N)
! DO I=N-1,1,-1
! F(I)=(Y(I)-C(I)*F(I+1))/(UT(I))
! ENDDO
!
!END SUBROUTINE



!SUBROUTINE TRI(N,FLT,UT,CT,F)
! INTEGER N
! REAL FLT(2:N),UT(1:N),CT(1:N-1)
! COMPLEX F(1:N)
! COMPLEX Y(1:N)
! INTEGER I
! 
! Y(1)=F(1)
! DO I=2,N
! Y(I)=F(I)-FLT(I)*Y(I-1)
! ENDDO
! 
! F(N)=Y(N)/UT(N)
! DO I=N-1,1,-1
! F(I)=(Y(I)-CT(I)*F(I+1))/(UT(I))
! ENDDO
!
!END SUBROUTINE









!************* SUBROUTINE OF GAUSSIAN ELIMINATION *************************
SUBROUTINE GAUSSIAN_ELI_COMPLEX(N,A,B)
INTEGER N,K,I,J
COMPLEX A(1:N,1:N),B(1:N),L(1:N,1:N),SOLU(1:N),SUM
 

DO K=1,N-1
  DO I=K+1,N
   L(I,K)=A(I,K)/A(K,K)
  ENDDO
  
 
  DO J=1,K
   DO I=K+1,N
   A(I,J)=0.0
   ENDDO
  ENDDO 
  
  DO I=K+1,N
   DO J=K+1,N
   A(I,J)=A(I,J)-L(I,K)*A(K,J)
   ENDDO
   B(I)=B(I)-L(I,K)*B(K)
  ENDDO 

ENDDO

SOLU(N)=B(N)/A(N,N)
DO I=N-1,1,-1
SUM=0.0
 DO J=I+1,N
  SUM=SUM+A(I,J)*SOLU(J)
 ENDDO
SOLU(I)=(B(I)-SUM)/A(I,I)
ENDDO
B=SOLU
END SUBROUTINE