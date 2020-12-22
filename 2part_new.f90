module parameters
  integer :: N 
  real(8) :: g, dt, dx, m1, m2, x0, vv
  real(8), parameter :: pi = 3.14159265359_8
end module parameters

program twoparticles
  use lapack95
  use blas95
  use parameters
  implicit none

  real(8) :: L, V, F, mi, x1, xm1
  real(8),dimension(:),allocatable :: en
  complex(8),dimension(:),allocatable  :: psi, psi1
  complex(8) :: im = (0._8,1._8), alpha
  complex(8),dimension(:,:),allocatable :: ham, h, expH, one, mat1, mat2, mat3
  real(8) :: anm1, an0, an1, a, aa, a1, disp, d, d2, gg
  integer :: i, j, k, info, ixm1, ix0, ix1
  integer, parameter :: np=4 ! order+1 of the interpolationg polynomial
  real(8) :: ac(np), vandermonde(np,np), y(np)
  integer :: ix(np)
  character(len=4) :: c


  namelist /params/ L, N, g, dt, m1, m2, x0, vv
  open(9,file='2part.inp')
  read(9,nml=params)
  close(9)
  allocate(en(-N:N), psi(-N:N), psi1(-N:N), ham(-N:N,-N:N), h(-N:N,-N:N), expH(-N:N,-N:N), one(-N:N,-N:N), mat1(-N:N,-N:N), mat2(-N:N,-N:N), mat3(-N:N,-N:N))

!
! T=2pi (okres)
dt = 2*pi*dt ! "j" w jednostkach 0.001 okresu
!

  alpha = cmplx(1._8, 0._8) ! initial value

  dx = L/N
  mi = 1._8/(24*m1*dx*dx)    ! 1/(2*m) * 1/(12*dx^2) 5-point c.d.f.
  x1 = x0
 
  one = (0._8, 0._8)
  do i=-N,N
     one(i,i) = (1._8, 0._8)
  enddo
gg=g    ! UWAGA: startuje ze stanu bez oddzialywania
g=0._8  ! UWAGA: startuje ze stanu bez oddzialywania
  ham = (0._8, 0._8)
! 5-point central difference formula
  do i = -N, N
     ham(i,i) = V(i*dx,x0,alpha) + 30*mi
     if (i < N) then
        ham(i,i+1) = -16*mi
        ham(i+1,i) = ham(i,i+1)
        if (i < N-1) then
           ham(i,i+2) = mi
           ham(i+2,i) = mi
        endif
     endif
  enddo
g=gg ! UWAGA: startuje ze stanu bez oddzialywania
  h = ham
  call heevd(h,en,jobz='V',info=info)
  psi = h(:,-N)/sqrt(dx) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! normalizacja f. falowej -- sprawdzić
  open(10,file='out0000.dat')
  do i=-N,N
     write(10,'(f7.2,4f15.9)') real(i*dx),abs(psi(i))**2,abs(psi(i))**2,real(x0),real(alpha)
  enddo
  close(10)
  open(11,file='position.dat')
  call gemv(ham, psi, psi1)
  !write(11,'(i,4f15.8)') 0,x0,abs(dot_product(psi,psi1)),0.5*m2*(x0*x0 + vv*vv)
        d = 0
        d2 = 0
        do i=-N,N
			d = d + abs(psi(i))**2*i
			d2 = d2 + abs(psi(i))**2*i**2
		enddo
		d = d * dx
		d2 = d2 * dx**2
		disp = sqrt(d2 - d*d)
  write(11,'(f6.3,5F15.8)') 0.,x0,real(alpha),imag(alpha),d,disp 
  mat1 = one - 0.5*im*ham*dt
  mat2 = one + 0.5*im*ham*dt
! First two steps in time to start the Verlet integration scheme
  xm1 = x0  ! x_(-1)
  vv = vv + F(psi,x0,alpha)*dt/m2
  x0 = xm1 + vv*dt ! Euler method, force only from the harmonic trap


  do j = 1,10000

     do i = -N, N
        ham(i,i) = V(i*dx,x0,alpha) + 30*mi  ! only these elements change
     enddo
  
     ! Crank-Nicolson
     do i = -N,N
        mat1(i,i) = (1.0_8,0.0_8) - (0.0_8,0.5_8)*ham(i,i)*dt  ! only these elements change
        mat2(i,i) = 2.0_8 - mat1(i,i)                          ! 1 + 0.5*i*H*dt = 2 - (1 - 0.5*i*H*dt)
     enddo
     mat3 = mat2
     call inv(mat3,N) ! in-place
     expH = matmul(mat1,mat3)
     call gemv(expH, psi, psi1)
     psi = psi1
!     vv = vv + F(psi,0.5*(x0+x1))*dt/m2    ! the Euler method
!     x1 = x0
!     x0 = x0 + vv*dt
     x1 = 2*x0 - xm1 + F(psi,x0,alpha)*dt*dt/m2   ! the Verlet method
     xm1 = x0
     x0 = x1
! alpha evolution:
	 !alpha = alpha + cmplx(0._8, 1._8)*(1-alpha**2+a*g/m2)
! d^2n/dx^2
	 ! polynomial interpolation
	 i = floor(x0/L*N)
	 if (i <= -N+1) then
		ix = (/-N,-N+1,-N+2,-N+3/)
	 else if (i >= N-1) then
		ix = (/-N,-N+1,-N+2,-N+3/)
	 else 
		ix = (/ i-1,i,i+1,i+2 /)
	 endif
	 do i=1,4
		do k=1,4
			vandermonde(i,k)=(L*ix(i)/N)**(k-1)
		enddo
		y(i)=abs(psi(ix(i)))**2
	 enddo
	 call gesv(vandermonde,y)
	 a = 6*y(4)*x0+2*y(3)
	 aa = sqrt(1+g/m2*a)
	 alpha = im*aa*tan(aa*dt + atan(alpha/(im*aa)))
	 !alpha = alpha + cmplx(0._8, 1._8)*(1-alpha**2+a*g/m2)*dt
	 print '(f6.2,6f17.7)',x0,alpha, sum(abs(psi)**2)*dx
!
!#############################  SAVE RESULTS #################################
!
     if (mod(j,10) == 0) then 
        call gemv(ham, psi, psi1)
        h=ham
        call heevd(h,en,jobz='V',info=info)
        d = 0
        d2 = 0
        do i=-N,N
           d = d + abs(psi(i))**2*i
           d2 = d2 + abs(psi(i))**2*i**2
	enddo
	d = d * dx**2
	d2 = d2 * dx**3
	disp = sqrt(d2 - d*d)
        !write(11,'(i,4F15.8)') j/50,x0,abs(dot_product(psi,psi1)),0.5*m2*(x0*x0 + vv*vv)
        write(11,'(f6.3,5F15.8)') j*dt,x0,real(alpha),imag(alpha),d,disp 
        write(c,'(I4.4)') j/10
        open(10,file='out'//c//'.dat')
        do i=-N,N
           write(10,'(f7.2,4f15.9)') i*dx,abs(psi(i))**2,abs(h(i,-N))**2,x0,real(alpha)
        enddo
        close(10)
     endif
!
!#############################################################################
!
  enddo
  close(11)
end program twoparticles
!
!#############################################################################
!
subroutine inv(a,N)
  use lapack95
  implicit none
  complex(8),intent(inout) :: a(-N:N,-N:N)
  integer :: N,info,ipiv(-N:N)
  call getrf(a,ipiv,info)
  if (info /= 0) then
     print *,'GETRF',info
     stop 'GETRF'
  endif
  call getri(a,ipiv,info)
  if (info /= 0) then
     print *,'GETRI',info
     stop 'GETRI'
  endif
end subroutine 
!
!#############################################################################
!
function V(x,x00,alpha)
use parameters
implicit none
real(8) :: V,x,x00
complex(8) :: alpha
  V = 0.5*x*x + g*sqrt(m2*real(alpha)/pi)*exp(-m2*real(alpha)*(x-x00)**2)
end function V
!
!#############################################################################
!
function F(psi,x00,alpha)
use parameters
implicit none
real(8) :: F,x,x00
complex(8) :: psi(-N:N), alpha
integer :: i
F = -m2*x00
do i = -N,N
   x = i*dx
   !F = F - m2*g*(x-x00)*exp(-0.5*m2*(x-x00)**2)*abs(psi(i))**2
   F = F - g*((m2*real(alpha))**1.5)/sqrt(pi)*(x-x00)*exp(-m2*real(alpha)*(x-x00)**2)*abs(psi(i))**2 * dx
enddo
end function F
