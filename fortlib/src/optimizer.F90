program gradient_descent_3d

implicit none

integer, parameter :: pr = selected_real_kind(15,3)

integer :: i

real(kind=16) :: f,f_prime
real(pr) :: fp,h
real(pr) :: x1_0,x2_0,x3_0,z0,tmp_z0
real(pr) :: tmp_x1_0,tmp_x2_0,tmp_x3_0
integer :: nb_iter,nb_max_iter
real(pr) :: alpha,eps,cond 

alpha = 0.1 ! learning rate
nb_max_iter = 100 ! Nb max d'iteration
eps = 0.0001 ! stop condition

x1_0 = 0.0 ! start point
x2_0 = 0.5 
x3_0 = 0.0 
z0 = f(x1_0,x2_0,x3_0)

cond = eps + 10.0 ! start with cond greater than eps (assumption)
nb_iter = 0 
tmp_z0 = z0
do while( cond > eps .and. nb_iter < nb_max_iter)
	tmp_x1_0 = x1_0 - alpha * f_prime(1,x1_0,x2_0,x3_0)
	tmp_x2_0 = x2_0 - alpha * f_prime(2,x1_0,x2_0,x3_0)
	tmp_x3_0 = x3_0 - alpha * f_prime(3,x1_0,x2_0,x3_0)
	x1_0 = tmp_x1_0
	x2_0 = tmp_x2_0
	x3_0 = tmp_x3_0
	z0 = f(x1_0,x2_0,x3_0)
	nb_iter = nb_iter + 1
	cond = abs( tmp_z0 - z0 )
	tmp_z0 = z0	
	write(6,*) x1_0,x2_0,x3_0,cond
end do
	
end program gradient_descent_3d

real(kind=16) function f(x1,x2,x3)
implicit none
integer, parameter :: pr = selected_real_kind(15,3)
real(pr), intent(in) :: x1,x2,x3
real(pr) :: x1_sphere,x2_sphere,x3_sphere
real(pr) :: r1,r2,r3
x1_sphere = 1.5
x2_sphere = 2.5
x3_sphere = 0.5
r1 = x1 - x1_sphere
r2 = x2 - x2_sphere
r3 = x3 - x3_sphere
f = r1**2 + r2**2 + r3**2
end function f

real(kind=16) function f_prime(dpi,x1,x2,x3)
implicit none
integer, parameter :: pr = selected_real_kind(15,3)
integer, intent(in) :: dpi
real(pr), intent(in) :: x1,x2,x3
real(pr) :: x1_p,x2_p,x3_p
real(kind=16) :: f
real(pr) :: h
h = 0.001
x1_p = x1
x2_p = x2
x3_p = x3
if( dpi == 1 ) x1_p = x1_p + h
if( dpi == 2 ) x2_p = x2_p + h
if( dpi == 3 ) x3_p = x3_p + h
f_prime = ( f(x1_p,x2_p,x3_p) - f(x1,x2,x3) ) / h
end function f_prime
