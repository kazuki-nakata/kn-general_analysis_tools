MODULE sensor_geometry
IMPLICIT NONE
CONTAINS

subroutine calc_cross_product(v1,v2,outer)
implicit none

real(8),intent(in),dimension(3) :: v1,v2
real(8),intent(out),dimension(3) :: outer

outer(1)=v1(2)*v2(3)-v1(3)*v2(2)
outer(2)=v1(3)*v2(1)-v1(1)*v2(3)
outer(3)=v1(1)*v2(2)-v1(2)*v2(1)

endsubroutine calc_cross_product

subroutine calc_local_az_el_angle(i,j,k,b,p0,p,az,el,nx)
    implicit none
    real(8):: x,y,z
    integer,intent(in) :: nx
    integer :: l
    real,intent(in),dimension(3) :: i,j,k
    real(8),intent(in),dimension(3) :: b,p0
    real(8),intent(in),dimension(nx,3) :: p
    real,intent(out),dimension(nx) :: az, el

    do l =1,nx
      x=dot_product((b+p(l,:)-p0),dble(i))
      y=dot_product((b+p(l,:)-p0),dble(j))
      z=dot_product((b+p(l,:)-p0),dble(k))
      az(l)=atan2(x,z)
      el(l)=atan2(y,z)
    enddo
endsubroutine calc_local_az_el_angle

subroutine calc_local_az_el_angle2(i,j,k,b,p0,p,az,el,nx)
    implicit none
    real(8):: x,y,z
    integer,intent(in) :: nx
    integer :: l
    real,intent(in),dimension(nx,3) :: i,j,k
    real(8),intent(in),dimension(nx,3) :: b,p0
    real(8),intent(in),dimension(nx,3) :: p
    real,intent(out),dimension(nx) :: az, el

    do l =1,nx
      x=dot_product((b(l,:)+p(l,:)-p0(l,:)),dble(i(l,:)))
      y=dot_product((b(l,:)+p(l,:)-p0(l,:)),dble(j(l,:)))
      z=dot_product((b(l,:)+p(l,:)-p0(l,:)),dble(k(l,:)))
      az(l)=atan2(x,z)
      el(l)=atan2(y,z)
    enddo
endsubroutine calc_local_az_el_angle2

subroutine calc_boresight_basis_vectors(p0,s,i,j,k,nx)
    implicit none
    integer(4) :: l
    real(8),dimension(3) :: a, b
    integer(4),intent(in)::nx
    real(8),intent(in),dimension(nx,3) :: p0,s
    real(4),intent(out),dimension(nx,3) :: i,j,k

    do l = 1, nx
      a=p0(l,:)-s(l,:)
      k(l,:)=a/norm2(a)
      call calc_cross_product(s(l,:),a,b)
      i(l,:)=b/norm2(b)
      call calc_cross_product(dble(k(l,:)),dble(i(l,:)),b)
      j(l,:)=b
    enddo
endsubroutine calc_boresight_basis_vectors


END MODULE sensor_geometry

