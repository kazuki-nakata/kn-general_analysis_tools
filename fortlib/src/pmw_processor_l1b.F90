MODULE PMW_Processor_L1B
USE sensor_geometry
IMPLICIT NONE
INTEGER(4),PARAMETER :: null = -32767
REAL(4),PARAMETER :: Undef = 9.9E33
REAL(4),PARAMETER :: mask_val = 999
CONTAINS

SUBROUTINE rSIR(grid_x,grid_y,vs,vb,vg,var,mask_grid,var_init,out_var,&
n_grid,nx,ny,wsize,ap,nx_ap,ny_ap,int_ap,window_radius,iterate)
!vs,vb,vg: ecef vector for s/c, obs_point(boresight), and obs_point()
IMPLICIT NONE
INTEGER :: i,ii,j,jj,k,kk,l,nz,grid_i,grid_j,isum
INTEGER :: i2,j2,i3,j3,p,dum
INTEGER(4),INTENT(IN) :: n_grid,nx,ny,wsize,nx_ap,ny_ap,iterate
! INTEGER(4),INTENT(IN) :: i_t,j_t
REAL(4),INTENT(IN) :: int_ap,window_radius
REAL(8),INTENT(IN),DIMENSION(1:nx,1:ny,3) :: vg
INTEGER(4),INTENT(IN),DIMENSION(1:nx,1:ny) :: mask_grid,var_init
REAL(4),INTENT(IN),DIMENSION(1:n_grid)  :: var, grid_x, grid_y
REAL(8),INTENT(IN),DIMENSION(1:n_grid,3) :: vs, vb
REAL(4),DIMENSION(nx,ny) :: pj_new
REAL(4),DIMENSION(n_grid) :: fi 
REAL(4),DIMENSION(n_grid,3) :: vi, vj, vk
REAL(4),DIMENSION(nx_ap,ny_ap) :: ap
REAL(8),DIMENSION(1,3) :: vdum
INTEGER(4) :: findex(1:nx,1:ny),pindex(1:n_grid),uindex(1:n_grid)
REAL(4) :: di,pj,pj_sim,uij
REAL(4),INTENT(OUT),DIMENSION(1:nx,1:ny):: out_var
! REAL(4),INTENT(OUT),DIMENSION(2000,6):: out_test

REAL(4) :: dist,adis
REAL(4),DIMENSION(1):: az,el
REAL(8) :: fr,asum,var_sum,pi

print *, "image nx=",nx,",","ny=",ny
print *, "antenna pattern nx=",nx_ap,",","ny=",ny_ap

pi=acos(-1.)
pindex(:)=-1
findex(:,:)=-1
uindex(:)=-1

do l = 1, n_grid
    i = nint(grid_x(l))
    j = nint(grid_y(l))
    if(mask_grid(i,j).eq.999) cycle
      pindex(l)=findex(i,j)
      uindex(l)=1
      findex(i,j)=l
enddo

!-------------------get basis vectors-------------------
call calc_boresight_basis_vectors(vb,vs,vi,vj,vk,n_grid)
out_var(:,:)=var_init(:,:)

!------------------rSIR--------------------------------

do k=1,iterate

fi(:)=0
do l = 1, n_grid
  if(uindex(l).eq.-1) cycle
    i = nint(grid_x(l))
    j = nint(grid_y(l))

    asum=0

    do jj =1,wsize
      do ii =1,wsize
        i2=i+ii-(wsize+1)/2
        j2=j+jj-(wsize+1)/2
        if((i2.lt.1).or.(i2.gt.nx)) cycle
        if((j2.lt.1).or.(j2.gt.ny)) cycle
        if(out_var(i2,j2).eq.9.9E33) cycle
          dist=sqrt((real(i2)-grid_x(l))**2+(real(j2)-grid_y(l))**2)
          if (dist.le.window_radius) then
            vdum(1,:)=vg(i2,j2,:)
            call calc_local_az_el_angle(vi(l,:),vj(l,:),vk(l,:),vb(l,:)-vs(l,:),vb(l,:),vdum,az,el,int(1))
            fr=ap(nint(az(1)*180.0/pi/int_ap)+(nx_ap-1)/2+1,nint(el(1)*180.0/pi/int_ap)+(ny_ap-1)/2+1)
            ! fr=exp(-dist**2./(sigma**2))
            asum=asum+fr
            fi(l)=fr*out_var(i2,j2)+fi(l)
          endif
      enddo
    enddo
    fi(l)=fi(l)/asum
enddo

pj_new(:,:)=9.9E33
do j = 1, ny
  do i = 1, nx
    if(mask_grid(i,j).eq.999) cycle
    if(out_var(i,j).eq.Undef) cycle
    pj_sim=0
    asum=0
    pj=out_var(i,j)
    vdum(1,:)=vg(i,j,:)

    do jj =1,wsize
    do ii =1, wsize
      i2=i+ii-(wsize+1)/2
      j2=j+jj-(wsize+1)/2
      if((i2.lt.1).or.(i2.gt.nx)) cycle
      if((j2.lt.1).or.(j2.gt.ny)) cycle
        p=findex(i2,j2)
        if(p.eq.-1) cycle
        do while (p>0)
          dist=sqrt((real(i)-grid_x(p))**2+(real(j)-grid_y(p))**2)
          if (dist.le.window_radius) then
            di=sqrt(var(p)/fi(p))
            if(di.ge.1) then
              uij=(1-1/di)/2./fi(p)+1/pj/di
              uij=1/uij
            else
              uij=(fi(p)/2.)*(1-di)+pj*di
            endif
            call calc_local_az_el_angle(vi(p,:),vj(p,:),vk(p,:),vb(p,:)-vs(p,:),vb(p,:),vdum,az,el,int(1))
            fr=ap(nint(az(1)*180.0/pi/int_ap)+(nx_ap-1)/2+1,nint(el(1)*180.0/pi/int_ap)+(ny_ap-1)/2+1)
            ! fr=exp(-dist**2./(sigma**2))
            asum=asum+fr
            pj_sim=fr*uij+pj_sim
          endif
          p=pindex(p)
        enddo
      enddo
    enddo
    pj_new(i,j)=pj_sim/asum
enddo
enddo
out_var=pj_new

enddo

END SUBROUTINE rSIR

SUBROUTINE Banach_gradient(grid_x,grid_y,vs,vb,vg,var,mask_grid,var_init,out_var,&
n_grid,nx,ny,wsize,ap,nx_ap,ny_ap,int_ap,window_radius,w,w2,iterate)
!vs,vb,vg: ecef vector for s/c, obs_point(boresight), and obs_point()
IMPLICIT NONE
INTEGER :: i,ii,j,jj,k,kk,l,nz,grid_i,grid_j,isum
INTEGER :: i2,j2,i3,j3,p,dum
INTEGER(4),INTENT(IN) :: n_grid,nx,ny,wsize,nx_ap,ny_ap,iterate
! INTEGER(4),INTENT(IN) :: i_t,j_t
REAL(4),INTENT(IN) :: int_ap,window_radius,w,w2
REAL(8),INTENT(IN),DIMENSION(1:nx,1:ny,3) :: vg
INTEGER(4),INTENT(IN),DIMENSION(1:nx,1:ny) :: mask_grid,var_init
REAL(4),INTENT(IN),DIMENSION(1:n_grid)  :: var, grid_x, grid_y
REAL(8),INTENT(IN),DIMENSION(1:n_grid,3) :: vs, vb
REAL(4),DIMENSION(nx,ny) :: pj_new
REAL(4),DIMENSION(n_grid) :: fi 
REAL(4),DIMENSION(n_grid,3) :: vi, vj, vk
REAL(4),DIMENSION(nx_ap,ny_ap) :: ap
REAL(8),DIMENSION(1,3) :: vdum
INTEGER(4) :: findex(1:nx,1:ny),pindex(1:n_grid),uindex(1:n_grid)
REAL(4) :: di,pj,pj_sim,uij
REAL(4),INTENT(OUT),DIMENSION(1:nx,1:ny):: out_var
! REAL(4),INTENT(OUT),DIMENSION(2000,6):: out_test

REAL(4) :: dist,adis
REAL(4),DIMENSION(1):: az,el
REAL(8) :: fr,asum,var_sum,pi

print *, "image nx=",nx,",","ny=",ny
print *, "antenna pattern nx=",nx_ap,",","ny=",ny_ap

pi=acos(-1.)
pindex(:)=-1
findex(:,:)=-1
uindex(:)=-1

do l = 1, n_grid
    i = nint(grid_x(l))
    j = nint(grid_y(l))
    if(mask_grid(i,j).eq.999) cycle
      pindex(l)=findex(i,j)
      uindex(l)=1
      findex(i,j)=l
enddo

!-------------------get basis vectors-------------------
call calc_boresight_basis_vectors(vb,vs,vi,vj,vk,n_grid)
out_var(:,:)=var_init(:,:)
!------------------rSIR--------------------------------

do k=1,iterate

fi(:)=0
do l = 1, n_grid
  if(uindex(l).eq.-1) cycle
    i = nint(grid_x(l))
    j = nint(grid_y(l))

    asum=0

    do jj =1,wsize
      do ii =1,wsize
        i2=i+ii-(wsize+1)/2
        j2=j+jj-(wsize+1)/2
        if((i2.lt.1).or.(i2.gt.nx)) cycle
        if((j2.lt.1).or.(j2.gt.ny)) cycle
        if(out_var(i2,j2).eq.9.9E33) cycle
          dist=sqrt((real(i2)-grid_x(l))**2+(real(j2)-grid_y(l))**2)
          if (dist.le.window_radius) then
            vdum(1,:)=vg(i2,j2,:)
            call calc_local_az_el_angle(vi(l,:),vj(l,:),vk(l,:),vb(l,:)-vs(l,:),vb(l,:),vdum,az,el,int(1))
            fr=ap(nint(az(1)*180.0/pi/int_ap)+(nx_ap-1)/2+1,nint(el(1)*180.0/pi/int_ap)+(ny_ap-1)/2+1)
            ! fr=exp(-dist**2./(sigma**2))
            asum=asum+fr
            fi(l)=fr*out_var(i2,j2)+fi(l)
          endif
      enddo
    enddo
    fi(l)=fi(l)/asum
enddo

pj_new(:,:)=9.9E33
do j = 1, ny
  do i = 1, nx
    if(mask_grid(i,j).eq.999) cycle
    if(out_var(i,j).eq.Undef) cycle
    pj_sim=0
    asum=0
    ! pj=out_var(i,j)
    vdum(1,:)=vg(i,j,:)
    call duality_map(out_var(i,j),pj,w)

    do jj =1,wsize
    do ii =1, wsize
      i2=i+ii-(wsize+1)/2
      j2=j+jj-(wsize+1)/2
      if((i2.lt.1).or.(i2.gt.nx)) cycle
      if((j2.lt.1).or.(j2.gt.ny)) cycle
        p=findex(i2,j2)
        if(p.eq.-1) cycle
        do while (p>0)
          dist=sqrt((real(i)-grid_x(p))**2+(real(j)-grid_y(p))**2)
          if (dist.le.window_radius) then
            call duality_map(fi(p)-var(p),di,w)
            ! di=fi(p)-var(p)

            call calc_local_az_el_angle(vi(p,:),vj(p,:),vk(p,:),vb(p,:)-vs(p,:),vb(p,:),vdum,az,el,int(1))
            fr=ap(nint(az(1)*180.0/pi/int_ap)+(nx_ap-1)/2+1,nint(el(1)*180.0/pi/int_ap)+(ny_ap-1)/2+1)
            asum=asum+fr
            pj_sim=fr*di+pj_sim
          endif
          p=pindex(p)
        enddo
      enddo
    enddo
    pj_sim=pj_sim/asum !Note: it is necessary to confirm this code (must introduce any antenna pattern scaling method?)
    ! pj_new(i,j)=pj-w2*pj_sim
    call duality_map(pj-w2*pj_sim,pj_new(i,j),real(w/(w-1)))!/asum

enddo
enddo
out_var=pj_new
enddo
END SUBROUTINE Banach_gradient

SUBROUTINE duality_map(x,y,p)
IMPLICIT NONE
REAL(4),INTENT(IN) :: x,p
REAL(4),INTENT(OUT) :: y
real(4) :: sgn
if(x>0) then
  sgn=1
else
  sgn=-1
endif
y=sgn*(abs(x)**(p-1))
end subroutine duality_map

SUBROUTINE get_antenna_pattern_matrix(grid_x,grid_y,vs,vb,vg,out_var,&
n_grid,nx,ny,nx_w,ny_w,ap,nx_ap,ny_ap,int_ap)
!vs,vb,vg: ecef vector for s/c, obs_point(boresight), and obs_point()
IMPLICIT NONE
INTEGER :: i,ii,j,jj,k,kk,l,nz,grid_i,grid_j,isum
INTEGER :: i2,j2,i3,j3,p,dum
INTEGER(4),INTENT(IN) :: n_grid,nx,ny,nx_w,ny_w,nx_ap,ny_ap
! INTEGER(4),INTENT(IN) :: i_t,j_t
REAL(4),INTENT(IN) :: int_ap
REAL(8),INTENT(IN),DIMENSION(1:nx,1:ny,3) :: vg
REAL(4),INTENT(IN),DIMENSION(1:n_grid)  :: grid_x, grid_y
REAL(8),INTENT(IN),DIMENSION(1:n_grid,3) :: vs, vb
REAL(4),DIMENSION(n_grid,3) :: vi, vj, vk
REAL(4),DIMENSION(nx_ap,ny_ap) :: ap
REAL(8),DIMENSION(1,3) :: vdum
INTEGER(4) :: findex(1:nx,1:ny),pindex(1:n_grid),uindex(1:n_grid)
REAL(4) :: di,pj,pj_sim,uij
REAL(4),INTENT(OUT),DIMENSION(1:nx_w,1:ny_w,n_grid):: out_var
! REAL(4),INTENT(OUT),DIMENSION(2000,6):: out_test

REAL(4) :: dist,adis
REAL(4),DIMENSION(1):: az,el
REAL(8) :: fr,asum,var_sum,pi

print *, "image nx=",nx,",","ny=",ny
print *, "antenna pattern nx=",nx_ap,",","ny=",ny_ap

pi=acos(-1.)
pindex(:)=-1
findex(:,:)=-1
uindex(:)=-1

do l = 1, n_grid
    i = nint(grid_x(l))
    j = nint(grid_y(l))
    pindex(l)=findex(i,j)
    uindex(l)=1
    findex(i,j)=l
enddo

!-------------------get basis vectors-------------------
call calc_boresight_basis_vectors(vb,vs,vi,vj,vk,n_grid)
out_var(:,:,:)=0

!------------------rSIR--------------------------------

do l = 1, n_grid
  if(uindex(l).eq.-1) cycle
    i = nint(grid_x(l))
    j = nint(grid_y(l))
    asum=0
    do jj =1,ny_w
      do ii =1,nx_w
        i2=i+ii-(nx_w+1)/2
        j2=j+jj-(ny_w+1)/2
        if((i2.lt.1).or.(i2.gt.nx)) cycle
        if((j2.lt.1).or.(j2.gt.ny)) cycle
        vdum(1,:)=vg(i2,j2,:)
        call calc_local_az_el_angle(vi(l,:),vj(l,:),vk(l,:),vb(l,:)-vs(l,:),vb(l,:),vdum,az,el,int(1))
        fr=ap(nint(az(1)*180.0/pi/int_ap)+(nx_ap-1)/2+1,nint(el(1)*180.0/pi/int_ap)+(ny_ap-1)/2+1)
        asum=asum+fr
        out_var(ii,jj,l)=fr
      enddo
    enddo
    out_var(:,:,l)=out_var(:,:,l)/asum
enddo

END SUBROUTINE get_antenna_pattern_matrix


SUBROUTINE get_antenna_pattern_matrix2(grid_x,grid_y,vs,vb,vg,out_var,&
n_grid,nx,ny,nx_w,ny_w,xmins,xmaxs,wn,ymin,ymax,ap,nx_ap,ny_ap,int_ap)
!vs,vb,vg: ecef vector for s/c, obs_point(boresight), and obs_point()
IMPLICIT NONE
INTEGER :: i,ii,j,jj,k,kk,l,nz,grid_i,grid_j,isum
INTEGER :: i2,j2,i3,j3,p,dum,xmin,xmax
INTEGER(4),INTENT(IN) :: n_grid,nx,ny,nx_w,ny_w,wn,ymin,ymax,nx_ap,ny_ap
INTEGER(4),INTENT(IN),DIMENSION(1:wn) :: xmins,xmaxs
REAL(4),INTENT(IN) :: int_ap
REAL(8),INTENT(IN),DIMENSION(1:nx,1:ny,3) :: vg
REAL(4),INTENT(IN),DIMENSION(1:n_grid)  :: grid_x, grid_y
REAL(8),INTENT(IN),DIMENSION(1:n_grid,3) :: vs, vb
REAL(4),DIMENSION(n_grid,3) :: vi, vj, vk
REAL(4),DIMENSION(nx_ap,ny_ap) :: ap
REAL(8),DIMENSION(1,3) :: vdum
INTEGER(4) :: findex(1:nx,1:ny),pindex(1:n_grid),uindex(1:n_grid)
REAL(4) :: di,pj,pj_sim,uij
REAL(4),INTENT(OUT),DIMENSION(1:nx_w, 1:ny_w, 1:nx_w, 1:ny_w, 1:wn):: out_var
REAL(4) :: dist,adis
REAL(4),DIMENSION(1):: az,el
REAL(8) :: fr,asum,var_sum,pi

print *, "image nx=",nx,",","ny=",ny
print *, "antenna pattern nx=",nx_ap,",","ny=",ny_ap

pi=acos(-1.)
pindex(:)=-1
findex(:,:)=-1
uindex(:)=-1

do l = 1, n_grid
    i = nint(grid_x(l))
    j = nint(grid_y(l))
    pindex(l)=findex(i,j)
    uindex(l)=1
    findex(i,j)=l
enddo

!-------------------get basis vectors-------------------
call calc_boresight_basis_vectors(vb,vs,vi,vj,vk,n_grid)
out_var(:,:,:,:,:)=0

do k = 1, wn
  xmin=xmins(k)
  xmax=xmaxs(k)
  do j = ymin, ymax
    do i =xmin, xmax
      asum=0
      if((i.lt.1).or.(i.gt.nx)) cycle
      if((j.lt.1).or.(j.gt.ny)) cycle

      vdum(1,:)=vg(i,j,:)
      do jj=1,ny_w
      do ii=1,nx_w
        i2=xmin+ii-1
        j2=ymin+jj-1
        if((i2.lt.1).or.(i2.gt.nx)) cycle
        if((j2.lt.1).or.(j2.gt.ny)) cycle
        p=findex(i2,j2)
        if(p.eq.-1) cycle
        call calc_local_az_el_angle(vi(p,:),vj(p,:),vk(p,:),vb(p,:)-vs(p,:),vb(p,:),vdum,az,el,int(1))
        i3=nint(az(1)*180.0/pi/int_ap)+(nx_ap-1)/2+1
        j3=nint(el(1)*180.0/pi/int_ap)+(ny_ap-1)/2+1
        if((i3.lt.1).or.(i3.gt.nx_ap)) cycle
        if((j3.lt.1).or.(j3.gt.ny_ap)) cycle
        asum=asum+ap(i3,j3)
        out_var(i-xmin+1,j-ymin+1,ii,jj,k)=ap(i3,j3)
      enddo
      enddo    
    enddo
    enddo
enddo

do k = 1, wn
  do jj = 1,ny_w
    do ii = 1,nx_w
    asum=0
     do j=1,ny_w
      do i=1,nx_w
        asum=out_var(i,j,ii,jj,k)+asum
      enddo
      enddo
        if(asum.ne.0) then
        out_var(:,:,ii,jj,k)=out_var(:,:,ii,jj,k)/asum
        endif
    enddo
    enddo
enddo

END SUBROUTINE get_antenna_pattern_matrix2

SUBROUTINE filter_type1(var,fil,out_var,nx,ny,nx_w,ny_w)
IMPLICIT NONE
INTEGER :: i,ii,j,jj,k,mx,my
INTEGER :: i2,j2,i3,j3
INTEGER(4),INTENT(IN):: nx,ny,nx_w,ny_w
REAL(4),INTENT(IN),DIMENSION(1:nx_w, 1:ny_w,1:nx):: fil
REAL(4),INTENT(IN),DIMENSION(1:nx, 1:ny):: var
REAL(4),INTENT(OUT),DIMENSION(1:nx, 1:ny):: out_var
REAL(4),ALLOCATABLE :: var2(:,:)
REAL(8) :: asum,var_sum

print *, "image nx=",nx,",","ny=",ny
mx=(nx_w-1)/2
my=(ny_w-1)/2
allocate(var2(-mx+1:nx+mx,-my+1:ny+my))
var2(:,:)=0
var2(1:nx,1:ny)=var(:,:)
out_var(:,:)=0

  do j = 1, ny
    do i =1, nx
    asum=0
      do jj=1,ny_w
      do ii=1,nx_w
        i2=i+ii-mx-1
        j2=j+jj-my-1
        asum=fil(ii,jj,i)+asum
        out_var(i,j)=var2(i2,j2)*fil(ii,jj,i)+out_var(i,j)
      enddo
      enddo
      out_var(i,j)=out_var(i,j)/asum
    enddo
    enddo

END SUBROUTINE filter_type1


SUBROUTINE init_conv(grid_x,grid_y,vs,vb,vg,var,mask_grid,out_var,n_grid,nx,ny,wsize,ap,nx_ap,ny_ap,int_ap)
!vs,vb,vg: ecef vector for s/c, obs_point(boresight), and obs_point()
IMPLICIT NONE
INTEGER :: i,ii,j,jj,k,kk,l,nz,grid_i,grid_j,isum
INTEGER :: i2,j2,i3,j3,p,dum
INTEGER(4),INTENT(IN) :: n_grid,nx,ny,wsize,nx_ap,ny_ap
! INTEGER(4),INTENT(IN) :: i_t,j_t
REAL(4),INTENT(IN) :: int_ap
REAL(8),INTENT(IN),DIMENSION(1:nx,1:ny,3) :: vg
INTEGER(4),INTENT(IN),DIMENSION(1:nx,1:ny) :: mask_grid
REAL(4),INTENT(IN),DIMENSION(1:n_grid)  :: var, grid_x, grid_y
REAL(8),INTENT(IN),DIMENSION(1:n_grid,3) :: vs, vb
REAL(4),DIMENSION(nx,ny) :: pj_new
REAL(4),DIMENSION(n_grid) :: fi 
REAL(4),DIMENSION(n_grid,3) :: vi, vj, vk
REAL(4),DIMENSION(nx_ap,ny_ap) :: ap
REAL(8),DIMENSION(1,3) :: vdum
INTEGER(4) :: findex(1:nx,1:ny),pindex(1:n_grid),uindex(1:n_grid)
REAL(4) :: di,pj,pj_sim,uij
REAL(4),INTENT(OUT),DIMENSION(1:nx,1:ny):: out_var
! REAL(4),INTENT(OUT),DIMENSION(2000,6):: out_test

REAL(4) :: sigma,scale_radius,window_radius,dist,adis
REAL(4),DIMENSION(1):: az,el
REAL(8) :: fr,asum,var_sum,pi

print *, "image nx=",nx,",","ny=",ny
print *, "antenna pattern nx=",nx_ap,",","ny=",ny_ap

pi=acos(-1.)
pindex(:)=-1
findex(:,:)=-1
uindex(:)=-1

do l = 1, n_grid
    i = nint(grid_x(l))
    j = nint(grid_y(l))
    if(mask_grid(i,j).eq.999) cycle
      pindex(l)=findex(i,j)
      uindex(l)=1
      findex(i,j)=l
enddo

!-------------------get basis vectors-------------------
call calc_boresight_basis_vectors(vb,vs,vi,vj,vk,n_grid)
!-------------------initialization----------------------
out_var(:,:)=0
k=0
! out_test(:,:)=0
do j =1,ny
  do i =1,nx
    var_sum=0
    asum=0
    isum=0

    if(mask_grid(i,j).eq.999) cycle
    
    vdum(1,:)=vg(i,j,:)

    do jj =1,wsize
      do ii =1,wsize
        i2=i+ii-(wsize+1)/2
        j2=j+jj-(wsize+1)/2
        if((i2.lt.1).or.(i2.gt.nx)) cycle
        if((j2.lt.1).or.(j2.gt.ny)) cycle
        
        p=findex(i2,j2)

        if(p.eq.-1) cycle
        do while (p>0)
         call calc_local_az_el_angle(vi(p,:),vj(p,:),vk(p,:),vb(p,:)-vs(p,:),vb(p,:),vdum,az,el,int(1))
         fr=ap(nint(az(1)*180.0/pi/int_ap)+(nx_ap-1)/2+1,nint(el(1)*180.0/pi/int_ap)+(ny_ap-1)/2+1)
         asum=asum+fr
         isum=isum+1
         var_sum=fr*var(p)+var_sum
         p=pindex(p)
        enddo
      enddo
    enddo

    if(isum.lt.4) then
      out_var(i,j)=Undef
    else
      out_var(i,j)=var_sum/asum
    endif
 enddo
enddo


END SUBROUTINE init_conv


SUBROUTINE test(grid_x,grid_y,val,out_var,n_grid,map,nx,ny)
!vs,vb,vg: ecef vector for s/c, obs_point(boresight), and obs_point()
IMPLICIT NONE
INTEGER :: i,j,k
INTEGER(4),INTENT(IN) :: n_grid,nx,ny
REAL(4),INTENT(IN),DIMENSION(1:n_grid)  :: grid_x, grid_y,val
REAL(4),INTENT(IN),DIMENSION(1:nx,1:ny)  :: map
REAL(4),INTENT(OUT),DIMENSION(1:nx,1:ny):: out_var

print *, "image nx=",nx,",","ny=",ny

do k = 1, n_grid
    i = nint(grid_x(k))
    j = nint(grid_y(k))
    out_var(i,j)=val(k)+map(i,j)
enddo

END SUBROUTINE test


SUBROUTINE test2(ap,nx_ap,ny_ap,out_var)
IMPLICIT NONE
INTEGER :: i,j
INTEGER(4),INTENT(IN) :: nx_ap,ny_ap
REAL(4),INTENT(IN),DIMENSION(nx_ap,ny_ap) :: ap
REAL(4),INTENT(OUT),DIMENSION(nx_ap,ny_ap):: out_var

print *, "antenna pattern nx=",nx_ap,",","ny=",ny_ap

do j=1,ny_ap
  do i =1,nx_ap
    print *,i,j,ap(i,j)
    out_var(i,j)=ap(i,j)
  enddo
enddo

END SUBROUTINE test2

SUBROUTINE test3(grid_x,grid_y,val,out_var,n_grid,map,map2,nx,ny,nx2,ny2)
!vs,vb,vg: ecef vector for s/c, obs_point(boresight), and obs_point()
IMPLICIT NONE
INTEGER :: i,j,k,l
INTEGER(4),INTENT(IN) :: n_grid,nx,ny,nx2,ny2
REAL(4),INTENT(IN),DIMENSION(1:n_grid)  :: grid_x, grid_y
REAL(4),INTENT(IN),DIMENSION(1:n_grid,3)  :: val
REAL(4),INTENT(IN),DIMENSION(1:nx,1:ny,3)  :: map
REAL(4),INTENT(IN),DIMENSION(1:nx,1:ny,1:nx2,1:ny2,3)  :: map2
REAL(4),INTENT(OUT),DIMENSION(1:nx,1:ny,3):: out_var

print *,nx,ny,nx2,ny2,n_grid

do l =1,3
  do k = 1, n_grid
      i = nint(grid_x(k))
      j = nint(grid_y(k))
      out_var(i,j,l)=val(k,l)+map(i,j,l)
  enddo
enddo

END SUBROUTINE test3


END MODULE PMW_Processor_L1B

