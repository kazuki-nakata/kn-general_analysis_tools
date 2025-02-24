MODULE grid_data
USE sensor_geometry
IMPLICIT NONE
INTEGER(4),PARAMETER :: null = -32767
REAL(4),PARAMETER :: Undef = 9.9E33
REAL(4),PARAMETER :: mask_val = 999
CONTAINS

SUBROUTINE weighted_mean_sigma(grid_x,grid_y,val,mask_grid,out_val,n_grid,nx,ny,wsize,sigma,res,rm_outer)
!WA: weighted averaging
IMPLICIT NONE
INTEGER :: i,ii,j,jj,k,kk,l,nz,grid_i,grid_j,isum
INTEGER :: i2,j2,i3,j3,p,dum
INTEGER(4),INTENT(IN) :: n_grid,nx,ny,wsize
REAL(4),INTENT(IN) :: res,sigma,rm_outer
INTEGER(4),INTENT(IN) :: mask_grid(1:nx,1:ny)
REAL(4),INTENT(IN) :: grid_x(1:n_grid),grid_y(1:n_grid)
REAL(4),INTENT(IN) :: val(1:n_grid)
INTEGER(4) :: findex(1:nx,1:ny),pindex(1:n_grid),uindex(1:n_grid)
REAL(4) :: di,pj,pj_sim,uij
REAL(4),INTENT(OUT):: out_val(1:nx,1:ny)
REAL(4) :: scale_radius,window_radius,dist,adis,max_fr
REAL(8) :: fr,asum,vsum

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

window_radius=(sigma*16)/2.

!-------------------initialization----------------------
out_val(:,:)=0

do j =1,ny
  do i =1,nx
    isum=0
    asum=0
    max_fr=0
    vsum=0
    if(mask_grid(i,j).eq.999) cycle
    do jj =1,wsize
      do ii =1,wsize
        i2=i+ii-(wsize+1)/2
        j2=j+jj-(wsize+1)/2
        if((i2.lt.1).or.(i2.gt.nx)) cycle
        if((j2.lt.1).or.(j2.gt.ny)) cycle
        p=findex(i2,j2)

        if(p.eq.-1) cycle
        do while (p>0)
          dist=sqrt((real(i)-grid_x(p))**2+(real(j)-grid_y(p))**2)*res          
          if (dist.le.window_radius) then
            fr=exp(-dist**2./(sigma**2))
            asum=asum+fr
            isum=isum+1
            vsum=fr*val(p)+vsum
            if(fr.ge.max_fr) max_fr=fr
          endif
            p=pindex(p)
        enddo
      enddo
    enddo

    ! if(isum.lt.4) then
    if(max_fr.lt.rm_outer) then
      out_val(i,j)=Undef
    else
      out_val(i,j)=real(vsum/asum)
    endif
 enddo
enddo

END SUBROUTINE weighted_mean_sigma

SUBROUTINE weighted_mean_sigma_warp(grid_x0,grid_y0,val,u_grid0,v_grid0,mask_grid0,out_val,n_grid,nx,ny,nz,wsize,sigma,rm_outer)
IMPLICIT NONE
INTEGER :: i,ii,j,jj,k,kk,l,grid_i,grid_j,isum
INTEGER :: i2,j2,i3,j3,p,dum
INTEGER(4),INTENT(IN) :: n_grid,nx,ny,nz,wsize
INTEGER(4),INTENT(IN) :: mask_grid0(1:nx,1:ny)
REAL(4),INTENT(IN) :: sigma,rm_outer
REAL(4),INTENT(IN) :: u_grid0(1:nx,1:ny),v_grid0(1:nx,1:ny)
REAL(4),INTENT(IN) :: grid_x0(1:n_grid),grid_y0(1:n_grid)
REAL(4),INTENT(IN) :: val(1:n_grid,nz)
REAL(4),INTENT(OUT):: out_val(1:nx,1:ny,1:nz)
REAL(4) :: grid_x(1:n_grid),grid_y(1:n_grid)
REAL(4) :: u_grid(0:nx+1,0:ny+1),v_grid(0:nx+1,0:ny+1),mask_grid(0:nx+1,0:ny+1)
REAL(4) :: var(4),u(4),v(4),mask(4),dx0,dy0,dx,dy
INTEGER(4) :: findex(1:nx,1:ny),pindex(1:n_grid),uindex(1:n_grid)
REAL(4) :: di,pj,pj_sim,uij
INTEGER(4),PARAMETER :: mask0=999
REAL(4) :: scale_radius,window_radius,dist,adis,max_fr
REAL(8) :: fr,asum,vsum(nz)

u_grid(:,:)=0
u_grid(1:nx,1:ny)=u_grid0
v_grid(:,:)=0
v_grid(1:nx,1:ny)=v_grid0
mask_grid(:,:)=mask0
mask_grid(1:nx,1:ny)=mask_grid0

grid_x(:)=grid_x0(:)
grid_y(:)=grid_y0(:)
pindex(:)=-1
findex(:,:)=-1
uindex(:)=-1

do l = 1, n_grid
    dx0=grid_x(l)
    dy0=grid_y(l)
    i = nint(dx0)
    j = nint(dy0)
    if(mask_grid(i,j).eq.999) cycle
      i2 = int(grid_x(l))
      j2 = int(grid_y(l))
      u(1)=u_grid(i2,j2)
      u(2)=u_grid(i2+1,j2+1)
      u(3)=u_grid(i2+1,j2)
      u(4)=u_grid(i2,j2+1)
      v(1)=v_grid(i2,j2)
      v(2)=v_grid(i2+1,j2+1)
      v(3)=v_grid(i2+1,j2)
      v(4)=v_grid(i2,j2+1)
      mask(:)=0
      if(mask_grid(i2,j2).eq.mask0) mask(1)=999
      if(mask_grid(i2+1,j2+1).eq.mask0) mask(2)=999
      if(mask_grid(i2+1,j2).eq.mask0) mask(3)=999
      if(mask_grid(i2,j2+1).eq.mask0) mask(4)=999
      dx0=dx0-real(i2)
      dy0=dy0-real(j2)
      dx=bilin(dx0,dy0,u,mask)
      dy=bilin(dx0,dy0,v,mask)
      grid_x(l)=grid_x(l)+dx
      grid_y(l)=grid_y(l)+dy
      i = nint(grid_x(l))
      j = nint(grid_y(l))
      if((i.lt.1).or.(i.gt.nx)) cycle
      if((j.lt.1).or.(j.gt.ny)) cycle      
      pindex(l)=findex(i,j)
      uindex(l)=1
      findex(i,j)=l
enddo

window_radius=(sigma*16)/2.

out_val(:,:,:)=0

do j =1,ny
  do i =1,nx
    isum=0
    asum=0
    max_fr=0
    vsum(:)=0
    if(mask_grid(i,j).eq.999) cycle
    do jj =1,wsize
      do ii =1,wsize
        i2=i+ii-(wsize+1)/2
        j2=j+jj-(wsize+1)/2
        if((i2.lt.1).or.(i2.gt.nx)) cycle
        if((j2.lt.1).or.(j2.gt.ny)) cycle
        p=findex(i2,j2)
        if(p.eq.-1) cycle
        do while (p>0)
          dist=sqrt((real(i)-grid_x(p))**2+(real(j)-grid_y(p))**2)      
          if (dist.le.window_radius) then
            fr=exp(-dist**2./(sigma**2))
            asum=asum+fr
            isum=isum+1
            do k =1,nz
              vsum(k)=fr*val(p,k)+vsum(k)
            enddo
            if(fr.ge.max_fr) max_fr=fr

          endif
            p=pindex(p)
        enddo
      enddo
    enddo

    ! if(isum.lt.4) then
    if(max_fr.lt.rm_outer) then
      out_val(i,j,:)=Undef
    else
      out_val(i,j,:)=real(vsum/asum)
    endif
 enddo
enddo

END SUBROUTINE weighted_mean_sigma_warp



! SUBROUTINE bilinear_warp(grid_x0,grid_y0,val,u_grid0,v_grid0,mask_grid0,out_val,n_grid,nx,ny,nz,wsize,sigma,rm_outer)
! IMPLICIT NONE
! INTEGER :: i,ii,j,jj,k,kk,l,grid_i,grid_j,isum
! INTEGER :: i2,j2,i3,j3,p,dum
! INTEGER(4),INTENT(IN) :: n_grid,nx,ny,nz,wsize
! INTEGER(4),INTENT(IN) :: mask_grid0(1:nx,1:ny)
! REAL(4),INTENT(IN) :: sigma,rm_outer
! REAL(4),INTENT(IN) :: u_grid0(1:nx,1:ny),v_grid0(1:nx,1:ny)
! REAL(4),INTENT(IN) :: grid_x0(1:n_grid),grid_y0(1:n_grid)
! REAL(4),INTENT(IN) :: val(1:n_grid,nz)
! REAL(4),INTENT(OUT):: out_val(1:nx,1:ny,1:nz)
! REAL(4) :: grid_x(1:n_grid),grid_y(1:n_grid)
! REAL(4) :: u_grid(0:nx+1,0:ny+1),v_grid(0:nx+1,0:ny+1),mask_grid(0:nx+1,0:ny+1)
! REAL(4) :: var(4),u(4),v(4),mask(4),dx0,dy0,dx,dy
! INTEGER(4) :: findex(1:nx,1:ny),pindex(1:n_grid),uindex(1:n_grid)
! REAL(4) :: di,pj,pj_sim,uij
! INTEGER(4),PARAMETER :: mask0=999
! REAL(4) :: scale_radius,window_radius,dist,adis,max_fr
! REAL(8) :: fr,asum,vsum(nz)

! u_grid(:,:)=0
! u_grid(1:nx,1:ny)=u_grid0
! v_grid(:,:)=0
! v_grid(1:nx,1:ny)=v_grid0
! mask_grid(:,:)=mask0
! mask_grid(1:nx,1:ny)=mask_grid0

! grid_x(:)=grid_x0(:)
! grid_y(:)=grid_y0(:)
! pindex(:)=-1
! findex(:,:)=-1
! uindex(:)=-1

! do l = 1, n_grid
!     dx0=grid_x(l)
!     dy0=grid_y(l)
!     i = nint(dx0)
!     j = nint(dy0)
!     if(mask_grid(i,j).eq.999) cycle
!       i2 = int(grid_x(l))
!       j2 = int(grid_y(l))
!       u(1)=u_grid(i2,j2)
!       u(2)=u_grid(i2+1,j2+1)
!       u(3)=u_grid(i2+1,j2)
!       u(4)=u_grid(i2,j2+1)
!       v(1)=v_grid(i2,j2)
!       v(2)=v_grid(i2+1,j2+1)
!       v(3)=v_grid(i2+1,j2)
!       v(4)=v_grid(i2,j2+1)
!       mask(:)=0
!       if(mask_grid(i2,j2).eq.mask0) mask(1)=999
!       if(mask_grid(i2+1,j2+1).eq.mask0) mask(2)=999
!       if(mask_grid(i2+1,j2).eq.mask0) mask(3)=999
!       if(mask_grid(i2,j2+1).eq.mask0) mask(4)=999
!       dx0=dx0-real(i2)
!       dy0=dy0-real(j2)
!       dx=bilin(dx0,dy0,u,mask)
!       dy=bilin(dx0,dy0,v,mask)
!       grid_x(l)=grid_x(l)+dx
!       grid_y(l)=grid_y(l)+dy
!       i = nint(grid_x(l))
!       j = nint(grid_y(l))
!       if((i.lt.1).or.(i.gt.nx)) cycle
!       if((j.lt.1).or.(j.gt.ny)) cycle      
!       pindex(l)=findex(i,j)
!       uindex(l)=1
!       findex(i,j)=l
! enddo


! out_val(:,:,:)=0

! do j =1,ny
!   do i =1,nx
!     isum=0
!     asum=0
!     max_fr=0
!     vsum(:)=0
!     if(mask_grid(i,j).eq.999) cycle
!     do jj =1,wsize
!       do ii =1,wsize
!         i2=i+ii-(wsize+1)/2
!         j2=j+jj-(wsize+1)/2
!         if((i2.lt.1).or.(i2.gt.nx)) cycle
!         if((j2.lt.1).or.(j2.gt.ny)) cycle
!         p=findex(i2,j2)
!         if(p.eq.-1) cycle
!         do while (p>0)
!           dist=sqrt((real(i)-grid_x(p))**2+(real(j)-grid_y(p))**2)      
!           if (dist.le.window_radius) then
!             fr=exp(-dist**2./(sigma**2))
!             asum=asum+fr
!             isum=isum+1
!             do k =1,nz
!               vsum(k)=fr*val(p,k)+vsum(k)
!             enddo
!             if(fr.ge.max_fr) max_fr=fr

!           endif
!             p=pindex(p)
!         enddo
!       enddo
!     enddo

!     ! if(isum.lt.4) then
!     if(max_fr.lt.rm_outer) then
!       out_val(i,j,:)=Undef
!     else
!       out_val(i,j,:)=real(vsum/asum)
!     endif
!  enddo
! enddo

! END SUBROUTINE bilinear_warp


function bilin(dx0,dy0,var,mask)
  implicit none
  REAL(4) :: bilin
  INTEGER(4),PARAMETER :: isize=4
  REAL(4),INTENT(in) :: var(isize),mask(isize),dx0,dy0
  REAL(4) ::dx(isize),dy(isize)
  REAL(4) :: num1,num2
  INTEGER(4) :: i
  
  dx(1)=dx0
  dy(1)=dy0
  dx(2)=1-dx0
  dy(2)=1-dy0
  dx(3)=1-dx0
  dy(3)=dy0
  dx(4)=dx0
  dy(4)=1-dy0

  num1=0
  num2=0
  do i = 1, isize
     if(mask(i).eq.999) cycle
         num1= num1 + (1-dx(i))*(1-dy(i))*var(i)
         num2= num2 + (1-dx(i))*(1-dy(i))
      enddo
      if(num2.ne.0) then
         bilin=num1/num2
      else
         bilin=120
      endif
end function bilin


SUBROUTINE weighted_mean_sat(grid_x,grid_y,vs,vb,vg,val,mask_grid,out_val,n_grid,nx,ny,&
wsize,ap,nx_ap,ny_ap,int_ap,res,sigma,rm_outer)
!vs,vb,vg: ecef vector for s/c, obs_point(boresight), and obs_point()
IMPLICIT NONE
INTEGER :: i,ii,j,jj,k,kk,l,nz,grid_i,grid_j,isum
INTEGER :: i2,j2,i3,j3,p,dum
INTEGER(4),INTENT(IN) :: n_grid,nx,ny,wsize,nx_ap,ny_ap
REAL(4),INTENT(IN) :: int_ap,res,sigma,rm_outer
INTEGER(4),INTENT(IN),DIMENSION(1:nx,1:ny) :: mask_grid
REAL(4),INTENT(IN),DIMENSION(1:n_grid)  :: val, grid_x, grid_y
REAL(8),INTENT(IN),DIMENSION(1:nx,1:ny,3) :: vg
REAL(8),INTENT(IN),DIMENSION(1:n_grid,3) :: vs, vb
REAL(4),DIMENSION(n_grid,3) :: vi, vj, vk
REAL(4),DIMENSION(nx_ap,ny_ap) :: ap
REAL(8),DIMENSION(1,3) :: vdum
INTEGER(4) :: findex(1:nx,1:ny),pindex(1:n_grid),uindex(1:n_grid)
REAL(4),INTENT(OUT),DIMENSION(1:nx,1:ny):: out_val

REAL(4) :: window_radius,dist,adis,max_fr
REAL(4),DIMENSION(1):: az,el
REAL(8) :: fr,asum,val_sum,time_sum,pi

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
window_radius=(sigma*8)/2.

out_val(:,:)=0
k=0
! out_test(:,:)=0
do j =1,ny
  do i =1,nx
    val_sum=0
    asum=0
    isum=0
    max_fr=0
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
        dist=sqrt((real(i)-grid_x(p))**2+(real(j)-grid_y(p))**2)*res
        if (dist.le.window_radius) then
         call calc_local_az_el_angle(vi(p,:),vj(p,:),vk(p,:),vb(p,:)-vs(p,:),vb(p,:),vdum,az,el,int(1))
          fr=ap(nint(az(1)*180.0/pi/int_ap)+(nx_ap-1)/2+1,nint(el(1)*180.0/pi/int_ap)+(ny_ap-1)/2+1)
          asum=asum+fr
          isum=isum+1
          val_sum=fr*val(p)+val_sum
          if(fr.ge.max_fr) max_fr=fr
          endif
          p=pindex(p)
        enddo
      enddo
    enddo

    if(max_fr.lt.rm_outer) then
      out_val(i,j)=Undef
    else
      out_val(i,j)=val_sum/asum
    endif
 enddo
enddo

END SUBROUTINE weighted_mean_sat

SUBROUTINE nearest_neighbor(grid_x,grid_y,val,mask_grid,out_val,n_grid,nx,ny,wsize,res,threshold)
!WA: weighted averaging
IMPLICIT NONE
INTEGER :: i,ii,j,jj,k,kk,l,nz,grid_i,grid_j,isum
INTEGER :: i2,j2,i3,j3,p,dum
INTEGER(4),INTENT(IN) :: n_grid,nx,ny,wsize
REAL(4),INTENT(IN) :: res,threshold
INTEGER(4),INTENT(IN) :: mask_grid(1:nx,1:ny)
REAL(4),INTENT(IN) :: grid_x(1:n_grid),grid_y(1:n_grid)
REAL(4),INTENT(IN) :: val(1:n_grid)
INTEGER(4) :: findex(1:nx,1:ny),pindex(1:n_grid),uindex(1:n_grid)
REAL(4) :: di,pj,pj_sim,uij
REAL(4),INTENT(OUT):: out_val(1:nx,1:ny)
REAL(4) :: scale_radius,window_radius,dist,adis,max_dist
REAL(8) :: fr,asum

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

!-------------------initialization----------------------
out_val(:,:)=0

do j =1,ny
  do i =1,nx
    isum=0
    asum=0
    max_dist=threshold
    if(mask_grid(i,j).eq.999) cycle
    
    do jj =1,wsize
      do ii =1,wsize
        i2=i+ii-(wsize+1)/2
        j2=j+jj-(wsize+1)/2
        if((i2.lt.1).or.(i2.gt.nx)) cycle
        if((j2.lt.1).or.(j2.gt.ny)) cycle
        p=findex(i2,j2)
        if(p.eq.-1) cycle
        do while (p>0)
          dist=sqrt((real(i)-grid_x(p))**2+(real(j)-grid_y(p))**2)*res          
          if ((dist.le.threshold).and.(dist.le.max_dist)) then 
                out_val(i,j)=val(p)
                max_dist=dist
          endif
            p=pindex(p)
        enddo
      enddo
    enddo

    if(max_dist.eq.threshold) then
      out_val(i,j)=Undef
    endif
 enddo
enddo

END SUBROUTINE nearest_neighbor


END MODULE grid_data

