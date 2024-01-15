MODULE Deblur
IMPLICIT NONE
INTEGER(4),PARAMETER :: null = 999
REAL(4),PARAMETER :: Undef = 999
CONTAINS

SUBROUTINE SIR(grid_x,grid_y,val,mask_grid,out_val,n_grid,nx,ny,wsize,ap,nx_ap,ny_ap,int_ap,res,fwhm,iterate)
IMPLICIT NONE
INTEGER :: i,ii,j,jj,k,kk,l,nz,grid_i,grid_j,isum
INTEGER :: i2,j2,i3,j3,p,dum,ip,jp
INTEGER(4),INTENT(IN) :: n_grid,nx,ny,wsize,nx_ap,ny_ap,iterate
REAL(4),INTENT(IN) :: int_ap,res,fwhm
INTEGER(4),INTENT(IN),DIMENSION(1:nx,1:ny) :: mask_grid
REAL(4),INTENT(IN),DIMENSION(1:n_grid)  :: val, grid_x, grid_y
REAL(4),DIMENSION(nx,ny) :: pj_new
REAL(4),DIMENSION(n_grid) :: fi 
REAL(4),INTENT(IN),DIMENSION(nx_ap,ny_ap) :: ap
INTEGER(4) :: findex(1:nx,1:ny),pindex(1:n_grid),uindex(1:n_grid)
REAL(4) :: di,pj,pj_sim,uij
REAL(4),INTENT(OUT),DIMENSION(1:nx,1:ny):: out_val

REAL(4) :: sigma,scale_radius,window_radius,dist,adis
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

scale_radius=fwhm/2
sigma=scale_radius*(2/2.35482) !2*sqrt(2ln2)=2.35482, sigma=sigma/2 for gaussian beam, fwhm/sqrt(2ln2) / 2
window_radius=(sigma*8)/2.

!-------------------initialization----------------------
out_val(:,:)=0
k=0
do j =1,ny
  do i =1,nx
    val_sum=0
    time_sum=0
    asum=0
    isum=0

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
          ip=nint((real(i)-grid_x(p))*res/int_ap)+(nx_ap-1)/2+1
          jp=nint((real(j)-grid_y(p))*res/int_ap)+(ny_ap-1)/2+1
          if((ip.gt.0).and.(ip.le.nx_ap)) then
          if((jp.gt.0).and.(jp.le.ny_ap)) then          
          fr=ap(ip,jp)
          asum=asum+fr
          isum=isum+1
          val_sum=fr*val(p)+val_sum
          endif
          endif
          endif
          p=pindex(p)
        enddo
      enddo
    enddo

    if(isum.lt.2) then
      out_val(i,j)=9.9E33
    else
      out_val(i,j)=val_sum/asum
    endif
 enddo
enddo


!------------------SIR--------------------------------

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
        if(out_val(i2,j2).eq.9.9E33) cycle
          dist=sqrt((real(i2)-grid_x(l))**2+(real(j2)-grid_y(l))**2)*res
          if (dist.le.window_radius) then
          ip=nint((real(i2)-grid_x(l))*res/int_ap)+(nx_ap-1)/2+1
          jp=nint((real(j2)-grid_y(l))*res/int_ap)+(ny_ap-1)/2+1
          
          if((ip.gt.0).and.(ip.le.nx_ap)) then
          if((jp.gt.0).and.(jp.le.ny_ap)) then          
            fr=ap(ip,jp)
            ! fr=exp(-dist**2./(sigma**2))
            asum=asum+fr
            fi(l)=fr*out_val(i2,j2)+fi(l)
          endif
          endif

          endif
      enddo
    enddo
    fi(l)=fi(l)/asum
enddo


pj_new(:,:)=9.9E33
do j = 1, ny
  do i = 1, nx
    if(mask_grid(i,j).eq.999) cycle
    if(out_val(i,j).eq.9.9E33) cycle
    pj_sim=0
    asum=0
    pj=out_val(i,j)

    do jj =1,wsize
    do ii =1, wsize
      i2=i+ii-(wsize+1)/2
      j2=j+jj-(wsize+1)/2
      if((i2.lt.1).or.(i2.gt.nx)) cycle
      if((j2.lt.1).or.(j2.gt.ny)) cycle
        p=findex(i2,j2)
        if(p.eq.-1) cycle
        do while (p>0)
          dist=sqrt((real(i)-grid_x(p))**2+(real(j)-grid_y(p))**2)*res
          if (dist.le.window_radius) then
            di=sqrt(val(p)/fi(p))
            if(di.ge.1) then
              uij=(1-1/di)/2./fi(p)+1/pj/di
              uij=1/uij
            else
              uij=(fi(p)/2.)*(1-di)+pj*di
            endif
            
            ip=nint((real(i)-grid_x(p))*res/int_ap)+(nx_ap-1)/2+1
            jp=nint((real(j)-grid_y(p))*res/int_ap)+(ny_ap-1)/2+1
            if((ip.gt.0).and.(ip.le.nx_ap)) then
            if((jp.gt.0).and.(jp.le.ny_ap)) then          
            fr=ap(ip,jp)
            ! fr=exp(-dist**2./(sigma**2))
            asum=asum+fr
            pj_sim=fr*uij+pj_sim
            endif
            endif

          endif
          p=pindex(p)
        enddo
      enddo
    enddo
    pj_new(i,j)=pj_sim/asum
enddo
enddo
out_val=pj_new

enddo

END SUBROUTINE SIR


SUBROUTINE SIR2(grid_x,grid_y,val,mask_grid,ap,out_val,n_grid,nx,ny,wsize,res,fwhm,w,iterate)
IMPLICIT NONE
INTEGER :: i,ii,j,jj,k,kk,l,nz,grid_i,grid_j,isum
INTEGER :: i2,j2,i3,j3,p,dum,ip,jp
INTEGER(4),INTENT(IN) :: n_grid,nx,ny,wsize,iterate
REAL(4),INTENT(IN) :: res,fwhm,w
REAL(4),INTENT(IN),DIMENSION(1:n_grid)  :: val, grid_x, grid_y
REAL(4),DIMENSION(nx,ny) :: pj_new
REAL(4),DIMENSION(n_grid) :: fi
REAL(4),INTENT(IN),DIMENSION(1:n_grid,nx,ny) :: ap
REAL(4),INTENT(IN),DIMENSION(nx,ny) :: mask_grid
INTEGER(4) :: findex(1:nx,1:ny),pindex(1:n_grid),uindex(1:n_grid)
REAL(4) :: di,pj,pj_sim,uij
REAL(4),INTENT(OUT),DIMENSION(1:nx,1:ny):: out_val

REAL(4) :: sigma,scale_radius,window_radius,dist,adis
REAL(8) :: fr,asum,val_sum,time_sum,pi

print *, "image nx=",nx,",","ny=",ny

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
      out_val(i,j)=val(l)
enddo

scale_radius=fwhm/2
sigma=scale_radius*(2/2.35482) !2*sqrt(2ln2)=2.35482, sigma=sigma/2 for gaussian beam, fwhm/sqrt(2ln2) / 2
window_radius=(sigma*8)/2.

!-------------------initialization----------------------
! out_val(:,:)=0
! k=0
! do j =1,ny
!   do i =1,nx
!     val_sum=0
!     time_sum=0
!     asum=0
!     isum=0

!     do jj =1,ny
!       do ii =1,nx
        
!         p=findex(ii,jj)
!         if(p.eq.-1) cycle
!         do while (p>0)
!           ip=int(grid_x(p))
!           jp=int(grid_y(p))
!           k=nx*(jp-1)+ip
!           fr=ap(k,i,j)
!           asum=asum+fr
!           val_sum=fr*val(p)+val_sum
!           p=pindex(p)
!         enddo
!       enddo
!     enddo

!     out_val(i,j)=val_sum/asum
!  enddo
! enddo

! ! !------------------SIR--------------------------------

do k=1,iterate

fi(:)=0
do l = 1, n_grid
  if(uindex(l).eq.-1) cycle
    i = nint(grid_x(l))
    j = nint(grid_y(l))
    kk=nx*(j-1)+i

    asum=0

    do jj =1,ny
      do ii =1,nx
            if(mask_grid(ii,jj).eq.999) cycle
            fr=ap(kk,ii,jj)
            asum=asum+fr
            fi(l)=fr*out_val(ii,jj)+fi(l)
      enddo
    enddo
    fi(l)=fi(l)/asum
enddo

pj_new(:,:)=9.9E33
do j = 1, ny
  do i = 1, nx
    pj_sim=0
    asum=0
    pj=out_val(i,j)
    if(mask_grid(i,j).eq.999) cycle

    do jj =1,ny
    do ii =1,nx
        p=findex(ii,jj)
        if(p.eq.-1) cycle
        do while (p>0)
            ! di=sqrt(val(p)/fi(p))
            di=(val(p)/fi(p))**w
            if(di.ge.1) then
              uij=(1-1/di)/2./fi(p)+1/pj/di
              uij=1/uij
            else
              uij=(fi(p)/2.)*(1-di)+pj*di
            endif
            ip=int(grid_x(p))
            jp=int(grid_y(p))
            kk=nx*(jp-1)+ip
            fr=ap(kk,i,j)
            asum=asum+fr
            pj_sim=fr*uij+pj_sim
            p=pindex(p)
        enddo
      enddo
    enddo
    pj_new(i,j)=pj_sim/asum
enddo
enddo
out_val=pj_new

enddo

END SUBROUTINE SIR2

SUBROUTINE MART(grid_x,grid_y,val,mask_grid,ap,out_val,n_grid,nx,ny,wsize,res,fwhm,w,iterate,frmin)
IMPLICIT NONE
INTEGER :: i,ii,j,jj,k,kk,l,nz,grid_i,grid_j,isum
INTEGER :: i2,j2,i3,j3,p,dum,ip,jp
INTEGER(4),INTENT(IN) :: n_grid,nx,ny,wsize,iterate
REAL(4),INTENT(IN) :: res,fwhm,w,frmin
REAL(4),INTENT(IN),DIMENSION(1:n_grid)  :: val, grid_x, grid_y
REAL(4),DIMENSION(nx,ny) :: pj_new
REAL(4),DIMENSION(n_grid) :: fi
REAL(4),INTENT(IN),DIMENSION(1:n_grid,nx,ny) :: ap
REAL(4),INTENT(IN),DIMENSION(nx,ny) :: mask_grid
INTEGER(4) :: findex(1:nx,1:ny),pindex(1:n_grid),uindex(1:n_grid)
REAL(4) :: di,pj,pj_sim,uij
REAL(4),INTENT(OUT),DIMENSION(1:nx,1:ny):: out_val

REAL(4) :: sigma,scale_radius,window_radius,dist,adis
REAL(8) :: fr,asum,val_sum,time_sum,pi

print *, "image nx=",nx,",","ny=",ny

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
      out_val(i,j)=val(l)
enddo

scale_radius=fwhm/2
sigma=scale_radius*(2/2.35482) !2*sqrt(2ln2)=2.35482, sigma=sigma/2 for gaussian beam, fwhm/sqrt(2ln2) / 2
window_radius=(sigma*8)/2.


do k=1,iterate

fi(:)=0
do l = 1, n_grid
  if(uindex(l).eq.-1) cycle
    i = nint(grid_x(l))
    j = nint(grid_y(l))
    kk=nx*(j-1)+i

    asum=0

    do jj =1,ny
      do ii =1,nx
            if(mask_grid(ii,jj).eq.999) cycle
            fr=ap(kk,ii,jj)
            asum=asum+fr
            fi(l)=fr*out_val(ii,jj)+fi(l)
      enddo
    enddo

    fi(l)=fi(l)/asum
    di=(val(l)/fi(l))

    do jj = 1, ny
      do ii = 1, nx
            fr=ap(kk,ii,jj)
            pj=out_val(ii,jj)
            if(fr.ge.frmin) then
              out_val(ii,jj)=out_val(ii,jj)*di**(w*fr)
            endif
       enddo
    enddo

enddo
enddo

END SUBROUTINE MART

SUBROUTINE MART2(grid_x,grid_y,val,ap,out_val,n_grid,nx,ny,wsize,res,fwhm,w,iterate)
IMPLICIT NONE
INTEGER :: i,ii,j,jj,k,kk,l,nz,grid_i,grid_j,isum
INTEGER :: i2,j2,i3,j3,p,dum,ip,jp
INTEGER(4),INTENT(IN) :: n_grid,nx,ny,wsize,iterate
REAL(4),INTENT(IN) :: res,fwhm,w
REAL(4),INTENT(IN),DIMENSION(1:n_grid)  :: val, grid_x, grid_y
REAL(4),DIMENSION(nx,ny) :: pj_new
REAL(4),DIMENSION(n_grid) :: fi 
REAL(4),INTENT(IN),DIMENSION(1:n_grid,nx,ny) :: ap
INTEGER(4) :: findex(1:nx,1:ny),pindex(1:n_grid),uindex(1:n_grid)
REAL(4) :: di,pj,pj_sim,uij
REAL(4),INTENT(OUT),DIMENSION(1:nx,1:ny):: out_val

REAL(4) :: sigma,scale_radius,window_radius,dist,adis
REAL(8) :: fr,asum,val_sum,time_sum,pi

print *, "image nx=",nx,",","ny=",ny

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
    out_val(i,j)=val(l)
enddo

scale_radius=fwhm/2
sigma=scale_radius*(2/2.35482) !2*sqrt(2ln2)=2.35482, sigma=sigma/2 for gaussian beam, fwhm/sqrt(2ln2) / 2
window_radius=(sigma*8)/2.

!-------------------initialization----------------------
! out_val(:,:)=0
! k=0
! do j =1,ny
!   do i =1,nx
!     val_sum=0
!     time_sum=0
!     asum=0
!     isum=0

!     do jj =1,ny
!       do ii =1,nx
        
!         p=findex(ii,jj)
!         if(p.eq.-1) cycle
!         do while (p>0)
!           ip=int(grid_x(p))
!           jp=int(grid_y(p))
!           k=nx*(jp-1)+ip
!           fr=ap(k,i,j)
!           asum=asum+fr
!           val_sum=fr*val(p)+val_sum
!           p=pindex(p)
!         enddo
!       enddo
!     enddo

!     out_val(i,j)=val_sum/asum
!  enddo
! enddo

! ! !------------------SIR--------------------------------

do k=1,iterate

fi(:)=0
do l = 1, n_grid
  if(uindex(l).eq.-1) cycle
    i = nint(grid_x(l))
    j = nint(grid_y(l))
    kk=nx*(j-1)+i

    asum=0

    do jj =1,ny
      do ii =1,nx
            fr=ap(kk,ii,jj)
            asum=asum+fr
            fi(l)=fr*out_val(ii,jj)+fi(l)
      enddo
    enddo
    fi(l)=fi(l)/asum
enddo


pj_new(:,:)=9.9E33
do j = 1, ny
  do i = 1, nx
    pj_sim=0
    asum=0
    pj=out_val(i,j)

    do jj =1,ny
    do ii =1,nx
        p=findex(ii,jj)
        if(p.eq.-1) cycle
        do while (p>0)
            ! di=sqrt(val(p)/fi(p))
            di=(val(p)/fi(p))**w
            uij=di!pj*di
            ip=int(grid_x(p))
            jp=int(grid_y(p))
            kk=nx*(jp-1)+ip
            fr=ap(kk,i,j)
            asum=asum+fr
            pj_sim=fr*uij+pj_sim
            p=pindex(p)
        enddo
      enddo
    enddo
    pj_new(i,j)=pj*pj_sim/asum
enddo
enddo
out_val=pj_new

enddo

END SUBROUTINE MART2


SUBROUTINE Banach_gradient(grid_x,grid_y,val,ap,out_val,n_grid,nx,ny,wsize,res,fwhm,w,w2,iterate)
IMPLICIT NONE
INTEGER :: i,ii,j,jj,k,kk,l,nz,grid_i,grid_j,isum
INTEGER :: i2,j2,i3,j3,p,dum,ip,jp
INTEGER(4),INTENT(IN) :: n_grid,nx,ny,wsize,iterate
REAL(4),INTENT(IN) :: res,fwhm,w,w2
REAL(4),INTENT(IN),DIMENSION(1:n_grid)  :: val, grid_x, grid_y
REAL(4),DIMENSION(nx,ny) :: pj_new
REAL(4),DIMENSION(n_grid) :: fi 
REAL(4),INTENT(IN),DIMENSION(1:n_grid,nx,ny) :: ap
INTEGER(4) :: findex(1:nx,1:ny),pindex(1:n_grid),uindex(1:n_grid)
REAL(4) :: di,pj,pj_sim,uij
REAL(4),INTENT(OUT),DIMENSION(1:nx,1:ny):: out_val

REAL(4) :: sigma,scale_radius,window_radius,dist,adis
REAL(8) :: fr,asum,val_sum,time_sum,pi

print *, "image nx=",nx,",","ny=",ny

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
    out_val(i,j)=val(l)
enddo

scale_radius=fwhm/2
sigma=scale_radius*(2/2.35482) !2*sqrt(2ln2)=2.35482, sigma=sigma/2 for gaussian beam, fwhm/sqrt(2ln2) / 2
window_radius=(sigma*8)/2.

do k=1,iterate

fi(:)=0
do l = 1, n_grid
  if(uindex(l).eq.-1) cycle
    i = nint(grid_x(l))
    j = nint(grid_y(l))
    kk=nx*(j-1)+i

    asum=0

    do jj =1,ny
      do ii =1,nx
            fr=ap(kk,ii,jj)
            asum=asum+fr
            fi(l)=fr*out_val(ii,jj)+fi(l)
      enddo
    enddo
    fi(l)=fi(l)/asum
enddo

pj_new(:,:)=9.9E33
do j = 1, ny
  do i = 1, nx
    pj_sim=0
    asum=0
    call duality_map(out_val(i,j),pj,w)
    
    do jj =1,ny
    do ii =1,nx
        p=findex(ii,jj)
        if(p.eq.-1) cycle
        do while (p>0)
            call duality_map(fi(p)-val(p),di,w)
            ip=int(grid_x(p))
            jp=int(grid_y(p))
            kk=nx*(jp-1)+ip
            fr=ap(kk,i,j)
            asum=asum+fr
            pj_sim=fr*di+pj_sim
            p=pindex(p)
        enddo
      enddo
    enddo
    pj_sim=pj_sim/asum !Note: it is necessary to confirm this code (must introduce any antenna pattern scaling method?)
    call duality_map(pj-w2*pj_sim,pj_new(i,j),real(w/(w-1)))!/asum
enddo
enddo
out_val=pj_new
enddo
! print *, "w=", w
! print *, "w2=",w2
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

SUBROUTINE init_conv(grid_x,grid_y,val,mask_grid,out_val,n_grid,nx,ny,wsize,ap,nx_ap,ny_ap,int_ap,res,fwhm)
IMPLICIT NONE
INTEGER :: i,ii,j,jj,k,kk,l,nz,grid_i,grid_j,isum
INTEGER :: i2,j2,i3,j3,p,dum,ip,jp
INTEGER(4),INTENT(IN) :: n_grid,nx,ny,wsize,nx_ap,ny_ap
REAL(4),INTENT(IN) :: int_ap,res,fwhm
INTEGER(4),INTENT(IN),DIMENSION(1:nx,1:ny) :: mask_grid
REAL(4),INTENT(IN),DIMENSION(1:n_grid)  :: val, grid_x, grid_y
REAL(4),INTENT(IN),DIMENSION(nx_ap,ny_ap) :: ap
INTEGER(4) :: findex(1:nx,1:ny),pindex(1:n_grid),uindex(1:n_grid)
REAL(4) :: di,pj,pj_sim,uij
REAL(4),INTENT(OUT),DIMENSION(1:nx,1:ny):: out_val

REAL(4) :: sigma,scale_radius,window_radius,dist,adis
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

scale_radius=fwhm/2
sigma=scale_radius*(2/2.35482) !2*sqrt(2ln2)=2.35482, sigma=sigma/2 for gaussian beam, fwhm/sqrt(2ln2) / 2
window_radius=(sigma*8)/2.

!-------------------initialization----------------------
out_val(:,:)=0
k=0
do j =1,ny
  do i =1,nx
    val_sum=0
    time_sum=0
    asum=0
    isum=0

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
          ip=nint((real(i)-grid_x(p))*res/int_ap)+(nx_ap-1)/2+1
          jp=nint((real(j)-grid_y(p))*res/int_ap)+(ny_ap-1)/2+1
          
          if((ip.gt.0).and.(ip.le.nx_ap)) then
          if((jp.gt.0).and.(jp.le.ny_ap)) then          
            fr=ap(ip,jp)
            asum=asum+fr
            isum=isum+1
            val_sum=fr*val(p)+val_sum
          endif
          endif

          endif
          p=pindex(p)
        enddo
      enddo
    enddo

    if(isum.lt.2) then
      out_val(i,j)=9.9E33
    else
      out_val(i,j)=val_sum/asum
    endif
 enddo
enddo



END SUBROUTINE init_conv

END MODULE Deblur