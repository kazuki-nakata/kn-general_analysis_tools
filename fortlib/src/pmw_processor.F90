MODULE PMW_Processor
IMPLICIT NONE
INTEGER(4),PARAMETER :: null = -32767
REAL(4),PARAMETER :: Undef = 9.9E33
REAL(4),PARAMETER :: mask_val = 999
CONTAINS

SUBROUTINE rSIR(grid_x,grid_y,vs,vb,vg,time,val,mask_grid,out_val,out_time,n_grid,nx,ny,&
wsize,ap,nx_ap,ny_ap,int_ap,res,fwhm,iterate)
!vs,vb,vg: ecef vector for s/c, obs_point(boresight), and obs_point()
IMPLICIT NONE
INTEGER :: i,ii,j,jj,k,kk,l,nz,grid_i,grid_j,isum
INTEGER :: i2,j2,i3,j3,p,dum
INTEGER(4),INTENT(IN) :: n_grid,nx,ny,wsize,nx_ap,ny_ap,iterate
! INTEGER(4),INTENT(IN) :: i_t,j_t
REAL(4),INTENT(IN) :: int_ap,res,fwhm
REAL(8),INTENT(IN),DIMENSION(1:nx,1:ny,3) :: vg
INTEGER(4),INTENT(IN),DIMENSION(1:nx,1:ny) :: mask_grid
REAL(4),INTENT(IN),DIMENSION(1:n_grid)  :: val, grid_x, grid_y, time
REAL(8),INTENT(IN),DIMENSION(1:n_grid,3) :: vs, vb
REAL(4),DIMENSION(nx,ny) :: pj_new
REAL(4),DIMENSION(n_grid) :: fi 
REAL(4),DIMENSION(n_grid,3) :: vi, vj, vk
REAL(4),DIMENSION(nx_ap,ny_ap) :: ap
REAL(8),DIMENSION(1,3) :: vdum
INTEGER(4) :: findex(1:nx,1:ny),pindex(1:n_grid),uindex(1:n_grid)
REAL(4) :: di,pj,pj_sim,uij
REAL(4),INTENT(OUT),DIMENSION(1:nx,1:ny):: out_val,out_time
! REAL(4),INTENT(OUT),DIMENSION(2000,6):: out_test

REAL(4) :: sigma,scale_radius,window_radius,dist,adis
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
scale_radius=fwhm/2
!sigma=scale_radius*sqrt(2/2.35482)
sigma=scale_radius*(2/2.35482) !2*sqrt(2ln2)=2.35482, sigma=sigma/2 for gaussian beam, fwhm/sqrt(2ln2) / 2
window_radius=(sigma*8)/2.

!-------------------initialization----------------------
out_val(:,:)=0
out_time(:,:)=0
k=0
! out_test(:,:)=0
do j =1,ny
  do i =1,nx
    val_sum=0
    time_sum=0
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
        dist=sqrt((real(i)-grid_x(p))**2+(real(j)-grid_y(p))**2)*res
        if (dist.le.window_radius) then
         call calc_local_az_el_angle(vi(p,:),vj(p,:),vk(p,:),vb(p,:)-vs(p,:),vb(p,:),vdum,az,el,int(1))
         fr=ap(nint(az(1)*180.0/pi/int_ap)+(nx_ap-1)/2+1,nint(el(1)*180.0/pi/int_ap)+(ny_ap-1)/2+1)
!          fr=exp(-dist**2./(sigma**2))     
        ! if((i.eq.i_t).and.(j.eq.j_t)) then
        !   k=k+1
        !   out_test(k,1)=grid_x(p)
        !   out_test(k,2)=grid_y(p)
        !   out_test(k,3)=exp(-dist**2./(sigma**2))
        !   out_test(k,4)=fr
        !   out_test(k,5)=az(1)*180.0/pi
        !   out_test(k,6)=el(1)*180.0/pi
          ! print *,sqrt((real(i)-grid_x(p))**2+(real(j)-grid_y(p))**2)*5000
          ! print *,i,j,ii,jj,az*180.0/3.14,el*180.0/3.14,fr
          ! print *,vi(p,:),vj(p,:),vk(p,:),vb(p,:)-vs(p,:),vb(p,:),vdum,az,el
        ! endif
          asum=asum+fr
          isum=isum+1
          val_sum=fr*val(p)+val_sum
          time_sum=time(p)+time_sum
          endif
          p=pindex(p)
        enddo
      enddo
    enddo

    if(isum.lt.4) then
      out_val(i,j)=Undef
      out_time(i,j)=Undef
    else
      out_val(i,j)=val_sum/asum
      out_time(i,j)=time_sum/real(isum)
    endif
 enddo
enddo


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
        if(out_val(i2,j2).eq.9.9E33) cycle
          dist=sqrt((real(i2)-grid_x(l))**2+(real(j2)-grid_y(l))**2)*res
          if (dist.le.window_radius) then
            vdum(1,:)=vg(i2,j2,:)
            call calc_local_az_el_angle(vi(l,:),vj(l,:),vk(l,:),vb(l,:)-vs(l,:),vb(l,:),vdum,az,el,int(1))
            fr=ap(nint(az(1)*180.0/pi/int_ap)+(nx_ap-1)/2+1,nint(el(1)*180.0/pi/int_ap)+(ny_ap-1)/2+1)
            ! fr=exp(-dist**2./(sigma**2))
            asum=asum+fr
            fi(l)=fr*out_val(i2,j2)+fi(l)
          endif
      enddo
    enddo
    fi(l)=fi(l)/asum
enddo


!------------------test----------------------
! out_val(:,:)=0
! out_time(:,:)=0

! do j =1,ny
!   do i =1,nx
!     isum=0
!     asum=0

!     if(mask(i,j).eq.999) cycle
    
!     do jj =1,wsize
!       do ii =1,wsize
!         i2=i+ii-(wsize+1)/2
!         j2=j+jj-(wsize+1)/2
!         if((i2.lt.1).or.(i2.gt.nx)) cycle
!         if((j2.lt.1).or.(j2.gt.ny)) cycle
!         p=findex(i2,j2)
!         if(p.eq.-1) cycle
!         do while (p>0)
!           dist=sqrt((real(i)-grid_x(p))**2+(real(j)-grid_y(p))**2)*res
!           if (dist.le.window_radius) then
!             fr=exp(-dist**2./(sigma**2))
!             asum=asum+fr
!             isum=isum+1
!             out_val(i,j)=fr*fi(p)+out_val(i,j)
!             out_time(i,j)=fr*time(p)+out_time(i,j)
!           endif
!             p=pindex(p)
!         enddo
!       enddo
!     enddo

!     if(isum.lt.4) then
!       out_val(i,j)=9.9E33
!       out_time(i,j)=9.9E33
!     else
!       out_val(i,j)=out_val(i,j)/asum
!       out_time(i,j)=out_time(i,j)/asum
!     endif
!  enddo
! enddo
!------------------------------kokomade-----------------------

pj_new(:,:)=9.9E33
do j = 1, ny
  do i = 1, nx
    if(mask_grid(i,j).eq.999) cycle
    if(out_val(i,j).eq.Undef) cycle
    pj_sim=0
    asum=0
    pj=out_val(i,j)
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
          dist=sqrt((real(i)-grid_x(p))**2+(real(j)-grid_y(p))**2)*res
          if (dist.le.window_radius) then
            di=sqrt(val(p)/fi(p))
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
out_val=pj_new

enddo

END SUBROUTINE rSIR


SUBROUTINE L1B_to_L3_WA_TYPE1(grid_x,grid_y,time,val,mask_grid,out_val,out_time,n_grid,nx,ny,wsize,fwhm,res,rm_outer)
!WA: weighted averaging
IMPLICIT NONE
INTEGER :: i,ii,j,jj,k,kk,l,nz,grid_i,grid_j,isum
INTEGER :: i2,j2,i3,j3,p,dum
INTEGER(4),INTENT(IN) :: n_grid,nx,ny,wsize
REAL(4),INTENT(IN) :: res,fwhm,rm_outer
INTEGER(4),INTENT(IN) :: mask_grid(1:nx,1:ny)
REAL(4),INTENT(IN) :: grid_x(1:n_grid),grid_y(1:n_grid),time(1:n_grid)
REAL(4),INTENT(IN) :: val(1:n_grid)
INTEGER(4) :: findex(1:nx,1:ny),pindex(1:n_grid),uindex(1:n_grid)
REAL(4) :: di,pj,pj_sim,uij
REAL(4),INTENT(OUT):: out_val(1:nx,1:ny),out_time(1:nx,1:ny)
REAL(4) :: sigma,scale_radius,window_radius,dist,adis,max_fr
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

scale_radius=fwhm/2
!sigma=scale_radius*sqrt(2/2.35482)
sigma=scale_radius*(2/2.35482) !2*sqrt(2ln2)=2.35482, sigma=sigma/2 for gaussian beam, fwhm/sqrt(2ln2) / 2
window_radius=(sigma*8)/2.

!-------------------initialization----------------------
out_val(:,:)=0
out_time(:,:)=0

do j =1,ny
  do i =1,nx
    isum=0
    asum=0
    max_fr=0
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
            out_val(i,j)=fr*val(p)+out_val(i,j)
            out_time(i,j)=fr*time(p)+out_time(i,j) 
            if(fr.ge.max_fr) max_fr=fr
          endif
            p=pindex(p)
        enddo
      enddo
    enddo

    ! if(isum.lt.4) then
    if(max_fr.lt.rm_outer) then
      out_val(i,j)=Undef
      out_time(i,j)=Undef
    else
      out_val(i,j)=out_val(i,j)/asum
      out_time(i,j)=out_time(i,j)/asum
    endif
 enddo
enddo

END SUBROUTINE L1B_to_L3_WA_TYPE1


SUBROUTINE l1B_to_l3_WA_TYPE2(grid_x,grid_y,vs,vb,vg,time,val,mask_grid,out_val,out_time,n_grid,nx,ny,&
wsize,ap,nx_ap,ny_ap,int_ap,res,fwhm)
!vs,vb,vg: ecef vector for s/c, obs_point(boresight), and obs_point()
IMPLICIT NONE
INTEGER :: i,ii,j,jj,k,kk,l,nz,grid_i,grid_j,isum
INTEGER :: i2,j2,i3,j3,p,dum
INTEGER(4),INTENT(IN) :: n_grid,nx,ny,wsize,nx_ap,ny_ap
REAL(4),INTENT(IN) :: int_ap,res,fwhm,rm_outer
INTEGER(4),INTENT(IN),DIMENSION(1:nx,1:ny) :: mask_grid
REAL(4),INTENT(IN),DIMENSION(1:n_grid)  :: val, grid_x, grid_y, time
REAL(8),INTENT(IN),DIMENSION(1:nx,1:ny,3) :: vg
REAL(8),INTENT(IN),DIMENSION(1:n_grid,3) :: vs, vb
REAL(4),DIMENSION(n_grid,3) :: vi, vj, vk
REAL(4),DIMENSION(nx_ap,ny_ap) :: ap
REAL(8),DIMENSION(1,3) :: vdum
INTEGER(4) :: findex(1:nx,1:ny),pindex(1:n_grid),uindex(1:n_grid)
REAL(4),INTENT(OUT),DIMENSION(1:nx,1:ny):: out_val,out_time

REAL(4) :: sigma,scale_radius,window_radius,dist,adis,max_fr
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
scale_radius=fwhm/2
!sigma=scale_radius*sqrt(2/2.35482)
sigma=scale_radius*(2/2.35482) !2*sqrt(2ln2)=2.35482, sigma=sigma/2 for gaussian beam, fwhm/sqrt(2ln2) / 2
window_radius=(sigma*8)/2.

out_val(:,:)=0
out_time(:,:)=0
k=0
! out_test(:,:)=0
do j =1,ny
  do i =1,nx
    val_sum=0
    time_sum=0
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
          time_sum=fr*time(p)+time_sum
          if(fr.ge.max_fr) max_fr=fr
          endif
          p=pindex(p)
        enddo
      enddo
    enddo

    if(max_fr.lt.rm_outer) then
      out_val(i,j)=Undef
      out_time(i,j)=Undef
    else
      out_val(i,j)=val_sum/asum
      out_time(i,j)=time_sum/asum
    endif
 enddo
enddo

END SUBROUTINE L1B_to_L3_WA_TYPE2



subroutine calc_latlon_from_L1B(lat,lon,olat,olon,imax,jmax,coregiA1,coregiA2)
 implicit none
 integer :: imax,jmax,oimax,ojmax
 real(4), intent(IN) :: lat(imax,jmax),lon(imax,jmax)
 real(4), intent(OUT) :: olat(imax/2,jmax),olon(imax/2,jmax)
 real(8),dimension(3) :: p1,p2,ex,ey,ez,pt
 real(8),dimension(3) :: eznume
 real :: theta,ezdeno,tmp

 integer :: i,j,o,e,k
 real*8 :: pi,degra,radeg

 real :: coregiA1,coregiA2

! imax=486; jmax=2039
!oimax=243; ojmax=2039
! coregiA1=1.40040; coregiA2=0.08020
oimax=imax/2
ojmax=jmax
 pi=acos(-1.)
 degra=pi/180. ; radeg=180./pi

 do j=1,ojmax
   do i=1,oimax
     o=2*i-1;e=2*i
     p1(:)=(/cos(lon(o,j)*degra)*cos(lat(o,j)*degra), &
           & sin(lon(o,j)*degra)*cos(lat(o,j)*degra), &
           & sin(lat(o,j)*degra)/)
     p2(:)=(/cos(lon(e,j)*degra)*cos(lat(e,j)*degra), &
           & sin(lon(e,j)*degra)*cos(lat(e,j)*degra), &
           & sin(lat(e,j)*degra)/)
     
     ex(:)=p1(:)
     call calc_cross_product(p1,p2,eznume)
     ezdeno=sqrt(dot_product(eznume,eznume))
     ez(:)=eznume(:)/ezdeno
     call calc_cross_product(ez,ex,ey)
     
     theta=acos(dot_product(p1,p2))
     pt(:)=cos(coregiA2*theta)*(cos(coregiA1*theta)*ex+sin(coregiA1*theta)*ey)&
          &+sin(coregiA2*theta)*ez
     
     olat(i,j)=asin(pt(3))*radeg
     olon(i,j)=atan2(pt(2),pt(1))*radeg
   enddo
 enddo

end subroutine calc_latlon_from_L1B

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

subroutine MAPLL(aii,ajj,LAT,LONG)
      implicit none
      REAL(4) :: X,Y,ALAT,ALONG,E,E2,CDR,PI,SLAT,MC,SGN
      REAL(4) :: RE,RHO,SL,T,TC,LAT,LONG
      REAL(4) :: aii,ajj

      SGN=1
      SLAT = 70.
      RE = 6378.273
      E2 = 0.006693883
      E =  0.081816153
!--------------------------------
!      CDR=57.29577951
      PI=3.141592654
      ALAT=abs(LAT)*PI/180.
      ALONG=(LONG+45.)*PI/180.

      IF (ABS(ALAT).LT.PI/2.) GOTO 250
      X=0.0
      Y=0.0
      GOTO 999
  250 CONTINUE
      T=TAN(PI/4.-ALAT/2.)/((1.-E*SIN(ALAT))/(1.+E*SIN(ALAT)))**(E/2.)
      IF (ABS(90.-SLAT).LT.1.E-5) THEN
      RHO=2.*RE*T/((1.+E)**(1.+E)*(1.-E)**(1.-E))**(1/2.)
      ELSE
      SL=SLAT*PI/180.
      TC=TAN(PI/4.-SL/2.)/((1.-E*SIN(SL))/(1.+E*SIN(SL)))**(E/2.)
      MC=COS(SL)/SQRT(1.0-E2*(SIN(SL)**2))
      RHO=RE*MC*T/TC
      END IF
      Y=-RHO*SGN*COS(SGN*ALONG)
      X= RHO*SGN*SIN(SGN*ALONG)
  999 CONTINUE

      aii=x !(x+3850.-8.0/2.)/8.0+1
      ajj=y !1400-(y+5350.-8.0/2.)/8.0

    end subroutine MAPLL

END MODULE PMW_Processor

