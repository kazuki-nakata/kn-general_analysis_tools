MODULE point_data
IMPLICIT NONE
INTEGER(4),PARAMETER :: null = -32767
REAL(4),PARAMETER :: Undef = 9.9E33
REAL(4),PARAMETER :: mask_val = 999
CONTAINS

SUBROUTINE weighted_mean_sigma(grid_x,grid_y,val,grid_xo,grid_yo,out_val,nx,ny,n_grid,n_grid_o,wsize,sigma_p)
!WA: weighted averaging
USE sensor_geometry
IMPLICIT NONE
INTEGER :: i,ii,j,jj,k,kk,l,nz,grid_i,grid_j,isum
INTEGER :: i2,j2,i3,j3,p,dum
INTEGER(4),INTENT(IN) :: n_grid,n_grid_o,nx,ny,wsize
REAL(4),INTENT(IN) :: sigma_p
REAL(4),INTENT(IN) :: grid_x(1:n_grid),grid_y(1:n_grid)
REAL(4),INTENT(IN) :: val(1:n_grid)
INTEGER(4) :: findex(1:nx,1:ny),pindex(1:n_grid),uindex(1:n_grid)
REAL(4) :: di,pj,pj_sim,uij
REAL(4),INTENT(IN) :: grid_xo(1:n_grid_o),grid_yo(1:n_grid_o)
REAL(4),INTENT(OUT):: out_val(1:n_grid_o)
REAL(4) :: scale_radius,window_radius,dist,adis,max_fr
REAL(8) :: fr,asum,out_val2

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

window_radius=(sigma_p*8)/2.

!-------------------initialization----------------------
out_val(:)=0

do k =1,n_grid_o
    i=nint(grid_xo(k))
    j=nint(grid_yo(k))
    isum=0
    asum=0
    out_val2=0
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
            fr=exp(-dist**2./(sigma_p**2))
            asum=asum+fr
            isum=isum+1
            out_val2=fr*val(p)+out_val2
          endif
            p=pindex(p)
        enddo
      enddo
    enddo

    if(isum.eq.0) then
      out_val(k)=Undef
    else
      out_val(k)=out_val2/asum
    endif
enddo

END SUBROUTINE weighted_mean_sigma


END MODULE point_data
