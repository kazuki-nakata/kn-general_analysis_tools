module coreg_tool
  use similarity
  implicit none
contains

  subroutine calc_offset_type1(img1,img2,nx,ny,band,ix,iy,sim,dx,dy,nz,sim_type,ndw,dst,nsw,sst,disdep)
    !img1: master,img2: slave
    !sim_type
    !ndw: window size for calculating similarity
    !dst,sst: stride 
    !nsw: search window size
    integer(4) :: i,j,k
    integer(4),intent(in) :: nx,ny,band,nz,sim_type,ndw,nsw,dst,sst
    integer(4) :: imin,imax,jmin,jmax
    integer(4) :: imin2,imax2,jmin2,jmax2
    integer(4),intent(IN) :: ix(nz),iy(nz)
    real(4),intent(IN) :: img1(nx,ny,band)
    real(4),intent(IN) :: img2(nx,ny,band)
    real(4) :: img_sub1(ndw,ndw,band),img_sub2(ndw,ndw,band)
    REAL(4),intent(OUT) :: sim(nz),dx(nz),dy(nz)
    REAL(4), INTENT(IN) :: disdep
    REAL(4) :: sim_max,dx0,dy0,sim0
    REAL(4) :: x(ndw,ndw)
    REAL(4) :: y(ndw,ndw)
   print *,"nx=",nx,"ny=",ny,"band=",band,"nz=",nz

   if(sim_type.eq.5) then
      call calc_rot_loc(x,y,ndw,ndw,real(0))
   endif

    do k =1,nz

       !the range of search image
      imin=ix(k)-(ndw-1)/2; imax=ix(k)+(ndw-1)/2
      jmin=iy(k)-(ndw-1)/2; jmax=iy(k)+(ndw-1)/2
      img_sub1(1:ndw,1:ndw,1:band)=img1(imin:imax,jmin:jmax,1:band)


      sim_max=-9.9E33
      do j =1, nsw,sst
         do i =1, nsw,sst
            imin2=imin+i-1-(nsw-1)/2
            imax2=imax+i-1-(nsw-1)/2
            jmin2=jmin+j-1-(nsw-1)/2
            jmax2=jmax+j-1-(nsw-1)/2

            img_sub2=img2(imin2:imax2,jmin2:jmax2,1:band)
            
            if(sim_type.eq.1) then
               call calculate_zncc(img_sub1,img_sub2,ndw,ndw,band,dst,sim0)
            elseif(sim_type.eq.2) then
               call calculate_ssd(img_sub1,img_sub2,ndw,ndw,band,dst,sim0)
               sim0=-sim0
            elseif(sim_type.eq.3) then
               call calculate_sad(img_sub1,img_sub2,ndw,ndw,band,dst,sim0)
               sim0=-sim0
            elseif(sim_type.eq.4) then
               call calculate_ncc(img_sub1,img_sub2,ndw,ndw,band,dst,sim0)
            elseif(sim_type.eq.5) then
               call calculate_bbs(img_sub1,x,y,img_sub2,x,y,ndw,ndw,band,dst,disdep,sim0)
            endif
            
            if(sim0.ge.sim_max) then
               dx0=real(i-1-(nsw-1)/2)
               dy0=real(j-1-(nsw-1)/2)
               sim_max=sim0
            endif
         enddo
      enddo

      sim(k)=sim_max
      dx(k)=dx0
      dy(k)=dy0

    enddo
  endsubroutine calc_offset_type1


endmodule coreg_tool
