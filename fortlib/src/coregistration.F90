module coreg_tool
  use similarity
  use Numerical_Module
  implicit none
contains

  subroutine calc_offset_type1(img1,img2,nx,ny,band,ix,iy,sim,dx,dy,nz,sim_type,ndw,dst,nsw,sst,disdep)
    !img1: master,img2: slave
    !sim_type
    !ndw: window size for calculating similarity
    !dst,sst: stride 
    !nsw: search window size
    implicit none
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

  subroutine lucas_kanade(gradx,grady,gradt,mask,weight,nx,ny,band,ix,iy,dx,dy,dc,nz,ndw)
    !img1: master,img2: slave
    !sim_type
    !ndw: window size for calculating similarity
    implicit none
    integer(4) :: i,j,k,kk,l
    integer(4),intent(in) :: nx,ny,band,nz,ndw
    integer(4) :: imin,imax,jmin,jmax,nsize,m_count
    integer(4) :: imin2,imax2,jmin2,jmax2
    integer(4),intent(IN) :: ix(nz),iy(nz)
    real(8),intent(IN) :: gradx(nx,ny,band)
    real(8),intent(IN) :: grady(nx,ny,band)
    real(8),intent(IN) :: gradt(nx,ny,band)
    real(8),intent(IN) :: weight(ndw,ndw,band)
    integer(4),intent(IN) :: mask(nx,ny)
    integer(4) :: info
    REAL(8),intent(OUT) :: dx(nz),dy(nz),dc(nz)
    REAL(8) ::dxy(2),AAt_inv(2,2),AAt(2,2)
    REAL(8), allocatable :: AAt_inv2(:,:),gradxy_sub(:,:),gradt_sub(:),gradxy_sub2(:,:)
         ! gradxy_sub(1,:)=reshape(gradx(imin:imax,jmin:jmax,1:band),[nsize])
         ! gradxy_sub(2,:)=reshape(grady(imin:imax,jmin:jmax,1:band),[nsize])
         ! gradt_sub=reshape(gradt(imin:imax,jmin:jmax,1:band),[nsize])
         ! mask_sub(:)=reshape(mask(imin:imax,jmin:jmax),[nsize])

    do k =1,nz
       !the range of search image
      if(mask(ix(k),iy(k)).eq.999) then
         dx(k)=999.
         dy(k)=999.
         dc(k)=999.
         cycle
      else
         imin=ix(k)-(ndw-1)/2; imax=ix(k)+(ndw-1)/2
         jmin=iy(k)-(ndw-1)/2; jmax=iy(k)+(ndw-1)/2

         if(imin<1) imin=1
         if(imax>nx) imax=nx
         if(jmin<1) jmin=1
         if(jmax>ny) jmax=ny

         ! if((imin<1).or.(imax>nx)) then
         !    dx(k)=999.
         !    dy(k)=999.
         !    dc(k)=999.
         !    cycle
         ! endif

         ! if((jmin<1).or.(jmax>ny)) then
         !    dx(k)=999.
         !    dy(k)=999.
         !    dc(k)=999.
         !    cycle
         ! endif


         nsize=0
         do j= jmin,jmax
            do i=imin,imax
               if(mask(i,j).ne.999) then
                  nsize=1+nsize
               endif
            enddo
         enddo    

         if(nsize.le.2) then
            dx(k)=999.
            dy(k)=999.
            dc(k)=999.
            cycle
         else
            nsize=nsize*band
            allocate(AAt_inv2(2,nsize),gradxy_sub(2,nsize),gradt_sub(nsize),gradxy_sub2(2,nsize))
            kk=0
            do l=1,band
            do j= jmin,jmax
               do i=imin,imax
                  if(mask(i,j).ne.999) then
                        kk=1+kk
                        gradxy_sub(1,kk)=gradx(i,j,l)
                        gradxy_sub(2,kk)=grady(i,j,l)
                        gradxy_sub2(1,kk)=gradx(i,j,l)*weight(i-imin+1,j-jmin+1,l)
                        gradxy_sub2(2,kk)=grady(i,j,l)*weight(i-imin+1,j-jmin+1,l)
                        gradt_sub(kk)=gradt(i,j,l)
                  endif
               enddo
            enddo
            enddo

            call matrix_product(gradxy_sub2,transpose(gradxy_sub), AAt, 2, 2, nsize)
            call get_inverse_matrix(AAt,AAt_inv,info,2)
            call matrix_product(AAt_inv, gradxy_sub2, AAt_inv2, 2, nsize, 2)
            call matrix_vector_product(AAt_inv2, gradt_sub, dxy, 2, nsize)
            dx(k)=real(dxy(1))
            dy(k)=real(dxy(2))
            dc(k)=nsize/band
            deallocate(AAt_inv2,gradxy_sub,gradt_sub,gradxy_sub2)
         endif
      endif
    enddo

    
  endsubroutine lucas_kanade


subroutine iterate_lucas_kanade(gradx0,grady0,val10,val20,mask0,weight,nx,ny,band, &
&uini,vini,ix,iy,dx,dy,dc,eval,nz,ndw,xs,ys,thres)
    !img1: master,img2: slave
    !sim_type
    !ndw: window size for calculating similarity
    implicit none
    integer(4) :: i,j,k,kk,l,ix2,iy2
    integer(4),intent(in) :: nx,ny,band,nz,ndw
    integer(4) :: imin,imax,jmin,jmax,nsize,m_count
    integer(4) :: imin2,imax2,jmin2,jmax2
    integer(4),intent(IN) :: ix(nz),iy(nz)
    REAL(4),intent(IN) :: thres
    real(4),intent(IN) :: uini(nz),vini(nz)
    real(8),intent(IN) :: gradx0(nx,ny,band)
    real(8),intent(IN) :: grady0(nx,ny,band)
    real(8),intent(IN) :: val10(nx,ny,band),val20(nx,ny,band)
    real(8) :: gradt1(nx,ny,band)
    real(8),intent(IN) :: weight(ndw,ndw,band)
    real(4),intent(IN) :: xs,ys
    integer(4),intent(IN) :: mask0(nx,ny)
    INTEGER(4),PARAMETER :: margin=5
    integer(4) :: info
    integer(4) :: mask(-margin:nx+margin,-margin:ny+margin)
    real(8) :: val1(-margin:nx+margin,-margin:ny+margin,band),val2(-margin:nx+margin,-margin:ny+margin,band)
    real(8) :: gradx(-margin:nx+margin,-margin:ny+margin,band),grady(-margin:nx+margin,-margin:ny+margin,band)
    REAL(4),intent(OUT) :: dx(nz),dy(nz),dc(nz),eval(nz)!,test(2,nz)
    REAL(8) ::dxy(2),AAt_inv(2,2),AAt(2,2)
    REAL(8), allocatable :: AAt_inv2(:,:),gradxy_sub(:,:),gradt_sub(:),gradxy_sub2(:,:)
    REAL(8) :: x_floor,y_floor,wx,wy,asum,dum,xs0,ys0
   mask(:,:)=0
   gradx(:,:,:)=0
   grady(:,:,:)=0
   val1(:,:,:)=0
   val2(:,:,:)=0
   mask(1:nx,1:ny)=mask0(1:nx,1:ny)
   gradx(1:nx,1:ny,:)=gradx0(1:nx,1:ny,:)
   grady(1:nx,1:ny,:)=grady0(1:nx,1:ny,:)
   val1(1:nx,1:ny,:)=val10(1:nx,1:ny,:)
   val2(1:nx,1:ny,:)=val20(1:nx,1:ny,:)
   ! test(:,:)=9.9E33

    do k =1,nz
       !the range of search image
      if((mask(ix(k),iy(k)).eq.0).or.(nint(uini(k)).eq.999).or.(nint(vini(k)).eq.999)) then
         dx(k)=999.
         dy(k)=999.
         dc(k)=999.
         eval(k)=999.
         ! test(:,k)=999.
         cycle
      else
         imin=ix(k)-(ndw-1)/2; imax=ix(k)+(ndw-1)/2
         jmin=iy(k)-(ndw-1)/2; jmax=iy(k)+(ndw-1)/2

         if(imin<1) imin=1
         if(imax>nx) imax=nx
         if(jmin<1) jmin=1
         if(jmax>ny) jmax=ny


         nsize=0
         do j= jmin,jmax
            do i=imin,imax
               if(mask(i,j).eq.1) then
                  nsize=1+nsize
               endif
            enddo
         enddo    

         if(nsize.le.2) then
            dx(k)=999.
            dy(k)=999.
            dc(k)=999.
            eval(k)=999.
            ! test(:,k)=999.
            cycle
         else
            xs0=-dble(uini(k))*dble(xs)
            ys0=-dble(vini(k))*dble(ys)
            x_floor = floor(xs0)
            y_floor = floor(ys0)
            ix2 = INT(x_floor)
            iy2 = INT(y_floor)
            wx = xs0 - x_floor
            wy = ys0 - y_floor

            nsize=nsize*band
            allocate(gradxy_sub(2,nsize),gradt_sub(nsize),gradxy_sub2(2,nsize))
            kk=0
            ! test(:,k)=0

            do l=1,band
            do j= jmin,jmax
               do i=imin,imax

                  if(mask(i,j).eq.1) then
                        kk=1+kk
                        asum = (1.0 - wx) * (1.0 - wy) * dble(mask(i+ix2, j+iy2)) + &
                              wx * (1.0 - wy) * dble(mask(i+ix2 + 1, j+iy2)) + &
                              (1.0 - wx) * wy * dble(mask(i+ix2, j+iy2+1)) + &
                              wx * wy * dble(mask(i+ix2+1, j+iy2+1))
                        if(asum.le.0.000001) then
                           kk=kk-1
                           cycle
                        endif

                        ! dum = (1.0 - wx) * (1.0 - wy) * gradx(i+ix2, j+iy2,l)* mask(i+ix2, j+iy2) + &
                        !       wx * (1.0 - wy) * gradx(i+ix2 + 1, j+iy2,l)* mask(i+ix2 + 1, j+iy2) + &
                        !       (1.0 - wx) * wy * gradx(i + ix2, j+ iy2 + 1,l)* mask(i + ix2, j+ iy2 + 1) + &
                        !       wx * wy * gradx(i + ix2 + 1, j + iy2 + 1,l)* mask(i + ix2+ 1, j+ iy2 + 1)
                        ! gradxy_sub(1,kk)=dum/asum

                        ! dum = (1.0 - wx) * (1.0 - wy) * grady(i+ix2, j+iy2,l)* mask(i+ix2, j+iy2) + &
                        !       wx * (1.0 - wy) * grady(i+ix2 + 1, j+iy2,l)* mask(i+ix2 + 1, j+iy2) + &
                        !       (1.0 - wx) * wy * grady(i + ix2, j+ iy2 + 1,l)* mask(i + ix2, j+ iy2 + 1) + &
                        !       wx * wy * grady(i + ix2 + 1, j + iy2 + 1,l)* mask(i + ix2+1, j+ iy2 + 1)
                        ! gradxy_sub(2,kk)=dum/asum

                        gradxy_sub(1,kk)=gradx(i,j,l)
                        gradxy_sub(2,kk)=grady(i,j,l)
                        gradxy_sub2(1,kk)=gradxy_sub(1,kk)*weight(i-imin+1,j-jmin+1,l)
                        gradxy_sub2(2,kk)=gradxy_sub(2,kk)*weight(i-imin+1,j-jmin+1,l)

                        dum = (1.0 - wx) * (1.0 - wy) * val2(i+ix2, j+iy2,l)*dble(mask(i+ix2, j+iy2)) + &
                              wx * (1.0 - wy) * val2(i+ix2+1, j+iy2,l)* dble(mask(i+ix2+1, j+iy2)) + &
                              (1.0 - wx) * wy * val2(i+ix2, j+iy2+1, l)* dble(mask(i+ix2, j+iy2+1)) + &
                              wx * wy * val2(i+ix2+1, j+iy2+1, l)* dble(mask(i+ix2+1, j+iy2+1))
                        gradt_sub(kk)=dum/asum-val1(i,j,l)
                        ! gradt_sub(kk)=val2(i,j,l)-dum/asum

                        ! if(val2(i+ix2, j+iy2,l)*0.0/=0.0) test(1,k)=i+ix2
                        ! if(val2(i+ix2+1, j+iy2,l)*0.0/=0.0) test(1,k)=i+ix2+1
                        ! if(val2(i+ix2, j+iy2+1, l)*0.0/=0.0) test(1,k)=i+ix2
                        ! if(val2(i+ix2+1, j+iy2+1, l)*0.0/=0.0) test(1,k)=i+ix2+1
                        ! if(val2(i+ix2, j+iy2,l)*0.0/=0.0) test(2,k)=j+iy2
                        ! if(val2(i+ix2+1, j+iy2,l)*0.0/=0.0) test(2,k)=j+iy2+1
                        ! if(val2(i+ix2, j+iy2+1, l)*0.0/=0.0) test(2,k)=j+iy2
                        ! if(val2(i+ix2+1, j+iy2+1, l)*0.0/=0.0) test(2,k)=j+iy2+1


                  endif
               enddo
            enddo
            enddo

            allocate(AAt_inv2(2,kk))
            call matrix_product(gradxy_sub2(:,1:kk),transpose(gradxy_sub(:,1:kk)), AAt, 2, 2, nsize)
            call get_inverse_matrix(AAt,AAt_inv,info,2)
            call matrix_product(AAt_inv, gradxy_sub2(:,1:kk), AAt_inv2, 2, nsize, 2)
            call matrix_vector_product(AAt_inv2, gradt_sub(1:kk), dxy, 2, nsize)

            dx(k)=real(dxy(1)+uini(k))
            dy(k)=real(dxy(2)+vini(k))
            if(dx(k)>thres) dx(k)=thres
            if(dx(k)<-thres) dx(k)=-thres
            if(dy(k)>thres) dy(k)=thres
            if(dy(k)<-thres) dy(k)=-thres
            dc(k)=nsize/band
            ! test(k)=info
            

            eval(k)=0.
            do i =1,kk
               eval(k)=eval(k)+gradt_sub(i)**2/real(kk)
            enddo
            eval(k)=sqrt(eval(k))
            
            deallocate(AAt_inv2,gradxy_sub,gradt_sub,gradxy_sub2)
            ! dx(k)=wx
            ! dy(k)=wy
            ! dc(k)=0
            ! eval(k)=0
         endif
      endif
    enddo

    
  endsubroutine iterate_lucas_kanade


subroutine calc_warp_image_errors(val10,val20,mask0,nx,ny,band,uini,vini,ix,iy,eval,nz,ndw,xs,ys)
    !(Ixu+Iyv-It+e)**2を計算するものではない。これは非線形項を考慮できないため、ワープをしてから評価するのが筋。
    !入力するuv0が0から推定されたuvでは、(It1-warp(It2))**2=(It+e+e_nonlinear)**2となる。イタレーションし、推定値uv=0に収束した場合、(It+e)**2となる。
    !上記のうち、uv=/0でない場合は、本来は(Ixu+Iyv-It+e+e_nonlinear)**2を計算すべき。
    !img1: master,img2: slave
    !sim_type
    !ndw: window size for calculating similarity
    implicit none
    integer(4) :: i,j,k,kk,l,ix2,iy2
    integer(4),intent(in) :: nx,ny,band,nz,ndw
    integer(4) :: imin,imax,jmin,jmax,nsize,m_count
    integer(4) :: imin2,imax2,jmin2,jmax2
    integer(4),intent(IN) :: ix(nz),iy(nz)
    real(4),intent(IN) :: uini(nz),vini(nz)
    real(8),intent(IN) :: val10(nx,ny,band),val20(nx,ny,band)
    real(4),intent(IN) :: xs,ys
    integer(4),intent(IN) :: mask0(nx,ny)
    INTEGER(4),PARAMETER :: margin=5
    integer(4) :: mask(-margin:nx+margin,-margin:ny+margin)
    real(8) :: val1(-margin:nx+margin,-margin:ny+margin,band),val2(-margin:nx+margin,-margin:ny+margin,band)
    REAL(4),intent(OUT) :: eval(nz)
    REAL(8), allocatable :: gradt_sub(:)
    REAL(8) :: x_floor,y_floor,wx,wy,asum,dum,xs0,ys0
   mask(:,:)=0
   val1(:,:,:)=0
   val2(:,:,:)=0
   mask(1:nx,1:ny)=mask0(1:nx,1:ny)
   val1(1:nx,1:ny,:)=val10(1:nx,1:ny,:)
   val2(1:nx,1:ny,:)=val20(1:nx,1:ny,:)

    do k =1,nz
       !the range of search image
      if((mask(ix(k),iy(k)).eq.0).or.(nint(uini(k)).eq.999).or.(nint(vini(k)).eq.999)) then
         eval(k)=999.
         ! test(:,k)=999.
         cycle
      else
         imin=ix(k)-(ndw-1)/2; imax=ix(k)+(ndw-1)/2
         jmin=iy(k)-(ndw-1)/2; jmax=iy(k)+(ndw-1)/2

         if(imin<1) imin=1
         if(imax>nx) imax=nx
         if(jmin<1) jmin=1
         if(jmax>ny) jmax=ny


         nsize=0
         do j= jmin,jmax
            do i=imin,imax
               if(mask(i,j).eq.1) then
                  nsize=1+nsize
               endif
            enddo
         enddo    

         if(nsize.le.2) then
            eval(k)=999.
            cycle
         else
            xs0=-dble(uini(k))*dble(xs)
            ys0=-dble(vini(k))*dble(ys)
            x_floor = floor(xs0)
            y_floor = floor(ys0)
            ix2 = INT(x_floor)
            iy2 = INT(y_floor)
            wx = xs0 - x_floor
            wy = ys0 - y_floor

            nsize=nsize*band
            allocate(gradt_sub(nsize))
            kk=0

            do l=1,band
            do j= jmin,jmax
               do i=imin,imax

                  if(mask(i,j).eq.1) then
                        kk=1+kk
                        asum = (1.0 - wx) * (1.0 - wy) * dble(mask(i+ix2, j+iy2)) + &
                              wx * (1.0 - wy) * dble(mask(i+ix2 + 1, j+iy2)) + &
                              (1.0 - wx) * wy * dble(mask(i+ix2, j+iy2+1)) + &
                              wx * wy * dble(mask(i+ix2+1, j+iy2+1))
                        if(asum.le.0.000001) then
                           kk=kk-1
                           cycle
                        endif
                        dum = (1.0 - wx) * (1.0 - wy) * val2(i+ix2, j+iy2,l)*dble(mask(i+ix2, j+iy2)) + &
                              wx * (1.0 - wy) * val2(i+ix2+1, j+iy2,l)* dble(mask(i+ix2+1, j+iy2)) + &
                              (1.0 - wx) * wy * val2(i+ix2, j+iy2+1, l)* dble(mask(i+ix2, j+iy2+1)) + &
                              wx * wy * val2(i+ix2+1, j+iy2+1, l)* dble(mask(i+ix2+1, j+iy2+1))
                        gradt_sub(kk)=(dum/asum-val1(i,j,l))!*weight(i-imin+1,j-jmin+1,l)
                  endif
               enddo
            enddo
            enddo

            eval(k)=0.
            do i =1,kk
               eval(k)=eval(k)+gradt_sub(i)**2/real(kk)
            enddo
            eval(k)=sqrt(eval(k))            
            deallocate(gradt_sub)
         endif
      endif
    enddo

    
  endsubroutine calc_warp_image_errors



  subroutine lucas_kanade_multi(gradx1,grady1,gradt1,gradx2,grady2,gradt2,mask1,mask2,nx1,ny1,nx2,ny2,band,ix,iy,dx,dy,dc,nz,ndw)
    !img1: master,img2: slave
    !sim_type
    !ndw: window size for calculating similarity
    implicit none
    integer(4) :: i,j,k,kk,l
    integer(4),intent(in) :: nx1,ny1,nx2,ny2,band,nz,ndw
    integer(4) :: imin,imax,jmin,jmax,nsize,m_count
    integer(4) :: imin2,imax2,jmin2,jmax2
    real(4),intent(IN) :: ix(nz),iy(nz)
    real(8),intent(IN) :: gradx1(nx1,ny1,band),gradx2(nx2,ny2,band)
    real(8),intent(IN) :: grady1(nx1,ny1,band),grady2(nx2,ny2,band)
    real(8),intent(IN) :: gradt1(nx1,ny1,band),gradt2(nx2,ny2,band)
    integer(4),intent(IN) :: mask1(nx1,ny1),mask2(nx2,ny2)
    integer(4) :: info
    REAL(8),intent(OUT) :: dx(nz),dy(nz),dc(nz)
    REAL(8) ::dxy(2),AAt_inv(2,2),AAt(2,2)
    REAL(8), allocatable :: AAt_inv2(:,:),gradxy_sub(:,:),gradt_sub(:)
         ! gradxy_sub(1,:)=reshape(gradx(imin:imax,jmin:jmax,1:band),[nsize])
         ! gradxy_sub(2,:)=reshape(grady(imin:imax,jmin:jmax,1:band),[nsize])
         ! gradt_sub=reshape(gradt(imin:imax,jmin:jmax,1:band),[nsize])
         ! mask_sub(:)=reshape(mask(imin:imax,jmin:jmax),[nsize])

    do k =1,nz
       !the range of search image
      if((mask1(int(ix(k)),int(iy(k))).eq.999).or.(mask2(int(ix(k)*2),int(iy(k)*2)).eq.999))  then
         dx(k)=999.
         dy(k)=999.
         dc(k)=999.
         cycle
      else
         imin=int(ix(k)-(ndw-1)/2); imax=int(ix(k)+(ndw-1)/2)
         jmin=int(iy(k)-(ndw-1)/2); jmax=int(iy(k)+(ndw-1)/2)
         imin2=imin*2;imax2=imax*2
         jmin2=jmin*2;jmax2=jmax*2

         if((imin<1).or.(imax>nx1)) then
            dx(k)=999.
            dy(k)=999.
            dc(k)=999.
            cycle
         endif

         if((jmin<1).or.(jmax>ny1)) then
            dx(k)=999.
            dy(k)=999.
            dc(k)=999.
            cycle
         endif


         nsize=0
         do j= jmin,jmax
            do i=imin,imax
               if(mask1(i,j).ne.999) then
                  nsize=1+nsize
               endif
            enddo
         enddo


         do j= jmin2,jmax2
            do i=imin2,imax2
               if(mask2(i,j).ne.999) then
                  nsize=1+nsize
               endif
            enddo
         enddo


         if(nsize.le.2) then
            dx(k)=999.
            dy(k)=999.
            dc(k)=999.
            cycle
         else
            nsize=nsize*band
            allocate(AAt_inv2(2,nsize),gradxy_sub(2,nsize),gradt_sub(nsize))
            kk=0
            do l=1,band
            do j= jmin,jmax
               do i=imin,imax
                  if(mask1(i,j).ne.999) then
                        kk=1+kk
                        gradxy_sub(1,kk)=2*gradx1(i,j,l)
                        gradxy_sub(2,kk)=2*grady1(i,j,l)
                        gradt_sub(kk)=gradt1(i,j,l)
                  endif
               enddo
            enddo
            enddo
            do l=1,band
            do j= jmin2,jmax2
               do i=imin2,imax2
                  if(mask2(i,j).ne.999) then
                        kk=1+kk
                        gradxy_sub(1,kk)=gradx2(i,j,l)
                        gradxy_sub(2,kk)=grady2(i,j,l)
                        gradt_sub(kk)=gradt2(i,j,l)
                  endif
               enddo
            enddo
            enddo

            call matrix_product(gradxy_sub,transpose(gradxy_sub), AAt, 2, 2, nsize)
            call get_inverse_matrix(AAt,AAt_inv,info,2)
            call matrix_product(AAt_inv, gradxy_sub, AAt_inv2, 2, nsize, 2)
            call matrix_vector_product(AAt_inv2, gradt_sub, dxy, 2, nsize)
            dx(k)=real(dxy(1))
            dy(k)=real(dxy(2))
            dc(k)=nsize/band
            deallocate(AAt_inv2,gradxy_sub,gradt_sub)
         endif
      endif
    enddo

    
  endsubroutine lucas_kanade_multi


endmodule coreg_tool
