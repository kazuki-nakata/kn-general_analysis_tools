MODULE segmentation
IMPLICIT NONE
INTEGER(4),PARAMETER :: null = -32767
REAL(4),PARAMETER :: Undef = 9.9E33
REAL(4),PARAMETER :: mask_val = 999
CONTAINS

SUBROUTINE region_merging(input,output,threshold,seed_x,seed_y,mask,isize,jsize,ksize,lsize)
   USE OBJECT_MODULE
   IMPLICIT NONE
   INTEGER(4),INTENT(IN) :: jsize,isize,lsize,ksize
   REAL(4),INTENT(OUT):: output(1:isize,1:jsize)
   REAL(4),INTENT(IN) :: input(1:isize,1:jsize)
   INTEGER(4),INTENT(IN) :: seed_x(lsize),seed_y(lsize)
   REAL(4),INTENT(IN) ::threshold(ksize)
   REAL(4),INTENT(IN) :: mask
   INTEGER(4) :: i, j, i2,j2,i3,j3,k,kk,ii,jj,l,irec,num
   integer(4), allocatable :: pairs1(:),pairs2(:)
   integer(4) :: ipos_8c(4),jpos_8c(4)
   integer(4) :: ipos_4c(2),jpos_4c(2)  
   REAL(4) :: img(1:isize,1:jsize,1)
   REAL(4) :: segment(1:isize,1:jsize)
   integer(4) :: pmax_uf,r1,r2,ui
   REAL(4) :: img1,img2
   INTEGER(4),PARAMETER::iNULL=0 

   irec=0
   img(1:isize,1:jsize,1)=input
   call init_unionfind(isize,jsize,img,1,mask,num)
   !------------connected system-------------------------------
   ipos_8c=(/-1,0,1,1/)
   jpos_8c=(/1,1,1,0/)
   ipos_4c=(/1,0/)
   jpos_4c=(/0,1/)
   allocate(pairs1(smax_uf*4),pairs2(smax_uf*4))
   pairs1(:)=-999;pairs2(:)=-999
   irec=0
   do j=1,jsize
   do i=1,isize
      !--------------------8-connected system-------------------------
      if(ptrt(i,j).eq.0) cycle
      do kk=1,4
         if((i+ipos_8c(kk).le.0).or.(i+ipos_8c(kk).gt.isize)) cycle
         if((j+jpos_8c(kk).le.0).or.(j+jpos_8c(kk).gt.jsize)) cycle           
         if(ptrt(i+ipos_8c(kk),j+jpos_8c(kk)).ne.0) then
            irec=1+irec
            pairs1(irec)=ptrt(i,j)
            pairs2(irec)=ptrt(i+ipos_8c(kk),j+jpos_8c(kk))
         endif
      enddo
   enddo
   enddo

   output(:,:)=9.9E33

   do k= 1,ksize

   pmax_uf=irec
   do i=1,pmax_uf
   r1=find_uf(pairs1(i))
   r2=find_uf(pairs2(i))
   img1=get_property_uf(r1,1)
   img2=get_property_uf(r2,1)
   if((img1.le.threshold(k)).and.(img2.le.threshold(k))) then
      call union_uf(pairs1(i),pairs2(i),ui)
   endif
   enddo

   do l =1,lsize
   segment(:,:)=9.9E33
   if((seed_x(l).le.0).and.(seed_x(l).gt.isize)) cycle
   if((seed_y(l).le.0).and.(seed_y(l).gt.jsize)) cycle
   if (ptrt(seed_x(l),seed_y(l)).eq.0) cycle
   jj=find_uf(ptrt(seed_x(l),seed_y(l)))
   ! print*,shape(ptrt),seed_x(l),seed_y(l)
   ! print *,ptrt(seed_x(l),seed_y(l)),jj,get_nchi_sll(jj)
   ! print *, get_ipos_uf(jj),get_jpos_uf(jj)


   if(jj.eq.get_nchi_sll(jj)) cycle
   do while (jj.ne.iNULL)
      i2=get_ipos_uf(jj);j2=get_jpos_uf(jj)
      segment(i2,j2)=threshold(k)
      jj=get_nchi_sll(jj)
   enddo

   do j=1,jsize
      do i=1,isize
         if(segment(i,j).ne.threshold(k)) cycle
         if(output(i,j).ne.9.9E33) cycle
         output(i,j)=threshold(k)
      enddo
   enddo


   enddo
   enddo

   call deallocate_uf
   deallocate(pairs1,pairs2)

   endsubroutine region_merging
endmodule segmentation
