MODULE feature_extraction
  USE Basic_processing
  USE util_shadow_module
  USE OBJECT_MODULE2
  USE UTIL_TEXTURE_MODULE
  IMPLICIT NONE
  INTEGER, PARAMETER :: MAX_IO = 64

CONTAINS

  
SUBROUTINE calculate_average_and_std(data_i,band,jsize,isize,ave,std)
  IMPLICIT NONE
  INTEGER :: i,j,k,ii,jj
  REAL(8):: xrate,yrate,xoffs,yoffs,num
  INTEGER(4),INTENT(IN) :: isize,jsize,band
  REAL(8), INTENT(OUT)  :: ave(2),std(2)
  REAL(4), INTENT(IN) :: data_i(band,jsize,isize)  
  
  ave(:)=0
  std(:)=0
  do k =1,2
     num = 0
     do j =1,jsize
        do i=1,isize
           if(data_i(k,j,i).eq.0) cycle
!              print *,data_i(k,j,i)
              ave(k)=255./20.*(data_i(k,j,i)-60)+ave(k)
              num = num + 1
!              print *,data_i(k,j,i),num,ave(k)
        enddo
     enddo
!     print *,ave(k),num
     ave(k)=ave(k)/num
  enddo
  
!  print *,ave(1),ave(2)
  
  do k =1,2
     num = 0
     do j =1,jsize
        do i=1,isize
           if(data_i(k,j,i).eq.0) cycle
!              print *,ave(k),data_i(k,j,i)
              std(k)=(255./20.*(data_i(k,j,i)-60)-ave(k))**2+std(k)
              num = num + 1
        enddo
     enddo
     std(k)=sqrt(std(k)/num)
  enddo
  
!  print *,std(1),std(2)
  
END SUBROUTINE calculate_average_and_std

  SUBROUTINE ndvi_and_shadow(data_i,data_o,band,jsize,isize,cbit)
  
  INTEGER :: i,j,k,l,ii,jj
  REAL, PARAMETER :: threshold_b= 0.4
  INTEGER,PARAMETER :: band1=3
  INTEGER,PARAMETER :: band2=4
  INTEGER,PARAMETER :: oband=6
  INTEGER, INTENT(IN) :: cbit
  INTEGER(4),INTENT(IN) :: isize,jsize,band
  REAL(4), INTENT(OUT) :: data_o(oband,jsize,isize)
  REAL(4), INTENT(IN) :: data_i(band,jsize,isize)
  REAL(4) :: data(isize,jsize,band)
  REAL(4) :: ndvi(isize,jsize),shadow(isize,jsize)
  REAL(4) :: alpha,beta,gamma,tau,norm
  REAL(4) :: dummy(isize,jsize,oband)
  INTEGER(4) :: mask(isize,jsize) 

  !---------------INITIALIZATION-----------------------------

  CALL reshape_backforward_3D(data_i(1:band,1:jsize,1:isize), &
       data(1:isize,1:jsize,1:band),band,jsize,isize)

  print *, band,jsize,isize,oband,band
  norm=2**(cbit)
  do j =1,jsize
     do i=1,isize
        if(data(i,j,2).eq.0) then
           mask(i,j)=999
        else
           data(i,j,:)=data(i,j,:)/norm
        endif
     enddo
  enddo
  print *,norm
  
  alpha=14
  beta=0.5
  gamma=2.2
  tau=10
  !--------------MAKE MASK DATA USING LAND COVER DATA-----------------------
  CALL Index(data(1:isize,1:jsize,1:band),ndvi(1:isize,1:jsize),isize,jsize,band,band1,band2,0)

  print *,'finish ndvi'
  CALL MAKE_SHADOW_PROB_MAP_1(data,isize,jsize,band,shadow,alpha,beta,gamma,tau)

  print *, 'hello'
  dummy(:,:,1:4)=data(:,:,:)
  dummy(:,:,5)=ndvi(:,:)
  dummy(:,:,6)=shadow(:,:)

!!$  do j =1,jsize
!!$     do i=1,isize
!!$        print *,data(i,j,:)
!!$     enddo
!!$  enddo
  
  
  CALL reshape_backforward_3D(dummy,data_o,isize,jsize,oband)

END SUBROUTINE ndvi_and_shadow

SUBROUTINE superpixel_features(data_o,data_i,sp_i,band,jsize,isize,cbit)
  USE util_texture_module   
  IMPLICIT NONE 
  INTEGER :: i,j,k,l,ii,jj,kk,i2,j2,i3,j3,i4,j4,irec,ui
  REAL, PARAMETER :: threshold_b= 0.4
  INTEGER,PARAMETER :: band1=3
  INTEGER,PARAMETER :: band2=4
  INTEGER,PARAMETER :: oband=6
  INTEGER, INTENT(IN) :: cbit
  INTEGER(4),INTENT(IN) :: isize,jsize,band
  REAL(4), INTENT(OUT)  :: data_o(oband,jsize,isize)
  REAL(4), INTENT(IN) :: data_i(band,jsize,isize)
  REAL(4), INTENT(IN) :: sp_i(jsize,isize)
  REAL(4) :: data(isize,jsize,band)
  REAL(4) :: sp(isize,jsize)
  REAL(4) :: ndvi(isize,jsize),shadow(isize,jsize)
  REAL(4) :: alpha,beta,gamma,tau,norm
  REAL(4) :: dummy(isize,jsize,oband)
  INTEGER(4) :: mask(isize,jsize) 
  integer(4) :: pmax_uf
  integer(4) :: ipos_8c(4),jpos_8c(4)
  integer(4) :: ipos_4c(2),jpos_4c(2)
  REAL(4) :: q  
  INTEGER(4) :: irec1,irec2
  INTEGER(4) :: dummy1(1000000),dummy2(1000000)
  real(4), allocatable::diffmax(:)
  INTEGER(4),ALLOCATABLE::edge(:,:,:)
  integer(4), allocatable :: pairs1(:),pairs2(:)
  REAL(4), ALLOCATABLE:: segment(:,:)
  INTEGER(4),PARAMETER::iNULL=0
  real(4) :: glcm_feature
  real(4) :: refdata(10000),adjdata(10000)
  real(4),allocatable :: refdata2(:),adjdata2(:)  
  real(4),allocatable :: homogeneity(:)
  integer(4) :: nbin,scale
  !---------------INITIALIZATION-----------------------------
  allocate(edge(isize,jsize,2))
  allocate(segment(isize,jsize))
  
  CALL reshape_backforward_3D(data_i(1:band,1:jsize,1:isize), &
       data(1:isize,1:jsize,1:band),band,jsize,isize)
  CALL reshape_backforward_2D_REAL(sp_i(1:jsize,1:isize), &
       sp(1:isize,1:jsize),jsize,isize)

  print *, band,jsize,isize,oband,band
  norm=2**(cbit)
  do j =1,jsize
     do i=1,isize
        if(data(i,j,2).eq.0) then
           mask(i,j)=999
           sp(i,j)=-999
        else
           data(i,j,:)=data(i,j,:)/norm
        endif
     enddo
  enddo
  print *,norm

  alpha=14
  beta=0.5
  gamma=2.2
  tau=10
  !--------------MAKE MASK DATA USING LAND COVER DATA-----------------------
  CALL Index(data(1:isize,1:jsize,1:band),ndvi(1:isize,1:jsize),isize,jsize,band,band1,band2,0)

  print *,'finish ndvi'
  CALL MAKE_SHADOW_PROB_MAP_1(data,isize,jsize,band,shadow,alpha,beta,gamma,tau)

  print *,'-----------statistical_region_merging---------------'
  !-------------initialization of uf----------------------------
  print *, 'initialize!'
  call init_unionfind(isize,jsize,sp(:,:),1,-999.,smax_uf)     
  !------------connected system-------------------------------
  ipos_8c=(/-1,0,1,1/)
  jpos_8c=(/1,1,1,0/)
  ipos_4c=(/1,0/)
  jpos_4c=(/0,1/)
  allocate(pairs1(smax_uf*4),pairs2(smax_uf*4))
  pairs1(:)=-999;pairs2(:)=-999
  print *,isize,jsize,smax_uf
  print *, 'make pixel pairs'
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
  !-------------------------------------------------------------------
  print *, "connect pixels"
  pmax_uf=irec
  irec=0
  edge(:,:,:)=0
  allocate(diffmax(1:pmax_uf))
  print *,pmax_uf
  do i=1,pmax_uf
     i2=get_ipos_uf(pairs1(i))
     j2=get_jpos_uf(pairs1(i))
     i3=get_ipos_uf(pairs2(i))
     j3=get_jpos_uf(pairs2(i))
     
     diffmax(i)=abs(get_property_uf(pairs1(i),1)-get_property_uf(pairs2(i),1))
     
     if(diffmax(i).eq.0) then
!        print *,i,pmax_uf
        call union_uf(pairs1(i),pairs2(i),ui)
        edge(i2,j2,1)=0
        edge(i2,j2,2)=0
     else
        if(edge(i3,j3,1).eq.0) then
           edge(i2,j2,1)=i3
           edge(i2,j2,2)=j3
        endif
     endif
  enddo

!!$  segment(:,:)=999.
!!$  do j=1,jsize
!!$     do i=1,isize
!!$        if(mask(i,j).eq.999) cycle
!!$!        print *,i,j,ii,sp(i,j)
!!$        ii=find_uf(ptrt(i,j))
!!$!        print *,i,j,ii,sp(i,j)
!!$           segment(i,j)=get_property_uf(ii,1)
!!$     enddo
!!$  enddo

  DEALLOCATE(diffmax,pairs1,pairs2)
  
  print *,'identification of first child node'
  !---------identification of first child node --------
  do i =1, smax_uf
     if(get_chi_dll(i).eq.0) then
        ii=i;exit
     endif
  enddo
  irec1=0
  irec2=0
  print *,'calc texture features'
  
  irec=0
  do while (ii.ne.smax_uf+1)
     irec1=1+irec1     
     dummy1(irec1)=find_uf(ii)
     ii=get_par_dll(ii)
  enddo

  print *,irec1
  !----------merging between residual similar color regions-------
  ALLOCATE(homogeneity(irec1))
  homogeneity(:)=-999
  scale=100
  nbin=100

  
  do i=1,irec1
     ii=dummy1(i)

     jj=get_nchi_sll(ii)
!        i2=get_ipos_uf(jj)
!        j2=get_jpos_uf(jj)     
!     print *,i,irec1,ii,jj,sp(i2,j2),i2,j2,mask(i2,j2),data(i2,j2,1),get_size_uf(jj)
     irec=0
     if(get_size_uf(jj).eq.1) cycle
     do while (jj.ne.iNULL)
        i2=get_ipos_uf(jj)
        j2=get_jpos_uf(jj)
        !        if(edge(i2,j2,1).ne.0) then
!      print *,i,irec1,ii,jj,sp(i2,j2),i2,j2,mask(i2,j2),data(i2,j2,1)       
        do j3=1,3
           do i3=1,3
              i4=i2+i3-2
              j4=j2+j3-2
              if((i4.lt.1).or.(i4.gt.isize)) cycle
              if((j4.lt.1).or.(j4.gt.jsize)) cycle
              if((i3.eq.2).and.(j3.eq.2)) cycle
              if(edge(i4,j4,1).eq.0) cycle
              irec=1+irec
              refdata(irec)=data(i2,j2,1)
              adjdata(irec)=data(i4,j4,1)
           enddo
        enddo
!        endif
        jj=get_nchi_sll(jj)
     enddo
     allocate(refdata2(irec),adjdata2(irec))  
     refdata2(:)=refdata(1:irec);adjdata2(:)=adjdata(1:irec)     
     call MAKE_GLCM_FEATURES_1D(refdata2,adjdata2,irec,glcm_feature,nbin,scale)      
     homogeneity(i)=glcm_feature      
     deallocate(refdata2,adjdata2)
!   print *,i,irec,homogeneity(i)
enddo

  print *,'make final image'
do i =1, irec1
   ii=dummy1(i)
   i2=get_ipos_uf(ii)
   j2=get_jpos_uf(ii)
     if(get_size_uf(ii).eq.1) cycle   
      do while(ii.ne.iNULL)
         dummy(get_ipos_uf(ii),get_jpos_uf(ii),4)=homogeneity(i)
         ii=get_nchi_sll(ii)
      enddo
enddo

  dummy(:,:,1:3)=data(:,:,1:3)
!  dummy(:,:,4)=(:,:)
  dummy(:,:,5)=ndvi(:,:)
  dummy(:,:,6)=shadow(:,:)

  CALL reshape_backforward_3D(dummy,data_o,isize,jsize,oband)

end subroutine SUPERPIXEL_FEATURES


END MODULE feature_extraction
