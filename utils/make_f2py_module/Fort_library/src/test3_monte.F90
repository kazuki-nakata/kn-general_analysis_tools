PROGRAM BeBuSi
  IMPLICIT NONE
  INTEGER :: i,j,ii,jj,k,kk,i2,i3,iii,iave,istd,n
  INTEGER(4) :: ii1,ii2,jj1,jj2,kx,ky,window
  INTEGER(4) :: inum
  INTEGER(4) :: mi,mj,isize2,jsize2,rad,nt,step      
  INTEGER(4) :: isize, jsize,band
  INTEGER(4), ALLOCATABLE :: minx2(:), miny2(:)  !tl 
  REAL(4):: bbp
  INTEGER(4),ALLOCATABLE :: minx(:), miny(:)
  REAL(4) :: ave,std,std2,x,y
  REAL(4),ALLOCATABLE :: z1(:),z2(:)
  REAL(8),PARAMETER :: disdep= 4!0.2  !dependency factor of pixels distance
  INTEGER(4),PARAMETER :: psize=1
  REAL(4) :: ai1,ai2,aj1,aj2,pi,dora,theta
  REAL(8) :: distance,maxdis,adis,dis_rgb,dis_xy
  real(4) :: bbp_max,asphalt,ave1,std1
  real(8),ALLOCATABLE :: dis(:) !TL
  real(4) :: sum,arad,bbp_ave
  integer(4) :: temp_sum  
  real(4) :: ndvi2,si2
  real(4) :: bbs(1000)
  real(4) :: bbs_ave,bbs_std
  !-------------for normalization----------------
  real(4) :: a,b
  INTEGER(4) :: inum2,imax,imin
  real(4),ALLOCATABLE :: list(:)
  open(70,file='lookuptable.txt',access='sequential',status='replace',form='formatted')
  
  pi=acos(-1.)
  dora=pi/180.

  std2=1
  
  do k=4,400,10
     n=k
     allocate(z1(n),z2(n))
     allocate(minx(n),minx2(n))
     allocate(dis(n))
     do istd=1,10
        std=real(istd)
        
        do iave=1,1
           ave=real(iave-1)

           std1=0
           ave1=0
           bbs_ave=0
           bbs(:)=0
           do inum=1,1000

              do i =1,n
                 call random_number(x)
                 call random_number(y)
                 z1(i)=std*(sqrt(-2*log(x))*cos(2*pi*y))+ave
                 !                 z2(i)=sqrt(-2*log(x))*sin(2*pi*y)
!                 print *,z1(i),x,ave,std
              enddo
              
              do i =1,n
                 call random_number(x)
                 call random_number(y)
                 z2(i)=std2*(sqrt(-2*log(x))*cos(2*pi*y))!-ave
              enddo

              bbp=0
              dis(:)=10000
              minx(:)=999
              minx2(:)=999             
              do ii=1,n
               distance=1000000.                
                 do iii=1,n!imax_tl
                    dis_rgb=(dble(z1(ii)-z2(iii)))**2
                    adis=dis_rgb
                    if(adis.lt.dis(iii)) then
                       minx2(iii)=ii
                       dis(iii)=adis
                    endif
                    if(adis.lt.distance) then
                       minx(ii) = iii
                       distance=adis
                    endif
                 enddo
              enddo

 
              do ii=1,n
                 if(minx(ii).eq.999) cycle
                 i2=minx(ii)
                 i3=minx2(i2)

                 if(ii.eq.i3) then
                    bbp=bbp+1
                 endif

              enddo
              bbs(inum)=bbp/real(n)
!              print *,bbs(inum)
              !              print *,bbp,bbp_ave
              bbs_ave=bbs(inum)+bbs_ave
!--------------------------------------------------------
!!$              do i=1,n
!!$                 ave1=ave1+z1(i)
!!$              enddo
!!$              ave1=ave1/real(n)
!!$              
!!$              do i=1,n
!!$                 std1=std1+(z1(i)-ave1)**2
!!$              enddo
!!$              std1=sqrt(std1/real(n))
!--------------------------------------------------------              
           enddo
           
           bbs_ave=bbs_ave/real(1000)
           bbs_std=0
           do inum=1,1000
              bbs_std=bbs_std+(bbs(inum)-bbs_ave)**2
           enddo
           bbs_std=sqrt(bbs_std/1000.)
           print *, ave,std,n,k,bbs_ave,bbs_std,bbs_ave+bbs_std*2
           write(70,*) k,int(std),bbs_ave+bbs_std*2
        enddo
     enddo
     deallocate(z1,z2)
     deallocate(minx,minx2)
     deallocate(dis)     
  enddo




END PROGRAM BEBUSI
