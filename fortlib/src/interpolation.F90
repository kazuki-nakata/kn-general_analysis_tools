MODULE interpolation
  IMPLICIT NONE
  SAVE !-----PRIVATE

  !----------type of interpolation_grid-------------------
  INTEGER :: NY_IP, NX_IP
  REAL, parameter :: lon_res=0.1
  REAL, parameter :: lat_res=0.05
  REAL, parameter :: lon_start=0.0
  REAL, parameter :: lat_start=-80.0
  REAL, parameter :: lon_end=360.0
  REAL, parameter :: lat_end=-55.0
  INTEGER, parameter :: num_variable=2 
  INTEGER, parameter :: NULL = 0
  TYPE interpolation_data
     real(4) :: lat,lon
     real(4) :: ptr,cnt
     real(4) :: property(num_variable)
     real(4) :: mask
  END TYPE interpolation_data
  TYPE(interpolation_data), ALLOCATABLE :: data_ip(:,:)

  !---------type of original data------------------------
  !----------gridded data----------------
  Integer, parameter :: NX_origin=1264
  integer, parameter :: NY_origin=1328
  integer, parameter :: NUM_origin=NX_origin*NY_origin
  !----------no-gridded data----------------
  !  integer, parameter :: NUM_origin=19321
  TYPE original_data
     INTEGER(4) :: id
     INTEGER(4) :: ipos, jpos, kpos
     REAL(4)    :: xpos, ypos, zpos
     REAL(4)    :: lat,lon
     REAL(4)    :: property(num_variable)
     REAL(4)    :: mask
     REAL(4)    :: area
     INTEGER(4)    :: ptr_next
     REAL(8), POINTER :: buffer(:,:)
  END TYPE original_data

  TYPE(original_data), ALLOCATABLE :: data_or(:)

  !------------------------------------------------------
contains

  subroutine set_interpolation_point_ll
    integer :: i,j,k
    NX_IP=int((lon_end-lon_start)/lon_res)
    NY_IP=int((lat_end-lat_start)/lat_res)+1
    print *,NX_IP,NY_IP
    ALLOCATE(data_ip(1:NX_IP,NY_IP))

    data_ip(:,:)%ptr = NULL
    data_ip(:,:)%cnt = NULL

    do k=1,num_variable
       data_ip(:,:)%property(k) = NULL
    enddo

    do j=1,NY_IP
       do i=1,NX_IP
          data_ip(i,j)%lon=real(i-1)*lon_res+lon_start
          data_ip(i,j)%lat=real(j-1)*lat_res+lat_start
       enddo
    enddo
  end subroutine set_interpolation_point_ll

  subroutine interpolation_gaussian
    integer :: i,j,k,ph,ii,jj,i3,j3
    real(8) :: dist
    real(8) :: pi,dora,alpha,sum,fr,dlon
    real(4) :: lat1,lon1,lon2,lat2,area
    real(4) :: mask
    real(8) :: dist_min
    real(4) :: mask_dist_min
    integer:: mask_num_all,mask_num
    real(8),parameter :: re=6378D0
    real(8),parameter :: e_folding=6D0
    real(8),parameter :: limit=18D0

    pi=acos(-1.)
    dora=pi/180.

    data_ip(:,:)%property(1) = 0.0D0
    data_ip(:,:)%property(2) = 0.0D0
    !    data_ip(:,:)%land = 0.0D0

    

    do j = 1,NY_IP
       do i = 1,NX_IP

          sum=0.
          dist_min=1000.
          mask_dist_min=1000.
          mask_num_all=0
          mask_num=0

          do jj=1,3
          do ii=1,3

             i3=i-2+ii
             j3=j-2+jj

             if(i3.lt.1) then
             i3=NX_IP
             elseif(i3.gt.NX_IP) then
             i3=1
             endif
             if(j3.lt.1) then
             j3=NY_IP
             elseif(j3.gt.NY_IP) then
             j3=1
             endif

          ph = data_handle_at(i3,j3)

          DO WHILE (ph /= NULL)
             lat1=data_ip(i3,j3)%lat
             lon1=data_ip(i3,j3)%lon
             lon2=get_data_lon(ph)
             lat2=get_data_lat(ph)
             mask=get_data_mask(ph)
             area=get_data_area(ph)
!             print *,area,lat2,lon2
             dlon=abs(lon2-lon1)
             if(dlon.gt.180.)then
                dlon=abs(360.-dlon)
             endif
             alpha=(cos(dora*lat1))*(cos(dora*lat2))*(cos(dora*dlon))+ (sin(dora*lat1))*(sin(dora*lat2))
             dist=re*acos(alpha)

             if (dist.le.limit) then

               if(dist < 4.0) then
                   mask_num_all=mask_num_all+1
                   if(mask.eq.0) mask_num=mask_num+1
                endif

                if(dist < dist_min) then
                   dist_min=dist
                   mask_dist_min=mask
                endif

                fr=exp(-dist**2./(e_folding**2))*area
                sum=sum+fr*mask

!                print *,fr,sum,mask,dist,dist_min,mask_dist_min
                DO k=1,num_variable
                   data_ip(i,j)%property(k)=data_ip(i,j)%property(k)+fr*get_data_property(ph,k)*mask
                end do
             endif

             ph = data_handle_next(ph)

          END DO

          end do
          end do

          if((sum.eq.0.).or.(mask_num_all.eq.mask_num)) then
!          if(sum==0) then
!          if((sum.eq.0.).or.(mask_dist_min.eq.0.0)) then

          DO k=1,num_variable
             data_ip(i,j)%property(k)=0.
          end do
          else
          DO k=1,num_variable
             data_ip(i,j)%property(k)=data_ip(i,j)%property(k)/sum
          end do

          endif
       end do
    end do

  end subroutine interpolation_gaussian



  subroutine set_pointer_on_llproj
    integer :: i,j,irec

    DO i=1,NUM_origin
       data_or(i)%id=i
       data_or(i)%ipos=nint((data_or(i)%lon-lon_start)/lon_res)+1
       data_or(i)%jpos=nint((data_or(i)%lat-lat_start)/lat_res)+1
       if((data_or(i)%ipos.lt.1).or.(data_or(i)%ipos.gt.NX_IP)) cycle
       if((data_or(i)%jpos.lt.1).or.(data_or(i)%jpos.gt.NY_IP)) cycle
       data_ip(data_or(i)%ipos,data_or(i)%jpos)%cnt= &
            1+data_ip(data_or(i)%ipos,data_or(i)%jpos)%cnt
       data_or(i)%ptr_next=data_ip(data_or(i)%ipos,data_or(i)%jpos)%ptr
       data_ip(data_or(i)%ipos,data_or(i)%jpos)%ptr=i
    ENDDO

    open(67,file='test.dat',  &
         status='replace',form='unformatted',access='direct',recl=4)
    irec=0
    do j=1,NY_IP
    do i=1,NX_IP
       irec=1+irec      
       write(67,rec=irec) real(data_ip(i,j)%cnt)
!       print *,data_or(i)%lat
    enddo
    enddo
  end subroutine set_pointer_on_llproj
!!$
!!$  SUBROUTINE set_data_property(ph, prop, value)
!!$    INTEGER, INTENT(IN) :: ph
!!$    INTEGER, INTENT(IN) :: prop
!!$    REAL(8), INTENT(IN) :: value
!!$    data_or(ph)%property(prop) = value
!!$  END SUBROUTINE set_data_property
!!$

  REAL(4) FUNCTION get_data_area(ph)
    INTEGER, INTENT(IN) :: ph
    get_data_area = data_or(ph)%area
  END FUNCTION get_data_area

  REAL(4) FUNCTION get_data_lon(ph)
    INTEGER, INTENT(IN) :: ph
    get_data_lon = data_or(ph)%lon
  END FUNCTION get_data_lon


  REAL(4) FUNCTION get_data_lat(ph)
    INTEGER, INTENT(IN) :: ph
    get_data_lat = data_or(ph)%lat
  END FUNCTION get_data_lat

  REAL(4) FUNCTION get_data_mask(ph)
    INTEGER, INTENT(IN) :: ph
    get_data_mask = data_or(ph)%mask
  END FUNCTION get_data_mask

  REAL(4) FUNCTION get_data_property(ph,prop)
    INTEGER, INTENT(IN) :: ph
    INTEGER, INTENT(IN) :: prop
    get_data_property = data_or(ph)%property(prop)
  END FUNCTION get_data_property

  !----------------------------------------------------------------------
  INTEGER PURE FUNCTION data_handle_at(i, j)
    INTEGER, INTENT(IN) :: i
    INTEGER, INTENT(IN) :: j
    data_handle_at = data_ip(i,j)%ptr
  END FUNCTION data_handle_at
  !--------------------------------------------------------------------------------------
  INTEGER PURE FUNCTION data_handle_next(ph)
    INTEGER, INTENT(IN) :: ph

    IF (ph == NULL) THEN
       data_handle_next = NULL
       RETURN
    END IF
    data_handle_next = data_or(ph)%ptr_next
  END FUNCTION data_handle_next


  subroutine read_original_data_point(lafile,lofile,arfile)
    implicit none
    integer :: i
    integer :: irec
    integer :: iglat,iglon,isize
    character:: lafile*66,lofile*66,arfile*66

    ALLOCATE(data_or(1:NUM_origin))

    open(57,file=lafile,  &
         status='old',form='unformatted',access='direct',recl=4)
    irec=0
    do i=1,NUM_origin
       irec=1+irec
       read(57,rec=irec) iglat
       data_or(i)%lat=iglat/1.e5
       !       print *,data_or(i)%lat
    enddo

    open(58,file=lofile &
         , status='old',form='unformatted',access='direct',recl=4)
    irec=0
    do i=1,NUM_origin
       irec=1+irec      
       read(58,rec=irec) iglon
       data_or(i)%lon=iglon/1.e5
       if(data_or(i)%lon.le.0.) then
          data_or(i)%lon=360.+data_or(i)%lon
       endif
       !       print *,data_or(i)%lon
    enddo

    open(59,file=arfile &
         , status='old',form='unformatted',access='direct',recl=4)
    irec=0
    do i=1,NUM_origin
       irec=1+irec      
       read(59,rec=irec) isize
       data_or(i)%area=real(isize)/1.e3
    enddo

  end subroutine read_original_data_point

  !--------------------------------------------------------
  subroutine read_original_data_binary(data_file)
    integer :: i,j
    character :: data_file*66
    open(55,file=data_file,form='unformatted',access='direct', &
         recl=4*num_origin,status='old')
    read(55,rec=1) data_or%mask
    read(55,rec=3) data_or%property(1)
    read(55,rec=4) data_or%property(2)

       do i=1,num_origin
          if((data_or(i)%mask.eq.9.9E33).or.(data_or(i)%mask.eq.-1)) then
!          print *,data_or(i)%mask,data_or(i)%property(1),data_or(i)%property(2)
          data_or(i)%mask=0
          else
          data_or(i)%mask=1
          endif
       enddo

  end subroutine read_original_data_binary
!------------------------------------------------------------------
  subroutine read_original_mask_binary(data_file)
    character :: data_file*66
    open(56,file=data_file,form='unformatted',access='direct', &
         recl=4*num_origin,status='old')
    read(56,rec=1) data_or%mask
!    data_or(:)%mask=1.
!    print *, data_or(:)%mask
  end subroutine read_original_mask_binary

!!$!------------------------------------------------------------
  subroutine write_interpolation_data_binary(data_file,n)
    implicit NONE
    integer :: i,j,k,irec,n
    character :: data_file*66
    SAVE
    open(61,file=data_file,form='unformatted',access='direct',recl=4,status='replace')
    irec=0
    do k =1,n
    do j=1,NY_IP
       do i=1,NX_IP
          irec=irec+1
          write(61,rec=irec) data_ip(i,j)%property(k)
       enddo
    enddo
    enddo
  end subroutine write_interpolation_data_binary

end module interpolation















