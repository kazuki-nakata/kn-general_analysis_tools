
  subroutine WRITE_TIF(infile,data,width,length,band)
    implicit none
    integer :: imax,jmax,jmin
    integer(4) :: width,length,band
    real(4) :: data(width,length,band)
    real(4) :: data2(width,length,band)    
    integer :: irec_entry_n,ifd  
    integer(2),parameter :: entry_n=16
    integer(4),parameter :: strip_n=2
    integer(1),parameter :: byte_order1=73
    integer(4),parameter :: compress=1
    integer(1),parameter :: version1=42
    integer(1),parameter :: version2=0
    integer(2),parameter :: dummy_bitpersample=32
    integer(2),parameter :: dummy_dataformart=3    !1=unsignint 2=signedint 3=float
    integer(2) :: bitpersample(band), dataformat(band)
    integer(2) :: dtty(entry_n),tag(entry_n)
    integer(4) :: ctfd(entry_n),dtfd(entry_n)
    integer(2) :: two_int
    integer(1) :: one_int
    integer(4) :: sin_int
    integer(4) :: wlb,hwlb,hwlbt1,hwlbt2,hwlbt3,hwlbt4
    real(4) :: fou_real
    real(8) :: eig_real
    integer(4) :: photometric
    integer(4) :: planarconfig
    integer :: NX,NY,NX2,NY2,start1,start2,finish1,finish2,iostat,kmax
    integer :: i, j, t,i2,j2,k,ii,jj,irec,irec2,irec3,irec4
    character(120) :: infile
    character(30) :: str0G30
    integer(4) :: bytenumber
    integer(4) :: stripoffset(strip_n),stripbytecount(strip_n)

    !      print *,'hi'

    do k =1,band
       do j =1,length
          do i =1,width
             data2(i,length-j+1,k)=data(i,j,k)
          enddo
       enddo
    enddo


    if(band.eq.1) then
       photometric=1         !0 = black mode
       planarconfig=1
    else
       photometric=2       !3 = color map
       planarconfig=1      !1 = pixel
    endif


    bitpersample(:)=dummy_bitpersample
    dataformat(:)=dummy_dataformart

    bytenumber=bitpersample(1)/8

    wlb=int(dble(width)*length*band*bytenumber)    
    hwlb=wlb+8      
    hwlbt1=hwlb
    hwlbt2=hwlbt1+band*2
    hwlbt3=hwlbt2+strip_n*4
    hwlbt4=hwlbt3+strip_n*4
    ifd=hwlbt4+band*2

    if(band.eq.1) then
       hwlbt1=bitpersample(1)
       hwlbt2=hwlb
       hwlbt3=hwlbt2+strip_n*4
       hwlbt4=1
       ifd=hwlbt3+strip_n*4
    endif

    do i =1,strip_n
       stripoffset(i)=(i-1)*width*length*band*bytenumber/strip_n+8
       stripbytecount(i)=width*length*band*bytenumber/strip_n
       !        print *,stripoffset(i),stripbytecount(i)
    enddo


    tag=(/254,  256,   257,  258,  259,  262,  273,  274,  277,  278,  279,  282, 283, 296, 284, 339/)
    dtty=(/4,    4,     4,    3,    3,    3,    4,    3,    3,    3,    4,     3,  3,   3,   3, 3/)
    ctfd=(/1,    1,     1,   band,    1,    1, strip_n, 1,    1,    1, strip_n,  1,  1,   1,   1, band/)
    dtfd(1:8)=(/0,  width,length,hwlbt1,compress, photometric,hwlbt2,-1441857535/)
    dtfd(9:16)= (/band,length/strip_n,hwlbt3,  1,  1,   2, planarconfig,hwlbt4/)


    open(21,file=infile,  &
         status='replace',form='unformatted',access='stream')


    write(21) byte_order1,byte_order1,version1,version2,ifd


    do i2 =1,strip_n
       imax=width
       jmax=stripbytecount(i2)/bytenumber* i2/width/band
       jmin=stripbytecount(i2)/bytenumber*(i2-1)/width/band+1
       irec=0
       !      do k =1, band
       do j =jmin, jmax
          !          print *,jmin,jmax,j
          do i =1, imax
             do k =1, band     
                !             two_int=data2(i,j,k)
                write(21) data2(i,j,k)
                !                        print *,two_int,data(i,j,k),i,j,k
             enddo
          enddo
       enddo
       !    print *, imax,jmin,jmax,band
    enddo

    if(band.gt.1) then
       do i =1, band
          write(21) bitpersample(i)
       enddo
    endif

    do i =1, strip_n
       write(21) stripoffset(i)
    enddo
    do i =1, strip_n
       write(21) stripbytecount(i)
    enddo

    if(band.gt.1) then
       do i =1, band
          write(21) dataformat(i)
          print *,dataformat(i)
       enddo
    endif

    write(21) entry_n

    do i =1, int(entry_n)
       write(21) tag(i),dtty(i),ctfd(i),dtfd(i)
    enddo


    print *,'File size is' , 12*entry_n+2+strip_n*2*4+width*length*band*2+8+band*2
  endsubroutine WRITE_TIF

  subroutine open_tif(infile,data,width,length)
    implicit none
    integer :: imax,jmax
    integer(4) :: width,length
    !      real, allocatable,intent(out) :: data(:,:)
    real :: data(width,length),dum(width,length)
    integer :: iglat,iglon,irec_entry_n,ifd
    integer :: strip_n
    integer(2) :: entry_n
    integer(2) :: dtty,atag
    integer(4) :: ctfd,dtfd,tag
    integer(2) :: dummy
    integer(4) :: photometric
    integer(4),allocatable :: stripoffset(:),stripbytecount(:)
    integer(4) :: planarconfig
    integer :: NX,NY,NX2,NY2,start1,start2,finish1,finish2,iostat
    integer :: i, j, t,i2,j2,k,ii,jj,irec,irec2,irec3,irec4
    character(120) :: infile
    INTEGER(1) :: kai(2)
    character(30) :: str0G30          



    open(21,file=infile,  &
         status='old',form='unformatted',access='stream')
    read(21,pos=5) ifd
    read(21,pos=ifd+1) entry_n
    do i =1, int(entry_n)
       read(21,pos=ifd+3+(i-1)*12) atag,dtty,ctfd,dtfd
       !       print *,atag,dtty,ctfd,dtfd
       if(atag.lt.0) then
          tag=int(256*256+atag)
       else
          tag=int(atag)
       endif

       if(tag.eq.256) then
          width=dtfd
       elseif(tag.eq.257) then
          length=dtfd
       elseif(tag.eq.262) then
          photometric=dtfd
       elseif(tag.eq.273) then
          ctfd=ctfd
          strip_n=ctfd
          allocate(stripoffset(ctfd))

          do i2 =1,ctfd
             read(21,pos=dtfd+1+(i2-1)*4) stripoffset(i2)
          enddo

       elseif(tag.eq.279) then

          allocate(stripbytecount(ctfd))
          do i2 =1,ctfd
             read(21,pos=dtfd+1+(i2-1)*4) stripbytecount(i2)
          enddo

       elseif(tag.eq.284) then
          planarconfig=dtfd

       elseif(tag.eq.50713) then
          read(21,pos=dtfd+1) dummy
          read(21,pos=dtfd+1) dummy
       endif

    enddo


    j2=0
    do i2 =1,strip_n
       imax=width
       jmax=stripbytecount(i2)/width/2
       irec=0
       do j =1, jmax
          j2=j2+1
          do i =1, imax
             irec=1+irec
             read(21,pos=stripoffset(i2)+1+(irec-1)*2) dummy
             if(dummy.lt.0) then
                dum(i,j2)=real(256*256+dummy)
             else
                dum(i,j2)=real(dummy)
             endif
          enddo
       enddo
    enddo

    do j =1,length
       do i =1,width
          data(i,j)=dum(i,length-j+1)
       enddo
    enddo


  endsubroutine open_tif



  subroutine read_size_tif(infile,width,length)
    implicit none
    integer :: imax,jmax
    integer :: iglat,iglon,irec_entry_n,ifd
    integer(2) :: entry_n
    integer(2) :: dtty,atag
    integer(4) :: ctfd,dtfd,tag
    integer(2) :: dummy
    integer(4) :: width,length
    integer(4) :: photometric
    integer(4),allocatable :: stripoffset(:),stripbytecount(:)
    integer(4) :: planarconfig
    integer :: NX,NY,NX2,NY2,start1,start2,finish1,finish2,iostat
    integer :: i, j, t,i2,j2,k,ii,jj,irec,irec2,irec3,irec4
    character(120) :: infile
    INTEGER(1) :: kai(2)
    character(30) :: str0G30          

    open(21,file=infile,  &
         status='old',form='unformatted',access='stream')
    read(21,pos=5) ifd
    read(21,pos=ifd+1) entry_n

    do i =1, int(entry_n)

       read(21,pos=ifd+3+(i-1)*12) atag,dtty,ctfd,dtfd
       if(atag.lt.0) then
          tag=int(256*256+atag)
       else
          tag=int(atag)
       endif

       if(tag.eq.256) width=dtfd
       if(tag.eq.257) length=dtfd

    enddo

  endsubroutine read_size_tif

  subroutine READ_EXIF_TIF(infile)
    implicit none
    integer :: imax,jmax
    integer :: iglat,iglon,irec_entry_n,ifd
    integer(2) :: entry_n
    integer(2) :: dtty,atag
    integer(4) :: ctfd,dtfd,tag
    integer(2) :: dummy
    integer(4) :: width,length
    integer(4) :: photometric
    integer(4),allocatable :: stripoffset(:),stripbytecount(:)
    integer(4) :: planarconfig
    integer :: NX,NY,NX2,NY2,start1,start2,finish1,finish2,iostat
    integer :: i, j, t,i2,j2,k,ii,jj,irec,irec2,irec3,irec4
    character(120) :: infile
    INTEGER(1) :: kai(2)
    character(30) :: str0G30          

    open(21,file=infile,  &
         status='old',form='unformatted',access='stream')
    read(21,pos=5) ifd
    read(21,pos=ifd+1) entry_n

    do i =1, int(entry_n)
       read(21,pos=ifd+3+(i-1)*12) atag,dtty,ctfd,dtfd

       if(atag.lt.0) then
          tag=int(256*256+atag)
       else
          tag=int(atag)
       endif

       !    print *,tag,dtty,ctfd,dtfd,entry_n

       if(tag.eq.34665) exifIFD=dtfd
       if(tag.eq.50714) then
          blacklevel=0
          do j =1,ctfd
             read(21,pos=dtfd+1+(j-1)*2) val
             blacklevel=real(val)/4+blacklevel
          enddo
       endif

    enddo

    !----------------------------------------------------    
    read(21,pos=exifIFD+1) entry_n
    !   print *,entry_n,exifIFD,ifd

    do i =1, int(entry_n)
       read(21,pos=exifIFD+3+(i-1)*12) atag,dtty,ctfd,dtfd
       if(atag.lt.0) then
          tag=int(256*256+atag)
       else
          tag=int(atag)
       endif
       !     print *,tag,dtty,ctfd,dtfd,entry_n
       if(tag.eq.33434) read(21,pos=dtfd+1) etime_u,etime_l
       if(tag.eq.33437) read(21,pos=dtfd+1) fnum_u,fnum_l
       if(tag.eq.37386) read(21,pos=dtfd+1) focal_length

    enddo
    !    print *,etime_u,etime_l,fnum_u,fnum_l,focal_length

  endsubroutine READ_EXIF_TIF
