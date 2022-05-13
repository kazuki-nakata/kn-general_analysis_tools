    subroutine open_tif(infile,data,width,length,band)          !8byte data
      implicit none
      integer :: imax,jmax
      integer(4) :: width,length,band
!      real, allocatable,intent(out) :: data(:,:)
      real(4) ::data(width,length,band)
!      real(8) :: dum(width,length)
      integer :: iglat,iglon,irec_entry_n,ifd
      integer :: strip_n,dataformat_true
      integer(2) :: entry_n
      integer(2) :: dtty,atag
      integer(4) :: ctfd,dtfd,tag
      integer(2) :: dummy
      real(4) :: dummy_dble
      integer(4) :: photometric
      integer(4),allocatable :: stripoffset(:),stripbytecount(:)
      integer(2),allocatable :: bitpersample(:),dataformat(:)      
      integer(4) :: planarconfig
      integer :: NX,NY,NX2,NY2,start1,start2,finish1,finish2,iostat
      integer :: i, j, t,i2,j2,k,ii,jj,irec,irec2,irec3,irec4,jt
      character(512) :: infile
      INTEGER(1) :: kai(2)
      character(30) :: str0G30
!      print *, 'hiiiiiiiiiiiiii'
      open(21,file=infile,  &
      status='old',form='unformatted',access='stream')
       read(21,pos=5) ifd
      read(21,pos=ifd+1) entry_n
      do i =1, int(entry_n)
!       print *,"hi"         
      read(21,pos=ifd+3+(i-1)*12) atag,dtty,ctfd,dtfd
      if(atag.lt.0) then
         tag=int(256*256+atag)
         else
         tag=int(atag)
         endif

      if(tag.eq.256) then
         width=dtfd
      elseif(tag.eq.257) then
         length=dtfd

      elseif(tag.eq.258) then
         ctfd=ctfd
      if(ctfd.eq.1) then
         allocate(bitpersample(1))
         bitpersample(1)=dtfd
      else
      allocate(bitpersample(ctfd))
      do i2 =1,ctfd
         read(21,pos=dtfd+1+(i2-1)*2) bitpersample(i2)
!         print *,bitpersample(i2),'bps'
      enddo
      endif
      elseif(tag.eq.262) then
      photometric=dtfd
      elseif(tag.eq.273) then
         ctfd=ctfd
         strip_n=ctfd
      allocate(stripoffset(ctfd))
      !      print *,strip_n,ctfd
    if(ctfd.ne.1) then
      do i2 =1,ctfd
         read(21,pos=dtfd+1+(i2-1)*4) stripoffset(i2)
!         print *,stripoffset(i2)
      enddo
   else
      stripoffset(ctfd)=dtfd
   endif
   
      elseif(tag.eq.279) then
         
         allocate(stripbytecount(ctfd))
    if(ctfd.ne.1) then        
      do i2 =1,ctfd
      read(21,pos=dtfd+1+(i2-1)*4) stripbytecount(i2)
      enddo
   else
      stripbytecount(ctfd)=dtfd
   endif
!      print *,stripbytecount(1)
!      print *,ctfd
      
      elseif(tag.eq.284) then
         planarconfig=dtfd

      elseif(tag.eq.50713) then
         read(21,pos=dtfd+1) dummy
         read(21,pos=dtfd+1) dummy
      
      elseif(tag.eq.339) then
         ctfd=ctfd
         dataformat_true=1
         if(ctfd.eq.1) then
            allocate(dataformat(1))
            dataformat(1)=dtfd
         else
      allocate(dataformat(ctfd))

      do i2 =1,ctfd
         read(21,pos=dtfd+1+(i2-1)*2) dataformat(i2)
         !         print *,dataformat(i2),'df'
      enddo
   endif
   
      endif

      enddo


      if(dataformat_true.ne.1) then
         allocate(dataformat(1))
         dataformat(1)=1
      endif

      print *,planarconfig,dataformat(1),"hi"
      
      if(planarconfig.eq.1) then     
      j2=0
      do i2 =1,strip_n
         imax=width
         jmax=stripbytecount(i2)/width/(bitpersample(1)/8)/band
         irec=0

         if(dataformat(1).eq.3) then
            do j =1, jmax
               j2=j2+1
               do i =1, imax
                  do k =1, band            
                     irec=1+irec   
                     read(21,pos=stripoffset(i2)+1+(irec-1)*bitpersample(1)/8) dummy_dble
                     data(i,j2,k)=real(dummy_dble)
!                                           print *,data(i,j2,k),stripoffset(i2)+1+(irec-1)*bitpersample(i2)/8
!                        if((i.eq.1404).and.(j2.eq.1457)) then
!                        print *,data(i,j2,k),dummy_dble
!                     endif                     
                  enddo

               enddo
            enddo

         else
!            print *,width,length,imax,jmax,stripoffset(1),stripbytecount(1)
!          print *,imax,jmax,stripoffset(i2)
            do j =1, jmax
               j2=j2+1
               do i =1, imax
                  do k =1, band            
                     irec=1+irec
!                     print *,i,j,k
                     print *,stripoffset(i2)+1+(irec-1)*bitpersample(1)/8
                     read(21,pos=stripoffset(i2)+1+(irec-1)*bitpersample(1)/8) dummy
!                           print *,dummy,i
                     if(dummy.lt.0) then
                        data(i,j2,k)=real(256+dummy)
                     else
                        data(i,j2,k)=real(dummy)
                     endif
!                     print *,data(i,j2,k),i,j2,j
                  enddo
                  
               enddo
            enddo

         endif


      enddo

   elseif(planarconfig.eq.2) then
      
      k=0
      jt=1
        j2=0       
      do i2 =1,strip_n
         imax=width
         jmax=stripbytecount(i2)/width/(bitpersample(1)/8)
         j2=j2+1
         if(j2.gt.length) then
            j2=1
            jt=1+jt
         endif
        k=jt
         irec=0
         if(dataformat(1).eq.3) then
            do j =1, jmax
               do i =1, imax
                     irec=1+irec   
                     read(21,pos=stripoffset(i2)+1+(irec-1)*bitpersample(1)/8) dummy_dble
                !     print *,dummy_dble
                     data(i,j2,k)=real(dummy_dble)                       
               enddo
            enddo

         else
            do j =1, jmax
               j2=j2+1
               do i =1, imax
                  do k =1, band            
                     irec=1+irec
                     read(21,pos=stripoffset(i2)+1+(irec-1)*bitpersample(i2)/8) dummy
                     !print *,dummy
                     if(dummy.lt.0) then
                        data(i,j2,k)=real(256+dummy)
                     else
                        data(i,j2,k)=real(dummy)
                     endif
                  enddo

               enddo
            enddo

         endif
      enddo
   endif
   

!!$   if(bitpersample(1).eq.8) then
!!$      data(:,:,:)=real((int(data(:,:,:))+1)/256-1)
!!$   endif
      
      
    endsubroutine open_tif

     subroutine write_tif(infile,data,width,length,band,strip_n)
      implicit none
      integer :: imax,jmax,jmin
      integer(4) :: width,length,band
!      real, allocatable,intent(out) :: data(:,:)
      real :: data(width,length,band)
      integer(4) :: strip_n 
      integer :: iglat,iglon,irec_entry_n,ifd
      integer(4) :: stripoffset(strip_n),stripbytecount(strip_n)      
      integer(2),parameter :: entry_n=16
      integer(1),parameter :: byte_order1=73            !!little 73 big 77
      integer(4),parameter :: compress=1
      integer(1),parameter :: version1=42               !little 42 big 0
      integer(1),parameter :: version2=0                !little 0 big 42
      integer(2),parameter :: dummy_bitpersample=32
      integer(2),parameter :: dummy_dataformart=3
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
      integer(4) :: i, j, t,i2,j2,k,ii,jj,irec,irec2,irec3,irec4
      character(512) :: infile
      character(30) :: str0G30
      integer(2) :: test(width,length,band)
      integer(4) :: bytenumber
!      test(:,:,1)=30
!      test(:,:,2)=60
!      test(:,:,3)=90
      
!      print *,'hi'

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
      
      print *,bytenumber
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
        stripoffset(i)=int(dble(i-1)*dble(width)*dble(length)*dble(band)*bytenumber/dble(strip_n))+8
        stripbytecount(i)=int(dble(width)*dble(length)*band*bytenumber/strip_n) 
!        print *,stripoffset(i), i,width, length, band,strip_n
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
         do j =jmin, jmax              
            do i =1, imax
          do k =1, band          
                  two_int=data(i,j,k) 
                  write(21) data(i,j,k)
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
   enddo
endif

      write(21) entry_n

      do i =1, int(entry_n)
         write(21) tag(i),dtty(i),ctfd(i),dtfd(i)
      enddo


      print *,'File size is' , 12*entry_n+2+strip_n*2*4+width*length*band*2+8+band*2
    endsubroutine write_tif


    
     subroutine read_size_tif(infile,width,length,band_num)
      implicit none
      integer :: imax,jmax
      integer :: iglat,iglon,irec_entry_n
      integer(4) :: ifd
!      integer(4) :: ifd
      integer(2) :: entry_n
      integer(2) :: dtty,atag
      integer(4) :: ctfd,dtfd,tag
      integer(1) :: dummy
      integer(4) :: width,length
      integer(4) :: photometric,band_num
      integer(4),allocatable :: stripoffset(:),stripbytecount(:)
      integer(4) :: planarconfig
      integer :: NX,NY,NX2,NY2,start1,start2,finish1,finish2,iostat
      integer :: i, j, t,i2,j2,k,ii,jj,irec,irec2,irec3,irec4
      character(512) :: infile
      INTEGER(1) :: kai(2)
      character(30) :: str0G30          


      
      open(21,file=infile,  &
           status='old',form='unformatted',access='stream')     
      read(21,pos=1) dummy
      read(21,pos=5) ifd
      print *,dummy
      print *,ifd
      read(21,pos=ifd+1) entry_n
 
      do i =1, int(entry_n)
         read(21,pos=ifd+3+(i-1)*12) atag,dtty,ctfd,dtfd
         print *,atag,dtty,ctfd,dtfd
      if(atag.lt.0) then
         tag=int(256*256+atag)
         else
         tag=int(atag)
         endif
 
      if(tag.eq.256) then
      width=dtfd
      elseif(tag.eq.257) then
         length=dtfd
      elseif(tag.eq.258) then
      band_num=ctfd         
      elseif(tag.eq.262) then
      photometric=dtfd
      elseif(tag.eq.273) then
         ctfd=ctfd
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
      endif
 
      enddo


    endsubroutine read_size_tif
