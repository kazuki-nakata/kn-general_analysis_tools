PROGRAM pstolalon
  USE interpolation
  IMPLICIT NONE
  integer :: mm,yy,dd,id,im,iy,ied
  character :: ms*2,ys*4,ds*2,wname*66,lon_file*66,lat_file*66,are_file*66,rname*66

  lat_file='../../pss06lats_v3.dat'
  lon_file='../../pss06lons_v3.dat'
  are_file='../../pss06area_v3.dat'
!--------内挿データの緯度経度設定(現在メルカトルのみ、随時更新)------------
  CALL set_interpolation_point_ll
! CALL set_interpolation_point_POLA
! CALL set_interpolation_point_EASE
!--------元データの緯度経度の読み込み------------------------
  CALL read_original_data_point(lat_file,lon_file,are_file)
! CALL read_original_data_point_1D
! CALL read_original_data_point_2D
!------内挿時の元データのポインター設定------------------------
  CALL set_pointer_on_llproj
! CALL set_pointer_on_POLA
! CALL set_pointer_on_EASE
  
  do iy=2003,2010
!2003,2011
     do im=2,3
        if ((im==1).or.(im==3).or.(im==5).or.(im==7).or.(im==8) &
             .or.(im==10).or.(im==12)) then 
           ied=31
        else if (im==2) then 
           if ((iy==1980).or.(iy==1984).or.(iy==1988).or.(iy==1992) &
                .or.(iy==1996).or.(iy==2000).or.(iy==2004).or.(iy==2008)) then 
              ied=29
           else 
              ied=28
           end if
        else 
           ied=30
        end if
        
        do id=1,ied
           yy=iy
           write(ys,403) yy
           mm=im
           if(mm.ge.10) write(ms,401) mm
           if(mm.lt.10) write(ms,402) 0,mm
           dd=id
           if(dd.ge.10) write(ds,401) dd
           if(dd.lt.10) write(ds,402) 0,dd
401        format(i2)
402        format(i1,i1)
403        format(i4)
           
           write(rname,501) ys,ms,ds
501        format(12H../version2/,a4,a2,a2,5H.data)
  !         rname='original_data/2003-2010_ave.icepro.data'
    open(42,file=rname,form='unformatted',access='direct', &
         recl=4*num_origin,status='old',err=100)
    close(42)
                  call read_original_data_binary(rname)

!           write(rname,502) ms
!502        format(5Hmask/,a2,10H.mask.data)
!           rname='mask/2003-2010.mask.data'
!                  call read_original_mask_binary(rname)

                  call interpolation_gaussian
!call box_averaging
!call interpolation_bilinear

           write(wname,503) ys,ms,ds
503        format(2H./,a4,a2,a2,5H.data)
!           wname='lldata/2003-2010_ave.ll.icepro.data'
           print *,wname
                 call write_interpolation_data_binary(wname,2)

100 cycle
        enddo
     enddo
  enddo

end program















