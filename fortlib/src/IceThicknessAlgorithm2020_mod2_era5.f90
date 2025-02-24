PROGRAM AMSR_ice_thickness_algorithm
  !-----created in 2018/02/13 by K. Nakata. The details are written on the README file.
  use netcdf_tool  
  IMPLICIT NONE
  INTEGER(4) :: NX,NY,NX2,NY2,start1,start2,finish1,finish2
  INTEGER(4) :: i, j, t,i2,j2,k,ii,jj,irec,irec2,irec3,irec4,l
  INTEGER(4) :: imax,jmax,kmax,inum,inum1,inum2
  
  INTEGER(4), PARAMETER:: imax2=1264
  INTEGER(4), PARAMETER:: jmax2=1328
  INTEGER(4), PARAMETER:: iams=1264
  INTEGER(4), PARAMETER:: jams=1328
  INTEGER(4), PARAMETER :: lamin=-85,lamax=-40,lomin=0,lomax=360
  INTEGER(4), PARAMETER :: xmin=1,xmax=1264,ymin=1,ymax=1328
  INTEGER(4) :: lxmin,lxmax,lymin,lymax,lxmin2,lxmax2,lymin2,lymax2,alx,aly
  
  CHARACTER, ALLOCATABLE::infile1(:)*55,infile2(:)*56,metaf(:)*57,orbit(:)*1
  CHARACTER :: rname*70,rname2*70,rname3*70, rname4*80,wname*70, maskfile*70

  !L2A DATA
  REAL(4), ALLOCATABLE :: lon(:,:), Hp(:,:), Vp(:,:),lat(:,:)
  REAL(4), ALLOCATABLE :: Hp2(:,:), Vp2(:,:)
  REAL(4), ALLOCATABLE :: lx1(:,:), ly1(:,:)

  !INTERPOLATED POLARSTEREO DATA
  REAL(4) :: V(imax2,jmax2),V2(imax2,jmax2),H(imax2,jmax2),H2(imax2,jmax2)

  !ANCILARY DATA
  REAL(4) :: lat2(imax2,jmax2),lon2(imax2,jmax2),land(imax2,jmax2)
  REAL(4) :: sir(imax2,jmax2),sir0(imax2,jmax2),seg(imax2,jmax2)
  REAL(4) :: sum1(imax2,jmax2),mdist(imax2,jmax2),sum2(imax2,jmax2)

  !STACK
  REAL(4) :: ipr_a(imax2,jmax2,20), hf_a(imax2,jmax2,20),jday_a(imax2,jmax2,20),pol_a(imax2,jmax2,20)
  REAL(4) :: ipr_d(imax2,jmax2,20), hf_d(imax2,jmax2,20),jday_d(imax2,jmax2,20),pol_d(imax2,jmax2,20)

  !FOR OUTPUT
  INTEGER(4) :: onum
  REAL(4) :: ip2_a, hf2_a, act_a, ina_a, hf2_af_a, hf2_ts_a, ip2_af_a, ip2_ts_a
  REAL(4) :: ip2_d, hf2_d, act_d, ina_d, hf2_af_d, hf2_ts_d, ip2_af_d, ip2_ts_d
  REAL(4) :: ip2(imax2,jmax2),hf2(imax2,jmax2)
  REAL(4) :: act(imax2,jmax2),ina(imax2,jmax2)
  REAL(4) :: hf2_af(imax2,jmax2),hf2_ts(imax2,jmax2),ip2_af(imax2,jmax2), ip2_ts(imax2,jmax2)
  
  INTEGER(4) :: lx2(imax2,jmax2),ly2(imax2,jmax2)
  INTEGER(4) :: dnum_a(imax2,jmax2),dnum_d(imax2,jmax2)
 
  INTEGER(4) :: yy,mm,mm2,dd,iy,im,im2,day,hour,minute,ied
  INTEGER(4) :: pflag
  INTEGER(4) :: jday5
  REAL(8) :: jday1,jday2,jday3,jday4,atime  
  CHARACTER :: ys*4,ms*2,ds*2,ms2*2

  INTEGER(4) :: iglat,iglon
  REAL(4) :: lonn1,lonn2,latn,theta,ati,anum,latl,lonl,alon,alat
  
  REAL(8) :: pi,re,dora,alpha,dist,dist_x,dist_y,adis,e,dlon,fr,sum  
  REAL(4) :: PR,GRV,TH,pol2,pol3,tii
  REAL(4) :: a1,b1,c1,d1,e1,f1,w1
  REAL(4) :: asir
  
  namelist /param/ inum,iy,im,day

  INTEGER(4),PARAMETER :: Imax_era=1440
  INTEGER(4),PARAMETER :: Jmax_era=241
  INTEGER(4),PARAMETER :: Imax_era2=1440
  INTEGER(4),PARAMETER :: Jmax_era2=721
  INTEGER(4),PARAMETER :: Nvars=6
  CHARACTER(len=3) :: var_list(Nvars)
  REAL(4),PARAMETER :: Xmin_era=0.0, Ymin_era=-90.0
  REAL(4),PARAMETER :: Xint_era=0.25, Yint_era=0.25
  REAL(4),PARAMETER::tf=-1.86,tk=273.15
  REAL(4)::hii,hs,ali,alw
  REAL(4)::hb,sw,lw,se,la,tw
  REAL(4)::dhi,gipr
  
  REAL(4) :: dat1(Imax_era,Jmax_era,Nvars)
  REAL(4),ALLOCATABLE :: dat01(:,:,:)
  REAL(4) :: era_land(Imax_era,Jmax_era)
  
  REAL(4) :: rlon_era(Imax_era), rlat_era(Jmax_era)
  REAL(4) :: slp,w10m,t2m,d2m,gc,ric     !from ERA
  REAL(4) :: d2m1,d2m2,d2m3,d2m4,factor_t,dx,dy
  REAL(4)::adat1(Nvars),adat2(Nvars),adat3(Nvars),adat4(Nvars),adat5(Nvars)
  REAL(4) :: bilin,dt,fc,hice
  INTEGER(4) :: nv,ipos,jpos,itime,ipos2,jpos2,iarea
  REAL(4),PARAMETER :: Undef=100. !9.9e33
  REAL(4),PARAMETER :: Undefm=-999.
  REAL(4) :: sigma,scale_radius,window_radius
  REAL(4),PARAMETER :: ic_thres=-0.5
  REAL(4) :: val

  var_list=(/'u10','v10','d2m','t2m','msl','tcc'/)
  pi=acos(-1.d0)
  re=6378.14d0 !赤道半径(km) 
  dora=pi/180.d0
  !-----------------------------------------------------------

!------------------------------------------------------
      open(21,file='latlon/pss06lats_v3.dat',  &
      status='old',form='unformatted',access='direct',recl=4)
       irec=0
       do 101 j=1,jams
          do 101 i=1,iams
       irec=1+irec     
       if((i.ge.xmin).and.(i.le.xmax)) then
       if((j.ge.ymin).and.(j.le.ymax)) then
       read(21,rec=irec) iglat
       lat2(i-xmin+1,j-ymin+1)=iglat/1.e5
       endif
       endif
 101   continue
!      close(57)

       open(22,file='latlon/pss06lons_v3.dat', &
      status='old',form='unformatted',access='direct',recl=4)
       irec=0
       do 102 j=1,jams
          do 102 i =1,iams
             irec=1+irec
      if((i.ge.xmin).and.(i.le.xmax)) then
       if((j.ge.ymin).and.(j.le.ymax)) then
       read(22,rec=irec) iglon
       lon2(i-xmin+1,j-ymin+1)=iglon/1.e5
           if (lon2(i-xmin+1,j-ymin+1).le.0.) then
           lon2(i-xmin+1,j-ymin+1)=360.+lon2(i-xmin+1,j-ymin+1)
           endif
       endif
       endif
!       print *,lat2(i,j),lon2(i,j)
 102   continue

       close(22)
       close(21)

!-----------------------------------------------------------------
      read(*,param)


        if ((im.eq.12).or.(im.eq.1).or.(im.eq.3).or.(im.eq.5) &
             .or.(im.eq.7).or.(im.eq.8).or.(im.eq.10)) then 
           ied=31
        else if (im.eq.2) then
           if ((iy.eq.2000).or.(iy.eq.2004).or.(iy.eq.2008) &
              .or.(iy.eq.2012).or.(iy.eq.2016).or.(iy.eq.2020)) then 
              ied=29
           else
              ied=28
           end if
        else 
           ied=30
        end if

        yy=iy
        write(ys,403) yy
        mm=im
        if(mm.ge.10) write(ms,401) mm
        if(mm.lt.10) write(ms,402) 0,mm
        dd=day
        if(dd.ge.10) write(ds,401) dd
        if(dd.lt.10) write(ds,402) 0,dd
        if(day.eq.ied) then
        mm2=im+1
        if(mm2.ge.10) write(ms2,401) mm2
        if(mm2.lt.10) write(ms2,402) 0,mm2
        endif
        
 401      format(i2)
 402      format(i1,i1)
 403      format(i4)


!!$       write(rname,502) ys,ms
!!$ 502    format(26H../iceproduction2/iprdata/, &
!!$        a4,a2,22H.icepro.antarctic.data)
!!$        open(23,file=rname,form='unformatted',access='direct',recl=4)
!!$        irec=3*iams*jams*(day-1)+2*iams*jams
!!$       do 103 j=1,jams
!!$          do 103 i=1,iams
!!$       irec=1+irec     
!!$       if((i.ge.xmin).and.(i.le.xmax)) then
!!$       if((j.ge.ymin).and.(j.le.ymax)) then
!!$       read(23,rec=irec) land(i-xmin+1,j-ymin+1)
!!$       endif
!!$       endif
!!$ 103   continue

       close(23)

       write(rname4,505) ys,ms,ds
 505    format(32H./iceconcentration/segmentation/, &
             a4,a2,a2,5H.data)
       open(23,file=rname4,form='unformatted',access='direct',recl=4*iams*jams, status='old')
       read(23,rec=1) sir0
       read(23,rec=2) seg
       sir(:,:)=0
       land(:,:)=0
       where((seg.eq.15).or.(seg.eq.30).or.(seg.eq.80).or.(seg.eq.100)) sir=100
       where(seg.eq.120) land=120
       
       close(23)
       
!---------------------------------------------------------
  rlon_era(1:Imax_era) = (/ (Xmin_era+Xint_era*(ii-1),ii=1,Imax_era) /)
  where(rlon_era > 360.)
    rlon_era=rlon_era-360.
  endwhere
  rlat_era(1:Jmax_era) = (/ (Ymin_era+Yint_era*(jj-1),jj=1,Jmax_era) /)
  era_land(:,:)=0
!  print *,rlat_era
      write(rname2,504) ys,ms
504   format(12H./ERA5/era5_,a4,a2,3H.nc)
      
  open(24,file=rname2,form='unformatted',access='direct', &
        recl=4*Imax_era2*Jmax_era2,status='old',action='read')
  call read_daily_netcdf(rname2,var_list,Nvars,dat01,day)
  do j=Jmax_era2-Jmax_era+1,Jmax_era2
     do i=1,Imax_era2
        dat1(i,Jmax_era2-j+1,:)=dat01(i,j,:)
     enddo
  enddo
!----------------------------------------------------------------
      Allocate(infile1(1:inum))
      Allocate(infile2(1:inum))
      Allocate(metaf(1:inum))
      Allocate(orbit(1:inum))
  open(12,file='file1.txt',status='old',action='read')
  open(13,file='file2.txt',status='old',action='read')
  open(14,file='file3.txt',status='old',action='read')
  
       write(wname,503) ys,ms,ds
503    format(11H./version1/,a4,a2,a2,5H.data)  
        print *, wname
        open(15,file=wname,form='unformatted',access='direct',recl=4,status='replace')

  do i =1,inum
  read(12,'(a55)') infile1(i)
  read(13,'(a57)') metaf(i)
  read(14,'(a1)') orbit(i)
enddo

            ipr_a(:,:,:)=0.
            hf_a(:,:,:)=0.
            pol_a(:,:,:)=0.
            ipr_d(:,:,:)=0.
            hf_d(:,:,:)=0.
            pol_d(:,:,:)=0.
            
            dnum_a(:,:)=0
            dnum_d(:,:)=0
!-------------------------------------------------------------------
!------------------------------------------------------------------


   do t= 1,inum
!      print *,metaf(t)
  print *,'-----status of the progress-----',t,'/',inum
  open(16,file=metaf(t),status='old',action='read')

  read(16,*) NX, NY
  read(16,*) jday1, jday2

  if(t.eq.1) then
     jday5=int(max(jday1,jday2))
  endif
  
     imax=NX
     jmax=NY

  ALLOCATE(Vp(1:imax,1:jmax))
  ALLOCATE(Hp(1:imax,1:jmax))
  ALLOCATE(Vp2(1:imax,1:jmax))
  ALLOCATE(Hp2(1:imax,1:jmax))
  ALLOCATE(lon(1:imax,1:jmax))
  ALLOCATE(lat(1:imax,1:jmax))
  ALLOCATE(lx1(1:imax,1:jmax))
  ALLOCATE(ly1(1:imax,1:jmax))

        open(17,file=infile1(t),form='unformatted',access='direct', &
     recl=4,status='old')


        print *,'       --data reading---'
!------------------------------------------------------------------
            irec=0
            lxmax=0
            lymax=0
            lxmin=2000
            lymin=2000
          lx2(:,:)=999
          ly2(:,:)=999

        do 87 j=1,jmax
        do 87 i=1,imax
           irec=1+irec  
           irec2=irec+imax*jmax
           irec3=irec+imax*jmax*2
           irec4=irec+imax*jmax*3
          read(17,rec=irec) lon(i,j)
          read(17,rec=irec2) lat(i,j)
          if(lon(i,j).lt.0) then
             lon(i,j)=lon(i,j)+360.
          endif
          CALL MAPLL(lx1(i,j),ly1(i,j),lat(i,j),lon(i,j))         
          read(17,rec=irec3) val
          Vp(i,j)=val-(-0.01019*val + 5.49799)
          read(17,rec=irec4) val
          Hp(i,j)=val-(-0.00561*val + 4.19181)
          irec4=irec+imax*jmax*4         
          read(17,rec=irec4) val
          Vp2(i,j)=val-(-0.01403*val + 5.32379)
          irec4=irec+imax*jmax*5
          read(17,rec=irec4) val
          Hp2(i,j)=val-(-0.00980*val + 3.75174)
!          print *, lon(i,j),lat(i,j),Vp(i,j),Hp(i,j),Vp2(i,j),Hp2(i,j)
          if((lat(i,j).le.lamin).or.(lat(i,j).ge.lamax)) goto 87
          if((lon(i,j).le.lomin).or.(lon(i,j).ge.lomax)) goto 87          
          if((nint(lx1(i,j)).lt.1).or.(nint(lx1(i,j)).gt.imax2)) goto 87
          if((nint(ly1(i,j)).lt.1).or.(nint(ly1(i,j)).gt.jmax2)) goto 87
          
!          print *,Vp(i,j)
          lx2(nint(lx1(i,j)),nint(ly1(i,j)))=i
          ly2(nint(lx1(i,j)),nint(ly1(i,j)))=j

          lxmax=max(nint(lx1(i,j)),lxmax)
          lymax=max(nint(ly1(i,j)),lymax)
          lxmin=min(nint(lx1(i,j)),lxmin)
          lymin=min(nint(ly1(i,j)),lymin)        
87        continue

        print *,'       --interpolation---'
!-------------------------------------------------------------------
       do 345 jj=lymin,lymax
       do 345 ii=lxmin,lxmax
          sum=0.
          V(ii,jj)=0.
          H(ii,jj)=0.
          V2(ii,jj)=0.
          H2(ii,jj)=0.
          if((sir(ii,jj).eq.0).or.(land(ii,jj).eq.120.)) goto 345    
          if((lx2(ii,jj).ne.999.).and.(ly2(ii,jj).ne.999)) then
             alx=lx2(ii,jj)
             aly=ly2(ii,jj)
          else
             alx=999
             aly=999
             adis=999

             do j =1,5
                do i =1,5
                   
                   if((ii+i-3.ge.lxmin).and.(ii+i-3.lt.lxmax)) then
                      if((jj+j-3.ge.lymin).and.(jj+j-3.lt.lymax)) then

                         if(lx2(ii+i-3,jj+j-3).ne.999) then
                            dist_x = real(ii)-lx1(lx2(ii+i-3,jj+j-3),ly2(ii+i-3,jj+j-3))
                            dist_y = real(jj)-ly1(lx2(ii+i-3,jj+j-3),ly2(ii+i-3,jj+j-3))
                            dist=sqrt(dist_x**2+dist_y**2)
                            if(dist.le.adis) then
                               adis=dist
                               alx=lx2(ii+i-3,jj+j-3)
                               aly=ly2(ii+i-3,jj+j-3)
                            endif
                         endif

                      endif
                   endif
                   
                enddo
             enddo
             
          endif


          if((alx.eq.999).or.(aly.eq.999)) goto 345

          lxmin2=alx-10
          lxmax2=alx+10
          lymin2=aly-10
          lymax2=aly+10
          if(lxmin2.le.0) lxmin2=1
          if(lymin2.le.0) lymin2=1
          if(lxmax2.ge.imax) lxmax2=imax
          if(lymax2.ge.jmax) lymax2=jmax

          scale_radius=10.0/2.
          sigma=scale_radius/1.177410
          window_radius=(sigma*6)/2.
          do 346 j =lymin2,lymax2
             do 346 i =lxmin2,lxmax2
                !-------
                dist=sqrt((real(ii)-lx1(i,j))**2+(real(jj)-ly1(i,j))**2)*6.25
                if((Vp(i,j).eq.0.).or.(Vp2(i,j).eq.0.)) goto 346
                if (dist.ge.window_radius) goto 346
                fr=exp(-dist**2./(sigma**2))
                sum=sum+fr
                V(ii,jj)=fr*Vp(i,j)+V(ii,jj)
                H(ii,jj)=fr*Hp(i,j)+H(ii,jj)
                V2(ii,jj)=fr*Vp2(i,j)+V2(ii,jj)
                H2(ii,jj)=fr*Hp2(i,j)+H2(ii,jj)
346             continue

                if(sum.gt.0.) then
                   V(ii,jj)=V(ii,jj)/sum
                   H(ii,jj)=H(ii,jj)/sum
                   V2(ii,jj)=V2(ii,jj)/sum
                   H2(ii,jj)=H2(ii,jj)/sum
                endif
345             continue



        print *,'       --polynya detection---'
        jday4=jday1+(jday2-jday1)
!        jday5=int(jday4)
        atime=jday4-jday5
        
        print *,jday4,jday5,atime
        
        do 88 j=lymin,lymax
        do 88 i=lxmin,lxmax
!           print *,i,j,lxmax,lymax
           if(sir(i,j).eq.0) then
           goto 88
           endif
           if(land(i,j).eq.120.) then
           goto 88
           endif

           if((V(i,j).eq.0.).or.(V2(i,j).eq.0.)) goto 88


           PR=(V(i,j)-H(i,j))/(V(i,j)+H(i,j))
           GRV=(V2(i,j)-V(i,j))/(V2(i,j)+V(i,j))
!           print *,PR,GRV
!           print *, PR,GRV,V2(i,j),V(i,j)
!           print *,PR,GRV,V2(i,j)
!--------------discrimination and estimation of sea ice thickness------- 

            w1=-192.86*PR+1002.23*GRV-0.700

            if((w1.ge.0).and.(PR.ge.0.05)) then
!active
            pflag=1
!            TH=exp(1/(PR*74))-1.06
            TH=exp(1/(PR*596-11.8))-1.008
!            TH=exp(1/(PR*72))-1.06
            if(TH.lt.0.001) then
               TH=0.001
            endif

            else

            TH=exp(1/(PR*72))-1.06
               
            if(TH.lt.0.001) then
               TH=0.001
            endif
            
            if(TH.gt.0.2) then
               TH=100.
               pflag=3
               else
             pflag=2
            endif

            endif

       ipos=-999
        do ii = 2, Imax_era
          if((lon2(i,j) > rlon_era(ii-1)) .and. &
                (lon2(i,j) <= rlon_era(ii))) then
            ipos=ii
            ipos2=ii-1
          endif
        enddo


        jpos=-999
        do jj = 2, Jmax_era
          if((lat2(i,j) > rlat_era(jj-1)) .and. &
                (lat2(i,j) <= rlat_era(jj))) then
            jpos=jj
            jpos2=jj-1
          endif
        enddo
 !       print *,jpos,ipos

        if(ipos == -999) then
           ipos2=1
           ipos=Imax_era
        endif

        if(jpos == -999) then
           cycle
        endif

        if(ipos == -999) then
        dx=(lon2(i,j)-rlon_era(ipos2)) &
             & / ((rlon_era(ipos)+360.)-rlon_era(ipos2))
        dy=(lat2(i,j)-rlat_era(jpos2)) &
             & / (rlat_era(jpos)-rlat_era(jpos2))
           else
        dx=(lon2(i,j)-rlon_era(ipos2)) &
             & / (rlon_era(ipos)-rlon_era(ipos2))
        dy=(lat2(i,j)-rlat_era(jpos2)) &
             & / (rlat_era(jpos)-rlat_era(jpos2))
        endif

        do nv = 1, Nvars
          adat1(nv) = bilin(ipos,jpos,dx,dy,Undef, &
                           dat1(ipos2,jpos2,nv), &
                           dat1(ipos  ,jpos2,nv), &
                           dat1(ipos  ,jpos  ,nv), &
                           dat1(ipos2,jpos  ,nv), &
                           era_land(ipos2,jpos2), &
                           era_land(ipos  ,jpos2), &
                           era_land(ipos  ,jpos  ), &
                           era_land(ipos2,jpos ))
!          print *, dat1(ipos2,jpos2,nv), dat1(ipos,jpos,nv),ipos,jpos,nv
        enddo

!------------
!  /'u10','v10','d2m','t2m','msl','tcc'/
!------------
        slp = adat1(5)
        w10m = sqrt((adat1(1))**2.0+(adat1(2))**2.0)
        t2m = adat1(4)
        d2m = adat1(3)
        gc =adat1(6)
        if(gc.le.0.) gc=0.
        if(gc.gt.1) gc=1.
!        print *,slp,w10m,t2m,d2m,gc
!!$        slp=100000
!!$        w10m=10
!!$        t2m=273.15-20
!!$        d2m=273.15-20
!!$        gc=1.0
!--------------------------------------------------------------           
           !*******
           tw=tf+tk
           !*******

           if(TH <= 0.1) then !new ice
              ric=100./100.
              hii=TH
              hs=0.
              ali=0.27; alw=0.06
!              print *,slp,t2m,d2m,w10m,ric,gc,tw,lat2(i,j),jday5
                 call hb_all(slp,t2m,d2m,w10m,ric,gc,tw,lat2(i,j), &
                             jday5,hii,hs,ali,alw,hb,sw,lw,se,la,tii)
 !             print *,'h2'
              if (hb < 0) then
                 call ice_production(hb,dhi)
              else
                 dhi =0.
              end if
!              print *,'hi3'
              gipr=dhi

           else if ((TH > 0.1).and.(TH <= 0.2)) then !young ice

              ric=100./100.
              hii=TH
              hs=0.
              ali=0.36; alw=0.06   

                call hb_all(slp,t2m,d2m,w10m,ric,gc,tw,lat2(i,j), &
                             jday5,hii,hs,ali,alw,hb,sw,lw,se,la,tii)

              if (hb < 0) then
                 call ice_production(hb,dhi)
              else
                 dhi =0.
              end if

              gipr=dhi
 
            else if (TH == 100.) then !First-year ice
               gipr=0.
               tii=t2m
               hb=0.
            end if
            !            print *,orbit(t)
            !
!            print *, dnum_d(352,1242),i,j,lymin,lymax,lxmin,lxmax
!            print *,orbit(t)
            if(orbit(t).eq."A") then
               
               dnum_a(i,j)=dnum_a(i,j)+1
               onum=dnum_a(i,j)            
               ipr_a(i,j,onum)=gipr
               hf_a(i,j,onum)=hb
               if(pflag.eq.1) then
                  pol_a(i,j,onum)=1
               elseif(pflag.eq.2) then
                  pol_a(i,j,onum)=10
               else
                  pol_a(i,j,onum)=100
               endif

               
               
            elseif(orbit(t).eq."D") then
               
               dnum_d(i,j) = dnum_d(i,j) + 1
!               if((t.eq.19).and.(i.eq.352).and.(j.eq.1242)) print *, dnum_d(352,1242),orbit(t)
               onum = dnum_d(i,j)
!               print *,onum
               ipr_d(i,j,onum)=gipr
               hf_d(i,j,onum)=hb
               if(pflag.eq.1) then
                  pol_d(i,j,onum)=1
               elseif(pflag.eq.2) then
                  pol_d(i,j,onum)=10
               else
                  pol_d(i,j,onum)=100
               endif
               
            endif
         
88              continue

            close(17)
            close(18)
           DEALLOCATE(Vp,Hp,Vp2,Hp2,lon,lat,lx1,ly1)
          enddo

!-----------------daily data ------------------------------------------t

       print *,'       --daily mean---'

       act(:,:)=0.
       ina(:,:)=0.
       ip2_af(:,:)=0.
       ip2_ts(:,:)=0.
       hf2_af(:,:)=0.
       hf2_ts(:,:)=0.
       hf2(:,:)=0.
       ip2(:,:)=0.          
       irec=0

       do j =1,jmax2
          do i =1,imax2

             if(land(i,j).eq.120.) then
                ip2(i,j)=9.9E33
                ip2_af(i,j)=9.9E33
                ip2_ts(i,j)=9.9E33                        
                act(i,j)=9.9E33
                ina(i,j)=9.9E33
                hf2_af(i,j)=9.9E33
                hf2_ts(i,j)=9.9E33
                hf2(i,j)=9.9E33                         
             elseif(sir(i,j).eq.0) then
                ip2(i,j)=0
                ip2_af(i,j)=0
                ip2_ts(i,j)=0                     
                act(i,j)=2
                ina(i,j)=2
                hf2_af(i,j)=0
                hf2_ts(i,j)=0
                hf2(i,j)=0         
             else
!                print *,'hi'
!                print *,i,j,dnum_a(i,j),dnum_d(i,j)
                if(dnum_a(i,j)+dnum_d(i,j).ge.1) then

                   ip2_a=0; hf2_a=0; act_a=0; ina_a=0; hf2_af_a=0; hf2_ts_a=0; ip2_af_a=0; ip2_ts_a=0
                   ip2_d=0; hf2_d=0; act_d=0; ina_d=0; hf2_af_d=0; hf2_ts_d=0; ip2_af_d=0; ip2_ts_d=0              
                   if(dnum_a(i,j).ge.1) then
 !                     print *,'a exist'
                      onum=dnum_a(i,j)!+dnum_d(i,j)
                      do t = 1,onum
                         hf2_a = hf2_a + hf_a(i,j,t)/real(onum)
                         ip2_a = ip2_a + ipr_a(i,j,t)/real(onum)
                         if(pol_a(i,j,t).eq.1) then
                            act_a = act_a + 1./real(onum)
                            hf2_af_a = hf2_af_a + hf_a(i,j,t)/real(onum)
                            ip2_af_a = ip2_af_a + ipr_a(i,j,t)/real(onum)                           
                         elseif(pol_a(i,j,t).eq.10) then
                            hf2_ts_a = hf2_ts_a  +  hf_a(i,j,t)/real(onum)
                            ip2_ts_a  = ip2_ts_a  +  ipr_a(i,j,t)/real(onum)
                            ina_a = ina_a + 1./real(onum)                         
                         endif
                      enddo
                      
                   endif

                   if(dnum_d(i,j).ge.1) then
!                      print *,'d exist'                      
                      onum=dnum_d(i,j) !dnum_a(i,j)+dnum_d(i,j)
                      do t = 1,onum
                         hf2_d = hf2_d + hf_d(i,j,t)/real(onum)
                         ip2_d = ip2_d + ipr_d(i,j,t)/real(onum)
                         if(pol_d(i,j,t).eq.1) then
                            act_d = act_d + 1./real(onum)
                            hf2_af_d = hf2_af_d + hf_d(i,j,t)/real(onum)
                            ip2_af_d = ip2_af_d + ipr_d(i,j,t)/real(onum)                           
                         elseif(pol_d(i,j,t).eq.10) then
                            hf2_ts_d = hf2_ts_d  +  hf_d(i,j,t)/real(onum)
                            ip2_ts_d  = ip2_ts_d  +  ipr_d(i,j,t)/real(onum)
                            ina_d = ina_d + 1./real(onum)                         
                         endif
                      enddo
                   endif
                   
                   act(i,j)=(act_a+act_d)/2
                   ina(i,j)=(ina_a+ina_d)/2
                   ip2_af(i,j)=(ip2_af_a+ip2_af_d)/2
                   ip2_ts(i,j)=(ip2_ts_a+ip2_ts_d)/2
                   hf2_af(i,j)=(hf2_af_a+hf2_af_d)/2
                   hf2_ts(i,j)=(hf2_ts_a+hf2_ts_d)/2
                   hf2(i,j)=(hf2_a+hf2_d)/2
                   ip2(i,j)=(ip2_a+ip2_d)/2

!                   if(dnum_a(i,j)+dnum_d(i,j).ge.7) then
!                   print *,act_a,act_d,act(i,j),dnum_a(i,j),dnum_d(i,j)
!                   endif
!                   print *,ip2(i,j),hf2(i,j),(ip2_af(i,j)+ip2_ts(i,j)),(hf2_af(i,j)+hf2_ts(i,j))
!!$                   act(i,j)=act_a+act_d
!!$                   ina(i,j)=ina_a+ina_d
!!$                   ip2_af(i,j)=ip2_af_a+ip2_af_d
!!$                   ip2_ts(i,j)=ip2_ts_a+ip2_ts_d
!!$                   hf2_af(i,j)=hf2_af_a+hf2_af_d
!!$                   hf2_ts(i,j)=hf2_ts_a+hf2_ts_d
!!$                   hf2(i,j)=hf2_a+hf2_d
!!$                   ip2(i,j)=ip2_a+ip2_d 


                else

                   ip2(i,j)=9.9E33
                   hf2(i,j)=9.9E33
                   act(i,j)=9.9E33
                   ina(i,j)=9.9E33
                   ip2_af(i,j)=9.9E33
                   ip2_ts(i,j)=9.9E33
                   hf2_af(i,j)=9.9E33
                   hf2_ts(i,j)=9.9E33                       
                endif

             endif

          enddo
       enddo

       do j =1,jmax2
          do i =1,imax2
             irec=1+irec
             write(15,rec=irec) act(i,j)
          enddo
       enddo

       do j =1,jmax2
          do i =1,imax2
             irec=1+irec
             write(15,rec=irec) ina(i,j)
          enddo
       enddo


       do j =1,jmax2
          do i =1,imax2
             irec=1+irec
             write(15,rec=irec) hf2_af(i,j)
          enddo
       enddo

       do j =1,jmax2
          do i =1,imax2
             irec=1+irec
             write(15,rec=irec) hf2_ts(i,j)
          enddo
       enddo

       do j =1,jmax2
          do i =1,imax2
             irec=1+irec
             write(15,rec=irec) ip2_af(i,j)
          enddo
       enddo

       do j =1,jmax2
          do i =1,imax2
             irec=1+irec
             write(15,rec=irec) ip2_ts(i,j)
          enddo
       enddo

       do j =1,jmax2
          do i =1,imax2
             irec=1+irec
             write(15,rec=irec) hf2(i,j)
          enddo
       enddo

       do j =1,jmax2
          do i =1,imax2
             irec=1+irec
             write(15,rec=irec) ip2(i,j)
          enddo
       enddo
       do j =1,jmax2
          do i =1,imax2
             irec=1+irec
             write(15,rec=irec) sir0(i,j)
          enddo
       enddo
       do j =1,jmax2
          do i =1,imax2
             irec=1+irec
             write(15,rec=irec) seg(i,j)
          enddo
       enddo
       print *,'-----finish process completely-----'
!-------------------------------------------------------------
end program
 
subroutine ice_production(hb,dhi)
  implicit none
!  REAL(4),PARAMETER::lf=3.02e5
!  REAL(4),PARAMETER::rhoi=900.
!  REAL(4),PARAMETER::lf=2.34e5
  REAL(4),PARAMETER::lf=3.34e5
  REAL(4),PARAMETER::rhoi=920.
  REAL(4),PARAMETER::dt=24.*60.*60.
  REAL(4)::hb,dhi
  dhi=abs(hb)*dt/rhoi/lf
end subroutine ice_production

subroutine hb_all(slp,t2m,td2m,w10m,ric,c,tw,lat,jday,hi,hs,ali,alw,hb,sw,lw,se,la,ti)
  implicit none
  INTEGER(4)::jday
  REAL(4),PARAMETER::tk=273.15,tf=-1.86,cki=2.03,cks=0.31
  REAL(4)::slp,t2m,td2m,w10m,ric,c,tw,ti,lat
  REAL(4)::hi,hs,ali,alw
  REAL(4)::ckb
  REAL(4)::res0,res,hres,dt
  REAL(4)::fc,qi,qw,ilw,olwi,olww,sehi,sehw,lahi,lahw
  REAL(4)::hbi,hbw,hb,sw,lw,se,la

  !--- first guess                    
  ti=t2m-10.
  if (ti >= tk) ti=tk-10. 

  ckb=cki*cks/(cki*hs+cks*hi)
  res0=1.; dt=1.               
  do 
     call swave(ali,alw,c,td2m,jday,lat,qi,qw)
     call lwave(c,ti,tw,t2m,ilw,olwi,olww)
     call theat(ti,tw,w10m,ric,t2m,td2m,slp,sehi,sehw,lahi,lahw)
     fc=ckb*(tf+tk-ti)

     res=qi+ilw+olwi+sehi+lahi+fc
     hres=res*res0

     if ((hres<=0.).or.(ti >= tk)) then     
        if (dt == 1.) then
           ti=ti-dt; dt=.1; cycle
        else if (dt == .1) then
           ti=ti-dt; dt=.01; cycle
        else if (dt == .01) then
           ti=ti-dt; dt=.001; cycle
        else if (dt == .001) then
           exit
        end if
     else
        ti=ti+dt; res0=res
        cycle
     end if
  end do

!  write(6,*) qi,qw
!  write(6,*) ilwi,ilww,olwi,olww

  !calc .heatbudjet
  if (ti.lt.tk) then
     hbi=-1.*fc*ric
  else
     hbi=(res-fc)*ric
  end if

  hbw=(qw+ilw+olww+sehw+lahw)*(1.-ric)
  if(ric.eq.0.) then
  hb=hbw
  elseif(ric.eq.1) then
  hb=hbi
  else
  hb=hbw+hbi
  endif
  sw=qi*ric+qw*(1.-ric)
  lw=ilw+olwi*ric+olww*(1.-ric)
  se=sehi*ric+sehw*(1.-ric)
  la=lahi*ric+lahw*(1.-ric)
end subroutine hb_all

subroutine heatTotemp(slp,t2m,td2m,w10m,ric,c,tw,lat,jday,hs,ali,alw,hb,ti)
  implicit none
  INTEGER(4)::jday
  REAL(4),PARAMETER::tk=273.15,tf=-1.86,cki=2.03,cks=0.31
  REAL(4)::slp,t2m,td2m,w10m,ric,c,tw,ti,lat
  REAL(4)::hi,hs,ali,alw
  REAL(4)::ckb
  REAL(4)::res0,res,hres,dt
  REAL(4)::fc,qi,qw,ilw,olwi,olww,sehi,sehw,lahi,lahw
  REAL(4)::hbi,hbw,hb,sw,lw,se,la

  ti=t2m-10.
  if (ti >= tk) ti=tk-10. 

!  ckb=cki*cks/(cki*hs+cks*hi)
  res0=1.; dt=1.               
  do 
     call swave(ali,alw,c,td2m,jday,lat,qi,qw)
     call lwave(c,ti,tw,t2m,ilw,olwi,olww)
     call theat(ti,tw,w10m,ric,t2m,td2m,slp,sehi,sehw,lahi,lahw)
     
     res=qi+ilw+olwi+sehi+lahi-hb
     hres=res*res0

     if ((hres<=0.).or.(ti >= tk)) then     
        if (dt == 1.) then
           ti=ti-dt; dt=.1; cycle
        else if (dt == .1) then
           ti=ti-dt; dt=.01; cycle
        else if (dt == .01) then
           ti=ti-dt; dt=.001; cycle
        else if (dt == .001) then
           exit
        end if
     else
        ti=ti+dt; res0=res
        cycle
     end if
  end do

!     print *,ti,qi,ilw,olwi,sehi,lahi,hb,hres

end subroutine heatTotemp

subroutine swave(ali,alw,clo,td2m,jday,lat,qi,qw)
  implicit none
  INTEGER(4)::i,jday
  REAL(4)::td2m,clo,lat,rjday
  REAL(4)::qi,qi1,qi2,qi3,qw,qw1,qw2,qw3
  REAL(4)::sz,szn,evi,evw,dec,decn,ha,han,st,an
  REAL(4)::ali,alw
  REAL(4),PARAMETER::s=1358.,pi=3.14159,b01=9.5,b02=265.3
! --- the Northern Hemisphere ---
  REAL(4),PARAMETER::a1=23.44, a2=172.
! --- the Southern Hemisphere ---
!  REAL(4),PARAMETER::a1=-23.44, a2=172.

  qi1=0.; qw1=0.; qi2=0.; qw2=0.; qi3=0.; qw3=0.
  st=0.; ha=0.; sz=0.; rjday=jday
  evi=0.; evw=0.
  
  ! --- calc. the vapor pressure -----------------------------
  evi=6.11*10.**((b01*(td2m-273.15))/(b02+td2m-273.15))
  evw=6.11*10.**((b01*(td2m-273.15))/(b02+td2m-273.15))
  ! --- calc. dec: declination -------------------------------
  dec=a1*cos((a2-rjday)*(pi/180.))
  
  do i=1,48              !every 30 minutes
     st=st+.5     
     ! --- calc. HA: hour angle ------------------------------
     ha=(12.-st)*pi/12.     
     ! --- calc. SZ: sine of zenith angle --------------------
     sz=sin(lat/180.*pi)*sin(dec/180*pi) &
          +cos(lat/180.*pi)*cos(dec/180.*pi)*cos(ha)
     if (sz <= 0.) sz=0.
     ! --- calc. Q: short wave radiation ---
     qi1=(s*(sz**2))/((sz+2.7)*evi*(1.e-3)+1.085*sz+0.10)
     qw1=(s*(sz**2))/((sz+2.7)*evw*(1.e-3)+1.085*sz+0.10)              
     ! --- sum ---
     qi2=qi2+qi1; qw2=qw2+qw1
  end do
  
  
  ! --- average ---
  qi3=qi2/48.; qw3=qw2/48.
  ! --- Qw&i is modified by a cloud factor -------------------
  decn=a1*cos((a2-rjday)*(pi/180.))
  han=(12.-12)*pi/12.
  szn=sin(lat/180.*pi)*sin(decn/180*pi) &
       +cos(lat/180.*pi)*cos(decn/180.*pi)*cos(han)
  an=asin(szn)*180./pi
  qi=qi3*(1-0.62*clo+0.0019*an) !(Andreas and Ackley,1982)
  qw=qw3*(1-0.62*clo+0.0019*an) !(Andreas and Ackley,1982)
  ! -- albedo ------------------------------------------------
  qi=(1-ali)*qi
  qw=(1-alw)*qw
end subroutine swave





subroutine lwave(clo,ti,sst,t2m,ila,olai,olaw)
  implicit none
  REAL(4)::clo,t2m,tw,ti,sst,ila,olai,olaw
  REAL(4),PARAMETER::sig=5.67E-8     !Stefan-Boltzmann constant
  REAL(4),PARAMETER::emw=.97,emi=.99 !emissivity
  tw=sst
!  write(6,*) sst
  ! --- incoming longwave radiation ---
!  ila=0.7855*(1.0+0.2232*clo**2.75)*sig*t2m**4.0 !(Maykut and Church, 1973)
  ! by Koenig-Langlo and Augstein, 1994
  ila=(0.765+0.22*clo**3.)*sig*t2m**4
 ! --- outgoing longwave radiation ---
  olai=-(emi*sig*ti**4)
  olaw=-(emw*sig*tw**4)
end subroutine lwave

subroutine theat(ti,sst,wg,ic,t2m,td2m,slp,sehi,sehw,lahi,lahw)
  implicit none 
  REAL(4)::sst,ic,hi,hw,lei,lew
  REAL(4)::ti,t2m,td2m,tsfc,wg,slp,tw
  REAL(4)::ch,ce,sehi,sehw,ea,esw,esi,lahi,lahw
  REAL(4)::a1,a2,b1,b2,c1,c2,p1,p2,cha,cea,s,s0
  REAL(4),PARAMETER::rhoa=1.3,cp=1004.,b01=9.5,b02=265.3
  hi=0.                     !sensible heat 
  hw=0.                     !sensible heat 
  ea=0.                     !vapor pressure of air  
  esw=0.                    !satulate vapor pressure of water  
  esi=0.                    !satulate vapor pressure of ice 
  lei=0.                    !latent heat 
  lew=0.                    !latent heat     
!  tw=272.15                 !water temp. -1C

  tw=sst

  if(wg == 0.0) wg=0.1
  ! ---- bulk transfer coefficients --- 
  if(ic < 0.15)then !no ice  
     tsfc=tw
  else
     tsfc=ic*ti+(1.-ic)*tw 
  endif
  ! --- set PARAMETER ---
  if(wg < 2.2)then
     a1=0.; a2=0.; b1=1.185; b2=1.23
     c1=0.; c2=0.; p1=-0.157; p2=-0.16
  else if((wg > 2.2).and.(wg < 5.))then
     a1=0.927; a2=0.969; b1=0.0546; b2=0.0521
     c1=0.; c2=0.; p1=1.; p2=1.
  else if((wg > 5.).and.(wg < 8.))then
     a1=1.15; a2=1.18; b1=0.01; b2=0.01
     c1=0.; c2=0.; p1=1.; p2=1.
  else if((wg > 8.).and.(wg < 25.))then

     a1=1.17; a2=1.196; b1=0.0075; b2=0.008
     c1=-0.00045; c2=-0.0004; p1=1.; p2=1.
  else if((wg > 25.).and.(wg < 50.))then
     a1=1.652; a2=1.68; b1=-0.017; b2=-0.016
     c1=0.; c2=0.; p1=1.; p2=1.
  endif
  cha=(a1+(b1*(wg**p1))+(c1*((wg-8)**2)))/1000.
  cea=(a2+(b2*(wg**p2))+(c2*((wg-8)**2)))/1000.
  ! --- 安定度 ---
  s0=(tsfc-t2m)*wg**(-2)
  s=s0*(abs(s0)/(abs(s0)+0.01))
  ! --- 中立の場合
  if (tsfc-t2m == 0.) then
     ch=cha
     ce=cea
     ! --- 安定の場合
  else if (tsfc-t2m < 0.) then
     if ((s > -3.3).and.(s < 0.)) then
        ch=cha*(0.1+0.03*s+0.9*exp(4.8*s))
        ce=cea*(0.1+0.03*s+0.9*exp(4.8*s))
     else if (s < -3.3) then
        ch=0.; ce=0.
     end if
     ! ------ 不安定の場合
  else if (tsfc-t2m > 0.) then
     ch=cha*(1.0+0.63*sqrt(s))
     ce=cea*(1.0+0.63*sqrt(s))
  endif
  ! --- Sensible heat ----------------------------------------
  sehi=rhoa*cp*ch*wg*(t2m-ti)
  sehw=rhoa*cp*ch*wg*(t2m-tw)
  ! --- Latent heat ------------------------------------------
  ea=6.11*10.**((b01*(td2m-273.15))/(b02+td2m-273.15))
  esw=6.11*10.**((b01*(tw-273.15))/(b02+tw-273.15))
  esi=6.11*10.**((b01*(ti-273.15))/(b02+ti-273.15))
  lahw=0.622*rhoa*2.52e6*ce/(slp/100.)*wg*(ea-esw)
  lahi=0.622*rhoa*2.86e6*ce/(slp/100.)*wg*(ea-esi)
end subroutine theat

function bilin(ipos,jpos,dx,dy,Undef,var1,var2,var3,var4,la1,la2,la3,la4)
  implicit none
  REAL(4) :: bilin
  REAL(4) :: w1,w2,w3,w4
  INTEGER(4),INTENT(in) :: ipos, jpos
  REAL(4),INTENT(in) :: dx, dy, Undef
  REAL(4),INTENT(in) :: var1,var2,var3,var4
  REAL(4),INTENT(in) :: la1,la2,la3,la4

  w1=1.
  w2=1.
  w3=1.
  w4=1.
!!$  if(la1.eq.1) w1=1/5.
!!$  if(la2.eq.1) w2=1/5.
!!$  if(la3.eq.1) w3=1/5.
!!$  if(la4.eq.1) w4=1/5.
  if((var1 /= Undef) .and. (var2 /= Undef) .and. &
       & (var3 /= Undef) .and. (var4 /= Undef)) then
    bilin = (1.-dx)*(1.-dy)*var1 &
         & +    dx *(1.-dy)*var2 &
         & +    dx *    dy *var3 &
         & +(1.-dx)*    dy *var4
!!$    bilin = ((1.-dx)*(1.-dy)*var1*w1 &
!!$          +    dx *(1.-dy)*var2*w2 &
!!$          +    dx *    dy *var3*w3 &
!!$          +(1.-dx)*    dy *var4*w4) &
!!$          /((1.-dx)*(1.-dy)*w1+dx*(1.-dy)*w2+dx*dy*w3+(1.-dx)*dy*w4)
  else
    bilin = Undef
  endif
end function bilin


subroutine MAPLL(aii,ajj,LAT,LONG)
      implicit none
      REAL(4) :: X,Y,ALAT,ALONG,E,E2,CDR,PI,SLAT,MC,SGN
      REAL(4) :: RE,RHO,SL,T,TC,LAT,LONG
      REAL(4) :: aii,ajj

      SGN=-1
      SLAT = 70.
      RE = 6378.273
      E2 = 0.006693883
      E =  0.081816153
!--------------------------------
!      CDR=57.29577951
      PI=3.141592654
      ALAT=abs(LAT)*PI/180.
      ALONG=LONG*PI/180.

      IF (ABS(ALAT).LT.PI/2.) GOTO 250
      X=0.0
      Y=0.0
      GOTO 999
  250 CONTINUE
      T=TAN(PI/4.-ALAT/2.)/((1.-E*SIN(ALAT))/(1.+E*SIN(ALAT)))**(E/2.)
      IF (ABS(90.-SLAT).LT.1.E-5) THEN
      RHO=2.*RE*T/((1.+E)**(1.+E)*(1.-E)**(1.-E))**(1/2.)
      ELSE
      SL=SLAT*PI/180.
      TC=TAN(PI/4.-SL/2.)/((1.-E*SIN(SL))/(1.+E*SIN(SL)))**(E/2.)
      MC=COS(SL)/SQRT(1.0-E2*(SIN(SL)**2))
      RHO=RE*MC*T/TC
      END IF
      Y=-RHO*SGN*COS(SGN*ALONG)
      X= RHO*SGN*SIN(SGN*ALONG)
  999 CONTINUE

      aii=(x+3950.-6.25/2.)/6.25+1
      ajj=1328-(y+3950.-6.25/2.)/6.25

      end subroutine MAPLL
