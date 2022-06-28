from matplotlib.style import reload_library
import numpy as np

def calc_ice_production(hb):
  dt=24.*60.*60.
  dhi=np.abs(hb)*dt/rhoi/lf
  return dhi

def calc_swave(al,clo,td2m,jday,lat,hour):
    s=1358.;b01=9.5;b02=265.3
    a1=-23.44; a2=172.
    st=0.; ha=0.; sz=0.; rjday=jday
    
#   --- calc the vapor pressure -----------------------------
    ev=6.11*10.**((b01*(td2m-273.15))/(b02+td2m-273.15))
#   --- calc dec: declination -------------------------------
    dec=np.radians(a1*np.cos(np.radians(a2-rjday)))    
    rad_la=np.radians(lat)
    
    st=hour
    ha=(12.-st)*np.pi/12.
    sz=np.sin(rad_la)*np.sin(dec)+np.cos(rad_la)*np.cos(dec)*np.cos(ha)
    if (sz <= 0.): sz=0.
    q1 = np.where(sz<0,0,(s*(sz**2))/((sz+2.7)*ev*(1.e-3)+1.085*sz+0.10))
    
    an=np.angle(np.arcsin(sz))
    q2=(1-al)*q1*(1-0.62*clo+0.0019*an) #(Andreas and Ackley,1982)
    return q2

def swave_seaice(clo,td2m,jday,lat,hour):
    return swave(ali,clo,td2m,jday,lat,hour)

def swave_water(clo,td2m,jday,lat,hour):
    return swave(alw,clo,td2m,jday,lat,hour)


def calc_lwave_MC1973(clo,st,t2m,em):
  # --- incoming longwave radiation ---
    ila=0.7855*(1.0+0.2232*clo**2.75)*sig*t2m**4.0 #(Maykut and Church, 1973)
    ola=-(em*sig*st**4)
    return ila,ola

def calc_lwave_KA1994(clo,st,t2m,em):
  # by Koenig-Langlo and Augstein, 1994
    ila=(0.765+0.22*clo**3.)*sig*t2m**4
    ola=-(em*sig*st**4)
    return ila,ola

def calc_theat_K1975(ti,tw,wg,ic,t2m,td2m,slp):
    b01=9.5; b02=265.3
    wg = np.where(wg==0,0.1,wg)
    tsfc = np.where(ic<0.15,tw,ic*ti+(1.0-ic)*tw)
  
  # --- set PARAMETER ---
   # ---- bulk transfer coefficients --- 
    calc_bulk_cha = lambda wg, p: (p["ah"]+(p["bh"]*(wg**p["ph"]))+(p["ch"]*((wg-8)**2)))/1000.
    calc_bulk_cea = lambda wg, p: (p["ae"]+(p["be"]*(wg**p["pe"]))+(p["ce"]*((wg-8)**2)))/1000.
    
    cha = np.where(wg<2.2, calc_bulk_cha(wg,params2[0]),
                np.where((wg>0.02)&(wg<5.0),calc_bulk_cha(wg,params2[1]),
                np.where((wg > 5.)&(wg < 8.),calc_bulk_cha(wg,params2[2]),
                np.where((wg > 8.)&(wg < 25.),calc_bulk_cha(wg,params2[3]),            
                np.where((wg > 25.)&(wg < 50.),calc_bulk_cha(wg,params2[4]),np.NaN)))))
  
    cea = np.where(wg<2.2, calc_bulk_cea(wg,params2[0]),
                np.where((wg>0.02)&(wg<5.0),calc_bulk_cea(wg,params2[1]),
                np.where((wg > 5.)&(wg < 8.),calc_bulk_cea(wg,params2[2]),
                np.where((wg > 8.)&(wg < 25.),calc_bulk_cea(wg,params2[3]),            
                np.where((wg > 25.)&(wg < 50.),calc_bulk_cea(wg,params2[4]),np.NaN)))))  
  
  
  # --- 安定度 ---
    s0=(tsfc-t2m)*wg**(-2)
    s=s0*(np.abs(s0)/(np.abs(s0)+0.01))
    dtsfc=tsfc-t2m
  # ---calc transfer coefficients----------------------
    ch=np.where(dtsfc==0, cha, # 中立
       np.where(dtsfc>0, cha*(1.0+0.63*np.sqrt(s)), #不安定
       np.where((s > -3.3)&(s < 0.),cha*(0.1+0.03*s+0.9*np.exp(4.8*s)),0))) #安定
     
    ch=np.where(dtsfc==0, cea, # 中立
       np.where(dtsfc>0, cea*(1.0+0.63*np.sqrt(s)), #不安定
       np.where((s > -3.3)&(s < 0.),cea*(0.1+0.03*s+0.9*np.exp(4.8*s)),0))) #安定
    
 # ! --- Sensible heat ----------------------------------------
    sehi=rhoa*cp*ch*wg*(t2m-ti)
    sehw=rhoa*cp*ch*wg*(t2m-tw)
 # ! --- Latent heat ------------------------------------------
    ea=6.11*10.**((b01*(td2m-273.15))/(b02+td2m-273.15))
    esw=6.11*10.**((b01*(tw-273.15))/(b02+tw-273.15))
    esi=6.11*10.**((b01*(ti-273.15))/(b02+ti-273.15))
    lahw=0.622*rhoa*2.52e6*ce/(slp/100.)*wg*(ea-esw)
    lahi=0.622*rhoa*2.86e6*ce/(slp/100.)*wg*(ea-esi)
    return sehi,sehw,lahw,lahi

def calc_heat_flux(slp,t2m,td2m,w10m,ric,c,tw,lat,jday,hi,hs,ali,alw,hb,sw,lw,se,la,ti):
  tk=273.15,tf=-1.86,cki=2.03,cks=0.31
       
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

