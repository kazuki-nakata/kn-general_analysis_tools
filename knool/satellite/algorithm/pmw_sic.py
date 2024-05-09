import math 
import numpy as np 
import os 
import glob 
import datetime

####D2.1 Sea Ice Concentration Algorithm Theoretical Basis Document (ATBD) 
#Doc Ref: SICCI-P2-ATBD(SIC) 
#Version: 1.0 
#Date: 22 September 2017 
#filename: sicci-p2-atbd-d2-1-sic-issue-1-0.pdf

def asi(tb85v, tb85h): 
 # ASI
    P0 = 47.0 #The coefficients are checked with Christian Melsheimer (see e-mail from him to Leif Friday, May 24, 2013 17:50). Same values are used for N and S. And for AMSR and SSM/I.
    P1 = 11.7
    A=np.matrix([[P1**3.0, P1**2.0, P1, 1.0],\ 
    [P0**3.0, P0**2.0, P0, 1.0],\ 
    [3.0*P1**3.0, 2.0*P1**2.0, P1, 0.0],\ 
    [3.0*P0**3.0, 2.0*P0**2.0, P0, 0.0]]) 
    b=np.matrix([1.0, 0.0, -0.14, -1.14]) 
    
    d=A.I * b.T 
    # d[0]=1.64/100000.0 d[1]=-0.0016 d[2]=0.0192 d[3]=0.971
    
    P = tb85v - tb85h 
    C = d[0] * P**3 + d[1] * P**2 + d[2] * P + d[3] 
    
    return C 

def bootstrap_f(tb18v, tb37v, tiepts): 
 
    tw18v = tiepts[6] 
    tw37v = tiepts[0] 
    tfy18v = tiepts[8] 
    tfy37v = tiepts[2] 
    tmy18v = tiepts[7] 
    tmy37v = tiepts[1] 


    if (tb18v-tw18v)==0: #open water
        cf=0.0
    else: 
        af = (tfy37v - tmy37v)/(tfy18v - tmy18v) 
        bf = (tmy37v - af*tmy18v) 
        qf = (tb37v - tw37v)/(tb18v - tw18v) 
        wf = (tw37v - qf*tw18v) 
        ti18vf = (bf - wf)/(qf - af) 
        cf = (tb18v - tw18v)/(ti18vf - tw18v) 
    return cf 
 
def bootstrap_p(tb37v, tb37h, tiepts):  
    tw37h = tiepts[3] 
    tw37v = tiepts[0] 
    tfy37h = tiepts[5] 
    tfy37v = tiepts[2] 
    tmy37h = tiepts[4] 
    tmy37v = tiepts[1] 
 
 
    if (tb37h-tw37h)==0: #open water
        cp=0.0
    else: 
        ap = (tfy37v - tmy37v) / (tfy37h - tmy37h) 
        bp = (tmy37v - ap * tmy37h) 
        qp = (tb37v - tw37v) / (tb37h - tw37h) 
        wp = (tw37v - qp * tw37h) 
        if (qp - ap)==0: 
            cp=-9.98
        else: 
            ti37hp = (bp - wp) / (qp - ap) 
            ti37vp = ap * ti37hp + bp 
            if (ti37vp - tw37v)==0: 
                cp=-9.98
            else: 
                cp = (tb37v - tw37v) / (ti37vp - tw37v) 
    return cp 

def bristol(tb18v, tb37v, tb37h, tiepts): 
 #Bristol ice concentration algorithm
    tw18v = tiepts[6] 
    tw37h = tiepts[3] 
    tw37v = tiepts[0] 
    tfy18v = tiepts[8] 
    tfy37h = tiepts[5] 
    tfy37v = tiepts[2] 
    tmy18v = tiepts[7] 
    tmy37h = tiepts[4] 
    tmy37v = tiepts[1] 

    xa = tmy37v + (1.045*tmy37h) + (0.525*tmy18v) 
    xd = tfy37v + (1.045*tfy37h) + (0.525*tfy18v) 
    xh = tw37v + (1.045*tw37h) + (0.525*tw18v) 
    xt = tb37v +(1.045*tb37h) + (0.525*tb18v) 
    ya = (0.9164*tmy18v) - tmy37v + (0.4965*tmy37h) 
    yd = (0.9164*tfy18v) - tfy37v + (0.4965*tfy37h) 
    yh = (0.9164*tw18v) - tw37v + (0.4965*tw37h) 
    yt = (0.9164*tb18v)- tb37v + (0.4965*tb37h) 
    a_ht = (yt - yh)/(xt - xh) 
    b_ht = yh - (a_ht*xh) 
    a_da = (ya - yd)/(xa - xd) 
    b_da = yd - (a_da*xd) 
    xi = (b_da - b_ht)/(a_ht - a_da) 
    cf = (xt - xh)/(xi - xh) 
    c = cf 
    return c

def calval(tb37v,tb18v,tiepts): 
    tw18v = tiepts[6] 
    tw37v = tiepts[0] 
    tfy18v = tiepts[8] 
    tfy37v = tiepts[2] 
    tmy18v = tiepts[7]

    tmy37v = tiepts[1] 


    A=np.matrix([[tw37v, tw18v, 1.0],\ 
    [tfy37v, tfy18v, 1.0],\ 
    [tmy37v, tmy18v, 1.0]]) 
    b=np.matrix([0.0, 1.0, 1.0]) 

    d=A.I * b.T 
    C=d[0]*tb37v+d[1]*tb18v+d[2] 
    return C 
 
def nasa(tb18v, tb18h, tb37v, tiepts): 
    #NASA-Team ice concetration algorithm
    ow18v = tiepts[6] 
    ow18h = tiepts[9] 
    ow37v = tiepts[0] 
    fy18v = tiepts[8] 
    fy18h = tiepts[11] 
    fy37v = tiepts[2] 
    my18v = tiepts[7] 
    my18h = tiepts[10] 
    my37v = tiepts[1] 
    a0 = - ow18v + ow18h 
    a1 = ow18v + ow18h 
    a2 = my18v - my18h - ow18v + ow18h 
    a3 = - my18v - my18h + ow18v + ow18h 
    a4 = fy18v - fy18h - ow18v + ow18h 
    a5 = - fy18v - fy18h + ow18v + ow18h 
    b0 = - ow37v + ow18v 
    b1 = ow37v + ow18v 
    b2 = my37v - my18v - ow37v + ow18v 
    b3 = - my37v - my18v + ow37v + ow18v 
    b4 = fy37v - fy18v - ow37v + ow18v 
    b5 = - fy37v - fy18v + ow37v + ow18v 
    gr = (tb37v - tb18v)/(tb37v + tb18v) 
    pr = (tb18v - tb18h)/(tb18v + tb18h) 
    d0 = (-a2*b4) + (a4*b2) 
    d1 = (-a3*b4) + (a5*b2) 
    d2 = (-a2*b5) + (a4*b3) 
    d3 = (-a3*b5) + (a5*b3) 
    dd = d0 + d1*pr + d2*gr + d3*pr*gr 

    f0 = (a0*b2) - (a2*b0) 
    f1 = (a1*b2) - (a3*b0) 
    f2 = (a0*b3) - (a2*b1) 
    f3 = (a1*b3) - (a3*b1) 
    m0 = (-a0*b4) + (a4*b0) 
    m1 = (-a1*b4) + (a5*b0) 
    m2 = (-a0*b5) + (a4*b1) 
    m3 = (-a1*b5) + (a5*b1) 
    cf = (f0 + f1*pr + f2*gr + f3*pr*gr)/dd 
    cm = (m0 + m1*pr + m2*gr + m3*pr*gr)/dd 
    cf = cf 
    cm = cm 
    ct = cm + cf 
    return ct, cm 
 
def near90(tb85v, tb85h, tiepts): 
    #tmy85v = tiepts[30]
    tfy85v = tiepts[31] 
    #tmy85h = tiepts[32]
    tfy85h = tiepts[33] 
    tw85v = tiepts[34] 
    tw85h = tiepts[35] 

    P = tb85v - tb85h 
    P0 = tw85v - tw85h 
    P1 = tfy85v - tfy85h

    A=np.matrix([[P1**3.0, P1**2.0, P1, 1.0],\ 
    [P0**3.0, P0**2.0, P0, 1.0],\ 
    [3.0*P1**3.0, 2.0*P1**2.0, P1, 0.0],\ 
    [3.0*P0**3.0, 2.0*P0**2.0, P0, 0.0]]) 
    b=np.matrix([1.0, 0.0, -0.14, -1.14]) 

    d=A.I * b.T 

    #d=np.linalg.solve(A,b.T)
    C = d[0] * P**3 + d[1] * P**2 + d[2] * P + d[3] 
    return np.float(C) 

def near90_linear_dyn(tb85v, tb85h, tiepts): 
    tmy85v = tiepts[30] 
    tfy85v = tiepts[31] 
    tmy85h = tiepts[32] 
    tfy85h = tiepts[33] 
    tw85v = tiepts[34] 
    tw85h = tiepts[35] 

    PFY = tfy85v - tfy85h 
    PMY = tmy85v - tmy85h 
    PW = tw85v - tw85h 

    PI = (PFY + PMY)/2
    P = tb85v - tb85h 

    c = (P - PW) / (PI-PW) 

    return c 

def norsex(tb18v,tb37v,sensor_name,lat): 
 
    SAT = 260.0
    T_sa = 270.0
    T_a = 250.0
    To = 272.0

    tau_sa19v = 0.0610
    tau_sa37v = 0.1000
    tau_a19v = 0.0440
    tau_a37v = 0.0700
    TB_w_19v, TB_w_37v, TB_fy_19v, TB_fy_37v, TB_my_19v, TB_my_37v = norsex_TPs(sensor_name,lat) 
    t_atm_surf=SAT #Initialize atmospheric surface temperature

    for i in range(0,2): 

    #interpolate opacity between arctic and subarctic values:
    tau19v = tau_a19v + (t_atm_surf - T_a) * (tau_sa19v - tau_a19v) / (T_sa - T_a) 
    tau37v = tau_a37v + (t_atm_surf - T_a) * (tau_sa37v - tau_a37v) / (T_sa - T_a) 

    #Constants to be used in computing ice concentrations:
    a11 = TB_fy_19v - TB_w_19v 
    a21 = TB_fy_37v - TB_w_37v 
    a12 = TB_my_19v - TB_w_19v 
    a22 = TB_my_37v - TB_w_37v 
    d_coef = a11 * a22 - a12 * a21 

    #find emitted brightness temperature at the surface by correcting for
    #atmospheric disturbances:
    TB_surf_19v = (tb18v - t_atm_surf * (2.0 * tau19v - tau19v**2.0 + 0.01)) / (1.0 - 2.0 * tau19v + 
    tau19v**2.0 - 0.01) 
    TB_surf_37v = (tb37v - t_atm_surf * (2.0 * tau37v - tau37v**2.0 + 0.01)) / (1.0 - 2.0 * tau37v + 
    tau37v**2.0 - 0.01) 
    #Find new atmospheric surface brightness temperature and mean surface
    #emissions by solving for first year and multi-year ice concentrations.
    c1 = TB_surf_19v - TB_w_19v 
    c2 = TB_surf_37v - TB_w_37v 
    Cmy = (a11 * c2 - a21 * c1) / d_coef 
    Cfy = (a22 * c1 - a12 * c2) / d_coef 
    CT = Cfy + Cmy 

    t_atm_surf = To + (SAT - To) * CT 
 
    return CT 
 
def P37(tb37v, tb37h,tiepts): 
    # instead of NRL
    tw37h = tiepts[3] 
    tw37v = tiepts[0] 
    tfy37h = tiepts[5] 
    tfy37v = tiepts[2] 
    tmy37h = tiepts[4] 
    tmy37v = tiepts[1] 

    PFY = tfy37v - tfy37h 
    PMY = tmy37v - tmy37h 
    PW = tw37v - tw37h 

    PI = (PFY + PMY)/2
    P = tb37v - tb37h 

    c = (P - PW) / (PI-PW) 
    return c 

def onechannel(tb6h, tiepts): 
    #Simple 1 channel algorithm

    fy6h = tiepts[17] 
    my6h = tiepts[16] 

    ow6h = 82.3 
    i6h = (fy6h+my6h)/2.0 
    ct = (tb6h - ow6h)/(i6h - ow6h) 
    return ct 
 
def osisaf(c0,c1,t): 
    wc = (abs(t - c0) + t - c0) / (2.0 * t) 
    c = c1 * (1.0 - wc) + wc * c0 
    if c0 < 0: 
    c = c0 

    return c 
 
def sicci(c0,c1): 
 
    if c0<0.7: 
        wCF=1.0

    if (c0 >= 0.7 and c0 < 0.9): 
        wCF=1.0-(c0-0.7)/(0.9-0.7) 

    if c0 >= 0.9: 
        wCF=0.0

    wBR = 1.0-wCF 

    c = c0 * wCF + c1 * wBR 

    return c 

def P90(tb85v, tb85h): 
    X=(tb85v-tb85h) 
    P=(X-2.63)/0.752

    d3=1.64/100000.0
    d2=-0.0016
    d1=0.0192
    d0=0.971

    c1 = d3 * P**3.0 + d2 * P**2.0 + d1 * P + d0 
    c = c1+(P-8)/700 #to adjust near SIC0
    if (P>48): 
    c=-0.026 #to prevent large P85 giving ice
    if (P<8.5): 
    c=1.03 #to prevent low P85 losing ice
    return c 
 
def pr(tb18v, tb18h, tb37v, tb37h, tiepts): 
 #Simple polarization ratio algorithm
 
    ow18v = tiepts[6] 
    ow18h = tiepts[9] 
    ow37v = tiepts[0] 
    ow37h = tiepts[3] 
    fy18v = tiepts[8] 
    fy18h = tiepts[11] 
    fy37v = tiepts[2] 
    fy37h = tiepts[5] 
    my18v = tiepts[7] 
    my18h = tiepts[10] 
    my37v = tiepts[1] 
    my37h = tiepts[4] 
    i18v = (fy18v+my18v)/2
    i18h = (fy18h+my18h)/2
    i37v = (fy37v+my37v)/2
    i37h = (fy37h+my37h)/2 
    PR18 = (tb18v - tb18h)/(tb18v + tb18h) 
    PR37 = (tb37v - tb37h)/(tb37v + tb37h) 
    c18 = (ow18v*(1 - PR18) - ow18h*(1 + PR18))/(PR18*(i18v + i18h - ow18v - ow18h) - (i18v - i18h - ow18v 
    + ow18h)) 
    c37 = (ow37v*(1 - PR37) - ow37h*(1 + PR37))/(PR37*(i37v + i37h - ow37v - ow37h) - (i37v - i37h - ow37v 
    + ow37h)) 
    c_old = (c18 + c37)/2
    c = c_old/(2-c_old) 
    return c, PR18, PR37 

def tud(c85, cf): 
    #TUD ice concentration alogorithm
    if ((c85>0) & (cf>10)): 
    c = np.sqrt(cf*c85) 
    else: 
    c = cf 

    return c 

def P10(tb10v, tb10h, tiepts): 
    #Simple 2 channel algorithm 10 GHz

    tw10h = tiepts[21] 
    tw10v = tiepts[18] 
    tfy10h = tiepts[23] 
    tfy10v = tiepts[20] 
    tmy10h = tiepts[22] 
    tmy10v = tiepts[19] 

    PFY = tfy10v - tfy10h 
    PMY = tmy10v - tmy10h 
    PW = tw10v - tw10h 

    PI = (PFY + PMY)/2
    P = tb10v - tb10h 

    c = (P - PW) / (PI-PW) 
    return c 
 
def P18(tb18v, tb18h,tiepts): 
 
    tw18v = tiepts[6] 
    tw18h = tiepts[9] 
    tfy18v = tiepts[8] 
    tfy18h = tiepts[11] 
    tmy18v = tiepts[7] 
    tmy18h = tiepts[10] 

    PFY = tfy18v - tfy18h 
    PMY = tmy18v - tmy18h 
    PW = tw18v - tw18h 

    PI = (PFY + PMY)/2
    P = tb18v - tb18h 

    c = (P - PW) / (PI-PW) 
    return c 
 
def UMass(tb18v,tb37v,tiepts): 
 #The code is based on C. T. Swift, L. S. Fedor, and R. O. Ramseier, ?An Algorithm 
 #to Measure Sea Ice Concentration With Microwave Radiometers,? Journal of Geophysical 
 #Research, vol. 90, no. C1, pages 1087 - 1099, 1985.
    tw19v = tiepts[6] 
    tw37v = tiepts[0] 
    tfy19v = tiepts[8] 
    tfy37v = tiepts[2] 
    tmy19v = tiepts[7] 
    tmy37v = tiepts[1] 
    #solution of the equations (11)-(12) in Swift et al 1985
    #here we use brightness temperatures instead!! (Rasmus 2012)
    #e19v=(TB19v-13)./(Ts-12);
    #e37v=(TB37v-26)./(Ts-26);
    a1 = (tfy19v - tb18v) / (tfy19v - tw19v) 
    a2 = (tfy19v - tmy19v) / (tfy19v - tw19v) 
    a3 = (tfy37v - tb37v) / (tfy37v - tmy37v) 
    a4 = (tfy37v - tw37v) / (tfy37v - tmy37v) 
    fw = (a1 - a2 * a3) / (1.0 - a2 * a4) 
    Cmy = a3 - fw * a4 
    Cfy = 1.0 - fw - Cmy 
    CT = Cmy + Cfy 
    return CT 

def norsex_TPs(sensor_name,lat): 
    # values produced by Leif, described in PVASR
    if lat >= 0: 
    # Northern hemisphere:
        if (sensor_name == 'AMSRE' or sensor_name == 'AMSR2'): 
            TB_w_19v = 170.01
            TB_w_37v = 193.19
            TB_fy_19v = 251.17
            TB_fy_37v = 244.47
            TB_my_19v = 222.11
            TB_my_37v = 184.02
        elif sensor_name == 'SSMI': 
            TB_w_19v = 171.56
            TB_w_37v = 191.87
            TB_fy_19v = 251.91
            TB_fy_37v = 241.53
            TB_my_19v = 219.20
            TB_my_37v = 175.93
        elif sensor_name == 'SMMR': 
            TB_w_19v = 162.61
            TB_w_37v = 190.80
            TB_fy_19v = 251.17
            TB_fy_37v = 244.47
            TB_my_19v = 222.11
            TB_my_37v = 184.02

        elif lat < 0: 
        # Southern hemisphere:
        if (sensor_name == 'AMSRE' or sensor_name == 'AMSR2'): 
            TB_w_19v = 171.86
            TB_w_37v = 196.65
            TB_fy_19v = 258.41
            TB_fy_37v = 252.57
            TB_my_19v = 244.39
            TB_my_37v = 219.62
        elif sensor_name == 'SSMI': 
            TB_w_19v = 171.52
            TB_w_37v = 192.94
            TB_fy_19v = 259.93
            TB_fy_37v = 253.25
            TB_my_19v = 244.59
            TB_my_37v = 219.59
        elif sensor_name == 'SMMR': 
            TB_w_19v = 160.77
            TB_w_37v = 190.92
            TB_fy_19v = 258.41
            TB_fy_37v = 252.57
            TB_my_19v = 244.39
            TB_my_37v = 219.62

    return TB_w_19v, TB_w_37v, TB_fy_19v, TB_fy_37v, TB_my_19v, TB_my_37v 
