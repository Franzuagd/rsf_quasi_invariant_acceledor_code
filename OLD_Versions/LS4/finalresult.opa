
Energy = 3.000000;

    BetaX   = 4.8032442; AlphaX  = 0.0000024;
    EtaX    = 0.0027806; EtaXP   = 0.0000000;
    BetaY   = 2.2903690; AlphaY  = 0.0000004;
    EtaY    = 0.0000000; EtaYP   = 0.0000000;
{-----Parametros-----}
KQF1 =3.742153457895532e+00;
KQD2 =-7.933030565890353e+00;
KQF3 =9.383357074512245e+00;
KQF1i =-1.961256869785227e+00;
KQD2i =-7.622311324145445e+00;
KQF3i =9.502339759326459e+00;
KQF1d =1.528008695772636e+00;
KQD2d =-6.400798167795065e+00;
KQF3d =1.049016046874180e+01;
KD1 =2.8791011093979669e+03;
KD4 =-4.6619314170935371e+02;
KD3D =1.0059424467319975e+02;
KD4D =-1.6765338409406527e+02;
KO2 =-2.0708328502237236e+12;
KO3 =2.4508597362125468e+08;
KO4 =-3.2161485700407800e+11;
KO3D =-8.7899899320609665e+10;
KO4D =8.0528105791989771e+11;
KD2 =-7.9630003934210026e+02;
KD3 =8.4643078766863107e+02;
KA1q =-3.291355350846638e+00;
KA1p =-4.138813010088813e+00;
{----- Fin-----------}
KO1 =-2.526041097299107e+12;
{----- Variables ----------------------------------------------------}
LD1       =   1.5037455548953357e-02 ;
LD2       =   1.1984645488606541e-01 ;
LD3       =   9.9224417359370098e-02 ;
LD4       =   1.2034242941591648e-01 ;
LD1i      =     2.0;                  
LD2i      =   1.0856575439288529e-02 ;
LD3i      =   1.6540270678201185e-02 ;
LD4i      =   1.4096214799570003e-01 ;
LD1d      =   1.5038380720489061e-02 ;
LD2d      =   4.3024817452685354e-02 ;
LD3d      =   1.0980500234939020e-01 ;
LD4d      =   1.1254100559506744e-01 ;
LQF1     =   2.7152587874047412e-01 ;
LQD2     =   2.8964226476521576e-01 ;
LQF3     =   2.7632435252776760e-01 ;
LQF1i    =   1.0335385712912697e-02 ;
LQD2i    =   1.5232578176198157e-02 ;
LQF3i    =   2.3313438322936453e-01 ;
LQF1d    =   2.9307154869324176e-01 ;
LQD2d    =   2.7535600743716460e-01 ;
LQF3d    =   2.9731377167539702e-01 ;                 

LA1      =   8.0113085336170697e-01 ;
LA1p     =   6.8888537112046666e-01 ;

{----- Table of elements ----------------------------------------------------}

{D1    : Drift, L = LD1 ,  Ax = 5.00, Ay = 5.00;}
{D2    : Drift, L = LD2 ,  Ax = 5.00, Ay = 5.00;}
{D3    : Drift, L = LD3 ,  Ax = 5.00, Ay = 5.00;}
{D4    : Drift, L = LD4 ,  Ax = 5.00, Ay = 5.00;}
D1i   : Drift, L = LD1i,  Ax = 5.00, Ay = 5.00;
D2i   : Drift, L = LD2i,  Ax = 5.00, Ay = 5.00;
D3i   : Drift, L = LD3i,  Ax = 5.00, Ay = 5.00;
D4i   : Drift, L = LD4i,  Ax = 5.00, Ay = 5.00;
D1d   : Drift, L = LD1d,  Ax = 5.00, Ay = 5.00;
D2d   : Drift, L = LD2d,  Ax = 5.00, Ay = 5.00;
{D3d   : Drift, L = LD3d,  Ax = 5.00, Ay = 5.00;  }                                                                         
{D4d   : Drift, L = LD4d,  Ax = 5.00, Ay = 5.00;}
{O1 no se usa}  
 
D1  : Sextupole, L = LD1,   K =  KD1*3, N = 4, Ax = 5.00, Ay = 5.00;
D4  : Sextupole, L = LD4,   K =  KD4*3, N = 4, Ax = 5.00, Ay = 5.00;
D3d : Sextupole, L = LD3d,  K =  KD3D*3, N = 4, Ax = 5.00, Ay = 5.00;
D4d : Sextupole, L = LD4d,  K =  KD4D*3, N = 4, Ax = 5.00, Ay = 5.00;
O1     : Octupole, L = 0,  K =  KO1*1e-8, N = 4;                   
O2     : Octupole, L = 0,  K =  KO2*1e-8, N = 4;                   
O3     : Octupole, L = 0,  K =  KO3*1e-8, N = 4;                   
O4     : Octupole, L = 0,  K =  KO4*1e-8, N = 4;                   
O3d    : Octupole, L = 0,  K =  KO3D*1e-8, N = 4;                   
O4d    : Octupole, L = 0,  K =  KO4D*1e-8, N = 4;                   
                                                                                      
D2  : Sextupole, L = LD2,   K =  KD2*3, N = 4, Ax = 5.00, Ay = 5.00;
D3  : Sextupole, L = LD3,   K =  KD3*3, N = 4, Ax = 5.00, Ay = 5.00;

QF1     : Quadrupole, L = LQF1 , K =  KQF1 , Ax = 5.00, Ay = 5.00;
QD2     : Quadrupole, L = LQD2 , K =  KQD2 , Ax = 5.00, Ay = 5.00; 
QF3     : Quadrupole, L = LQF3 , K =  KQF3 , Ax = 5.00, Ay = 5.00;
QF1i    : Quadrupole, L = LQF1i, K =  KQF1i, Ax = 5.00, Ay = 5.00;
QD2i    : Quadrupole, L = LQD2i, K =  KQD2i, Ax = 5.00, Ay = 5.00;
QF3i    : Quadrupole, L = LQF3i, K =  KQF3i, Ax = 5.00, Ay = 5.00;
QF1d    : Quadrupole, L = LQF1d, K =  KQF1d, Ax = 5.00, Ay = 5.00;
QD2d    : Quadrupole, L = LQD2d, K =  KQD2d, Ax = 5.00, Ay = 5.00;
QF3d    : Quadrupole, L = LQF3d, K =  KQF3d, Ax = 5.00, Ay = 5.00;

A1   : Bending, L = LA1 , T = 2.57142857, K = KA1q, T1 = 0, T2 = 0, Ax = 5.00, Ay = 5.00;
A1p  : Bending, L = LA1p, T = 2.57142857, K = KA1p, T1 = 0, T2 = 0, Ax = 5.00, Ay = 5.00;

                       
{----- Table of segments ----------------------------------------------------}

 MCELL :   D1,  QF1,  O2,   D2,  QD2,  O3,   D3,  QF3,  O4,   D4,   A1,   D4,  O4,  QF3,   D3,  O3,  QD2,   D2,  O2,  QF1,   D1;              
 ECELL :  D1i, QF1i,  D2i, QD2i,  D3i, QF3i,  D4i,  A1p,  O4d,  D4d, QF3d, O3d,  D3d, QD2d,  D2d, QF1d,  D1d;
                                                                                      
 Cell  : ECELL,  5*MCELL,  -ECELL; 
{Cell  : ECELL,   -ECELL;}
RING  : 20*CELL;

{..rs\Antillon\Dropbox\Harim\ESRFconSLS\ESRF_SLS_21_mayo_2020_20C.opa.txt}
