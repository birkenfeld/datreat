#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
//=================================================================================
//  generalized local reptation expression along the lines of deGennes, 
//  but with finite summation of integrals and lenthscale, 
//  timescale and fluctuation ratio as parameters
//  Journal de Physique, 1981, 42 (5), pp.735-740. <10.1051/jphys:01981004205073500>
//  and modified to account for the Gaussian chain condition (Debye)
//==================================================================================
//
// Author: Michael Monkenbusch, JCNS-1, Forschungszentrum Juelich
//         m.monkenbusch@fz-juelich.de
//
//==================================================================================
//  GNU GENERAL PUBLIC LICENSE 
//  Version 3, 29 June 2007
//==================================================================================
//
//  complile with 
//  gfortran -c reptation_module -fopenmp -O2
//
//==================================================================================

#define pi (3.141592654)

typedef struct {  
                 double sqt;
                 double sq0;
               } sqt_type;

 




  double local_reptationdr( 
    double q      , 
    double t      ,   
    double lseg   ,  // lsegment;
    double w      ,  // Rouse rate;
    double n      ,  // No Segments;
    double ne       // No Segments per entanglement strand;
                          ) {

    double xlimser = 15.0;  // where to switch erf expression to series expansion;

    double T1, T2, dec;
    double x, xa, xa2, xa3, xa4;
    double Z;

    double val;

    double lseg2, lseg4, x2, x3, x4, x5, q2,n2;
    
    lseg2 = lseg  * lseg;
    lseg4 = lseg2 *lseg2;
    q2    = q*q;
    n2    = n*n;
    


    Z    = n/ne;

    T1 = 72.0 * (exp(-q*q * n * lseg*lseg / 6.0) + q*q * n * lseg*lseg / 6.0 - 1.0) / (q*q*q*q) / (lseg*lseg*lseg*lseg);
////>> limit(T1(q=0)) is  n**2  ;
////>> here we tweak the normalisation to the new interpretation;
    T1 = T1 / (n*n);   //// Normalisation to 1 at q=0 ;
//// now T2 the local reptation part;
    x  = sqrt(lseg*lseg*lseg*lseg*q*q*q*q*w*t/36.0);

    x2=x*x;
    x3=x2*x;
    x4=x3*x;
    x5=x4*x;
    

    xa = n/(2*sqrt(w*t));
    xa2 = xa*xa;
    xa3 = xa2*xa;
    xa4 = xa3*xa;


    if( abs(x) < xlimser ) {;
      dec = exp(x2) *(erfc(x)-erfc(x+xa));
    } else {;
      dec =  1.0/(sqrt(pi)*x)-1.0/(2*sqrt(pi)*x3)+3.0/(4*sqrt(pi)*x5) 
          - (x4-xa*x3+(xa2-1.0/2.0)*x2+(-xa3+3.0/2.0*xa)*x+xa4-3*xa2+3.0/4.0) 
           * exp(-xa2)*exp(-2*xa*x)/(sqrt(pi)*x5);
    };

    T2  = ((exp(-1/6*lseg2*n*q2-n2/(4*w*t))-1)*2/3*ne*sqrt(w*t/pi)+ 
        ne*(n/3+lseg2*q2*w*t/9) 
        * dec );

    T2  = T2*3.0 / (n*ne);   //// Normalisation to T2(t=0) = 1;

    val = 1.0/(1+1/Z) *  (T1 + (2.0/3.0 + (1.0/3.0) * T2) / Z );

    return (val);

  } // endfunction local_reptationdr;



       sqt_type nrouse_ngcorr(
//      =========================================================================;
//;
// Rouse expression for a chain of finite length:;
// Input parameters:;
//    q             ----> momentum transfer in A**-1;
//    t             ----> time in nano-sec;
//    N             ----> effective number of segments;
//    R             ----> Re of "blob", l will be computed from N and R;
//    W             ----> Rouse rate;
//    alpha0        ----> strength of Non-Gaussianity (NG) correction;
//    talphamax     ----> peak of the NG-alpha function;
//    talphawd      ----> width of the (log) distribution alpha(t);
//    sqt(1:2)      <--- [S(Q), S(Q,t)];
// ------------------------------------------------------------------------------;
//;
       double q  ,
       double t  ,
       int N  ,
       double R  ,
       double Wl4  ,
       double alpha0  ,
       double talphamax  ,
       double talphawd  
                             ) {

       sqt_type Sq0t;



       int nn,mm,ip;

       double Sq, Sqt;
       double Dr, l, ff2, W;

       double cosarray[N][N], ewfac[N];
       double rmm, fqq, fqq0;

       int  i;

       double sum;



   double alpha(double t) { //// see Guenza PHYSICAL REVIEW E 89, 052603 (2014) Eq(4);
      double a;
   // since "contained//  in th_nrosueaplha, the decribing parameters if not declared explictitly here;
   // the model;
      a = ((log(t+1e-3)-log(talphamax))/talphawd);
      a = alpha0 * exp(-(a*a) / 2.0) ;
      return (a);
   } // endfunction alpha;



//       integer iout;

       if(N <= 0) {
         Sq0t.sq0 = 0;
         Sq0t.sqt = 0;
         printf("\nERROR: number of segments=%d (nrouse_ngcorr) is <=0\n",N);
         return(Sq0t);
       };

// ---- determine the segment length l ----;
       l = sqrt(R*R/N)       ;
       Dr  = 0;               // in this module: no com-diffusion of the "blobs";
       W = Wl4/(l*l*l*l);

//$OMP PARALLEL DO     ;
        for(nn=0;nn<N;nn++){
          for(ip=0;ip<N;ip++) {
              cosarray[nn][ip] = cos((pi*(ip+1)*(nn+1))/(double)(N)) / (ip+1) ;
          }
        }
//$OMP } // endPARALLEL DO   ;

//$OMP PARALLEL DO    ;
       for(i=0;i<N;i++) {
         ewfac[i] = (1.0-exp(-2*W*(1-cos((pi*(i+1))/(double)(N)))*t)) ;
       }
//$OMP } // endPARALLEL DO    ;
// ---- init sums ----;
       Sq  = 0;
       Sqt = 0;
       ff2  = -2*N*(l*q)*(l*q)/(3*pi*pi);
// ---- Do the sums -----;
       rmm = 0;
//$OMP PARALLEL DO REDUCTION(+:rmm);  ip SUMME explicit AUFLÖSEN HIER !!
       for(mm=0;mm<N;mm++){
             sum = 0.0;
             for(ip=0;ip<N;ip++) {
               sum = sum + cosarray[mm][ip] * cosarray[mm][ip] *  ewfac[ip];
             }
             rmm = rmm + 4.0*N*l*l/(pi*pi) * sum;
       }
//$OMP } // endPARALLEL DO;
       rmm = rmm/N;
       fqq  = 1.0 -q*q * rmm/12.0 * alpha(t);       //// see Guenza PHYSICAL REVIEW E 89, 052603 (2014) Eq(4);
       fqq0 = 1.0;                                   //// since rmm(0)=0;

//$OMP PARALLEL DO REDUCTION(+:Sq,Sqt); ip SUMME explicit AUFLÖSEN !!
       for(nn=0;nn<N;nn++) {
         for(mm=0;mm<N;mm++) {
             sum = 0;
             for(ip=0;ip<N;ip++) {
               sum = sum + cosarray[nn][ip] * cosarray[mm][ip] *  ewfac[ip];
             }
          Sq  = Sq  + exp(- fqq0*(q*q)*(       abs(nn-mm)*(l*l)/6.0));
          Sqt = Sqt + exp(- fqq* (q*q)*(Dr*t + abs(nn-mm)*(l*l)/6.0) + 
                fqq*ff2* sum);
         };
       };
//$OMP } // endPARALLEL DO;

       Sq  = Sq /N;
       Sqt = Sqt/N;
       Sq0t.sq0 = Sq;
       Sq0t.sqt = Sqt;

       return (Sq0t);

} // endfunction nrouse_ngcorr;

  sqt_type reptation_sqt( 
    double q         ,// the Q-value ;
    double t         ,// the time;
    double n         ,// total number of segments in chain;
    double lseg      ,// segment length;
    double ne        ,// average number of segments per entanglement strand;
    double Re        ,// entanglement "blob" size (R-end-end);
    double Wl4       ,// Rouse rate;
    double alpha0    ,// Non-Gaussianity correction amplitude;
    double talphamax ,//     "  distribution max-location;
    double talphawd   //     "        "      width;
  ) {

    sqt_type sqr0t;

    double teps = 1e-6;
    double W, Wr, lr;
    double tmax, twidth;
    sqt_type sq_rouse, sq_locrep;
    int Nrsum;

    double b=1;

    if(alpha0 == 0.0) {;
       tmax   = 1.0;
       twidth = 1.0;
    } else {;
       tmax   = talphamax;
       twidth = talphawd;
    };

    W         = Wl4 / (lseg*lseg*lseg*lseg);
    Nrsum     = (int)(ne);
    if(Nrsum > 300) Nrsum=300;

    sq_rouse = nrouse_ngcorr( q, t, Nrsum, Re, Wl4, alpha0, tmax, twidth) ;
   
    if(t < teps) t=teps;

    sq_locrep.sq0   =  local_reptationdr( q, teps          ,  lseg, W, n, ne);
    sq_locrep.sqt   =  local_reptationdr( q, t             ,  lseg, W, n, ne);

    sqr0t.sq0 = sq_rouse.sq0 * sq_locrep.sq0;
    sqr0t.sqt = sq_rouse.sqt * sq_locrep.sqt;

    return (sqr0t);

  } // endfunction reptation_sqt;





 
int  main( int argc , char **argv ) {

    double q          ; //  the Q-value ;
    double t          ; //  the time;
    double n          ; //  total number of segments in chain;
    double lseg       ; //  segment length;
    double ne         ; //  average number of segments per entanglement strand;
    double Re         ; //  entanglement "blob" size (R-end-end);
    double wl4        ; //  Rouse rate;
    double alpha0     ; //  Non-Gaussianity correction amplitude;
    double talphamax  ; //      "  distribution max-location;
    double talphawd   ; //      "        "      width;


    double teps = 1e-6;
    double W, Wr, lr;
    double tmax, twidth;
    sqt_type sqrep, sqrouse;
    int Nrsum, i;

    double b=1;
    double w;

    wl4    = 30000;
    Re     = 50.0;
    Nrsum  = 100 ;
    n      = 5000;
    ne     = 100;
    lseg   = 5.0 ;
    alpha0 = 0.30;
    tmax   = 1;
    twidth = 1;
    q      = 0.150;

    for(i=0;i<51;i++){
      t = 0.01 * pow(1.30,(double)i);
      sqrouse = nrouse_ngcorr( q, t, Nrsum, Re, wl4, alpha0, tmax, twidth) ;
      w  = wl4/(lseg*lseg*lseg*lseg);
      lr = local_reptationdr( q, t, lseg, w,  n ,ne );
  
      sqrep = reptation_sqt(q,t, n, lseg, ne, Re, wl4, alpha0, tmax,twidth);
      printf("%18.9f %18.9f %18.9f  %18.9g  %18.9g\n",t,sqrep.sqt/sqrep.sq0, lr, sqrouse.sq0, sqrouse.sqt);
    }

 

} // end main


 
