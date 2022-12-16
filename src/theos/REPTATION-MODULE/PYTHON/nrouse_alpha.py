import numpy as np
from scipy import special
import time


# local rpetation part
 
def local_reptationdr(q, t, lseg, W, n, ne):   
    Z    = n/ne
    T1 = 72.0 * (np.exp(-q**2 * n * lseg**2 / 6.0) + q**2 * n * lseg**2 / 6.0 - 1.0) / q**4 / lseg**4
    T1 = T1 / n**2  

    x  = np.sqrt(lseg**4*q**4*W*t/36)
    xa = n/(2*np.sqrt(W*t))
    
    if np.abs(x) < 15.0 : 
        dec = np.exp(x**2) *(special.erfc(x)-special.erfc(x+xa))
    else:
        dec =  1.0/(np.sqrt(np.pi)*x)-1.0/(2*np.sqrt(np.pi)*x**3)+3.0/(4*np.sqrt(np.pi)*x**5) \
          - (x**4-xa*x**3+(xa**2-1.0/2.0)*x**2+(-xa**3+3.0/2.0*xa)*x+xa**4-3*xa**2+3.0/4.0) \
           * np.exp(-xa**2)*np.exp(-2*xa*x)/(np.sqrt(np.pi)*x**5)
    

    T2  = ((np.exp(-1/6*lseg**2*n*q**2-n**2/(4*W*t))-1)*2/3*ne*np.sqrt(W*t/np.pi)+ \
        ne*(n/3+lseg**2*q**2*W*t/9) \
        * dec )

    T2  = T2*3.0 / (n*ne)   

    val = 1.0/(1+1/Z) *  (T1 + (2.0/3.0 + (1.0/3.0) * T2) / Z )

    return val

 

# Rouse S(q,t) and S(q) via direct summation (including Non-Gaussianity correction)

def nrouse_a(q, t, wl4, Re, N, alpha0, talphamax, talphawd): 



    def alpha_ng(t):
        a = alpha0 * np.exp(-(((np.log(t+1e-3)-np.log(talphamax))/talphawd)**2) / 2.0) 
        return a

     
    l = np.sqrt(Re**2/N)       
    Dr  = 0.0                
    W = wl4/l**4
          
    nr    = range(N)
    nn    = np.zeros(N)
    nn    = np.array(nr)+1
    cosarray  = np.zeros((N,N))
    for i in nr:
         Ip=i+1.0   
         cosarray[:,i]  = np.cos(nn*Ip*np.pi/N) / Ip
   
    ewfac = np.zeros(N)
    ewfac = (1.0-np.exp(-2.0*W*(1.0-np.cos((np.pi*nn)/N))*t))     
    ff2   = -2*N*(l*q)**2/(3*np.pi**2)
   
    rmm = 0    
    for mm in nr:
         rmm = rmm + 4.0*N*l**2/(np.pi**2)*np.sum(cosarray[mm,:]*cosarray[mm,:]*ewfac[:])    

    rmm = rmm/N
   
    fqq  = 1.0 -q**2 * rmm/12.0 * alpha_ng(t)    
    fqq0 = 1.0                 
    
# i guess most of the time is spnt in this double loop...
    Sq  = 0.0
    Sqt = 0.0
    for nnn in nr:
         for mmm in nr:
              Sq  = Sq  + np.exp(-fqq0*(q**2)*(np.abs(nnn-mmm)*(l**2)/6.0))
              Sqt = Sqt + np.exp(-fqq* (q**2)*(Dr*t + np.abs(nnn-mmm)*(l**2)/6.0) \
                                 +fqq*ff2*np.sum(cosarray[nnn,:]*cosarray[mmm,:]*ewfac[:]))

    Sq  = Sq /N
    Sqt = Sqt/N
    
    return (Sqt, Sq)




    
# the normaliised scattering function S(q,t)/S(q) for a long reptation chain
# n    = total number of segments in the long chain
# ne   = number of segments per (average) entanglement strand
# lseg = segment length
# re   = lateral extension of the tube (i.e. entanglement blob size)
# wl4  = the Rouse rate
# alpha0          = amplitude of non-Gaussianty (NG) correction fi the Rouse blon
# talphamax       = center of alpha(t) NG distribution
# talphawd        = width(log) of alpha(t) NG distribution

def reptation_fqt(q, t, n, ne, lseg, re, wl4, alpha0, talphamax, talphawd):

    Nro =  min(int(ne),300) 
    W   =  wl4/lseg**4
    sr  =  nrouse_a(q, t, wl4, re, Nro, alpha0, talphamax, talphawd) 
    lr0 =  local_reptationdr(q, 1e-6, lseg, W, n, ne)
    lr  =  local_reptationdr(q, t   , lseg, W, n, ne)
    
    return sr[0]*lr/(sr[1]*lr0)







####################################################################################################


seconds0 = time.time()


wl4    = 30000
Re     = 50.0
N      = 100  
alpha0     = 0.1
talphamax  = 1.0
talphawd   = 1.0
ne     = 100.0
W      = 10.0
lseg   = 5.0
n      = 5000.0

q      = 0.150
t      = 0.0

    
for it in range(21):
    t =  0.01 * 1.3**it
    s = reptation_fqt(q, t, n, ne, lseg, Re, wl4, alpha0, talphamax, talphawd)
    print(t,s)

seconds1 = time.time()

print(" ")
print("net execution time/sec:")
print(seconds1-seconds0)

### BUT:
### fortran (openmp, O2): 
# real	0m0.051s
# user	0m0.094s
# sys	0m0.026s

# THIS:
# real	0m7.810s
# user	0m7.832s
# sys	0m1.342s

 
