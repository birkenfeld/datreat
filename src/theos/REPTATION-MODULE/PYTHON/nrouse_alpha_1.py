import numpy as np
import time








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
    
    



####################################################################################################


seconds0 = time.time()


wl4    = 30000
Re     = 50.0
N      = 100  
alpha0     = 0.1
talphamax  = 1.0
talphawd   = 1.0


q      = 0.150
t      = 0.0

    
for it in range(51):
    t = 1.1**it
    sx =  nrouse_a(q, t, wl4, Re, N, alpha0, talphamax, talphawd) 
    print(t,sx)

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

 
