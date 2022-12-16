import numpy as np

a=np.array([11,22,33,44,55])
print(a)
b=range(5)
print(b)
range(0, 5)
for i in b:
     print(i)

def alpha_ng(t):
    a = alpha0 * np.exp(-(((np.log(t+1e-3)-np.log(talphamax))/talphawd)**2) / 2.0) 
    return a
    


#       do nn=1,N
#        do ip=1,N
#         cosarray(nn,ip) = cos((pi*ip*nn)/dfloat(N)) / ip 
#        enddo
#       enddo
Ip       = 2.0
N        = 10
W        = 1.0
t        = 1
l        = 1
q        = 0.1
alpha    = 0.0
Dr       = 0.0

alpha0     = 0
talphamax  = 1.0
talphawd   = 1.0




nr        = range(N)
nn    = np.zeros(N)
nn    = np.array(nr)+1
cosarray  = np.zeros((N,N))
for i in nr:
     Ip=i+1.0   
     cosarray[:,i]  = np.cos(nn*Ip*np.pi/N) / Ip
#     cospinn[:,i]  = nn*Ip

print(nn)
print(cosarray)


ewfac = np.zeros(N)
ewfac = (1.0-np.exp(-2.0*W*(1.0-np.cos((np.pi*nn)/N))*t)) 
      

print(np.pi*nn/N)
print(np.cos((np.pi*nn)/N))
print(1.0-np.cos((np.pi*nn)/N))
print(ewfac) 

ff2  = -2*N*(l*q)**2/(3*np.pi**2)

print(ff2)



rmm = 0

for mm in nr:
     rmm = rmm + 4.0*N*l**2/(np.pi**2)*np.sum(cosarray[mm,:]*cosarray[mm,:]*ewfac[:])

rmm = rmm/N

print(rmm)

#fqq  = 1.0 -q**2 * rmm/12.0 * alpha     
fqq  = 1.0 -q**2 * rmm/12.0 * alpha_ng(t)    
fqq0 = 1.0                 


# !$OMP PARALLEL DO REDUCTION(+:Sq,Sqt)
#        do nn = 1,N
#         do mm = 1,N
# 
#           Sq  = Sq  + exp(- fqq0*(q**2)*(       abs(nn-mm)*(l**2)/6.0d0))
#           Sqt = Sqt + exp(- fqq* (q**2)*(Dr*t + abs(nn-mm)*(l**2)/6.0d0) + &
#                 fqq*ff2* sum(cosarray(nn,ipmin:ipmax) * cosarray(mm,ipmin:ipmax) *  ewfac(ipmin:ipmax) ))
# 
#         enddo
#        enddo
# !$OMP END PARALLEL DO
Sq  = 0.0
Sqt = 0.0
for nnn in nr:
     for mmm in nr:
          Sq  = Sq  + np.exp(-fqq0*(q**2)*(np.abs(nnn-mmm)*(l**2)/6.0))
          Sqt = Sqt + np.exp(-fqq* (q**2)*(Dr*t + np.abs(nnn-mmm)*(l**2)/6.0)+fqq*ff2*np.sum(cosarray[nnn,:]*cosarray[mm,:]*ewfac[:]))



Sq  = Sq /N
Sqt = Sqt/N

print(Sq)
print(Sqt)



# cosarray = np.zeroes((N,N))
 
