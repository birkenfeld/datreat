       implicit none
      
       double precision pi
       parameter (pi=3.141592654d0)
  
       double precision kb
       parameter (kb=1.380662d-18) ! in g A**2 / ns**2
      

       double precision p,L,zeta,l0,Cinf,Temp,q,tim,diffsca,gamx
       common /cstiffz1/p,L,zeta,l0,Cinf,Temp,q,tim,diffsca,gamx
 
       integer nsegment
       common /cstiffz1/nsegment

       integer maxn
       parameter(maxn=40)
       double precision zeros1(maxn), zeros2(maxn)
       double precision exzeros1(maxn), exzeros2(maxn)
       integer          infozer1(maxn), infozer2(maxn)
       integer iz1, iz2, ip,iq, maxmode

       common /cszero/exzeros1,exzeros2,iz1,iz2,ip,iq,maxmode

       double precision norm1(maxn),norm2(maxn)
       common /cnrm/norm1,norm2


