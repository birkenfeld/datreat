#!/usr/bin/env python

"""

==================================
NSE scattering evaluator for a Rouse polymer chain
==================================

-----------------------------------------------------------------------
Original author:
M. Monkenbusch
Juelich Center of Neutron Science
Forschungszentrum Juelich
m.monkenbusch@fz-juelich.de

Fortran to Python translation:
Piotr Zolnierczuk
zolnierczukp@ornl.gov

-----------------------------------------------------------------------
This is a little toolbox to estimate the effect of contrast and
H/D content on the dynamics data of a polymer melt in the Rouse
regime as obtained by NSE spectroscopy.

It is not intended to scrutinize Rouse behaviour and deviations from that,
for that the time functions should be computed with the methods
given e.g. in the book of Doi-Edwards.

But it will give a first impression on scattering contributions etc.

In particular:
consideration concerning the roles of coherent an incoherent scattering
in polymer samples

please observe:
-----------------------------------------------------------
units for
Q                                        --> 1/Angstroem
lambda                               --> Angstroem

scattering length                 --> cm
scattering length density    --> cm**-2
scattering cross section      --> cm**2
neutron flux                       --> neutrons/cm**2/s

molecular weights              --> g/mol
density                              --> g/cm**3

sample dimensions:

thichkness,d                      --> cm
area                                   --> cm**2

Rouse rate parameter

Wl4                                    --> Angstroem**4/nanoseconds
t                                      --> nanoseconds
------------------------------------------------------------
"""

import numpy as np
from numpy import pi, exp, cos, expm1, sqrt, mgrid, where
import time
#pylint: disable=no-name-in-module
from scipy.special   import erfc
from scipy.integrate import quad
from scipy.constants import Boltzmann as k_B


_debye_cut  = 1e-5  # cut-off value for f_debye
_rouse_cut  = 1e-33 # cut-off value for skernel_rouse
_rouse_eps  = 1e-6  # relative eps for Rouse kernel integration
def f_debye(x):
    """Debye function f(y)
    where y=x^2=(kR)^2"""
    y = x**2
    y = where(y<_debye_cut,_debye_cut,y)
    return 2/y**2*(expm1(-y)+y)  #pylint: disable=invalid-unary-operand-type


def rouse_parameters(Wl4, Re, N, T=300.0):
    """
    Calculate Rouse parameters using Wl4, Re and N
    """
    l    = Re/sqrt(N)          # segment length
    kT   = k_B * T * 100.0     # kg*m^2/s^2 -> kg * A^2/ns^2
    zeta = 3*kT*l**2/Wl4       # friction factor
    W    = 3*kT/(zeta*l**2)    # Rouse rate in 1/ns
    Dr   = kT/(N*zeta)         # diffusion constant
    tR   = (N/pi)**2/W         # Rouse time
    Rg   = Re/sqrt(6)          # radius of gyration
    return dict(W=W, l=l, zeta_0=zeta, Rg=Rg, D=Dr, tau=tR,
                Wl4=Wl4, Re=Re, N=N, T=T)

def rouse_sqt(t, q, W, Re, Dr=0, N=100, **kwargs):
    """
    Rouse expression for a chain of finite length:
    args::
        t    - time in ns
        q    - momentum transfer in A^(-1)
        W    - Rouse rate in 1/ns
        Re   - end-to-end distance of the polymer molecule
        Dr   - center of mass diffusion constant in A**2/ns
        N    - number of chain segments
        pmin - minimum p (default 1)
        pmax - maximum p (default N)
    returns:
        Sq   - S(Q)
        Sqt  - S(Q,t)

    Based on a Fortran subroutine nrousep by M. Monkenbusch (JCNS)
    Fortran -> Python translation by P. Zolnierczuk (JCNS)
    """
    if N <= 0:
        raise RuntimeError('Number of chain segments must be positive')
    pmin = max(kwargs.get('pmin', 1), 1)
    pmax = min(kwargs.get('pmax', N), N) + 1

    l   = Re/sqrt(N)       # calculate segment length
    fac = 2/3*(Re*q/pi)**2 #

    # ---- Do the sums -----
    PP, NN, MM  = mgrid[pmin:pmax, 0:N, 0:N]
    # exponent for S(Q) and S(Q,t)
    anm  = exp(-(q**2*abs(NN[0]-MM[0])/6*l**2))
    # exponent for S(Q,t)
    cosf = exp(-fac* np.sum(cos(pi*(NN+1)*PP/N)*cos(pi*(MM+1)*PP/N)/PP**2
                           *(1-exp(-2*W*(1-cos(pi*PP/N))*t)), axis=(0,))
                            )
#    print(cosf)
    Sq  = np.sum(anm)
    Sqt = exp(-q**2*Dr*t) * np.sum(anm*cosf)

    return Sq/N, Sqt/N

# ###################################################################################
# WARNING THIS IS OBSOLETE APPROXIMATION FOR QR>>1
# Rouse kernel $h(u) = \frac{2}{\pi} \int_0^\infty \cos(ux)/x^2 (1-e^{-x^2}) dx$
#              $h(u) = \frac{2}{\sqrt{pi}} \exp{-(u/2)^2} - u erfc(u/2) $
def h_rouse_hiq(u):
    "Rouse kernel h(u)"
    u2 = u/2
    return 2/sqrt(pi)*exp(-(u2)**2) - u*erfc(u2)

# exp(-(u+sqrt(omegat)*h_rouse(u/sqrt(omegat))))
def skernel_rouse_hiq(u, xot):
    "Rouse S(Q,t) inegrand"
    #xot  = sqrt(omegat)
    if xot<_rouse_cut:
        return exp(-u)  # h(oo) = 0
    return exp(-(u+xot*h_rouse_hiq(u/xot)))

#  $ S(Q,t)/S(Q,0) = \int_0^\infty dx exp(-u + \sqrt(omegat)
def sqt_rouse_hiq(xomegat):
    "Rouse dynamic structure factor"
    if np.isscalar(xomegat):
        return quad(skernel_rouse_hiq, 0, np.inf, epsrel=_rouse_eps,args=(xomegat,))
    res = [    quad(skernel_rouse_hiq, 0, np.inf, epsrel=_rouse_eps,args=(_xot,   ))
               for _xot in xomegat ]
    return np.asarray(res)


def test():
    # results from datreat
    #                  tau      S(Q)      S(Q,t)
    dt = np.asarray( [[ 0.00000, 48.835216, 48.835216],
                      [10.00000, 48.835216, 42.918585],
                      [20.00000, 48.835216, 42.702855],
                      [30.00000, 48.835216, 42.692559],
                      [40.00000, 48.835216, 42.692053],
                      [60.00000, 48.835216, 42.692025],
                      [70.00000, 48.835216, 42.692025],
                      [80.00000, 48.835216, 42.692025],
                      [90.00000, 48.835216, 42.692025],
                      [100.00000, 48.835216, 42.692025],
                      [110.00000, 48.835216, 42.692025],
                      [120.00000, 48.835216, 42.692025],
                      [130.00000, 48.835216, 42.692025],
                      [140.00000, 48.835216, 42.692025],
                      [150.00000, 48.835216, 42.692025],
                      [160.00000, 48.835216, 42.692025],
                      [170.00000, 48.835216, 42.692025],
                      [180.00000, 48.835216, 42.692025],
                      [190.00000, 48.835216, 42.692025],
                      [200.00000, 48.835216, 42.692025] ])
    Wl4 = 39500.0
    Re  =    40.0
    T   =   300.0
    N   =   100
    seconds0 = time.time()
    rouse_pars =rouse_parameters(Wl4=Wl4, Re=Re, T=T, N=N)
    W = rouse_pars.get('W')
    q = 0.1
    for tau, Sq_dt, Sqt_dt in dt:
        Sq, Sqt = rouse_sqt(tau, q , W, Re, N=N, pmin=1, pmax=N)
        print("q=%.1fA^{-1} tau=%4.1fns S(q)=%.3f/%.3f S(q,t)=%.3f/%.3f" % (q, tau, Sq, Sq_dt, Sqt, Sqt_dt))
    seconds1 = time.time()

    print(" ")
    print("net execution time/sec:")
    print(seconds1-seconds0)


if __name__ == "__main__":
    test()
