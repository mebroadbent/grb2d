import numpy as np
cimport numpy as np

import mpmath

from cython.view cimport array as cvarray

import cython
cimport cython

cdef extern from "math.h":
    double exp(double)
    double log(double)
    double sqrt(double)
    double erf(double)


@cython.cdivision(True)
def k2d(np.ndarray[double, ndim = 1] times not None, np.ndarray[double, ndim = 1] energies not None, np.ndarray[double, ndim = 2]  theta not None):

#double t0, double A, double gamma, double m_E, double s_E, double lam):

    cdef int J = theta.shape[0]
    cdef int L = times.shape[0]-1
    
    cdef int M =energies.shape[0]
    
    cdef np.ndarray[double,ndim = 2]y = np.zeros((L, M))
    
    cdef float time1
    cdef float time2
    cdef float temp
    cdef float t
    cdef float e
    
    cdef float t0
    cdef float A
    cdef float gamma
    cdef float lam    
    cdef float alpha
    cdef float beta
    
    cdef int i
    cdef int iplus1
    cdef int k
    cdef int j
    cdef int current_loc
    
    cdef double current_val
    cdef float Ec
    
    cdef float constant 
    for k in range(J):
        
        t0  = theta[k][0]
        A  = theta[k][1]
        gamma = theta[k][2]
        lam = theta[k][3]
        alpha = theta[k][4]
        beta = theta[k][5]
        Ec= theta[k][6]

        
        
        
       
        constant = A 
        for i in range(L):
            iplus1 = i + 1
            time1 = times[i]
            time2 = times[iplus1]
            
            time_int =exp(-lam*time1) - exp(-lam*time2)    #exp(-lam*(time2+time1)/2)*(time2-time1)#
            t = 0.5*(time1 + time2)
            for j in range(M):
            
                e = energies[j]
                energypart = ((e/100.)**alpha) * exp(-e/Ec)

                if t > (t0-gamma*(log(e/100) )) :
                      time_part = exp(lam*(t0 - gamma*(log(e) -log(100)) ))
                      current_val = time_int*energypart*time_part*constant
                      y[i,j] += current_val         
                 

    return y