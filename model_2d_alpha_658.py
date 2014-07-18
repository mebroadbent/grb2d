# -*- coding: utf-8 -*-
## this script controls the 2d modeling
import numpy as np
import pickle
import os
import matplotlib as mpl
mpl.use('Agg')#
import matplotlib.pylab as pylab
import matplotlib.pyplot as plt
#import pandas as p
import batse5bp
import alpha_2d_simple as code
import scipy as sp
import scipy.special
import cProfile
from mb2d_plots_functions import plot_rj


grb_cat = batse5bp.catalog.load_catalog(root_dir='BATSE_5Bp_Data')
class bunch:
    def __init__(self, **kwds):
        self.__dict__.update(kwds)
import colorsys

import datetime
import pdb
#from plots_2d_simple import *
import string


from utils import *
#restore_rng()
GRBno = 658
directory = str(datetime.date.today())[-4::]+str(GRBno)+'_a'

while os.path.exists(directory):
   print directory, ' already in use'
   index = directory[-1]
   directory = directory[0:-1]
 
   directory = directory + string.lowercase[ord(index)-97+1]
print directory
if not os.path.exists(directory):
    os.mkdir(directory)
os.system('cp ' + __file__ + ' ' + directory)    


#directory = 'sim-data_'+ str(datetime.date.today())
grb = grb_cat[GRBno]
ascii = grb.ascii64

ascii.counts = np.random.binomial(ascii.counts, p = 1.)# round(ascii.counts*0.05)

B = [np.mean(ascii.counts[i][0:50])/1.024 for i in range(4)]


drms = grb.discsc_drms
drms.set_sum_prodquad() 
config = bunch(grb = grb,
            drms = drms,
            ascii = ascii,
            tmin = -5., #tmin and tmax indicate truncation
            tmax = 20., 
            M = 100,  #number of iterations
            B = B, #bsackground
            burn =0, # number of burnin iterations
            thin = 5, #thinning 
            eps = 1., #min pulse height
            prob = [0.15,0.15,0.7,0,0], #probs for birth, death, wal
            d_t = .15, #d_x is the MH proposal standard deviation
            d_A = 1.,
            d_ab = .25,
            d_lambda = 0.075,
            d_gamma = 0.075,
            d_Epeak = 75., 
            psi = 1./(1.),
            N =1, #starting # of pulses
            p_thin =[1.], #chain thinnings
            redraw = 10000000000, #how often to redraw thinned chains (you see that right now i am not redrawing ever..."
            swap = 100,          
            zeta = 1.,
            xsi = 0.1,
            overdisp = True,
            prior = False)

prior = bunch(alpha =1.5, #np.exp(0.4)/37,
              #beta = 0.001,#3/2.,
                    
             # gamma = 1/(37*37.), #Amplitudes
              scale_alpha = 1.,
              scale_beta = .2,
              mu_Epeak = 300.,#np.log(500.),
              s_Epeak = 100.,

              sigma_gamma = 2.2,# 0.000001,#0.3,#1.,#1.4,#prior variance of linear lag
             # shape_scale_sigma_E = (.5,2.), #parameters for sigma_E
              shape_scale_lambda_t = (2.,1.) #parameters for lambda
            )
            
log_c =  ( np.log( 1./np.pi) + sp.special.gammaln(prior.alpha) + np.log(np.sin(np.pi*prior.alpha /2) ) )
num_pulses = 3
gamma = num_pulses/(2*(config.tmax-config.tmin)*np.exp(log_c)*config.eps**(-prior.alpha) )
prior.gamma = gamma
log_gamTalph2 = np.log( 2*prior.gamma*(config.tmax-config.tmin)*prior.alpha)  

prior.log_gamTalpha2 =  log_gamTalph2
prior.log_c = log_c

starting_value = np.array([
 [ -4.85727046e-01,   1.03108509e+00,  -4.21827918e-01,   1.14380161e+00, -1.87134646e-01,  -2.00000000e+00,   8.20594236e+01],
 [  1.97586223e+00,   3.06228855e+00,   7.53773189e-01,   7.26595812e-01, -6.51899120e-01,  -2.00000000e+00,   6.36506118e+01]] )

OUT = code.grb_mcmc(prior,config,starting_value)
os.chdir(directory)
   
myfile = open('datafile','w')
pickle.dump(OUT,myfile)
myfile.close()

import scipy.stats.distributions as dist
J = [len(th) for th in OUT[0][-1] ] 
plt.hist(J, normed = True)
plt.plot(range(max(J)),dist.poisson.pmf(range(max(J)),3) )
plt.show()


TH,A,LLH,Xt,Yt,swap, params_accepted,r_thinned,STEPS,CHISQ = OUT
print STEPS[0][0]
plot_rj(TH[-1], config, prior, Xt,Yt[-1],None, None,LLH, CHISQ,r_thinned,p_thin = 1,weights = None)
#plot_norj(TH[-1],params_band,'.',config,prior,Xt,Yt[-1],STEPS[0],r_thinned)