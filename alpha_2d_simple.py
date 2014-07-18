# -*- coding: utf-8 -*-
import numpy as np
import batse5bp
import scipy as sp
import scipy.stats
import scipy.misc as misc
import os
import pickle
import random

from scipy.interpolate import interp1d
import sys
import scipy.stats.distributions as dist
import scipy.special 
import pdb
import collections
import functools
##########
## I re-compile the kernel function using setup.py each time I load this module
sys.argv=['','build_ext','--inplace']
execfile('setup.py')
from k_2d  import k2d as K

batse5bp.catalog.load_catalog(root_dir='BATSE_5Bp_Data')


## Memoization for the likelihood; the cache is a bounded-ordered-dictionary, which 
#  has a limit on the number of items and discards the oldest ones. 
class BoundedOrderedDict(collections.OrderedDict):
    def __init__(self, *args, **kwds):
        self.maxlen = kwds.pop("maxlen", None)
        collections.OrderedDict.__init__(self, *args, **kwds)
        self._checklen()

    def __setitem__(self, key, value):
        collections.OrderedDict.__setitem__(self, key, value)
        self._checklen()

    def _checklen(self):
        if self.maxlen is not None:
            while len(self) > self.maxlen:
                self.popitem(last=False)

def memoize(func=None, maxlen=None):
    if func:
        cache = BoundedOrderedDict(maxlen=maxlen)
        @functools.wraps(func)
        def memo_target(*args):
            lookup_value = args
            if lookup_value not in cache:
                cache[lookup_value] = func(*args)
            return cache[lookup_value]
        return memo_target
    else:
        def memoize_factory(func):
            return memoize(func, maxlen=maxlen)
        return memoize_factory

##############################                
#
@memoize(maxlen = 100000)
def mbs_quad_nodes_wts(l,u,config):
    '''this is a wrapper for Tom's function drms.quad_nodes_wts, 
      so I can memoize it.
      l, u: floats, bounds for the integration
      config: bunch object, contains the drms object for the current grb_mcmc
      
      returns:
      a: nodes (list of floats)
      b: weights (list of floats)
      
      '''
    a,b = config.drms.quad_nodes_wts(l,u)
    return a, b
    
## transdim kernel steps
def birth(th,config,prior,Xt,Yt,llh_old,p_thin, r, m = 0):
    '''
    conducts a birth step (draws new pulse, adds it, and does MH)
    th: np array with J rows for the J pulses
    config, prior: dictionaries of parameters
    Xt: list of times (truncated according to the parameters in config)
    Yt: List of nchan (typically 4) lists of photon counts.
    llh_old: scalar, loglikelihood of th
    p_thin: scalar, thinning parameter
    r: scalar, current value of the negative binomial parameter.
    m: mcmc iteration, not used currently
    
    
    RETURNS:
    th or th_star, depending on if move was accepted
    a list [0,0,0,0,0] or [1,0,0,0,0], to indicate whether or not a walk move was accepted
    llh_old or llh_star, depending on if move was accepted
    stepsizes: None 
    chisq: None if there has been no change in chisquare statistic; or the new value if the move was accepted.
    '''
    
    
    J = len(th)

    th_star = np.copy(th)
    
    new_row = draw_new_state(config,prior)
    

    if J > 0:
        th_star = np.vstack([th,[ new_row ]] )
    else: 
        th_star = np.array([new_row])
    llh_star = 0
    chisq = None
    if config.prior == False:
        llh_star,trash, chisq = llh(th_star,config,prior,Xt,Yt,r,p_thin = p_thin)

    log_trans_bwd=np.log(config.prob[1]/(J+1))
    log_trans_fwd=log_b(new_row,config,prior)+np.log(config.prob[0]/(J+1))
    
    gamma_proc = prior.log_gamTalpha2 + prior.log_c - np.log( J +1) -(prior.alpha + 1)*np.log(new_row[1])

    LH = (gamma_proc#- prior.beta*new_row[1]  + np.log(prior.alpha)  - np.log( J+1) - np.log(new_row[1])
          + log_trans_bwd
          -log_trans_fwd
          + eval_prior(new_row,config,prior)
          +llh_star - llh_old)
   
    if np.random.exponential(1) + LH > 0 :
    
        return [th_star,[1,0, 0,0,0], llh_star,None,chisq]        
    else: return [th,[0,0,0,0,0], llh_old,None,None]    
    
def death(th,config,prior,Xt,Yt,llh_old,p_thin, r, m = 0): 
    '''
    conducts a death step (proposes th_star by removing one pulse from th, and does MH)
    th: np array with J rows for the J pulses
    config, prior: dictionaries of parameters
    Xt: list of times (truncated according to the parameters in config)
    Yt: List of nchan (typically 4) lists of photon counts.
    llh_old: scalar, loglikelihood of th
    p_thin: scalar, thinning parameter
    r: scalar, current value of the negative binomial parameter.
    m: mcmc iteration, not used currently 
    
    RETURNS
    th or th_star, depending on if move was accepted
    a list [0,0,0,0,0] or [0,1,0,0,0], to indicate whether or not a walk move was accepted
    llh_old or llh_star (depending on if move was accepted
    stepsizes: None 
    chisq: None if there has been no change in chisquare statistic; or the new value if the move was accepted.
    '''
    
    J = len(th)
    if J>0:
        J_star     = np.random.randint(0,J)
        th_jstar   = th[J_star]

        log_q_forward = np.log(config.prob[1]/J)
        log_q_backwards = log_b(th_jstar,config,prior) + np.log(config.prob[0]/J)
        
        
        th_star=np.delete(np.copy(th),J_star,axis=0)
        
        llh_star = 0
        chisq = None
        if config.prior == False:
            llh_star,trash, chisq = llh(th_star,config,prior,Xt,Yt,r,p_thin = p_thin)
         

        log_prior = (np.log(J) + (prior.alpha + 1)*np.log(th_jstar[1]) - prior.log_gamTalpha2 - prior.log_c
                     -eval_prior(th_jstar,config,prior)
                     )
                                  
        LH = (log_prior- log_q_forward
              +log_q_backwards  + llh_star - llh_old
              )    
         
        if (LH+np.random.exponential(1)>0 ):
            
            return [th_star,[0,1,0,0,0], llh_star,None,chisq]   
  
    return(th,[0,0,0,0,0],llh_old,None,None)


def walk(th,config, prior,Xt,Yt, llh_old, p_thin,r,m =0 ):
    '''
    conducts a walk step (proposes th_star by wiggling one pulse from th, and does MH)
    th: np array with J rows for the J pulses
    config, prior: dictionaries of parameters
    Xt: list of times (truncated according to the parameters in config)
    Yt: List of nchan (typically 4) lists of photon counts.
    llh_old: scalar, loglikelihood of th
    p_thin: scalar, thinning parameter
    r: scalar, current value of the negative binomial parameter.
    m: mcmc iteration, not used currently
    
    RETURNS:
    th or th_star, depending on if move was accepted
    a list [0,0,0,0,0] or [0,0,1,0,0], to indicate whether or not a walk move was accepted
    llh_old or llh_star (depending on if move was accepted
     stepsizes: list of 2 elements: first element: seven innovations. second element: whether or not accepted
    chisq: None if there has been no change in chisquare statistic; or the new value if the move was accepted.
   '''
    
    J = len(th)
    which_params = [0 for i in range(7)]
    stepsize = [None for i in range(7)]
    
    if J>0:
        J_star = np.random.randint(0,J)
        th_star = np.copy(th)
        row = th[J_star]
        A_old = row[1]
        step_adj = A_old**0.5
        ##A  = np.sort(th[:,1])
        A_star = np.random.normal(row[1],config.d_A/step_adj)
        
        if A_star <config.eps:
            A_star = config.eps+(config.eps-A_star)
            
        
        t_star =  np.random.normal(row[0],config.d_t/step_adj)
        g_star =  np.random.normal(row[2],config.d_gamma/step_adj)
        lam_star = np.random.normal(row[3],config.d_lambda/step_adj)
        
        if lam_star < 0.1:
            lam_star = 2*0.1 -lam_star
       
        b_star = row[5]#np.random.normal(row[5],config.d_ab/step_adj)
        a_star=  -abs(np.random.normal(row[4],config.d_ab/step_adj))
        
        if b_star >0:
            b_star = - b_star
        if a_star >0:
            a_star = - a_star
        E_star =np.abs( np.random.normal(row[6],config.d_Epeak/step_adj) )
        how_many_params = np.random.binomial(7, 1.)
        
        if E_star < 0:
            E_star = -E_star
        
        param_index = random.sample(range(0,7),how_many_params)
      
        new_row = [t_star,
                            A_star,
                            g_star ,
                            lam_star,
                            a_star,
                            b_star,
                            E_star]

        th_star[J_star] = np.copy(new_row)
        if  th_star[J_star][0] < config.tmin:
            th_star[J_star][0] = 2*config.tmin -th_star[J_star][0]
        if th_star[J_star][0] > config.tmax:
            th_star[J_star][0] = 2*config.tmax - th_star[J_star][0]

        A_old = row[1]
        # t, A, gamma, mu_e, sigma_e, lambda
        log_prior_A = (prior.alpha+1)*(np.log(A_old) - np.log(A_star))
        
        log_prior_other_new = eval_prior(new_row,config,prior)
        log_prior_other_old  =eval_prior(row,config,prior)
      
        log_trans =( np.log( (A_star/A_old)**(6/2) ) 
                   -(1./(2*config.d_A**2) * (A_star - A_old)**3) 
                   -(1./(2*config.d_t**2) *  (A_star - A_old) *(t_star - row[0])**2) 
                 -(1./(2*config.d_gamma**2) *  (A_star - A_old) *(g_star - row[2])**2) 
                  -(1./(2*config.d_lambda**2) *  (A_star - A_old) *(lam_star - row[3])**2) 
                   -(1./(2*config.d_ab**2) *  (A_star - A_old) *(a_star - row[4])**2) 
                   -(1./(2*config.d_Epeak**2) *  (A_star - A_old) *(E_star - row[6])**2) 
                  )

        llh_new = 0
        chisq = None
        if config.prior == False:
            llh_new,trash, chisq = llh(th_star,config,prior,Xt,Yt,r,p_thin = p_thin)
        log_hastings = log_prior_A +  log_trans+llh_new - llh_old + log_prior_other_new - log_prior_other_old
        stepsize = [(new_row- row,0)] #not accepted

        which_params = [ 1 if ind in param_index else 0 for ind in range(7)]    
        if np.random.exponential(1) +  log_hastings > 0 :
            stepsize = [(new_row- row,1)]
            return [np.copy(th_star),[0,0,1,0,0 ], llh_new,stepsize,chisq]
            
        else: return [np.copy(th),[0,0,0,0,0], llh_old,stepsize,None]   
    else: return [np.copy(th),[0,0,0,0,0], llh_old,stepsize,None]     


def split_row(row,config):
    ''' This function splits "row" into two rows, "row1" and "row2". each row represents a pulse
    row: list of 7 floats.
    config: bunch object of parameters
    
    returns: 
    row1, row2: each a list of seven parameters: T, A, gamma, lam, alpha, beta = 2, E_c
    log_fwd_split: float. gives the likelihood of drawing these degree of freedom parameters
    reject: boolean. 1 means that some of the parameters of row1 or row2 are not allowed, and the move must be rejected
    '''

    T, A, gamma, lam,alpha, beta, E = row

    dof     = np.random.beta(config.zeta,config.zeta,size=1)
    delta   = np.random.normal(0,config.xsi,size=6)

    A1 = dof[0]*A
    A2 = A - A1
    
    T1      = (A2*delta[0]+A*T)/A
    T2      =  T1-delta[0]

    gamma1 = (A2*delta[1]+A*gamma)/A
    gamma2 = gamma1 - delta[1]
    
    lam1 = (A2*delta[2]+ A*lam)/A
    lam2 = lam1 - delta[2]
    
    alpha1 = (A2*delta[3] +A*alpha)/A
    alpha2 = alpha1 - delta[3]
    
    beta1 = -2.#(A2*delta[4]+A*beta)/A
    beta2 = -2.#beta1 - delta[4]
    
    E1 = (A2*delta[5] +A*E)/A
    E2 = E1 - delta[5]
    
    reject = 0
    if min(E1,E2)<0:
        reject = 1
    if max(alpha1,alpha2)>0:
        reject = 1
    if max(beta1,beta2)>0:
        reject = 1
    if min(lam1,lam2) <0:
        reject = 1
    if min(A1,A2)< config.eps:
        reject = 1
    if min(lam1,lam2) < 0.1:
        reject = 1
    
    
    row1 = [T1,A1,gamma1,lam1,alpha1,beta1,E1]
    row2 = [T2,A2,gamma2,lam2,alpha2,beta2,E2]
    
    log_fwd_split = np.log(np.sum(dist.norm.pdf(delta,0,config.xsi))) + np.log(np.sum(dist.beta.pdf(dof,config.zeta,config.zeta)))
    
    #if reject == 1:
     #   log_fwd_split = -np.inf
    
    return [row1,row2,log_fwd_split,reject]

def merge_row(row1,row2,config):
    ''' This function performs the merging of two rows into one.
    row1, row2: each a list of seven parameters: T, A, gamma, lam, alpha, beta = 2, E_c
    config: gives the beta probabilities for degree of freedom calculation 
   
   
    returns:
    row_star: merged row. list of 7 floats
    log_bwd_merge. probability that row_star would be split into row1 and row2. float
    
    '''
   
    T1, A1, gamma1, lam1,alpha1, beta1, E1 = row1
    T2, A2, gamma2, lam2,alpha2, beta2, E2 = row2
    
    A_star = A1 + A2
    T_star = (A1*T1 + A2*T2)/A_star
    gamma_star = (A1*gamma1 + A2*gamma2)/A_star
    lam_star = (A1*lam1 + A2*lam2)/A_star
    alpha_star = (A1*alpha1 + A2*alpha2)/A_star
    beta_star =-2.# (A1*beta1 + A2*beta2)/A_star
    E_star = (A1*E1 + A2*E2)/A_star

    row_star = [T_star, A_star, gamma_star,lam_star,alpha_star, beta_star,E_star]
    
    delta = [T1 - T2, 
            gamma1 - gamma2,
            lam1 - lam2,
            alpha1 - alpha2,
            beta1 - beta2,
            E1 - E2]
            
    dof = A1/(A1+A2) 

    log_bwd_merge = (np.log(np.sum(dist.norm.pdf(delta,0,config.xsi))) 
                     + np.log(np.sum(dist.beta.pdf(dof,config.zeta,config.zeta))))
    
    return [row_star,log_bwd_merge]
    

def split(th,config, prior,Xt,Yt, llh_old, p_thin,r,m =0 ):
    '''
    conducts a split step 
    th: np array with J rows for the J pulses
    config, prior: dictionaries of parameters
    Xt: list of times (truncated according to the parameters in config)
    Yt: List of nchan (typically 4) lists of photon counts.
    llh_old: scalar, loglikelihood of th
    p_thin: scalar, thinning parameter
    r: scalar, current value of the negative binomial parameter.
    RETURNS:
    th or th_star, depending on if move was accepted
    a list [0,0,0,0,0] or [0,0,0,1,0], to indicate whether or not a walk move was accepted
    llh_old or llh_star (depending on if move was accepted
     stepsizes: None
    chisq: None if there has been no change in chisquare statistic; or the new value if the move was accepted.
   '''

    J = len(th)
    if J > 0:
        J_star = np.random.randint(0,J)
        th_old = th[J_star,]
        
        row1, row2, log_fwd_split,reject = split_row(th_old,config)
        theta_hold = np.delete(np.copy(th),J_star, axis = 0)
        theta_star = np.vstack([np.copy(theta_hold),row1,row2] )
        
        llh_star = 0
        chisq = None
        if config.prior == False:
            llh_star,trash, chisq = llh(theta_star,config,prior,Xt,Yt,r,p_thin = p_thin)
        
        prior_old  = eval_prior(th_old,config,prior) 
        prior_new  = eval_prior(row1,config,prior) + eval_prior(row2,config,prior)
    
        
        prior_gamma_proc = (prior.log_gamTalpha2+prior.log_c 
                            - np.log(J+1) + (prior.alpha+1)*(np.log(th_old[1]) - np.log(row1[1]) - np.log(row1[1])) )


        log_trans_bwd=(np.log(config.prob[3]) #should we split?
                      -np.log((J+1)*J/2) #where to put new pulses
                      -np.log(J)#which to split
                      )
                          
        log_trans_fwd=(log_fwd_split
                       +np.log(config.prob[4]/J) #should we merge
                       -np.log( J*(J+1)/2 ) #have to put the new peaks somewhere
                       )

        LH = llh_star - llh_old + prior_new - prior_old +prior_gamma_proc + log_trans_bwd - log_trans_fwd
        if reject !=1:
            if LH + np.random.exponential(1) > 0:
              
                return [np.copy(theta_star),[0,0,0,1,0 ], llh_star,None, chisq]
    return [np.copy(th),[0,0,0,0,0], llh_old,None, None]
            

def merge(th,config, prior,Xt,Yt, llh_old, p_thin,r,m =0 ):  
    '''
    conducts a merge step 
    th: np array with J rows for the J pulses
    config, prior: dictionaries of parameters
    Xt: list of times (truncated according to the parameters in config)
    Yt: List of nchan (typically 4) lists of photon counts.
    llh_old: scalar, loglikelihood of th
    p_thin: scalar, thinning parameter
    r: scalar, current value of the negative binomial parameter.
    RETURNS:
    th or th_star, depending on if move was accepted
    a list [0,0,0,0,0] or [0,0,0,1,0], to indicate whether or not a walk move was accepted
    llh_old or llh_star (depending on if move was accepted
     stepsizes: None
    chisq: None if there has been no change in chisquare statistic; or the new value if the move was accepted.
   '''


    J = len(th)

    if J > 1:
        mylist=range(0,J)
        np.random.shuffle(mylist)
        J_star=mylist[0:2]
        row1 = th[J_star[0]]
        row2 = th[J_star[1]]
        
        row_star,  log_bwd_merge = merge_row(row1,row2,config)

        theta_temp=np.delete(np.copy(th),J_star,axis=0)
        theta_star=np.vstack([np.copy(theta_temp), row_star] )
        del theta_temp
        
        llh_star = 0
        chisq = None
        if config.prior == False:
            llh_star,trash, chisq = llh(theta_star,config,prior,Xt,Yt,r,p_thin = p_thin)
        prior_old  = eval_prior(row1,config,prior) + eval_prior(row2,config,prior)
        prior_new  = eval_prior(row_star,config,prior)
        

        
        prior_gamma_proc = (np.log(J) - prior.log_gamTalpha2 -prior.log_c + (prior.alpha + 1)*(np.log(row1[1]) + np.log(row2[1]) - np.log(row_star[1])))#np.log((J*row1[1]*row2[2] ) /(prior.alpha*row_star[1]) ) 
        
        
        log_trans_fwd = np.log(2*config.prob[3]) -np.log(J*(J-1) )-np.log(J-1)
        log_trans_bwd = (np.log(config.prob[4]) #do we split
                        +log_bwd_merge 
                        -np.log( (J-1)*J/2) #where do we put new ones
                        -np.log(J-1) #which do we split

                        ) 
        LH = llh_star - llh_old + log_trans_bwd - log_trans_fwd + prior_gamma_proc + prior_new - prior_old
        if LH + np.random.exponential(1) > 0:
            
            return [np.copy(theta_star),[0,0,0,0,1], llh_star,None, chisq]
    return [np.copy(th),[0,0,0,0,0], llh_old,None, None]

def E1(x):
    return dist.gamma.logsf(x,0.00000001,0,1) - np.log(0.00000001)
    
def starting_val(config, prior):
    ## determines a starting value for TH
    starting = []
    for j in range(config.N):
        starting.append(draw_new_state(config,prior))
    starting = np.array(starting)
    
    return starting

    
def draw_new_state(config,prior):
    ## draws a new pulse
    ## returns the new pulse and the log of its probability of being drawn
    T = config.tmax-config.tmin

    T_star   =T*np.random.beta(1,1)+config.tmin
    ups_star   =config.eps + np.random.exponential(scale=1./config.psi)
    gamma =np.random.normal(0,prior.sigma_gamma)
    lam =np.random.gamma(prior.shape_scale_lambda_t[0],prior.shape_scale_lambda_t[1])
    alpha =-abs( np.random.normal(-1,prior.scale_alpha) )
    beta = - abs(2.) #np.random.normal(-2,prior.scale_beta)
 
    if lam < 0.1:
        lam = 0.1*2 - lam

    Epeak= abs(np.random.normal(prior.mu_Epeak, prior.s_Epeak))
    
    row    = [T_star] + [ups_star] + [gamma] +  [lam] +[alpha]+[beta]+[Epeak]
    return row

def log_b(row, config, prior):
    ''' evaluates birth probability of a single pulse (row - list of parameters). 
    depends on dictionaries config and prior
    
     inputs: list row, dictionaries config, prior
      outputs: scalar evaluation of prior'''
    T, ups,gamma,lam,alpha,beta, Epeak = row
    x = (
         np.log(dist.expon.pdf(ups-config.eps, scale = 1./config.psi )) 
         + np.log(dist.norm.pdf(gamma,0,scale = prior.sigma_gamma))
         + np.log(dist.gamma.pdf(lam,prior.shape_scale_lambda_t[0], 0, prior.shape_scale_lambda_t[1]))
         + np.log(dist.norm.pdf(alpha,-1,prior.scale_alpha))
         + np.log(dist.norm.pdf(beta,-2,prior.scale_beta))
         + np.log(dist.norm.pdf(Epeak,prior.mu_Epeak,prior.s_Epeak))
         )
         
    return x
         

def eval_prior(row, config,prior):
   '''evals the prior for the adaptive parameters) for a single row of th
      inputs: list row, dictionaries config, prior
     outputs: scalar evaluation of prior'''
   x= (np.log(dist.norm.pdf(row[2],0,prior.sigma_gamma)) #gamma 
            +np.log(sp.stats.gamma.pdf(1./row[3], prior.shape_scale_lambda_t[0], 0, scale = prior.shape_scale_lambda_t[1])) #lambda
            +np.log(dist.expon.pdf(-row[4], scale = prior.scale_alpha))
          #  + np.log(dist.norm.pdf(row[5],-2, scale = prior.scale_beta))           
            + np.log( dist.norm.pdf(row[6],prior.mu_Epeak,prior.s_Epeak))
               )
  
   return x

   
########################################    
import mpmath

def exp_cts(th, config,prior,Xt):
    '''calc's the expected counts in bins with endpoints defined by Xt, with pulses def'd by th
    th: J x 7 np.array for J pulses, 7 parameters per pulse.
    config, prior: bunchs of parameters
    Xt: list of times (truncated according to the parameters in config)
    
    returns:
    exp_counts:  4 x (len(Xt)-1) nparray with the expected counts in the len(Xt)-1 bins for the 4 channels. note that 
    nchan is hardcoded here.

'''
    int_length=np.array([float(Xt[ind+1]-Xt[ind]) for ind in xrange(len(Xt)-1)])
    k = len(th) # number  of pulses

    
    drms = config.drms
    out  = 0
    quad_nodes = drms.quad_nodes


    n_intervals = len(int_length)
    n_chan = len(config.B)

    exp_counts = np.zeros((4,n_intervals))
    thinned_counts = np.zeros((4,n_intervals))
    
    energy_int = []
    for j in range(k):
       T,A,gam,lam,alpha,bet,Ec = th[j]
       temp =   (Ec)*((Ec/100.)**alpha*mpmath.gammainc(alpha+1,a = 10.*(1./Ec), )  )/lam

       energy_int.append(temp  )

    E_l, E_u = drms.E_l, drms.E_u 
    exp_cts_0 = np.zeros((4,n_intervals))
    
    for i in range(n_intervals):
       for j in range(k):
            gam = th[j,2] 
            if gam >= 0:
                lower = 100*np.exp( ((Xt[i+1]+Xt[i])/2 - th[j,0] )/(-th[j,2]) )
                l = round(max(lower, E_l),2)
                u = E_u
            if gam < 0:
                upper = 100*np.exp( ((Xt[i+1]+Xt[i])/2 - th[j,0] )/(-th[j,2]) )
                l = E_l
                u = round(min(E_u,upper),2)
                
            Range = [l,u]
            if l>=u:
                Range = None
            if Range != None:
                nodes, wts = mbs_quad_nodes_wts(l, u,config)
                
                exp_cts_0[:,i] += np.copy(np.sum(K(np.array([Xt[i],Xt[i+1]]),nodes,np.array([th[j]]))*wts, axis = 1))/energy_int[j]
                

       exp_counts[:,i] =  np.array(config.B)*int_length[i]+ np.array(exp_cts_0[:,i])
    return exp_counts

def llh(th,config,prior,Xt,Yt, r, p_thin =1.):
    ''' evaluates likelihood of array of pulses th.
    th: np array with J rows for the J pulses
    config, prior: dictionaries of parameters
    Xt: list of times (truncated according to the parameters in config)
    Yt: List of nchan (typically 4) lists of photon counts.
    llh_old: scalar, loglikelihood of th
    p_thin: scalar, thinning parameter
    
    returns:
    [likelihood value (scalar), expected counts (nchan x len(Xt)-1) ] 
    '''


    exp_counts = exp_cts(th,config,prior,Xt)
    shape  = [r*p_thin*ind for ind in exp_counts]
    p =r/(r + 1)  
    
    chisq = np.sum((Yt - exp_counts)**2/exp_counts)
    out = np.sum(dist.nbinom.logpmf(Yt,shape,p))
    return (out,exp_counts,chisq)


#########################################            
def pois_thin(Y_in,p):
    '''Thins the count data at level p
    
    Inputs:
    Y_in: list of photon counts (unthinned)
    p: scalar. This is the fraction of the info we want to retain
    
    Returns:
    Y_out: list of photon counts (thinned). same length as Y_in
    '''
    
    Y_out=np.zeros((4,len(Y_in[0])))#[[0]*len(Y_in[0])]*len(Y_in)
    for chan in range(len(Y_in)):
       for index in range(len(Y_in[0])):
          #  if Y_in[chan,index]>0:
            Y_out[chan,index]=np.random.binomial(Y_in[chan,index],p)
    return Y_out
    

##########################################    
def grb_mcmc(prior,config,starting_value = None):
    ''' This is the function that handles the parallel thinning and swaps the chains_swapped 
'''
    N = config.N
    
    X  = config.ascii.times
    Xt = X[X>config.tmin]
    Xt = Xt[Xt<config.tmax]
    Y=config.ascii.counts

    Y_full = Y[:,X[0:-1]>config.tmin]
    Y_full = Y_full[:,Xt[0:-1]<config.tmax]
    

    tmin = config.tmin
    tmax = config.tmax
    T    =  tmax-tmin
    
    X_mid=np.array([float(X[ind]-X[ind-1])/2+ X[ind-1]for ind in range(1,len(X))])
    int_length=[float(X[ind]-X[ind-1]) for ind in range(1,len(X))]
    
    X_data = [Xt, int_length]



    p_thin = config.p_thin
    #initializing lists of data, Thetas
    Y  = []
    TH = [[] for index in range(len(p_thin))]
    theta_hold =[]
    R_STORE =[]
    
    CHISQ = [[] for i in range(len(p_thin))]
    
    n = len(p_thin)

    ## create thinned data
    for t in p_thin:
        Y=Y+[pois_thin(Y_full,t)]
    llh_store   =[]#[]]
    for t in p_thin:
        llh_store.append([])
        if config.overdisp == True:
            R_STORE.append([10.])
        else:
            R_STORE.append([10000.])
    
    params_tried = [[0 for i in range(7) ] for j in range(n)]
    params_accepted = [[0 for i in range(7) ] for j in range(n)]
    STEPS = [ [] for j in range(n)]
    ## create initial starting positions
    for j in range(n):
        if starting_value == None: 
            theta_start = starting_val(config,prior)
        else: 
            theta_start=starting_value
        print theta_start, 'starting value'    
        index = 0
        llh_start = 0
        
        if config.prior == False:
            llh_start = llh(theta_start,config,prior,Xt,Y[j], R_STORE[j][-1], p_thin =p_thin[j])[0]
        while  llh_start == -np.inf:
            theta_start = starting_val(config,prior)
            index+=1
            print index, 'failed to find starting position'
        theta_hold.append(theta_start)
        llh_store[j].append( llh_start)
        del theta_start
        print llh_store, 'beginning llh',config.prior, llh_start
    tried = np.array([np.array([0,0,0,0,0])]*n  )
    accept = np.array([np.array([0,0,0,0,0])]*n )
    swap_tries  =[0]*n
    swap_success=[0]*n
    testing = []
    chains_swapped=[]
   
    M = config.M
    for m in range(config.M):
         if m%config.swap == 0:
             chains_swapped=chains_swapped+[None]
         if m%1000 == 0:
             print 'iteration ', m, 'llh',llh_store[0][-1]
         if m%config.redraw==0:
            Y  = []
            for t in range(len(p_thin)):
                 Y=Y+[pois_thin(Y_full,p_thin[t])]
                 llh_store[t][-1] = llh(theta_hold[t],config,prior,Xt,Y[t], R_STORE[t][-1],p_thin =p_thin[t])[0]
                 
        ## here we take a regular step for each chain
         for j in range(n):
   
             ### gibbs step for "r", the negative binomial parameter
             r_in = R_STORE[j][-1]

             invover_disp = r_in/(1+r_in)
             inv_star = np.random.normal(invover_disp,0.1)
             while inv_star > 1:
               inv_star = 2 - inv_star
             while inv_star < 0:
                inv_star = - inv_star
             r_star = inv_star/(1-inv_star)
             if config.overdisp == True: # only need to innovate if we are allowing overdispersion
                 llh_old = 0
                 llh_new = 0
                 if config.prior == False:
                     llh_old =llh_store[j][-1]
                     llh_new = llh(theta_hold[j],config,prior,Xt,Y[j],r_star,p_thin[j])[0]
                 LH = llh_new - llh_old 
                 if (np.random.exponential(scale=1)+LH< 0): #move is not accepted
                    r_star = r_in
             else: 
                r_star = 10000.
             R_STORE[j].append(r_star)
             ###

             step = np.random.beta(1,1)
             if step < config.prob[0]:
                 tried[j][0]+=1
             #   'birth'                 
                 temp_out, accept_temp, llh_temp,stepsizes,chisq = birth(np.copy(theta_hold[j]),config,prior,Xt,Y[j],llh_store[j][-1],p_thin[j],R_STORE[j][-1],m=m)                    
                 theta_hold[j]= np.copy(temp_out)
                 
             elif step < config.prob[0] + config.prob[1] and step> config.prob[0]: 
              #  'death'
                 tried[j][1]+=1
                 temp_out, accept_temp, llh_temp, stepsizes,chisq = death(np.copy(theta_hold[j]),config,prior,Xt,Y[j],llh_store[j][-1],p_thin[j],R_STORE[j][-1],m=m)               
                 theta_hold[j]= np.copy(temp_out)
                     
             elif step < config.prob[0] + config.prob[1] + config.prob[2] and step> config.prob[0] + config.prob[1]:
              #  'walk'
                 tried[j][2]+=1
                 temp_out, accept_temp, llh_temp,stepsizes,chisq = walk(np.copy(theta_hold[j]),config,prior,Xt,Y[j],llh_store[j][-1],p_thin[j],R_STORE[j][-1],m=m)   
                 theta_hold[j]= np.copy(temp_out)
                 STEPS[j].append(stepsizes)
 
             elif step < config.prob[0] + config.prob[1] + config.prob[2] + config.prob[3] and step> config.prob[0] + config.prob[1] + config.prob[2]:
              #  'split'
                 tried[j][3]+=1
                 temp_out, accept_temp, llh_temp,stepsizes,chisq= split(np.copy(theta_hold[j]),config,prior,Xt,Y[j],llh_store[j][-1],p_thin[j],R_STORE[j][-1],m=m)   
                 theta_hold[j]= np.copy(temp_out)    

             elif step> config.prob[0] + config.prob[1] + config.prob[2] + config.prob[3]:
              #  'merge'
                 tried[j][4]+=1
                 temp_out, accept_temp, llh_temp,stepsizes,chisq = merge(np.copy(theta_hold[j]),config,prior,Xt,Y[j],llh_store[j][-1],p_thin[j],R_STORE[j][-1],m=m)   
                 theta_hold[j]= np.copy(temp_out)    

                 
             CHISQ[j].append(chisq)
             llh_store[j]=llh_store[j]+[llh_temp]                 
                 
             if m > config.burn:
                 if m % config.thin == 0:
                    TH[j].append(np.copy(temp_out))
             
             accept[j]+=accept_temp

             del temp_out
        
         rej=1
        #here we are trying to switch between two chains
         if m%20==0:
            if n>1:
                rej = 0
                n1  = np.random.randint(0,n-1)
                n2  = n1 + 1

            if rej==0:

                r2 = np.copy(R_STORE[n2][-1])
                r1 = np.copy(R_STORE[n1][-1])
         
                swap_tries[min(n1,n2)]=swap_tries[min(n1,n2)]+1
                th1=np.copy(theta_hold[n1] )
                #likelihood of state 1 under data 1
                llh_11=llh(th1,config,prior,Xt,Y[n1],r1, p_thin[n1])[0]
                #prior of state 1 under prior 1
                th2=np.copy(theta_hold[n2] )
                
                #likelihood of state 2 under data 2
                llh_22=llh(th2,config,prior,Xt,Y[n2],r2, p_thin[n2])[0]
  
     
                #the likelihood of state one under data #2
                llh_12=llh(th1,config,prior,Xt,Y[n2],r1,  p_thin[n2])[0]
                #the likelihood of state two under data #1
                llh_21=llh(th2,config,prior,Xt,Y[n1],r2, p_thin[n1])[0]

                LH=(llh_12+llh_21-llh_11-llh_22)
                if (np.random.exponential(scale=1)+LH>0):

                    swap_success[min(n1,n2)]=swap_success[min(n1,n2)]+1
                    
                    theta_hold[n1] = np.copy(th2)
                    
                    R_STORE[n2][-1] = np.copy(r1)
                    R_STORE[n1][-1] = np.copy(r2)

                    theta_hold[n2]= np.copy(th1)
                    llh_store[n2][-1] = np.copy(llh_12)
                    llh_store[n1][-1] = np.copy(llh_21)
                    
                    chains_swapped[-1]=min(n1,n2)

    print tried, accept, 'bdw accepts'
    print swap_tries, swap_success, 'swap accepts'
    print params_tried, params_accepted, 'params accepted'
    
    r_thinned = R_STORE[-1][config.burn::config.thin]
    print TH[-1][-1], 'final'
    return [TH, accept, llh_store,Xt,Y,swap_success, [tried,accept],r_thinned, STEPS, CHISQ] 
              
##########################











