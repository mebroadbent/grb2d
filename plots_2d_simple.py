import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import colorsys
import alpha_2d_simple as code3
import matplotlib.pylab as pylab
def rgb_range(color_num):
    h_rng = [ix / float(color_num) for ix in range(color_num)]
    return [colorsys.hsv_to_rgb(h, 1.0, 1.0) for h in h_rng]

import ColorScale as c

def find_rgb(color_hex):
     return [float(int(color_hex[i:i+2], 16))/255. for i in range(0, len(color_hex), 2)]


import matplotlib.animation as animation




def summed_pulses(TH, prior,config,Xt,Yt,p_thin):
    plt.figure()
    N = len(TH)

    int_length = [Xt[i+1] - Xt[i] for i in range(0, len(Xt)-1)]   
    Y_star = [[Yt[j][i]/int_length[i] for i in range(len(int_length)) ] for j in range(4)]
  #  print Y_star[0][0:5]
    plt.subplot(221)
    plt.scatter(Xt[1::],Y_star[0],alpha = 0.2)    
   # plt.xlim(0,max(Xt))
    
    plt.xlabel('time since trigger (s)')
  #  plt.ylim(np.min(Yt[0])*0.9, 400)
    plt.ylabel('photon rate (/s)')
    plt.axis('tight')
    plt.title('Channel 1')
    
    plt.subplot(222)
    plt.scatter(Xt[1::],Y_star[1],alpha = 0.2)

    plt.xlabel('time since trigger (s)')
    plt.ylabel('photon rate (/s)')
    plt.axis('tight') 
    plt.title('Channel 2')
    
   
    plt.subplot(223)
    plt.scatter(Xt[1::],Y_star[2],alpha = 0.2)
   # plt.xlim(0,max(Xt))
   # plt.ylim(np.min(Yt[2])*0.9, 400)
    plt.xlabel('time since trigger (s)')
    plt.ylabel('photon rate (/s)')
    plt.axis('tight')
    plt.title('Channel 3')    

    plt.subplot(224)
    #plt.xlim(0,max(Xt))
    #plt.ylim(np.min(Yt[3])*0.9, 400)
    plt.scatter(Xt[1::],Y_star[3],alpha = 0.2)
    plt.xlabel('time since trigger (s)')
    plt.ylabel('photon rate (/s)')
    plt.axis('tight')
    plt.title('Channel 4')
   

    num = 100
    scale = c.ColorScale(0,num)
    colors = scale.range('red','green')
    #rgb_range(num)
    for i in range(num):
       n_i = int(i*N/num)
       print i, colors(i)
       like, count, val = code3.llh(TH[n_i],config, prior, Xt, Yt, 10000., p_thin)

       for j in range(4):
           count[j] = count[j]/int_length
       plt.subplot(221)
       plt.plot(Xt[1::],count[0], color = find_rgb(colors(i)),alpha = 0.1)
       plt.axis('tight')        
       plt.subplot(222)
       plt.plot(Xt[1::],count[1],color = find_rgb(colors(i)),alpha = 0.1)
       plt.axis('tight')  
       plt.subplot(223)
       plt.plot(Xt[1::],count[2], color = find_rgb(colors(i)),alpha = 0.1)
       plt.axis('tight')   
       plt.subplot(224)
       plt.plot(Xt[1::],count[3],color = find_rgb(colors(i)),alpha = 0.1)
       plt.axis('tight')   
    plt.tight_layout()
    plt.show()
    
    
def individual_pulses(th,prior,config,Xt,Yt,r,p_thin):
    ''' plots the individual pulses at a single iteration of the MCMC'''
    plt.figure() 
    

    int_length = [Xt[i+1] - Xt[i] for i in range(0, len(Xt)-1)]   
    Y_star = [[Yt[j][i]/int_length[i] for i in range(len(int_length)) ] for j in range(4)]
    
    K = ('b','g','r','c','m','y')
    J = len(th)
    ALL_COUNTS = []
    raw_counts = []
    for j in range(J):
          counts = code3.exp_cts(np.array([th[j]]), config,prior,Xt)
          ALL_COUNTS.append(np.array(counts)) #/int_length)

    legend_facts = []      
    # the shape of ALL_COUNTS IS : J x 4 x N_times
    ALL_COUNTS =  np.array(ALL_COUNTS)/int_length
    plt.subplot(221)

    X_star = Xt[1::]

    
    plt.scatter(X_star,Y_star[0],alpha = 0.2)
    counts_raw= code3.exp_cts(th,config, prior, Xt)
    counts_total = np.array(counts_raw/int_length)
    plt.plot(X_star,counts_total[0],color = 'black') 
    for j in range(J):
        x = plt.plot(X_star,ALL_COUNTS[j,0,:],color = K[j%(len(K))])
 #   plt.plot(X_star,counts_total[0],color = 'black')


    plt.xlabel('time since trigger (s)')
    plt.ylabel('photon rate (/s)')
    plt.axis('tight')
    plt.title('Channel 1')   
    print np.shape(Yt)
    chisq = np.sum( (counts_raw[0] - Yt[0])**2/counts_raw[0])
    dof = str((len(Y_star[0])*4 - J*7 )/4)
    plt.text(0,max(Yt[0])*0.7,r'$\chi^2$ ' + str(round(chisq)) + ' on ' + dof +' dof')

    
    plt.subplot(222)
    plt.plot(X_star,counts_total[1],color = 'black')
    plt.scatter(X_star,Y_star[1],alpha = 0.2)
    for j in range(J):
        plt.plot(X_star,ALL_COUNTS[j,1,:])
    plt.axis('tight')

    plt.xlabel('time since trigger (s)')
    plt.ylabel('photon rate (/s)')
    #plt.xlim(0,max(X_star))
    #plt.ylim(np.min(Y_star[1])*0.9, 400)
    plt.title('Channel 2')    
    
    chisq = np.sum( (counts_raw[1] - Yt[1])**2/counts_raw[1])
    plt.text(0,max(Yt[1])*0.7,r'$\chi^2$ ' + str(round(chisq)) + ' on ' + dof + ' dof')

        
    plt.subplot(223)
    plt.plot(X_star,counts_total[2],color = 'black')
    plt.scatter(X_star,Y_star[2],alpha = 0.2)
    for j in range(J):
       plt.plot(X_star,ALL_COUNTS[j,2,:])
    plt.axis('tight')
    plt.xlabel('time since trigger (s)')
    plt.ylabel('photon rate (/s)')
    plt.title('Channel 3')
    
    chisq = np.sum( (counts_raw[2] - Yt[2])**2/counts_raw[2])
    plt.text(0,max(Yt[2])*0.7,r'$\chi^2$ ' + str(round(chisq)) + ' on ' + dof + ' dof')
    print chisq


    plt.subplot(224)
    plt.plot(X_star,counts_total[3],color = 'black')
    plt.scatter(X_star,Y_star[3],alpha = 0.2)
    for j in range(J):
        plt.plot(X_star,ALL_COUNTS[j,3,:])
    plt.axis('tight')
    plt.xlabel('time since trigger (s)')
    plt.ylabel('photon rate (/s)')
    plt.title('Channel 4')   
    chisq = np.sum( (counts_raw[3] - Yt[3])**2/counts_raw[3])
    plt.text(0,max(Yt[3])*0.7,r'$\chi^2$ ' + str(round(chisq)) + ' on ' + dof + ' dof')

    
    plt.tight_layout()      
    plt.show()
    
def hist_r(r_full,weights = None):
    print weights
    if weights is None:
      plt.figure()
      plt.hist(r_full, color = '0.5',normed = True)
      plt.title(r'$r$')
      plt.savefig('rhist.png')
      
      plt.figure()
      plt.hist( (1+np.array(r_full))/np.array(r_full), color = '0.5',normed = True)
      plt.title(r'$(1+r)/r$')
      plt.savefig('rhist2.png')
      
      plt.figure()
      plt.hist( np.array(r_full)/(1+np.array(r_full)), color = '0.5',normed = True)
      plt.title(r'$r/(1+r)$')
      plt.savefig('rhist3.png')

    if weights is not None:
      plt.figure()
      plt.hist(r_full, color = 'white',normed = True,weights = weights)
      plt.title(r'$r$')
      plt.savefig('rhist.png')
      
      plt.figure()
      plt.hist( (1+np.array(r_full))/np.array(r_full), color = 'white',normed = True,weights = weights)
      plt.title(r'$(1+\beta)/\beta$')
      plt.savefig('rhist2.png')
      
      plt.figure()
      plt.hist( np.array(r_full)/(1+np.array(r_full)), color = 'white',normed = True,weights = weights)
      plt.title(r'$\beta/(1+\beta)$')
      plt.savefig('rhist3.png')

    
def plot_r(r_out,param):
    mb =  np.arange(param.burn-1, param.M, param.thin) 
    print len(mb), len(r_out)
    plt.figure()
    plt.plot(mb,np.array((1+np.array(r_out))/r_out))
    plt.xlabel('MCMC iterations')
    plt.ylabel(r'$(1+\beta)/\beta$')
    plt.title(r'Inflation parameter $(1+\beta)/\beta$ ')
    plt.savefig('r_out.png')

    plt.figure()
    plt.plot(mb,np.array(r_out))
    plt.xlabel('MCMC iterations')
    plt.ylabel(r'$\beta$')
    plt.title(r'Inflation parameter $\beta$ ')
    plt.savefig('r_out2.png')

    plt.figure()
    plt.plot(mb, np.array(r_out)/(1+np.array(r_out)))
    plt.title(r'$\beta/(1+\beta)$')
    plt.savefig('rout3.png')
    
def amplitude_traceplots(TH,prior,config,true=None,weights = None):
    plt.figure()
    J = config.N
    mb =  np.arange(config.burn, config.M, config.thin) 
    

    Amps = []
    for i in range(len(TH)):
        Amps.append( TH[i][:,1]) 
        
    Amps = np.array(Amps)
    if len(mb)!=len(Amps):
        mb = range(len(Amps))
    
    plt.plot(mb, Amps)
    plt.xlabel('MCMC iterations (thinned)')
    plt.ylabel('Volumes')
    plt.title('Pulse Volume traceplot')
    if true is not None:
        plt.axhline(y = true,xmin = 0,xmax = 100,color = 'k')
    plt.axis('tight')
    plt.tight_layout()
    pylab.savefig('V.png')

    
    plt.figure()
    plt.hist(Amps, color = 'white',normed = True,weights = weights)
    plt.title('Pulse Volume Histogram')
    plt.axvline(x = true, ymin = 0, ymax =1,color ='k')
    plt.tight_layout()
    plt.show()
    plt.savefig('Vhist.png')
    
def e_peak_traceplots(TH,prior,config,true=0,weights = None):
    plt.figure()
    mb =  np.arange(config.burn, config.M, config.thin) 
    E_peak = []
    for i in range(len(TH)):
        E_peak.append(TH[i][:,-1])
    E_peak = np.array(E_peak)
    if len(mb)!=len(E_peak):
        mb = range(len(E_peak))    
    
    plt.plot(mb, E_peak)
    plt.xlabel('MCMC iterations (thinned)')
    plt.ylabel('Peak Energy (keV)')
    plt.title(r'$E_{c}$ traceplot')
    plt.axhline(y = true, xmin = 0, xmax =100,color = 'k')
    plt.axis('tight')
    plt.tight_layout()
    plt.show()
    pylab.savefig('epeaktrace.png')
    
    plt.figure()
    plt.hist(E_peak, color = 'white',normed = True,weights = weights)
    plt.title(r'$E_c$ histogram')
    plt.axvline(x=true, ymin =0 , ymax = 1000, color = 'k')
    plt.tight_layout()
    plt.show()
    pylab.savefig('ec-hist.png')

    #  plt.show()
    
def plot_gamma(TH,prior,config,true=None,weights = None):            
    plt.figure()
    mb =  np.arange(config.burn, config.M, config.thin) 
    J = config.N

    gamma = []
    for i in range(len(TH)):
        gamma.append(TH[i][:,2])
    gamma = np.array(gamma)
    
   # gamma = [g[0] for g in gamma]
    if len(mb)!=len(gamma):
        mb = range(len(gamma))    
    print np.shape(gamma), np.shape(mb)
    print gamma[0:5], mb[0:5]

    plt.plot(mb, gamma)
    if true is not None:
        plt.axhline(y = true, xmin = 0, xmax = 50000, color = 'k')
    plt.xlabel('MCMC iterations (thinned)')
    plt.title('Lag traceplot')
    plt.axis('tight') 
    plt.tight_layout()
    pylab.savefig('gamma-trace.png')
    plt.figure()
    plt.hist(gamma, color = 'white',normed = True,weights = weights)
    plt.title(r'$\gamma$ histogram')
    if true is not None:
    
        plt.axvline(x = true,ymin = 0,ymax = 100, color = 'k')
    plt.tight_layout()
    plt.savefig('gamma-hist.png')
    #  plt.show()
    
def plot_lambda(TH,prior,config, true = None,weights = None):            
    mb =  np.arange(config.burn, config.M, config.thin) 

    plt.figure()
    J = config.N
    lam = []
    for i in range(len(TH)):
        lam.append(TH[i][:,3])
    lam = np.array(lam)
    print np.shape(lam)
    if len(mb)!=len(lam):
        mb = range(len(lam))    
    plt.plot(mb,lam)
    plt.ylim(np.min(lam), np.max(lam))
    plt.xlabel('MCMC iterations (thinned)')
    plt.title(r'Time decay parameter $\lambda$ traceplot')
    plt.axis('tight')
    if true is not None:
        plt.axhline(y =true, xmin = 0, xmax = 55000, color = 'k')
    plt.tight_layout()
    plt.show()
    pylab.savefig('lambda.png')
    
    plt.figure()
    plt.hist(lam, color = 'white',normed = True,weights = weights)
    if true is not None:
        plt.axvline(x = true, ymin = 0, ymax = 100, color = 'k')
    plt.title(r'$\lambda$ histogram')
    plt.tight_layout()
    plt.savefig('lam-hist.png')

   # plt.show()
        
    
def plot_ab(TH,prior,config,Xt,Yt,a_compare = None, b_compare = None,weights = None):            
    plt.figure()
    J = config.N
    a = []
    mb =  np.arange(config.burn, config.M, config.thin) 
    b = []
    for i in range(len(TH)):
        a.append(TH[i][:,4])
        b.append(TH[i][:,5])
        
    if len(mb)!=len(a):
        mb = range(len(a))        
    a = np.array(a)
    b = np.array(b)
    plt.subplot(111)
    plt.plot(mb,a)
    if a_compare is not None:
        print np.max(mb), config.burn
        print config.M
        plt.axhline(y =a_compare, xmin = 0, xmax = np.max(mb), color = 'k')
    plt.xlabel('MCMC iterations (thinned)')
    plt.title(r'$\alpha$ traceplot')
    plt.axis('tight')
    plt.show()
    
    plt.figure()
    plt.hist(a, color = 'white',normed = True,weights = weights)
    plt.title(r'$\alpha$')
    if a_compare is not None:
        plt.axvline(x = a_compare, ymin = 0, ymax = 1)
    plt.savefig('alpha.png')

def hist_J(theta_final,weights = None):
    J = [len(th) for th in theta_final]
    minJ = min(J)
    maxJ = max(J)
    N =max(4,maxJ)
    plt.tight_layout()
    bins = np.arange(-0.5,N+0.5,1)
    plt.figure()
    
    plt.hist(J, bins = bins, color = '0.5',normed = True,weights = weights)
    plt.title('Histogram of number of pulses')
    plt.tight_layout()
    plt.ylim([0,1])
    
    plt.xticks([i for i in range(N)],[str(i)[0] for i in range(N)])
    plt.savefig('Jhist.png')
        
def J_traceplot(J,config):
    plt.figure()
    mb =  np.arange(config.burn, config.M, config.thin) 
    plt.plot(J)
    plt.title('J traceplot')

    
def llh_traceplots(config,LLH):
    mb =  np.arange(config.burn, config.M, config.thin) 
    print len(mb)
    print np.shape(LLH[:][config.burn::config.thin] )
    #if len(mb) !=len(LLH[0]):
    #    mb = range(len(LLH[0]))
    plt.figure()
    plt.plot(mb,LLH[config.burn ::config.thin][1::] )
    plt.xlabel('MCMC iterations (thinned)')
    plt.ylabel('LLH')
    plt.tight_layout()
    plt.title('Log-likelihood for All Chains')

from helpers import effectiveSampleSize as eff
def chisq_traceplots(TH,config, prior,Xt =None,Yt =None):
    mb =  np.arange(config.burn, config.M, config.thin) 

    plt.figure()
    print len(mb)
    print len(TH)
    #plt.plot(mb,CHISQ[:][config.burn::config.thin] )
    if Xt == None:
        Xt,Yt = code3.get_XY(config.ascii.times,config.ascii.counts,config.tmin,config.tmax)
    dof = len(Xt)*4 
    print dof
    plt.axhline(y = dof,color = 'red')
    
    chisq = []
    z = 0
    
    for th in TH:
        z +=1
        if z%1000 == 0:
            print z, 'chisq calc',chisq[-1]
            print [np.min(chisq),max(np.max(chisq),dof+5)]
        chisq.append(code3.llh(th, config,prior,Xt,Yt, 10000., p_thin =1.)[2])
    plt.plot(mb[1::],chisq )
    plt.xlim(mb[1],mb[-1])
    plt.xlabel('MCMC iterations (thinned)')
    plt.ylabel(r'$\chi^2$')
    plt.tight_layout()
    plt.ylim(np.min(chisq),max(np.max(chisq),dof+5))
    plt.title(r'$\chi^2$ for All Chains')
    plt.show()
    E = eff(chisq)
    print 'Effective Sample Size',E
from scipy.misc import factorial as fac
def redo_traceplots(config,prior,Xt,Yt,TH,LLH):
    ''' plots a new chisq and a new prior + likelihood plot '''
    chisq = []
    llh   = LLH[:][config.burn+1::config.thin]
    prior_list_1 = []
    prior_list_2 = []
    prior_rlw = []
    for th in TH:
        J = len(th)
        #llh_temp,exp_counts,chisq_temp = code3.llh(th,config,prior,Xt,Yt, 10000., p_thin =1.)
        
      #  llh.append(llh_temp)
        #chisq.append( chisq_temp)       

        prior_list_1.append( J*np.log(prior.alpha) - prior.beta*np.sum(th[:,1]) - np.log(fac(J)) - np.log(np.product(th[:,1])) )       
        prior_temp = 0
        prior_rlw_temp = 0
        for th_j in th:
            mb = code3.eval_prior(th_j, config,prior)
            prior_temp +=  mb
            
            prior_rlw_temp += ( -np.log(2*3.14*prior.sigma_gamma*prior.s_Epeak/prior.scale_alpha) 
                                + th_j[2]**2/(2*prior.sigma_gamma**2) + th_j[3] + prior.scale_alpha*th_j[4] + (th_j[6] - prior.mu_Epeak)**2/(2*prior.s_Epeak**2)
                                + J * np.log(prior.alpha) )
            
        prior_list_2.append(prior_temp)
        prior_rlw.append(prior_rlw_temp)
        
        
        
        
    llh = np.array(llh)
    prior_list_1 = np.array(prior_list_1)
    prior_list_2 = np.array(prior_list_2)
    prior_mb = prior_list_1 + prior_list_2 - 50
    
    
    
  #  prior_list =np.array(prior_list)
    llh = llh - np.max(llh)
   # prior_list = prior_list - np.max(prior_list)
 #   posterior = llh[:,0] + prior_list
        
    mb = np.arange(config.burn,config.M,config.thin)
    #plt.figure()
    #plt.plot(mb,chisq)
    #plt.xlabel('MCMC iterations (thinned)')
    #plt.ylabel(r'$\chi^2$')
    #plt.title(r'$\chi^2$ for All Chains')
    #plt.tight_layout()
    #plt.show()
    #plt.savefig('chisq_new.png')
   # print np.shape(mb), np.shape(prior_list), np.shape(llh), np.shape(posterior)
    plt.figure()
    plt.plot(mb,llh)
    plt.plot(mb,prior_mb)
    plt.plot(mb,prior_rlw)
    plt.legend(['LLH',r'Log-prior MB',r'Log-prior RLW'], loc = 3)
   # plt.legend(['LLH',r'Log-prior $\nu_x \times \nu_u$',r'Log-prior $\nu_a$'], loc = 3)
    plt.xlabel('MCMC iterations (thinned)')
    plt.title('Log-densities')
    plt.tight_layout()
    plt.ylim([-200,0])
    plt.show()
    plt.savefig('densities.png')
    
    
        
    

def t_traceplots(TH,config,true=None,weights = None):

    mb =  np.arange(config.burn, config.M, config.thin) 

    plt.figure()
    J = config.N
    lam = []
    for i in range(len(TH)):
        lam.append(TH[i][:,0])
    lam = np.array(lam)
    if len(mb) !=len(lam):
        mb = range(len(lam))
    plt.plot(mb,lam)
    plt.xlabel('MCMC iterations (thinned)')

    plt.ylabel('Pulse Start Times')
    plt.title('Pulse Start Times')
    plt.axis('tight')
    plt.tight_layout()
    if true is not None:
        plt.axhline(y = true,xmin = 0,xmax = 100,color = 'k')
    plt.show()
    plt.savefig('t-traceplot.png')
    
    plt.figure()
    plt.hist(lam, color = 'white',normed = True,weights = weights)
    plt.title('Pulse Start Times')
    if true is not None:
        plt.axvline(x = true, ymin = 0, ymax = 1, color = 'k')
    plt.show()
    plt.tight_layout()
    pylab.savefig('t-hist.png')
  #  plt.show()

from scipy.stats.mstats import mquantiles as quantile
def ppi(TH,prior,config, Xt, Yt,p_thin,r_thinned):
    plt.figure()
    #############
    int_length = [Xt[i+1] - Xt[i] for i in range(0, len(Xt)-1)]   
    Y_star = [[Yt[j][i]/int_length[i] for i in range(len(int_length)) ] for j in range(4)]
    

    
  #  plt.scatter(X_star,Y_star[0],alpha = 0.2)
  #  for j in range(J):
  #      x = plt.plot(X_star,ALL_COUNTS[j,0,:],color = K[j%(len(K))])
    #########33

    n = 250
    TH = TH[(len(TH)/2)::]
    SHAPE = []
    P =[]
    for i in range(n):
         n_temp = i*len(TH)/n
         likelihood1, counts,vals1 = code3.llh(TH[n_temp],config,prior,Xt,Yt,p_thin)
         SHAPE.append(counts*r_thinned[n_temp]*p_thin)
    #the shape of counts is n x 4 x n_times
         P.append(r_thinned[n_temp]/(1+r_thinned[n_temp]))
    SHAPE = np.array(SHAPE)
   # P = r_thinned/(1+r_thinned)
    for chan in range(4):
        lower =[]
        upper =[]
        for t in range(len(Xt) -1):
          #  print np.shape(SHAPE[:,chan,t]), np.shape(P)
            relevant_counts = np.random.negative_binomial(SHAPE[:,chan,t],P)/int_length[t]
            #relevant_counts = np.random.poisson(COUNTS[:,chan,t])/int_length[t]
            quants = quantile(relevant_counts,[0.025,0.975]) 
            lower.append(quants[0])
            upper.append(quants[1])
        index = '22'+str(chan+1)
        plt.subplot(int(index))
        
        print len(Xt), len(Yt[0]), np.shape(Yt)
        plt.fill_between(Xt[1::],lower,upper,color = '0.5')
        plt.scatter(Xt[1::],Y_star[chan],alpha = 0.2)
        plt.xlabel('time (s) since trigger')
        plt.ylabel('PPI')
        plt.title('Channel ' + str(chan + 1))
        plt.axis('tight')
    plt.tight_layout()        
    plt.show()

def ci(TH,prior,config, Xt, Yt,p_thin,r_thinned,theta):
    int_length = [Xt[i+1] - Xt[i] for i in range(0, len(Xt)-1)]      
    Y_star = [[Yt[j][i]/int_length[i] for i in range(len(int_length)) ] for j in range(4)]
    print theta, 'true model'
    plt.figure()

    n = 100
    TH = TH[(len(TH)/2)::]
    COUNTS = []
    P = []
    for i in range(n):
         n_temp = i*len(TH)/n
         likelihood1, counts,vals1 = code3.llh(TH[n_temp],config,prior,Xt,Yt,p_thin)
         P.append(r_thinned[n_temp]/(1+r_thinned[n_temp]))
         COUNTS.append(counts*p_thin)
    #the shape of counts is n x 4 x n_times
      #   P.append(r_thinned[n_temp]/(1+r_thinned[n_temp]))
    SHAPE = np.array(COUNTS)
   # P = r_thinned/(1+r_thinned)
    for chan in range(4):

        lower_ci = []
        median_ci = []
        upper_ci = []
        for t in range(len(Xt) -1):
          #  print np.shape(SHAPE[:,chan,t]), np.shape(P)
            relevant_counts = SHAPE[:,chan,t]/int_length[t]
            #relevant_counts = np.random.poisson(COUNTS[:,chan,t])/int_length[t]
            quants = quantile(relevant_counts,[0.025,0.5,0.975]) 
      #      quant_ci = quantile(SHAPE[:,chan,t],[0.025,0.5,0.975])
            lower_ci.append(quants[0])
            median_ci.append(quants[1])
            upper_ci.append(quants[2])
        index = '22'+str(chan+1)
        plt.subplot(int(index))
        
        if theta is not None:
            l, c,vals = code3.llh(theta,config,prior,Xt,Yt,1. )
                        
                        
            rel_c_temp = [c[chan][i]/int_length[i] for i in range(len(c[chan]))]
        

        plt.fill_between(Xt[1::],lower_ci,upper_ci,color = '0.5')
        if theta is not None: 
            plt.plot(Xt[1::],rel_c_temp,color = 'black')
        else:
            plt.plot(Xt[1::],median_ci, color = 'black')
        plt.scatter(Xt[1::],Y_star[chan],alpha = 0.2)
        plt.xlabel('time (s) since trigger')
        plt.ylabel('95% Interval')
        plt.title('Channel ' + str(chan + 1))
        plt.axis('tight')
    plt.tight_layout()   
    plt.show()        
    
    
    
def convert_STEPS(STEPS):
    ''' STEPS is a list of the form [ (row,accepted) for i in range(len(TH))]
        where row is a row [T,A,gam,lam,alpha, beta(deprecated),E_peak] and accepted is a boolean
        
        we want to return a list of 7 arrays. Each array has two lists of #s, which represents the proposed and accepted steps
        '''
    NEW_STEPS = [ [ [],[] ] for i in range(7)]
    N = len(STEPS)
    for i in range(N):
        step, accept = STEPS[i][0]
        for j in range(7):
            NEW_STEPS[j][0].append(step[j])
            if accept == 1:
                NEW_STEPS[j][1].append(step[j])
    return NEW_STEPS

def plot_steps(STEPS_new):
    plt.figure()
    T_prop,T_steps = STEPS_new[0]
    plt.hist(T_prop,alpha = 0.5, normed = True)
    plt.hist(T_steps,alpha =0.5,normed = True)
    plt.title('T')
    plt.legend(['Proposed','Accepted'])
    plt.show()
    plt.savefig('compare.png')
#######################################################
    A_prop, A_steps = STEPS_new[1]
    plt.figure()
    plt.hist(A_prop,alpha = 0.5, normed = True)
    plt.hist(A_steps,alpha =0.5,normed = True)
    plt.title('V')
    plt.legend(['Proposed','Accepted'])
    plt.show()
    plt.savefig('compare_V.png')

#######################################################
    g_prop, g_steps = STEPS_new[2]
    plt.figure()
    plt.hist(g_prop,alpha = 0.5, normed = True)
    plt.hist(g_steps,alpha =0.5,normed = True)
    plt.title(r'$\gamma$')
    plt.legend(['Proposed','Accepted'])
    plt.show()
    plt.savefig('compare_g.png')
#######################################################
    l_prop,l_steps = STEPS_new[3]
    plt.figure()
    plt.hist(l_prop,alpha = 0.5, normed = True)
    plt.hist(l_steps,alpha =0.5,normed = True)
    plt.title(r'$\lambda$')
    plt.legend(['Proposed','Accepted'])
    plt.show()
    plt.savefig('compare_l.png')
########################################################
    alph_prop, alph_steps = STEPS_new[4]
    plt.figure()
    plt.hist(alph_prop,alpha = 0.5, normed = True)
    plt.hist(alph_steps,alpha =0.5,normed = True)
    plt.title(r'$\alpha$')
    plt.legend(['Proposed','Accepted'])
    plt.show()
    plt.savefig('compare_alph.png')
########################################################
    E_prop, E_steps = STEPS_new[6]
    plt.figure()
    plt.hist(E_prop,alpha = 0.5, normed = True)
    plt.hist(E_steps,alpha =0.5,normed = True)
    plt.title(r'$E$')
    plt.legend(['Proposed','Accepted'])
    plt.show()
    plt.savefig('compare_E.png')
######################################################

def sort_params(TH, n = 3,weights = None):
    ''' this function sorts the posterior samples of TH into those of len 1...n. for J>1, it also orders the pulses in a sample wrt T'''
    TH_new = [ [] for i in range(n)]
    
    w_new = None
    if weights is not None:
        w_new = [ [] for i in range(n)]
    
    for i in range(len(TH)):
        th = TH[i]
        J = len(th)
        #print J
        if J < n+1:
            if J > 0:
                TH_new[J-1].append( th[np.argsort(th[:,1])])
                if w_new is not None:
                     w_new[J-1].append(weights[i])
    return TH_new,w_new

 ########################3   
def plot_J_frac(theta_final):
    '''plts and saves the fraction of amplitudes vs #J included'''
    N = len(theta_final)
    
    J = [len(th) for th in theta_final]
    max_J = max(J)
    print max_J,'maxJ'
    frac = []
    
    for i in range(N):
        # fraction for this single draw
        
        current_frac = [1. for index in range(max_J)]
        th = theta_final[i]
      #  print np.shape(th)
        A_sorted = th[:,1][th[:,1].argsort()][::-1]
        A_normed = A_sorted/np.sum(A_sorted)
        #print A_normed
        J_temp = len(A_normed)
        current_frac[:J_temp] = np.cumsum(A_normed)
        current_frac = np.array(current_frac)

        #print current_frac
        frac.append(current_frac)
        
    frac = np.array(frac)

    frac_robert = np.mean(frac,axis = 0)
    
    print np.shape(frac_robert)
    plt.figure()
    plt.plot(range(1,max_J+1),frac_robert,color = 'k')
    plt.xlabel('Number of Pulses Used (J)')
    plt.ylabel('Volume Fraction')
    plt.ylim(0,1.03)
    plt.axhline(1,color='blue')
    
    plt.title('Fraction of Pulse Volumes using Top J Pulses')
    plt.show()
    plt.savefig('robert.png')
    
    plt.figure()
    bp = plt.boxplot(frac)
    plt.setp(bp['boxes'], color='black')
    plt.setp(bp['whiskers'], color='black')
    plt.setp(bp['fliers'], color='red', marker='+')
    pylab.xticks(range(max_J+1),['']+range(1,max_J+1))
    
    plt.xlabel('Number of Pulses Used (J)')
    plt.ylabel('Volume Fraction')
    plt.ylim(0,1.03)
    plt.axhline(1,color = 'blue')
    
    plt.title('Fraction of Photons in Largest J Pulses')
    plt.tight_layout()
    plt.show()
    plt.savefig('robert2.png')
    
    
    
    
    
    
    
    
    
    
    