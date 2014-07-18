import numpy as np
import pickle
import os
import matplotlib as mpl
mpl.use('Agg')#
import matplotlib.pylab as pylab
import matplotlib.pyplot as plt
#import pandas as p
import batse5bp
import myplot
import cProfile
import alpha_2d_simple as code
from plots_2d_simple import *
######################## load data
#grb_cat = batse5bp.catalog.load_catalog(root_dir='BATSE_5Bp_Data')
#grb = grb_cat[537]
#ascii = grb.ascii64
#ascii.counts = np.random.binomial(ascii.counts, p = 1.)# round(ascii.counts*0.05)
#B = [np.mean(ascii.counts[i][0:50])/1.024 for i in range(4)]
#drms = grb.discsc_drms
#drms.set_sum_prodquad() 
########################
#class bunch:
    #def __init__(self, **kwds):
        #self.__dict__.update(kwds)
#import colorsys

#import datetime
#import pdb

#directory = '/home/grad/meb67/Dropbox/2d/2014-05-18537_b'
#os.chdir(directory)


###myfile = open('single_pulse_smc_21','r')
###draws,weights, llh=pickle.load(myfile)
###myfile.close()
###TH= [draws[i][1] for i in range(len(draws))]
###r_thinned = [draws[i][0] for i in range(len(draws))]


#os.chdir(directory)
#myfile = open('datafile','r')
#data = pickle.load(myfile)
#TH,A,LLH,Xt,Yt,swap, params_accepted,r_thinned,STEPS,CHISQ= data
#print params_accepted
##,chisq = data
#myfile.close()

#print params_accepted
import matplotlib
def hscatter(values1, values2, true=None, nburn=None, member=None):
    """
    Produce a "historical scatterplot," with points color-coded by their
    time in the time series (early=cool=blue, late=hot=red).
    """
    plt.figure()
    if member is not None:
        values1, values2 = values1[:,member], values2[:,member]
    if nburn:
        values1, values2 = values1[nburn:], values2[nburn:]
    times = np.arange(len(values1), dtype=float)/len(values1)
    cmap = matplotlib.cm.jet  # the color map
    plt.scatter(values1, values2, c=times, cmap=cmap, edgecolors='none',
               alpha=.3)
    
    #hscatter(pts[:,0],pts[:,1],nburn = 0)
   # plt.xlim(x_l,x_u)
    #plt.ylim(y_l,y_u)
    plt.xlabel(r'$V_1$')
    plt.ylabel(r'$V_2$')
    plt.scatter(true[:,0],true[:,1],marker = '+',s = 120, color = 'black')
    #plt.legend(['Posterior Draws','Source Locations'])
    plt.title(r"Historical Scatterplot $V_2$ vs $V_1$")
    plt.savefig('historical.png')




#os.chdir(newplots)

#config = bunch(grb = grb,
            #drms = drms,
            #ascii = ascii,
            #tmin = -12., #tmin and tmax indicate truncation
            #tmax = 25., 
            #M = 60000,  #number of iterations
            #B = B, #bsackground
            #burn =50, # number of burnin iterations
            #thin = 5, #thinning 
            #eps = .1, #min pulse height
            #prob = [0.2,0.2,0.5,0.05,0.05], #probs for birth, death, walk
            #d_t = .15, #d_x is the MH proposal standard deviation
            #d_A = 8.,
            #d_ab = .15,
            #d_lambda = 0.1,
            #d_gamma = 0.05,
            #d_Epeak = 8.3, 
            #psi = 1./(200.),
            #N =1, #starting # of pulses
            #p_thin =[1.], #chain thinnings
            #redraw = 10000000000, #how often to redraw thinned chains (you see that right now i am not redrawing ever..."
            #swap = 100,          
            #zeta = 1.,
            #xsi = 0.1,
            #overdisp = True,
            #prior = False)

#prior = bunch(alpha =np.exp(0.4)/37,
              #beta = 0.001, #Amplitudes
              #scale_alpha = .2,
              #scale_beta = .2,
              #mu_Epeak = 300.,#np.log(500.),
              #s_Epeak = 150.,

              #sigma_gamma = 5.,# 0.000001,#0.3,#1.,#1.4,#prior variance of linear lag
             ## shape_scale_sigma_E = (.5,2.), #parameters for sigma_E
              #shape_scale_lambda_t = (1.,1.) #parameters for lambda
            #)

#TH[-1] = TH[-1][1::]
#r_thinned = r_thinned[1::]
#LLH=  np.array(LLH)[1::]
#CHISQ = np.array(CHISQ)[1::] 



## Plots for a 1D MODEL ONLY
from plots_2d_simple import *
def plot_norj(TH, params_true,directory, config,prior,X,Y,STEPS, r_thinned,p_thin = 1):
    ''' TH is only the chain of interest, with thinning p_thin 
    Xt, Yt are the corresponding data for that chain only '''
    os.chdir(directory)
    if len(params_true) !=1:
       print 'params_true fis the wrong dimension for this function!'
    
    ## summed pulses and individual iterations
    N = len(TH[-1])/2
    n1 = N
    n2 = N + N/5
    n3 = N + N/5*2
    n4 = N + N/5*3
    n5 = N + N/5*4 
    summed_pulses(TH,prior,config,X,Y,config.p_thin[-1])
    pylab.savefig('severaldraws.png')
    individual_pulses(TH[n1],prior,config,X,Y,1000.,config.p_thin[-1])
    print n1, TH[n1]
    pylab.savefig('iter-n1.png')
    individual_pulses(TH[n2],prior,config,X,Y,10000.,config.p_thin[-1])
    print n2, TH[n2]
    pylab.savefig('iter-n2.png')
    individual_pulses(TH[n3],prior,config,X,Y,1000.,config.p_thin[-1])
    print n3, TH[n3]
    pylab.savefig('iter-n3.png')
    individual_pulses(TH[n4],prior,config,X,Y,1000.,config.p_thin[-1])
    print n4, TH[n4]
    pylab.savefig('iter-n4.png')
    individual_pulses(TH[n5],prior,config,X,Y,1000.,config.p_thin[-1])
    pylab.savefig('iter-n5.png')
    print n5, TH[n5]
    

    print params_true[0]
    plot_gamma(TH,prior,config,params_true[0][2])
    plot_lambda(TH,prior,config,params_true[0][3])
    e_peak_traceplots(TH,prior,config,params_true[0][6])
    plot_ab(TH,prior,config,X,Y,a_compare = params_true[0][4], b_compare= params_true[0][5])
    amplitude_traceplots(TH,prior,config,params_true[0][1])
    t_traceplots(TH,config, params_true[0][0])

    
    ppi(TH,prior,config,X,Y,config.p_thin[-1],np.array(r_thinned))       
    pylab.savefig('ppi.png')

    ci(TH,prior,config,X,Y,config.p_thin[-1],np.array(r_thinned),params_true )       
    pylab.savefig('ci2.png')
    
    STEPS_new = convert_STEPS(STEPS)
    plot_steps(STEPS_new)

    
def plot_rj(TH, config, prior, X,Y,steps, params_true, LLH,CHISQ, r_thinned, p_thin = 1,weights = None):
      ''' TH is only the chain of interest, with thinning p_thin '''
      #summed pulses and individual iterations
      
      newplots = 'plots'
      if not os.path.exists(newplots):
        os.makedirs(newplots)      
      os.chdir(newplots)
      N = len(TH)/2
      n1 = N
      n2 = N + N/5
      n3 = N + N/5*2
      n4 = N + N/5*3
      n5 = N + N/5*4 
      summed_pulses(TH,prior,config,X,Y,config.p_thin[-1])
      pylab.savefig('severaldraws.png')
      individual_pulses(TH[n1],prior,config,X,Y,1000.,config.p_thin[-1])
      print n1, TH[n1]
      pylab.savefig('iter-n1.png')
      individual_pulses(TH[n2],prior,config,X,Y,10000.,config.p_thin[-1])
      print n2, TH[n2]
      pylab.savefig('iter-n2.png')
      individual_pulses(TH[n3],prior,config,X,Y,1000.,config.p_thin[-1])
      print n3, TH[n3]
      pylab.savefig('iter-n3.png')
      individual_pulses(TH[n4],prior,config,X,Y,1000.,config.p_thin[-1])
      print n4, TH[n4]
      pylab.savefig('iter-n4.png')
      individual_pulses(TH[n5],prior,config,X,Y,1000.,config.p_thin[-1])
      pylab.savefig('iter-n5.png')
      print n5, TH[n5]
      
      llh_traceplots(config, np.transpose(LLH)) 
      pylab.savefig('llh.png')      
      chisq_traceplots(config, np.transpose(CHISQ ))
      pylab.savefig('chisq.png')

      
 #     redo_traceplots(config,prior,Xt,Yt,TH,np.transpose(LLH ))    
      
      J = [len(th) for th in TH] 
      J_traceplot(J,config)
      pylab.savefig('JTrace.png')
      hist_J(TH,weights = weights)
      
      
      ppi(TH,prior,config,X,Y,config.p_thin[-1],np.array(r_thinned))       
      pylab.savefig('ppi.png')

      ci(TH,prior,config,X,Y,config.p_thin[-1],np.array(r_thinned),params_true )       
      pylab.savefig('ci2.png')

      n = 4
      w_i = None #initializing weights for unweighted case
      sorted_params,sorted_weights = sort_params(TH,n = n,weights = weights)
      
      V1 = []
      V2 = []
      if params_true is not None:
	if len(params_true) > 1:
	  for th in TH:
	    if len(th) > 0:
		vtemp = np.sort(th[:,1])[::-1]
		if len(vtemp) >1:
		    V1.append(vtemp[0])
		    V2.append(vtemp[1])
	  hscatter(V1,V2,true = np.array([[10.,params_true[1][1]]]) )
        #plt.figure()
	#plt.scatter(V1,V2,color = 'purple',alpha = 0.3)
	#plt.xlabel('V1')
	#plt.ylabel('V2')
	#plt.scatter([10.],[params_true[1][1]],marker = '+', s = 200)
	#plt.show()
	#plt.savefig('v1v2.png')
             
      plot_J_frac(TH)
      
      for i in range(n): #this is how many pulses there are
	  if not os.path.exists(str(i+1)+'pulses'):
	      os.mkdir(str(i+1)+'pulses')
	  os.chdir(str(i+1)+'pulses')
	  th_i = list(sorted_params[i])
	  if sorted_weights is not None:
	      w_i = sorted_weights[i]
	  print np.shape(th_i), 'th_i'
	  if len(th_i) > 1:
	    for j in range(i+1): #now  we look at a single pulse
	      print i,j, np.shape(th_i)
	      if not os.path.exists('pulse'+str(j+1)):
		  os.mkdir('pulse'+str(j+1))
	      os.chdir('pulse'+str(j+1))

	     
                    
	      th_i_j = [ np.array([th_i[ind][j]]) for ind in range(len(th_i))]
	      if params_true is not None:
	          params_true_temp = params_true[min(j,len(params_true)-1),:]
	      else:
	          params_true_temp = [None for index in range(7)]
	      plot_gamma(th_i_j,prior,config, true = params_true_temp[2],weights =w_i)
	      plot_lambda(th_i_j,prior,config, true = params_true_temp[3],weights = w_i)
	      e_peak_traceplots(th_i_j,prior,config, true = params_true_temp[-1],weights =  w_i)
	      plot_ab(th_i_j,prior,config,X,Y, a_compare = params_true_temp[4],weights =  w_i)
	      amplitude_traceplots(th_i_j,prior,config, true = params_true_temp[1],weights =  w_i)
	      t_traceplots(th_i_j,config, true = params_true_temp[0],weights =  w_i)
	      
	      os.chdir('..')
	  os.chdir('..')
      hist_r(r_thinned)
      plot_r(r_thinned,config)
        
        
#print np.shape(Yt)        , np.shape(Yt[-1])

##plot_norj(TH[-1], params_band,directory, config,prior,Xt,Yt[-1],STEPS, p_thin = 1)
#print B, 'B'
#plot_rj(TH[-1], config, prior, Xt,Yt[-1],None, None,r_thinned,p_thin = 1,weights = None)

#for smc
#X = ascii.times
#Y = ascii.counts

#Xt,Yt = code.get_XY(X,Y,-12,25)
#plot_rj(TH, config, prior,Xt,Yt,None, None,r_thinned,p_thin = 1,weights = weights)
