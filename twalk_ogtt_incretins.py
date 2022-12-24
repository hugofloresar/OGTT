import scipy.stats as ss
from xlrd import open_workbook
from ogtt_Erlang_ode_incretins import rhs
import numpy as np
import scipy as sp
import pylab as pl
import pytwalk
import sys
import os
from tempfile import TemporaryFile
from scipy import integrate, optimize


if not os.path.exists('chain_pat'+str(p_number)+'_incretins'):
    os.makedirs('chain_pat'+str(p_number)+'_incretins')


book_ogtt = open_workbook('Datos_OGTT.xlsx')
sheet_ogtt = book_ogtt.sheet_by_index(0)


if p_number == 0:
    noisy_data = np.array([90.0,125.0,65.0,115.0,130.0])
else:
    noisy_data = sheet_ogtt.row_values(p_number,4,9)

t_data = np.array([0.0,0.5,1.0,1.5,2.0])
time2 = np.linspace(0.0,2.0,1000)

p5 = 60.0/31.0
p7 = 60.0/31.0

m = 2  ## glucagon
n = 2  ## insulin
s = 2  ## digestive

def soln(p,m,n,s):
    x0 = np.zeros(m+n+s+1)
    x0[0] = noisy_data[0]
    x0[n + m + 1] = 400
    return integrate.odeint(rhs,x0,t_data,args=(m,n,s,p,))

def solen(p,m,n,s):
    x0 = np.zeros(m+n+s+1)
    x0[0] = noisy_data[0]
    x0[n + m + 1] = 400
    return integrate.odeint(rhs,x0,time2,args=(m,n,s,p,))


# gamma distribution for parameters
k0 = 2.0
theta0 = 1.0
k1 = 10.0
theta1 = 1.0
k2 = 10.0
theta2 = 1.0
k3 = 90**2/20
theta3 = 20.0/90
var=5.0**2
k4 = 10.0
theta4 = 1.0


def energy(p): # -log of the posterior
    my_soln = soln(p, m, n, s)

    # uncomment the desired likelihood
    log_likelihood = -0.5*(np.linalg.norm(noisy_data - my_soln[:,0]))**2/var # Gaussian


    a0 = (k0-1)*np.log(p[0])- (p[0]/theta0)
    a1 = (k1-1)*np.log(p[1])- (p[1]/theta1)
    a2 = (k2-1)*np.log(p[2])- (p[2]/theta2)
    a3 = (k3-1)*np.log(p[3])- (p[3]/theta3)
    a4 = (k4-1)*np.log(p[4])- (p[4]/theta4)

    log_prior = a3 + a1 + a0 + a2 + a4

    return -log_likelihood - log_prior

def support(p):
    rt = True
    rt &= 0.5<p[0]<4.0
    rt &= 0.0<p[1]<100.0
    rt &= 0.0<p[2]<100.0
    rt &= 60.0<p[3]<200.0
    rt &= 0.0<p[4]<100.0
    return rt

def init():
    p = np.zeros(5)
    p[0] = np.random.uniform(0.5,2.0)
    p[1] = np.random.gamma(k1, theta1)
    p[2] = np.random.gamma(k2, theta2)
    p[3] = np.random.gamma(k3, theta3)
    p[4] = np.random.gamma(k4, theta4)
    return p

burnin = 2000

T=10000
ogtt = pytwalk.pytwalk(n=5,U=energy,Supp=support)
y0=init()
yp0=init()
ogtt.Run(T,y0,yp0)

chain = ogtt.Output
energy = chain[:,-1]

pl.figure()
pl.plot(range(T+1),energy)
pl.savefig('chain_pat'+str(p_number)+'_incretins/energy.png')
pl.close()

pl.figure()
pl.plot(range(T+1-burnin),energy[burnin:])
pl.savefig('chain_pat'+str(p_number)+'_incretins/energy_burn.png')
pl.close()

cadena=TemporaryFile()
np.save('chain_pat'+str(p_number)+'_incretins/cadena',ogtt.Output)

pl.figure()
pl.hist(ogtt.Output[burnin:,0],bins=30)
#pl.title('Theta0 Patient '+str(p_number))
pl.savefig('chain_pat'+str(p_number)+'_incretins/Hist0.png')
pl.close()

pl.figure()
pl.hist(ogtt.Output[burnin:,1],bins=30)
#pl.title('Theta1 Patient '+str(p_number))
pl.savefig('chain_pat'+str(p_number)+'_incretins/Hist1.png')
pl.close()

pl.figure()
pl.hist(ogtt.Output[burnin:,2],bins=30)
#pl.title('Theta2 Patient '+str(p_number))
pl.savefig('chain_pat'+str(p_number)+'_incretins/Hist2.png')
pl.close()

pl.figure()
pl.hist(ogtt.Output[burnin:,3],bins=30)
#pl.title('G0 Patient '+str(p_number))
pl.savefig('chain_pat'+str(p_number)+'_incretins/Hist3.png')
pl.close()

pl.figure()
pl.hist(ogtt.Output[burnin:,4],bins=30)
#pl.title('G0 Patient '+str(p_number))
pl.savefig('chain_pat'+str(p_number)+'_incretins/Hist4.png')
pl.close()

###############################################
#### Computing the MAP estimate
energy_MAP = min(energy)
loc_MAP = np.where(energy==energy_MAP)[0]
MAP = chain[loc_MAP[-1]]
MAP = MAP[:-1]

np.save('chain_pat'+str(p_number)+'_incretins/MAP',MAP)

my_solnMAP = solen(MAP,m,n,s)
glucose_MAP = my_solnMAP[:,0]
insulin_MAP = my_solnMAP[:,n]
glucagon_MAP = my_solnMAP[:,n + m]
digestive_MAP = my_solnMAP[:,n + m + s]

#### Computing the posterior mean
Post_mean = np.ones(5)
Post_mean[0] = np.mean(chain[burnin:,0])
Post_mean[1] = np.mean(chain[burnin:,1])
Post_mean[2] = np.mean(chain[burnin:,2])
Post_mean[3] = np.mean(chain[burnin:,3])
Post_mean[4] = np.mean(chain[burnin:,4])

np.save('chain_pat'+str(p_number)+'_incretins/CM',Post_mean)

my_solnCM = solen(Post_mean,m,n,s)
glucose_CM = my_solnCM[:,0]
insulin_CM = my_solnCM[:,n]
glucagon_CM = my_solnCM[:,n + m]
digestive_CM = my_solnCM[:,n + m + s]

median = np.median(chain[burnin:],axis=0)[:-1]

np.save('chain_pat'+str(p_number)+'_incretins/median',median)

my_solnMedian = solen(median,m,n,s)
glucose_Median = my_solnMedian[:,0]
insulin_Median = my_solnMedian[:,n]
glucagon_Median = my_solnMedian[:,n + m]
digestive_Median = my_solnMedian[:,n + m + s]

##############################################

pl.figure()
pl.plot(t_data,noisy_data,'ro',label='data')
pl.plot(time2,glucose_CM,'g',label='CM')
pl.plot(time2,glucose_MAP,'b',label='MAP')
pl.plot(time2,glucose_Median,'m',label='Median')
pl.xlabel('Time', size =14)
pl.ylabel('Blood Glucose level', size =14)
pl.legend(loc=0, shadow=True)
pl.savefig('chain_pat'+str(p_number)+'_incretins/glucose.png')
pl.close()

pl.figure()
pl.plot(time2,digestive_CM,'g',label='CM')
pl.plot(time2,digestive_MAP,'b',label='MAP')
pl.plot(time2,digestive_Median,'m',label='Median')
pl.xlabel('Time', size =14)
pl.legend(loc=0, shadow=True)
pl.savefig('chain_pat'+str(p_number)+'_incretins/digestive.png')
pl.close()

pl.figure()
pl.plot(time2,insulin_CM,'g',label='CM')
pl.plot(time2,insulin_MAP,'b',label='MAP')
pl.plot(time2,insulin_Median,'m',label='Median')
pl.xlabel('Time', size =14)
pl.legend(loc=0, shadow=True)
pl.savefig('chain_pat'+str(p_number)+'_incretins/insulin.png')
pl.close()

pl.figure()
pl.plot(time2,glucagon_CM,'g',label='CM')
pl.plot(time2,glucagon_MAP,'b',label='MAP')
pl.plot(time2,glucagon_Median,'m',label='Median')
pl.xlabel('Time', size =14)
pl.legend(loc=0, shadow=True)
pl.savefig('chain_pat'+str(p_number)+'_incretins/glucagon.png')
pl.close()

pl.figure()
for k in np.arange(200):
    qq = chain[-1-10*k][:-1]
    my_soln = solen(qq,m,n,s)
    pl.plot(time2,my_soln[:,0],'k',color='0.75')
pl.plot(time2,glucose_Median,'m',label='Median')
pl.plot(t_data,noisy_data,'ro',label='data')
pl.xlabel('Time', size =14)
pl.ylabel('Blood Glucose level', size =14)
pl.legend(loc=0, shadow=True)
pl.savefig('chain_pat'+str(p_number)+'_incretins/Fit_glucoseUQ.png')
pl.close()

iat = ogtt.Ana()

IAT = np.array(iat)[0]
np.save('chain_pat'+str(p_number)+'_incretins/IAT',IAT)
