import numpy as np
import scipy as sp
import pylab as pl
import sys
import os
from scipy import integrate, optimize


t_data = np.array([0.0,0.5,1.0,1.5,2.0])
time2 = np.linspace(0.0,2.0,1000)

p5 = 60.0/31.0
p7 = 60.0/31.0

### m is the box number for the glucagon comp
### n is the box number for the insulin comp
### s is the box number for the digestic systems

def rhs(x, t, m, n, s, p):
    fx = np.zeros(m + n + s + 1)
    fx[0] = x[n + m] - x[n] + p[0] * x[n + m + s]

    ### insulin box
    fx[1] = p[1] * np.max([x[0]-p[3],0]) + p[4] * x[n + m + s] - n * p5 * x[1]
    if n > 1:
        for i in range(n-1):
            fx[2 + i] = n * p5 * x[2 + i -1] -n * p5 * x[2 + i]

    ### glucagon box
    fx[n+1] = p[2] * np.max([p[3]-x[0],0]) - m * p7 * x[n+1]
    if m > 1:
        for i in range(m-1):
            fx[n + 2 + i] = m * p7 * x[n + 2 + i -1] - m * p7 * x[n + 2 + i]

    fx[n + m + 1] = - s * p[0]* x[n + m + 1]
    if s > 1:
        for i in range(s-1):
            fx[n + m + 2 + i] = s * p[0]* x[n + m + 2 + i -1] - s * p[0]* x[n + m + 2 + i]
    return fx

#x0 = np.array([90,0.0,0.0,300.0])
p1 = np.array([0.8,8.0,10.0,90.0,1.0])  ### p1
p2 = np.array([0.8,20.0,35.0,90.0,3.0])  ### p1
# p = np.array([30.0,30.0,100.0,100.0,30.0])  ### p2
#p = np.array([0.77,1.87,0.01,100.0])  ### p3

def soln(p,m,x0):
    x0 = np.concatenate((x0,np.zeros(m-1)))
    return integrate.odeint(rhs,x0,t_data,args=(m,p,))

def solen(p,m,n,s):
    x0 = np.zeros(m+n+s+1)
    x0[0] = 90.0
    x0[n + m + 1] = 300.0
    return integrate.odeint(rhs,x0,time2,args=(m,n,s,p,))

# m = 2
# n = 2
# s = 2
# my_soln = solen(p1,m,n,s)
# glucose1 = my_soln[:,0]
# insulin1 = my_soln[:,m]
# glucagon1 = my_soln[:,m+n]
# digestive1 = my_soln[:,m+n+s]

# m = 2
# n = 2
# s = 2
# my_soln = solen(p2,m,n,s)
# glucose2 = my_soln[:,0]
# insulin2 = my_soln[:,m]
# glucagon2 = my_soln[:,m+n]
# digestive2 = my_soln[:,m+n+s]

# m = 3
# n = 3
# s = 3
# my_soln = solen(p,m,n,s)
# glucose3 = my_soln[:,0]
# insulin3 = my_soln[:,m]
# glucagon3 = my_soln[:,m+n]
# digestive3 = my_soln[:,m+n+s]

# pl.figure()
# pl.plot(time2,glucose1, label='p1')
# pl.plot(time2,glucose2,'g', label='p2')
# pl.xlabel('Time (minutes)', size =14)
# pl.ylabel('Blood Glucose Level (mg/dL)', size =14)
# pl.xticks([0.0,0.5,1.0,1.5,2.0],[0,30,60,90,120])
# # pl.plot(time2,glucose3, label='m=3')
# #pl.title('Theta0 Patient '+str(p_number))
# pl.legend()
# pl.savefig('2scenario.png')
# pl.show()

# pl.figure()
# pl.plot(time2,insulin1, label='m=1')
# pl.plot(time2,insulin2, label='m=2')
# pl.plot(time2,insulin3, label='m=3')
# #pl.title('Theta0 Patient '+str(p_number))
# pl.legend()
# #pl.savefig('insulin_all_p2.png')
# pl.close()

# pl.figure()
# pl.plot(time2,glucagon1, label='m=1')
# pl.plot(time2,glucagon2, label='m=2')
# pl.plot(time2,glucagon3, label='m=3')
# #pl.title('Theta0 Patient '+str(p_number))
# pl.legend()
# #pl.savefig('glucagon_all_p2.png')
# pl.close()

# pl.figure()
# pl.plot(time2,digestive1, label='m=1')
# pl.plot(time2,digestive2, label='m=2')
# pl.plot(time2,digestive3, label='m=3')
# #pl.title('Theta0 Patient '+str(p_number))
# pl.legend()
# #pl.savefig('digestive_all_p2.png')
# pl.close()
