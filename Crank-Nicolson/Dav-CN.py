# 1D - Diffusion with source term

import math
import numpy as np
from scipy.linalg import solve_banded
import matplotlib.pyplot as plt
import time

start = time.time()

def analy_sol(x):
    return math.cos(math.pi * 0.5 * x)

def sourc_term(x,diff):
    return math.cos(math.pi * 0.5 * x) * pow(math.pi,2) * 0.25 * diff

H = 1
xMin = -H
xMax = H
D = 0.1
dt = pow(10,-3)

exps=[8,9,10]
N1 = pow(2,exps[0])
N2 = pow(2,exps[1])
N3 = pow(2,exps[2])
N=[N1,N2,N3]

rms= []

tmax= 100
time_plot=range(1,tmax+1)
time_evolution=[99,98,97,96,95,94]


fig, axs = plt.subplots(2,3)
fig.set_figheight(10)
fig.set_figwidth(15)
fig.suptitle('Evolution of $u(x)$ and RMS error for different sampling', fontsize = 16)
for k in range(len(N)):

    x = np.linspace(xMin,xMax,N[k])
    dx = (xMax-xMin)/(N[k]-1)
    #dx = x[1]-x[0]

    Q = [sourc_term(x[i],D) for i in range(N[k])]
    Q[0]=0
    Q[-1]=0

    S = [analy_sol(x[i]) for i in range(N[k])]
    S[0]=0
    S[-1]=0

    u = np.zeros(N[k])
    rhs = np.zeros(N[k]-2)

    alpha = (dt * D) / (2 * pow(dx,2))
    central_diag = (1+2*alpha)
    upper_diag = -alpha
    lower_diag = -alpha

    error = 0
    difference = []
    
    evolution_u= np.zeros((len(time_evolution),N[k]))

    error_counter=0
    plot_counter=0

    t=tmax
    while t > dt:
        for i in range(1,N[k]-1):
            rhs[i-1] = alpha * u[i+1] + (1-2*alpha) * u[i] + alpha * u[i-1] + Q[i] * dt 

        Ab = np.zeros((3,N[k]-2))

        Ab[0,1:] = upper_diag
        Ab[1,0:] = central_diag
        Ab[2,:-1] = lower_diag

        sol = solve_banded((1,1),Ab,rhs)

        for i in range(1,N[k]-1):
            u[i]=sol[i-1]

        t -= dt

        if error_counter == 1000-1:
            for i in range(N[k]):
                error += pow(u[i] - S[i],2)
            difference.append(math.sqrt(error/(N[k])))
            error = 0
            error_counter=0

        error_counter+=1

        point_plot = round(t,3)
        if point_plot in time_evolution:
            for i in range(N[k]-1):
                evolution_u[plot_counter,i]=u[i]
            plot_counter+=1    


    rms.append(difference)


    axs[0,k].set_title('N = 2$^{%i}$' %(exps[k]))
    axs[0,k].plot(evolution_u.T)
    axs[0,k].set_xlabel('x')
    axs[0,k].set_ylabel('$u(x)$')

    axs[1,k].set_ylabel('$||f(x) - f_{a}(x)||$')
    axs[1,k].set_yscale('log')
    axs[1,k].set_xlabel('t')
    axs[1,k].plot(time_plot,difference,'bo')


end = time.time()

print('Time elapsed:', end - start)

plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.3, hspace=None)
plt.savefig("Results.pdf")

plot2= plt.figure(2)
plot2.set_figheight(5)
plot2.set_figwidth(10)
plt.yscale('log')
plt.ylabel('$||f(x) - f_{a}(x)||$')
plt.xlabel('t')
plt.title('RMS comparison for different sampling')
for i in range(len(rms)):
    plt.plot(time_plot,[pt for pt in rms[i]],'o',label = '2$^{%i}$' %(exps[i]))
plt.legend()
plt.savefig("RMS_error_comparison.pdf")
#plt.show()

diff_highN_medN = rms[2][99]/rms[1][99]
diff_medN_lowN = rms[1][99]/rms[0][99]

print('When t goes to infinite, error(N=2^9)/error(N=2^8) = {:.2E}, \n' 
'and error(N=2^8)/error(N=2^7) = {:.2E}'.format(diff_highN_medN,diff_medN_lowN))