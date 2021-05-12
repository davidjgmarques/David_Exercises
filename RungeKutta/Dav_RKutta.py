## Runge - Kutta
import numpy as np
import math
import matplotlib.pyplot as plt
import sys

#--------------------/ Functions \--------------------

def d_theta(omega):
    return omega

def d_omega(a,s,theta):
    return -(a/s)*math.sin(theta)

def warp(x):
    if x > 2*math.pi:
        x = x - 2*math.pi
    elif x < (-2*math.pi):
        x = x + 2*math.pi
    return x
#This can be substituted with the function 'fmod'


#--------------------/ Constants \--------------------

g = 9.8
L = 0.1
m = 1


#--------------------/ Step \-------------------------

exp = int(sys.argv[1])                       #steps order of magnitude
periods = 10                                #half periods (1 swing)
h = (math.pi)*periods/(10**exp)
t = np.arange(0,periods,h)
print('Number of steps:', len(t))


#--------------------/ Variables \--------------------

theta = []
omega = []
init_theta = math.radians(179)
init_omega = 0
theta.append(init_theta)
omega.append(init_omega)


#--------------------/ RungeKutta \--------------------

for i in range(len(t)-1):

    k1_t = h*d_theta(omega[i])
    k1_w = h*d_omega(g,L,theta[i])
    
    k2_t = h*d_theta(omega[i] + 0.5*k1_w)
    k2_w = h*d_omega(g,L,theta[i] + 0.5*k1_t)

    loop = warp(theta[i])
    theta[i]=loop

    theta.append(theta[i]+k2_t)
    omega.append(omega[i]+k2_w)


#--------------------/ System's Energy \--------------------

Energy = []
Energy0 = []
for i in range(len(t)):
    Energy.append(0.5*m*((omega[i]**2)*(L**2)) + m*g*L*(1-math.cos(theta[i]))) 
    Energy0.append(0.5*m*((omega[0]**2)*(L**2)) + m*g*L*(1-math.cos(theta[0]))) 
    #Line above can be substituted with:' plt.plot(t,np.full(t.size,Energy[0]))


#--------------------/ Plots \--------------------------------

fig, axs = plt.subplots(3)
fig.tight_layout(pad=3.0)
fig.set_figheight(5)
fig.set_figwidth(10)
axs[0].set_title('Theta (t)')
axs[0].plot(t,theta,'b')
axs[1].set_title('Omega (t)')
axs[1].plot(t,omega, 'r')
axs[2].set_title('System\'s energy')
axs[2].plot(t,Energy,'g', label= 'RK Energy')
axs[2].plot(t,Energy0,'k', label='Constant Energy')
plt.legend()
#plt.show()
plt.savefig("plots.pdf")

