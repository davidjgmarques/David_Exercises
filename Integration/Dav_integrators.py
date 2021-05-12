import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.integrate import simps

# ---------------/ Function to integrate \------------------------
def f1(x):
    return x**4 - 2*x + 1

def f1_1p(x):
    return 4*x**3 - 2    

def f1_3p(x):
    return 24*x


# ---------------/ Trapezium integrator \------------------------
def Trapezium(func,low_lim,high_lim,steps):
    a=low_lim
    b=high_lim
    N=steps
    h = (b-a)/N

    I_trap = (0.5*func(a) + 0.5*func(b))

    for k in range(1,N):
        I_trap += func(a+k*h) 
    I_trap *= h
    return(I_trap)

# ---------------/ Simpson integrator \------------------------
def Simpson(func,low_lim,high_lim,steps):
    a=low_lim
    b=high_lim
    N=steps
    h = (b-a)/N

    I_simp = func(a)+ func(b)

    for k in range(1,steps):
        if k % 2 == 0:               #even
            I_simp += 2*func(a+k*h) 
        else:                       #odd
            I_simp += 4*func(a+k*h) 
    I_simp *= (h/3)
    return(I_simp)

# ---------------/ Simpson with sums (no loops) \-------------------------

# x = np.linspace(a,b,n_steps+1)
# func= f1(x)
# my_A_simp1 = h/3 * (func[0] + func[n_steps] + 4*sum(func[1:n_steps:2]) + 2*sum(func[2:n_steps-1:2]))
# print('\nMy simpson integral_1: \n',my_A_simp1)


#----------------/ -------------- \--------------------------


#----------------/ Quick approach \--------------------------
a=0
b=2
exponent = 7
n_steps = 10**exponent
Integral_teo = 4.4

print('Trapezium method for %2.0e samples: ' % n_steps, Trapezium(f1,a,b,n_steps))
print('Simpson method for %2.0e samples: ' % n_steps, Simpson(f1,a,b,n_steps))  

x = np.linspace(a,b,n_steps+1)
python_A_simp= simps(f1(x),x)
print('Python\'s simpson integral for %2.0e samples: ' % n_steps, python_A_simp)

print(
'\nAnalytic method: \n'
' Int(x^4 - 2x + 1) in x{0,2} (=) \n'
' (=) x^5/5 - x^2 + x in x{0,2} (=) \n'
' (=) 32/5 - 4 + 2 - (0) (=) \n'
' (=) 6.4 - 4 + 2 (=) \n'
' (=) 4.4')

plot10= plt.figure(7)
plt.plot(x,f1(x),'-r', label=r'$y = x^4 - 2x + 1$')
plt.grid()
plt.axvline(c='k')
plt.axhline(c='k')
plt.title('Polynomial Curve')
plt.legend(loc='upper left')


#----------------/ Errors \--------------------------


S_error_approx = []
S_error_diff = []
S_error_teo = []

T_error_approx = []
T_error_diff = []
T_error_teo = []

N_iterative=[]

for g in range(1,exponent+1):
    exp=10**g
    h = (b-a)/exp

    N_iterative.append(exp)

    Simp_1n = Simpson(f1,a,b,exp)
    Simp_2n = Simpson(f1,a,b,exp*2)

    S_error_approx.append((1/15)*(abs(Simp_2n-Simp_1n)))
    S_error_diff.append(abs(Integral_teo - Simp_1n))
    S_error_teo.append((1/90)*h**4*(abs(f1_3p(a) - f1_3p(b))))

    Trap_1n = Trapezium(f1,a,b,exp)
    Trap_2n = Trapezium(f1,a,b,exp*2)

    T_error_approx.append((1/3)*(abs(Trap_2n-Trap_1n)))
    T_error_diff.append(abs(Integral_teo - Trap_1n))
    T_error_teo.append((1/12)*h**2*(abs(f1_1p(a) - f1_1p(b))))

plot1= plt.figure(1)
plt.xlabel('N')
plt.grid()
plt.xscale('log')
plt.yscale('log')
plt.ylabel('(1/15) * (I$_{2N}$ - I$_{N}$)')
plt.title('Simpson - Approximated error')
plt.plot(N_iterative,S_error_approx,'o')

plt.savefig("Simpson - Approximated error.pdf")

plot2= plt.figure(2)
plt.xlabel('N')
plt.grid()
plt.xscale('log')
plt.yscale('log')
plt.ylabel('I$_{teo}$ - I$_{N}$')
plt.title('Simpson - Difference between theoretical and computed value')
plt.plot(N_iterative,S_error_diff,'o')

plt.savefig("Simpson - Difference between theoretical and computed value.pdf")

plot3= plt.figure(3)
plt.xlabel('N')
plt.grid()
plt.xscale('log')
plt.yscale('log')
plt.ylabel('(1/90) * h$^{4} * (f\'\'\'(a) - f\'\'\'(b))$')
plt.title('Simpson - Theoretical error')
plt.plot(N_iterative,S_error_teo,'o')

plt.savefig("Simpson - Theoretical error.pdf")

plot4= plt.figure(4)
plt.xlabel('N')
plt.grid()
plt.xscale('log')
plt.yscale('log')
plt.ylabel('(1/3) * (I$_{2N}$ - I$_{N}$)')
plt.title('Trapezium - Approximated error')
plt.plot(N_iterative,T_error_approx,'o')

plt.savefig("Trapezium - Approximated error.pdf")

plot5= plt.figure(5)
plt.xlabel('N')
plt.grid()
plt.xscale('log')
plt.yscale('log')
plt.ylabel('I$_{teo}$ - I$_{N}$')
plt.title('Trapezium - Difference between theoretical and computed value')
plt.plot(N_iterative,T_error_diff,'o')

plt.savefig("Trapezium - Difference between theoretical and computed value.pdf")

plot6= plt.figure(6)
plt.xlabel('N')
plt.grid()
plt.xscale('log')
plt.yscale('log')
plt.ylabel('(1/12) * h$^{2} * (f\'(a) - f\'(b))$')
plt.title('Trapezium - Theoretical error')
plt.plot(N_iterative,T_error_teo,'o')

plt.savefig("Trapezium - Theoretical error.pdf")

#plt.show()