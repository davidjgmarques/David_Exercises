# Numerical Methods - Exercise 2

Implementation of Runge Kutta method to solve differential equations.

Problem: Simple Pendullum Oscillation

## Language

Python3

## Libraries necessary

numpy\ 
matplotlib\
math

## Usage

The script requires an argument: the exponential(**exp**) in **steps = 10^{exp}**.  
This number is used to calculate the **step size** used in the Runge Kutta.  
The *size* and *number* of steps is calculated for a given number of half-periods, here set to 10.

The user is recommended to use **exp=3** and **exp=4**.

Windows 10

```bash
py ./Dav_RKutta.py [exp]
```

Unix/MacOS

```bash
python ./Dav_RKutta.py [exp]
```

### Python (conflicting) versions:

If several versions of python are installed, the initialization is made with 'python3' instead of 'python'.

## Contributing

Pull requests are welcome, as well as verbal comments & tips.

# Solutions

When using **exp=3**:

The user will note that the pendullum's theta angle ($\theta$) will be greater than 2$\pi$.  
This happens due to the low number of steps used to sample $\theta (t)$, which generates greater variations in $\theta (t)$.
At this point, $\theta$ is 'warped' back. This results in a $\theta (t)$ with several sharp rises.  
The energy of the system is slowly increasing thus proving the *non energy conservation* of the Runge Kutta method. 

When using **exp=4**:
The sampling is good enough to avoid great variations of ($\theta$).  
This results in a correct behaviour of the pendullum, *i.e.*, a variation of ($\theta$) between -90 and 90 degrees.  
The energy system also increases in this case, but much slower, due to the larger number of samples.  
Using **exp>4** has the same result, with the difference of the energy increase being even smaller.