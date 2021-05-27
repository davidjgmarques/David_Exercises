# Numerical Methods - Exercise 3

Problem: Solve the 1-D diffusion equation with a source term.

Method: Crank-Nicolson

Result: Difference between function evolution calculated and analytic solution.

## Language

Python3

## Libraries necessary

math  
matplotlib  
time  
numpy  
scipy.linalg

## Usage

Windows 10

```bash
py ./Dav-CN.py
```

Unix/MacOS

```bash
python ./Dav-CN.py
```

### Python (conflicting) versions:

If several versions of python are installed, the initialization is made with 'python3' instead of 'python'.

## Contributing

Pull requests are welcome, as well as verbal comments & tips.

# Solutions

Once the simulation is done, two output files are generated.  
The first, "Results", presents the evolution of the function in the first seconds and the difference between the result obtained and analytical solution, for different sampling.  
The second, "RMS_error_comparison", shows the difference referred above, for different sampling, in the same plot for easier comparison.

In terminal, the compilation time is shown as well as the limit ratio between the differences referred, between different samplings.

## Solving tridiagonal matrix equations

The package used to solve these matrixes was *solve_banded* from the *scipy.linalg* library.  
The reference guide to use this package can be found in [here](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.solve_banded.html).  

