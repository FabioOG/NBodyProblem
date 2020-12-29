# NBodyProblem

This program contains a few integrators, able to compute more or less precisely the solution of an N-Body problem. 

### Program arguments

The program arguments should be :

```nBodyProblem [options] initFile [outputFile]```

The `initFile` contains the presets of the simultation, it has to be specified. The `outputFile` will contains the trajectories of the bodies, its default value is `out.csv`

### Program options

The program options are the following :

* `-euler` Uses the explicit Euler integrator
* `-taylor` Uses a 4-order Taylor serie
* `-symplectic4` Uses a 4-order symplectic integrator
* `-symplectic6` Uses a 6-order symplectic integrator
* `-energy` Corrects the energy of the system

The default integrator running is the classic 4-order Runge-Kutta method

### Initialisation file format

The file containing the simulation informations should have the following format :

```
bodyNumber,computationStep,computationTime[,outputStepRatio,G]
x1,y1,z1,vx1,vy1,vz1,m1
...
```

Where `outputStepRatio` indicates the ratio between the outputed step number and the computed step number. One step out of this much will be recorded in the file. The default value of this is 100. `G` is the gravitational constant used, its default value is 1, adjust depending on the units used. These two arguments are then optional. Try to choose the units so the system dimensions are as close as possible to 10^0.

### Trajectories file format

The output file will have the following format :

```
bodyNumber,stepNumber,maxPosition,maxSpeed
t,x1,y1,z1,vx1,vy1,vz1 ... ,xn,yn,zn,vxn,vyn,vzn
...
```

`maxPosition` and `maxSpeed` gives the maximum intensity of each vectors, it can be used to quickly determine the system dimensions

_______________

## Authors :

* Fabio d'ORTOLI-GALERNEAU
