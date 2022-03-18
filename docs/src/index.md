# Cuas.jl Documentation

*The anti aerial package.*

## Table of contents


```@contents
Pages = ["index.md"]
```
## Introduction
Cuas is a anti aerial simulation package designed to compute efficiency of a weapon system against an aerial target.
It offers the following scenarios: stationary vs stationary, stationary vs moving, moving vs stationary and moving vs moving.
For one on one, one on many, many on one and many on many configurations.

```@repl
using Cuas
```


## Package Features
- Hit probability
- Fire Control System

### Fire Control System
```
acquisition(sim,target, tank, proj, aero,fixedError,varibleErrorIn,variableErrorOut,randomError,zones, shapes; optional arguments)
```
    Is a simJulia process which simulate the firing chain process for an ABM projectile. From target acquisition to damage assessment.
    Input parameters are :
    * sim is the simulation environment
    * target is a target object which contains target information
    * tank is a tank object (WeaponSystems) 
    * proj contains the projectile informaiton
    * aero contains the aerodynamic coefficient
    * fixedError is a dictionary which contains the fixed bias (in mils) measured in a plane perpendicular to the LOS
    * varibleErrorIn is a dictionary which contains the standard deviation of measured quantities that are propagated via Monte Carlo
    * variableErrorOut is a dictionary which contains variable bias (in mils) measured in a plane perpendicular to the LOS 
    * randomError is a dictionary which contains random errors (in mils) measured in a plane perpendicular to the LOS
    * zones is an object of ZoneList (ExternalBallistics) which contains informaition on generated fragments
    * shapes is a vector which contains the shapes of the fragments

    Optional arguments are:
    * Δt_fireCommand is the time delay for issuing the fire command default value is 1.0
    * Δ_slew is the time delay to lay the canon in a desired position. The default value is 2.0
    * Δ_aim is the time delay to aim at the target. The default value is  2.0
    * Δ_range is the time delay to get the target range. The default value is 0
    * Δ_fireControl is the time delay to get the ballistic correction. The default value is 1.0
    * Δ_fireDemand is the time delay for the fire demand. The default value is 0.0
    * Δ_fireDelay is the time delay for the fire. The default value is  0.0
    * Δ_shotExit is the time delay for the shot to exit. The default value is 1.0
    * w is the wind (cross and range). The default value is  Wind(0.0,0.0)
    * deltaR defines the detonation distance. It expressed in projectile revolution. To the number of revolutions achieved 
    to hit the target, deltaR is soustracted. The default value is 5
    * N is the number of Monte Carlo simulations. The defaults is 100
    * p defines the number in % of fragment that hit target to be considered a hit
### Function documentation

```@autodocs
Modules = [Cuas]
Order = [:function]
```

### Type documentation

```@autodocs
Modules = [Cuas]
Order = [:type]
```


```@index
```
