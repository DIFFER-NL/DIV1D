# DIV 1 D 
A 1D time-dependent simulation environment for the scrape-off layer from the x-point to the target.
## Introduction
The heat exhaust in ITER is about 100 MW/m^2 [Chen 2019] while solid materials can cope with heat loads < 15 MW/m^2 [Pitts 2019].
This code allows time-dependent investigations such as the effect from edge-localized modes [Frankemolle 2021].
It can be used to investigate dynamics and compare to measurements such as [Ravensbergen 2020].
## Technologies
* Preferably on the rekenserver at DIFFER
* Fortran90 and cmake compilers
## Setup
clone to your folder: git clone hash

Go into div1d: cd div1d

Make obj folder: mkdir obj

Use makefile: make 

Copy or make input.txt and div1d_restart_old.txt 

Run a simulation, go to the directory of input.txt and div1d_restart_old.txt and type: ./<dir>/obj/div1d.exe < input.txt

## Status
A place to add ideas for development
### Plasma
- [x] Evaluation of plasma particle, momentum and energy balance including atomic modeling using the AMJUEL database.
- [x] Flux expansion for magnetics and transport.
- []  Core particle inventory (e.g. Blanken 2018).
### Neutrals
- [x] Evaluation of neutral particle balance including charge-exchange and ionisation from AMJUEL database.
- []  neutral model for particles outside plasma
- []  molecular-proton elastic collisions and dissociation
### Wall
- [x] Recycling of ions.
- []  Wall particle inventory (e.g. Kallenbach 2010).
### Inputs
- [x] Input gas puff possible.
- [x] Input model for edge-localized modes.
### Outputs
- [x] Output consists of plasma and neutral parameters along the "flux tube" of DIV1D
- [] Synthetic diagnostic output model for impurity emission front based on ADAS.

## Version History Short
Although a history is kept by git, this readme also presents an overview of changes.
### Semantics
Version X.Y.Z of DIV1D:
* X is raised when backwards compatibility is lost.
* Y is raised when new features are added with compatibility.
* Z is raised when a new version only contains small changes.
Tag: git tag -a v1.0.0 -m "Version 1.0.0 of DIV1D".
Show: git show v1.0.0
Push: git push origin v1.0.0
### Version History
* Version 1.0.0 merged to master by Westerhof and Frankemolle in 2020
* Version 2.0.0 added flux expansion and dynamic inputs in 2021

