

![pipeline status](https://git.differ.nl/imm/div1d/badges/main/pipeline.svg)
![Static Badge](https://img.shields.io/badge/license%20-%20LGPLv3%20-%20blue)


# DIV 1 D 
A 1D time-dependent simulation environment for the scrape-off layer from the x-point to the target, stagnation point to target, double target, extended with reservoirs for atoms, molecules, and core plasma.

#### Copyright notice
CopyRight (c) 2018-2024, Egbert Westerhof, Gijs Derks, Jens Peter Frankemolle, Stijn Kobussen under contract at the Dutch Institute For Fundamental Energy Research.

#### CopyRight Disclaimer
The Dutch Institute For Fundamental Energy Research (DIFFER) hereby disclaims all copyright interest in the program “DIV1D” (which models the scrape-off layer plasma of magnetically confined fusion reactors) written by Egbert Westerhof, Gijs Derks, Jens-Peter Frankemölle, and Stijn Kobussen being under contract at DIFFER.

/signature/ of Martin van Breukelen March 2024

Martin van Breukelen, Managing Director of DIFFER

#### License
This program is free software: you can redistribute it and/or modify it under the terms of the GNU  Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

## Introduction
The heat exhaust in ITER is about 100 MW/m^2 [Chen 2019] while solid materials can cope with heat loads < 15 MW/m^2 [Pitts 2019].
This code allows time-dependent investigations such as the effect from edge-localized modes [Frankemolle 2021].
It can be used to investigate dynamics and compare to measurements such as [Ravensbergen 2020].

## Technologies
* Preferably on the rekenserver at DIFFER
* Fortran90 and cmake compilers
* tickets helpdesk@differ mail
### Eurofusion gateway
* Module load itm-intel/17.0 (for intel compiler)
* Module load imasenv/3.42.0/intel/2023b/1.0 (for IMAS)
* go into makefile -> turn on IMAS compiler flag
* Module show imas ( to see what is loaded )
* tickets: https://jira.eufus.psnc.pl/projects/ACH04SUPP/issues/ACH04SUPP-199?filter=allopenissues
### ITER cluster 
* tbd

## Setup DIV1D
clone to your folder: git clone hash

Go into div1d: cd div1d

Make obj folder: mkdir obj

Use makefile: make 

Copy or make input.txt and div1d_restart_old.txt 

Run a simulation, go to the directory of input.txt and div1d_restart_old.txt and type: ./<dir>/obj/div1d.exe < input.txt

Run a full code integration test: ./test_div1d.sh ->when prompted type 12 and press enter. (this will run integration test 12)

To run unit and partial integration tests  tests: cd div1d, make test

Tip when tests dont compile: rm test/build/CMakeCache.txt

## Status
A place to add ideas for development
### Plasma
- [x] Evaluation of plasma particle, momentum and energy balance including atomic modeling using the AMJUEL database.
- [x] Flux expansion for magnetics and transport.
- [x] Core Scrape-off layer up to flow stagnation point.
- [x] Core particle inventory (e.g. Blanken 2018).
### Neutrals
- [x] Evaluation of neutral particle balance including charge-exchange and ionisation from AMJUEL database.
- [x] neutral model for particles outside plasma
- [x] molecular-proton elastic collisions and dissociation
### Wall
- [x] Recycling of ions.
- [x] Wall association of atoms into molecules
- []  Wall particle inventory (e.g. Kallenbach 2010, Schmidt 2024).
### Inputs
- [x] Input gas puff possible.
- [x] Input model for edge-localized modes.
- [] IMAS I/O with Musscle3
### Outputs
- [x] Output consists of plasma and neutral parameters along the "flux tube" of DIV1D
- [] Synthetic diagnostic output model for impurity emission front based on ADAS.
- [] solutions written to IMAS dictionary (a 1D data structure has been setup)

### Library
- [x] call through memory, using internal DIV1D initialization
- [] call through memory and push grid array per time-step rescaling solutions
- [] accelerate non-particle balances by solving not per cell, but projecting on dominant POD modes

## Version History Short
Although a history is kept by git, this readme also presents an overview of changes.
### Semantics
Version X.Y.Z of DIV1D:
* X is raised when backwards compatibility is lost.
* Y is raised when new features are added with compatibility.
* Z is raised when a new version only contains small changes.
Tag: git tag -a v1.0.0 -m "Version 1.0.0 of DIV1D".
Show: git show v1.0.0
Push: git push origin v1.0.0 (you have to specifically use the tag)
Alternatively use git push origin (upstream) --follow-tags
### Version History
* Version 1.0.0 merged to master by Westerhof and Frankemolle in 2020
* Version 2.0.0 added flux expansion and dynamic inputs in 2021
* Version 2.1.2 neutral background (important corrections for recombination and ionization friction)
* Version 3.0.0 time-dependent impurity concentration with homogeneous profile + time dependent neutral background
* Version 4.0.2 double leg, 5 external interacting neutral reservoirs, realistic geometry (todo time-dependent + pumping)
* Version 6.0.0 time-dependent atom and molecule reservoirs, core, library functionality with time-dependent geometry, wall association

