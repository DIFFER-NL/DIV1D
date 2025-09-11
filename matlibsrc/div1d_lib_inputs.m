function [out] = div1d_lib_inputs()

out=struct;
floatinphys(1:141) = 0.0;
intinphys(1:31) = 0;
floatinnum(1:90) = 0.0;
intinnum(1:30)  = 0;


% set values
intinnum(0+1) = 200; % Nx
intinnum(1+1) = 10; % nout
intinnum(2+1) = 227; % method
intinnum(3+1) = 10000; % itate_mod
intinnum(4+1) = 10000; % max_step
intinnum(5+1) = 100; % max attempts
intinnum(6+1) = -100; % nzswag
intinnum(7+1) = 1; % evolve density
intinnum(8+1) = 1; % evolve momentum
intinnum(9+1) = 1; % evolve energy
intinnum(10+1) = 1; % evolve neutral
intinnum(11+1) = 1; % evolve neutral momentum
intinnum(12+1) = 1; % evolve molecule
intinnum(13+1) = 1; % evolve background
intinnum(14+1) = 1; % evolve core
intinnum(15+1) = 0; % evolve core neutral
intinnum(16+1) = 1; % mol dens model
intinnum(17+1) = 1; % D_new harmonic average of temperature for diffusion
% norms
floatinnum(0+1) = 1.0*10^19 ; % density norm
floatinnum(1+1) = 1.0; % temperature norm
floatinnum(2+1) = 1.0*10^4 ; % velocity norm
floatinnum(3+1) = 1.0; % neutral norm
% equation switch
floatinnum(9+1) = 1.0; %  density source
floatinnum(10+1) = 1.0; % momentum source
floatinnum(11+1) = 1.0; % energy source
floatinnum(12+1) = 1.0; % neutral source
floatinnum(13+1) = 1.0; % neutral momentum source
floatinnum(14+1) = 1.0; % molecule source
% switches energy terms
floatinnum(19+1) = 1.0; % conv heat
floatinnum(20+1) = 1.0; % energy flux
floatinnum(21+1) = 1.0; % energy compr
floatinnum(22+1) = 0.0; % only core source Q
% switches atomic reaction rates + impurity
floatinnum(29+1) = 1.0; % cx
floatinnum(30+1) = 1.0; % rec
floatinnum(31+1) = 1.0; % rec ene
floatinnum(32+1) = 1.0; % ion 
floatinnum(33+1) = 1.0; % exc 
floatinnum(34+1) = 1.0; % impurity  
floatinnum(35+1) = 1.0; % mom atoms
floatinnum(36+1) = 1.0; % mom molecules
% switches molecular rates
floatinnum(39+1) = 1.0; %  diss mol
floatinnum(40+1) = 1.0; %  ion mol
floatinnum(41+1) = 1.0; %  cx mol
floatinnum(42+1) = 1.0; %  diss h2pl 
floatinnum(43+1) = 1.0; %  dis rec h2pl
floatinnum(44+1) = 1.0; %  dis ion h2pl
floatinnum(45+1) = 1.0; %  energy molecules
floatinnum(46+1) = 1.0; %  diss att
floatinnum(47+1) = 1.0; %  cx hmin
floatinnum(48+1) = 1.0; % ion hmin
% numerical settings
floatinnum(59+1) = 0.05; % dxmin
floatinnum(60+1) =10^-5 ; %  delta t
floatinnum(61+1) = 10^-6 ; %  abstol
floatinnum(62+1) = 10^-6 ; %  reltol
floatinnum(63+1) = 5; % vis
floatinnum(64+1) = 0.5; % cent diff
floatinnum(65+1) = 1.0; % lax switch


 % physics settings
% phys integers
intinphys(0+1) = 5  ;% num imp
intinphys(1+1) = 6  ;% imp1
intinphys(2+1) = 7 ;% imp2
intinphys(3+1) = 0 ;% imp3
intinphys(4+1) = 0 ;% imp4
intinphys(5+1) = 0 ;% imp5
intinphys(6+1) = 0  ;% switch imp dist
intinphys(7+1) = 0 ;% widt pfr

 % phys floats
% geometry, grid distributions
floatinphys(0+1) = 20.0 ;% L
floatinphys(1+1) = 10.0 ;% L cor sol
floatinphys(2+1) = 0.0 ;% x cor sol
floatinphys(3+1) = 0.1 ;% sintheta
floatinphys(4+1) = 90.0 ;% pol tar ang 1
floatinphys(5+1) = 90.0 ;% pol tar ang 2
floatinphys(6+1) = 1.1*10^-2  ;% sol wid omp 
floatinphys(7+1) = 0.0 ;% loc omp 
floatinphys(8+1) = 1.4 ;% maj rad

floatinphys(9+1) = 1.0 ;% alpha core prof q
floatinphys(10+1) = 1.0*10^-3  ;% alpha core prof n
floatinphys(11+1) = 2.0 ;% flux expansion
floatinphys(12+1) = 3.0 ;% trans expansion
floatinphys(13+1) = 0.0 ;% gas puf loc
floatinphys(14+1) = 10^20  ;% gas puf wid
floatinphys(15+1) = 0.4500 ;% sig nb backgrounds
floatinphys(16+1) = 0.0; % L_baffle

% sol influx settings
floatinphys(19+1) = 0.0 ;% gamma_x
floatinphys(20+1) = 0.0 ;% qparx
floatinphys(21+1) = 10^20  ;% gamma core
floatinphys(22+1) = 2.2*10^5  ;% q core
floatinphys(23+1) = 0.0 ;% density ramp rate
floatinphys(24+1) = 0.0 ;% gas puf src

% Q_CORE =  2.21735701291665e+05
%  GAMMA_CORE = 1.0000e+20

% initial values
floatinphys(29+1) = 10^19  ;% initial n
floatinphys(30+1) = 0.0 ;% initial v 
floatinphys(31+1) = 100.0 ;% initial T 
floatinphys(32+1) = 6.2*10^17  ;% initial a
floatinphys(33+1) = 0.0 ;% initial vn
floatinphys(34+1) = 5.5*10^17  ;% initial m
floatinphys(35+1) = 2.7*10^18  ;% initial ncor 
floatinphys(36+1) = 1.0*10^14  ;% initial core neutral

floatinphys(39+1) = 6.28*10^17 ;
floatinphys(40+1) = 6.28*10^17 ;
floatinphys(41+1) = 5.4*10^16 ; 
floatinphys(42+1) = 6.28*10^17 ;
floatinphys(43+1) = 6.28*10^17 ; % initial nb(5 

floatinphys(44+1) = 3.89*10^18 ;
floatinphys(45+1) = 6.44*10^18 ;
floatinphys(46+1) = 1.87*10^17 ;
floatinphys(47+1) = 6.44*10^18 ;
floatinphys(48+1) = 3.89*10^18 ;% initial mb(5 

% limits
floatinphys(49+1) = 10^10 ;% min den
floatinphys(50+1) = 10^25 ;% max den
floatinphys(51+1) = 0.1000;% min tem

% sol behavior
floatinphys(59+1) = 6 ;% gamma
floatinphys(60+1) = 3.3436*10^-27  ;% mass
floatinphys(61+1) = 30 ;% energy loss ion
floatinphys(62+1) = 0.4 ;% recycling
floatinphys(63+1) = 1.0;% mol_rec
floatinphys(64+1) = 5.0*10^-4  ;% neutral residence time
floatinphys(65+1) = 8.0*10^-3  ;% molecule residence time
floatinphys(66+1) = 0.0000 ; % far sol ion losses 
% new array
floatinphys(69+1) =2.0*10^-3 ; 
floatinphys(70+1) = 0.0*10^-12 ;
floatinphys(71+1) = 0.0*10^-12 ;
floatinphys(72+1) = 0.0*10^-12 ;
floatinphys(73+1) = 0.0*10^-12 ;  % imp con 1:5

% reservoirs
floatinphys(79+1) = 10^-22 ;
floatinphys(80+1) = 10^-22 ;
floatinphys(81+1) = 2.05*10^2 ;
floatinphys(82+1) = 10^-22 ; 
floatinphys(83+1) = 10^-22 ;% cor res ato pmp(1:5 

floatinphys(84+1) = 10^-22 ;
floatinphys(85+1) = 10^-22 ; 
floatinphys(86+1) = 2.05*10^2 ; 
floatinphys(87+1) = 10^-22 ;
floatinphys(88+1) = 10^-22 ; % cor res mol pmp(1:5 

floatinphys(89+1) = 1.0599; 
floatinphys(90+1) = 1.826200;
floatinphys(91+1) = 30.000;
floatinphys(92+1) = 1.8268;
floatinphys(93+1) = 1.0599;% ext vol 1:5

floatinphys(99+1) = 0.000;
floatinphys(100+1) = 0.000;
floatinphys(101+1) = 6.67*10^2 ;% ext ato ext 1:3
floatinphys(102+1) = 0.000;
floatinphys(103+1) = 0.000; 
floatinphys(104+1) = 4.717*10^2  ;% ext mol ext 1:3

floatinphys(109+1) = 0.000;
floatinphys(110+1) = 0.000;
floatinphys(111+1) = 1.25*10^4 ;
floatinphys(112+1) = 6.25*10^2 ;
floatinphys(113+1) = 1.875*10^3 ;% pmp rat n(1:5 

floatinphys(114+1) = 0.000;
floatinphys(115+1) = 0.000; 
floatinphys(116+1) = 1.25*10^4 ;
floatinphys(117+1) = 6.25*10^2 ; 
floatinphys(118+1) = 1.875*10^3 ;% pmp rat m(1:5 

% core
floatinphys(129+1) = 0.2000 ;% core conf tim
floatinphys(130+1) = 0.0*10^-22 ;% cor sol ato exchange
floatinphys(131+1) = 0.0*10^-22 ;% cor sol mol exchange
floatinphys(132+1) = 4.877; % cor volumes




out.floatinphys = floatinphys;
out.intinphys = intinphys;

out.floatinnum = floatinnum;
out.intinnum  = intinnum;
out.loginnum(1:10)  = false;
out.loginnum(3) = true;



end