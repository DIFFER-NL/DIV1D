\documentclass[amsmath,amssymb,a4]{revtex4-2}
% updated from revtex4 tot revtex4-2 because revtex4 would cause compiler malfunctioning in the most recent Miktex package. It is superseded by revtex4-2 (https://journals.aps.org/revtex/revtex-faq) so that was used as of 12-2-2020.
\usepackage{graphicx}
\begin{document}

\title[DIV1D manual]{DIV1D manual: a 1D code for divertor plasma simulation}

\author{E. Westerhof}

\address{DIFFER -- Dutch Institute for Fundamental Energy Research, PO Box 6336, 5600HH Eindhoven, The Netherlands, www.differ.nl}

\email[E-mail: ]{e.westerhof@differ.nl}
\vspace{10pt}
\date\today


\begin{abstract}
Abstract goes here.
\end{abstract}


\maketitle

\section{Introduction}

This manual describes the {\tt DIV1D} code for the quick simulation of the 1D plasma and neutral behaviour along a flux tube in a tokamak divertor from the X-point up to the target. The code is inspired by similar works described in the papers by Nakazawa et al.~\cite{nakazawa2000} and Dudson et al.~\cite{dudson2019, SD1D}. In fact, the numerical implementation borrows heavily from the methods used in the SD1D code. The code is maintained on a SVN repository at DIFFER (https://www.differ.nl/svn/shared/div1d). This document describes the basic model equations in Section~\ref{basic_equations} including the boundary conditions and the sources and sinks. Section~\ref{rates} contains a full description of the different choices that are implemented to calculate the various atomic processes like charge exchange, ionization, recombination, and impurity radiation. Some notes on the discretization and the numerical implementation are given in Section~\ref{numerics}.


\section{Equations}\label{basic_equations}

The equations solved are: the particle balance equation, the plasma momentum balance, the plasma energy balance, and an equation for the evolution of the neutral particle density. As in Nakazawa et al.~\cite{nakazawa2000} the neutral momentum and energy are ignored.

\noindent The particle balance is given by
\begin{equation}\label{particle_balance}
    {\partial n \over \partial t} = - B{\partial \over \partial x} \left({\Gamma_n \over B}\right) + S_n,
\end{equation}
where $n$ is the plasma (electron) density, $\Gamma_n = n v_\parallel$ with $v_\parallel$ the parallel velocity is the convective particle flux (a possible effect of diffusion is ignored), and $S_n$ represent the sum of all particle sources and sinks. $B$ represents the magnitude of the total magnetic field. Here and in the equations below, the inclusion of $B$ accounts for the effect of flux expansion due to a varying total magnetic field~\cite{dudson2019,havlickova2013}.

\noindent The momentum balance is given by
\begin{equation}\label{momentum_balance}
    {\partial n m v_\parallel \over \partial t} = - B{\partial \over \partial x} \left({n m v_\parallel^2 \over B}\right) - {\partial \over \partial x} p + S_{\rm mom},
\end{equation}
where $m$ is the mass of the dominant ion species (default Deuterium $m = 3.3436 \times 10^{-27}$~kg) , $p = 2 n e T$ is the total plasma pressure (we will use T in units eV such that the Boltzmann constant can be equated with the elementary charge $e$ to obtain the pressure in Pascal), and $S_{\rm mom}$ represent the sum of all momentum sources and sinks.

\noindent Ion and electron temperatures are considered equal such that only one energy balance needs to be solved. The internal energy ($3 n k T$) balance of the plasma is then given by
\begin{equation}\label{energy_balance}
    {\partial 3 n e T \over \partial t} = - B{\partial \over \partial x} \left({q_\parallel \over B}\right) + v_\parallel {\partial \over \partial x} p + Q,
\end{equation}
where the heat flux $q_\parallel$ is given by the equation
\begin{equation}\label{heat_flux}
    q_\parallel = 5 n e T v_\parallel - \kappa_\parallel {\partial \over \partial x} T,
\end{equation}
where the first term on the right hand side represents the total enthalpy flux and the second term represents the parallel heat conduction with (for T in eV) the parallel conductivity being given by $\kappa_\parallel = 2 \times 10^3 T^{5/2}$ J/eVms (see Chapter 4.10.1 of \cite{stangeby}). $Q$ represent the sum of all internal energy sources and sinks.

\noindent The neutral particle balance is given by
\begin{equation}\label{neutral_particle_balance}
    {\partial n_n \over \partial t} = {\partial \over \partial x} D {\partial \over \partial x} n_n + S_{\rm neutral},
\end{equation}
where $n_n$ is the neutral particle density and D is the neutral particle diffusion coefficient which is given by \cite{nakazawa2000}
\begin{equation}\label{neutral_diffusion_coefficient}
    D = { eT \over m \nu_{\rm cx} \sin^2\theta }
\end{equation}
where $\theta$ is the angle at which the magnetic field hits the target, and $\nu_{\rm cx} = n \langle\sigma_{\rm cx} v\rangle$ is the charge exchange collision frequency for the neutrals and the average charge exchange rate $\langle\sigma_{\rm cx} v\rangle$ is specified below. The factor $\sin^2\theta$ appears in the equation for the neutral diffusion because the neutrals are free to move across magnetic field lines, such that their motion perpendicular to the target surface results in an effective parallel velocity that is increased by a factor $1/\sin\theta$.


\subsection{Boundary Conditions}

Each of these equations requires boundary conditions at the X-point $x=0$ and at the target $x=L$, where $L$ is the given length of the flux tube. At the X-point we use the following boundary conditions:
the plasma density at the X-point is given: $n(x=0) = n_{\rm X}$; the plasma momentum flux at the X-point assumed to be constant: the parallel heat flux at the X-point is given:  $q_\parallel(x=0) = q_{\parallel,\rm X}$; and finally, the neutral particle density is assumed to have zero gradient.

The boundary conditions at the target are given by the usual sheath boundary conditions assuming that density and temperature are constant across the sheath, while the plasma particle flux and momentum are given by the Bohm condition, $\Gamma_n(x=L) \ge n c_{\rm s}$ and  where $c_{\rm s} = \sqrt{2eT/m}$ is the plasma sound speed, and the the heat flux on the target is given by the sheath heat transmission factor $\gamma$ which must be specified at input:
\begin{equation}\label{sheath_heat_transmission}
    q_\parallel(x=L) = \gamma n e T c_{\rm s}.
\end{equation}
The neutral particle flux coming from the target is determined by the recycling coefficient $R$ and a redistribution fraction $f_r$ given no input, i.e. $\Gamma_{\rm neutral}(x=L) = - R (1 - f_r) \Gamma_n(x=L)$. Note that when $R=1$ the total number of particles should be conserved requiring a zero plasma inflow velocity at the X-point in case of a steady state solution.


\subsection{Sources and Sinks}

The various sources and sinks are determined mostly by atomic processes like charge exchange, ionization, excitation, and recombination. In addition a neutral particle loss term can be specified in terms of an average residence time for the neutral in the flux tube to account for cross field neutral particle transport losses. For the plasma density the sources and sinks are given by ionization, recombination and target recycling, respectively:
\begin{equation}\label{particle_source}
    S_n = + n n_n \langle\sigma_{\rm ion} v\rangle - n^2 \langle\sigma_{\rm rec} v\rangle.
\end{equation}

Because the neutral particle momentum is neglected there is no momentum source, while the momentum sinks are induced by charge exchange and recombination
\begin{equation}\label{momentum_source}
    S_{\rm mom} = - n m v_\parallel \left(n_n \langle\sigma_{\rm cx} v\rangle + n \langle\sigma_{\rm rec} v\rangle \right).
\end{equation}

The energy balance contains both heat sinks and sources associated with charge exchange, hydrogenic ionization, and excitation, radiative and three-body recombination as well as impurity radiation
\begin{equation}\label{energy_source}
\begin{split}
    Q = &- 1.5 n e T \left(n_n \langle\sigma_{\rm cx} v\rangle + 2 n \langle\sigma_{\rm rec} v\rangle \right) \\
        &+ \frac{1}{2} n m v^2_\parallel n_n \langle\sigma_{\rm ion} v\rangle\\
        &- n n_n E_{\rm ion} \langle \sigma_{\rm ion} v\rangle - n n_n \langle E_{\rm exc} \sigma_{\rm exc} v\rangle \\
        &+ n^2 \;E_{\rm ion} \langle \sigma_{\rm rec} v\rangle - n^2 \;\langle E_{\rm rad} \sigma_{\rm rec} v\rangle \\
        &- n^2 \xi_Z L_Z(T)
\end{split}
\end{equation}
where the first line represents the losses to the neutrals due to charge exchange and recombination and the second line is a source term from friction with the neutrals that comes from ionization of (assumed momentum less) neutrals. The third line comes from ionization and the associated excitation, where $E_{\rm ion} = 13.6$~eV is the ionization energy and $E_{\rm exc}$ is the energy loss per excitation. In practice these are typically combined into an effective energy loss per ionization (see the paragraph below on the implemented reaction rates). The fourth line represents the effective energy balance from three-body and radiative recombination. Finally, the last line $L_Z(T)$ is the radiative cooling rate of the impurity with concentration $\xi$ for up to five specified impurities. 

What is a sink for the plasma density is a source for the neutral particles, while additional neutral particle sources are associated with gas puff, redistribution of recycled neutrals or a possible exchange with a neutral background, such that
\begin{equation}\label{neutral_source}
    S_{\rm neutral} = - S_n + S_{\rm puff} - \frac{n_n - n_b}{\tau_n} + {n(L)v_{\|}(L)R f_r \over L}
\end{equation}
where $S_{\rm puff}$ is a neutral particle source from an additional gas puff, and $\tau_n$ accounts for a finite residence time of the neutral particles. When a finite background neutral density $n_b$ is specified this term represents the exchange of neutrals in the leg with these background neutrals. A fraction $f_r$ of the recycling flux is redistributed homogeneously over the entire leg. The atomic rates for charge exchange, recombination, ionization, and excitation as used in the code as given in a later paragraph where also the impurity radiative cooling rates are be given.


\subsection{Radial power losses}
The power losses in the radial direction can be modelled through a radial power loss factor. Although one should be cautious when using this together with flux expansion as  they have similar effects on the power balance.
To simulate this effect, DIV1D includes several versions of an effective radial loss term. The simplest version is
\begin{equation}
\label{constant_rad_losses}
S_{\perp}=-\eta_{\perp}\frac{q_\parallel}{L}
\end{equation}

where $S_{\perp}$ denotes the radial power sink, and $\eta_{\perp}$ is the fraction of the upstream heat flux that is lost over the length of the flux tube due to radial losses. A slightly more advanced version of these losses is

\begin{equation}
S_\perp(x)=-\frac{\eta_{\perp}q_{\parallel}}{\sqrt{2\pi(\mathrm{\Delta}x)^2}}
\mathrm{exp}{-\frac{1}{2}\left(\frac{x-x_0}{\mathrm{\Delta}x}\right)^2}
\end{equation}

where losses are no longer incurred equally over the flux tube, but are instead peaked at $x_0$ with a typical width of $\mathrm{\Delta}x$ around it. Note that the prefactor is calculated analytically in the code to make sure that the total losses always amount to $-\eta_{\perp}\times q_{\parallel}$. A comparison of either loss profile show that the target parameters are rather indifferent to the choice of loss profile.\\

Finally, a heat-flux-dependent variant is also included as

\begin{equation}
S_{\perp}=-\eta_{\perp}\frac{q_\parallel(x)}{L}
\end{equation}

which is crucially different from~\ref{constant_rad_losses} in that it takes into account the local heat flux in each grid cell rather than the upstream heat flux. As a result $\eta_{\perp}$ no longer reflects the total fraction of the upstream heat flux that is lost.\\

Four input parameters exist which the user can employ to shape the radial losses to his or her desire. These are:
radial\_loss\_factor ($\eta_{\perp}$),
radial\_loss\_gaussian (choice of loss profile),
radial\_loss\_width ($\mathrm{\Delta}x$), and
radial\_loss\_location ($x_0$). See table~2. Although this feature was translated into the flux expansion factor, it might be usefull again when extending to the core scrape-off layer. 

\section{Atomic Rates}\label{rates}

A number of different options is available to calculate the atomic rates. For each of the rates a number of options is provided for different expressions that have been used in the literature. These include approximate formulas which are used amongts others in the SD1D code \cite{SD1D} of which the origin is however not always fully clear. The to our knowledge more accurate rates are obtained from the AMJUEL data base that is also being used for the EIRENE neutral particle Monte Carlo code and which can be found at the web site of the code \cite{EIRENE}. By default the rates from the AMJUEL data base are selected.


\subsection{Charge Exchange}

The default option uses the charge exchange rate as specified in section 2.19 reaction 3.1.8 of the AMJUEL data base for the total charge exchange rate of Hydrogen~\cite{EIRENE}. This uses a fit function of the form
\begin{equation}\label{charge_exchange_AMJUEL}
    \langle\sigma_{\rm cx} v\rangle = 10^{-6} \exp\left( \sum_{i=0}^8 b_i (\ln T)^i \right)  {\rm m}^3/{\rm s}
\end{equation}
with the fitting constants $b_i$ defined by
\begin{small}\begin{verbatim}
   b0 -1.850280000000E+01    b1  3.708409000000E-01    b2  7.949876000000E-03
   b3 -6.143769000000E-04    b4 -4.698969000000E-04    b5 -4.096807000000E-04
   b6  1.440382000000E-04    b7 -1.514243000000E-05    b8  5.122435000000E-07
\end{verbatim}\end{small}
Note that the factor $10^{-6}$ stems from the use of the units ${\rm cm}^3/{\rm s}$ for the reaction rates in AMJUEL.

In order to allow comparisons with models from other authors a number of alternative options are provided for the caclculations of the reaction rates. This includes the charge exchange rate as implemented in SD1D originating from the work by Havlickova \cite{havlickova2013} which is given by \cite{SD1D}
\begin{equation}\label{charge_exchange_SD1D}
    \langle\sigma_{\rm cx} v\rangle = \begin{cases} 1.0 \times 10^{-14} {\rm m}^3/{\rm s}             & \mbox{for } T \le 1 {\rm eV} \\
                                        1.0 \times 10^{-14} T^{1/3} {\rm m}^3/{\rm s} & \mbox{for } T >   1 {\rm eV}. \end{cases}
\end{equation}
This option is selected by setting {\tt case\_cx = "Havlickova"}.

Another expression that has been used for example in the work by Nakazawa et al. \cite{nakazawa2000}, comes from Table 3 of the report by Freeman and Jones \cite{freeman1974} in which case the charge exchange rate is given by a fit function of the same form as above
\begin{equation}\label{charge_exchange_FJ}
    \langle\sigma_{\rm cx} v\rangle = 10^{-6} \exp\left( \sum_{i=0}^8 b_i (\ln T)^i \right)  {\rm m}^3/{\rm s}
\end{equation}
with fitting coefficients $b_i$ now given as
\begin{small}\begin{verbatim}
   b0 -1.841757E+01  b1  5.282950E-01  b2 -2.200477E-01
   b3  9.750192E-02  b4 -1.749183E-02  b5  4.954296E-04
   b6  2.174910E-04  b7 -2.530206E-05  b8  8.230751E-07
\end{verbatim}\end{small}
This option is selected by setting {\tt case\_cx = "Freeman"}. The validity range of this fit is indicated as 1 to $10^5$~eV.

Charge exchange rates for the Hydrogen isotopes like Deuterium and Tritium are obtained by the same expressions given above, but using a rescaled temperature multiplied with the factor $m_p/m_i$, i.e. the ratio of the proton mass over the mass of the Deuterium or Tritium ion, respectively.


\subsection{Ionization}

The default option for calculation of the ionization rate is also obtained from the AMJUEL data base, which provides an effective ionization rate as calculated using a double fit function
\begin{equation}\label{ionization_AMJUEL}
    \langle\sigma_{\rm ion} v\rangle = 10^{-6} \exp\left( \sum_{i=0}^8\sum_{j=0}^8 \alpha_{ij} (\ln \bar n)^j (\ln T)^i \right)  {\rm m}^3/{\rm s}
\end{equation}
where the density is normalized as $\bar n \equiv n / 10^{14}$ and the fitting coefficients $\alpha_{ij}$ are given in section 4.3 reaction 2.1.5 of the AMJUEL document for the case of the total ionization rate (including all excited states of the neutral hydrogen atoms)
\begin{small}\begin{verbatim}
        n-Index:     0                     1                     2
  T-Index:
        0   -3.248025330340D+01   -5.440669186583D-02    9.048888225109D-02
        1    1.425332391510D+01   -3.594347160760D-02   -2.014729121556D-02
        2   -6.632235026785D+00    9.255558353174D-02   -5.580210154625D-03
        3    2.059544135448D+00   -7.562462086943D-02    1.519595967433D-02
        4   -4.425370331410D-01    2.882634019199D-02   -7.285771485050D-03
        5    6.309381861496D-02   -5.788686535780D-03    1.507382955250D-03
        6   -5.620091829261D-03    6.329105568040D-04   -1.527777697951D-04
        7    2.812016578355D-04   -3.564132950345D-05    7.222726811078D-06
        8   -6.011143453374D-06    8.089651265488D-07   -1.186212683668D-07

        n-Index:     3                     4                     5
  T-Index:
        0   -4.054078993576D-02    8.976513750477D-03   -1.060334011186D-03
        1    1.039773615730D-02   -1.771792153042D-03    1.237467264294D-04
        2   -5.902218748238D-03    1.295609806553D-03   -1.056721622588D-04
        3    5.803498098354D-04   -3.527285012725D-04    3.201533740322D-05
        4    4.643389885987D-04    1.145700685235D-06    8.493662724988D-07
        5   -1.201550548662D-04    6.574487543511D-06   -9.678782818849D-07
        6    8.270124691336D-06    3.224101773605D-08    4.377402649057D-08
        7    1.433018694347D-07   -1.097431215601D-07    7.789031791949D-09
        8   -2.381080756307D-08    6.271173694534D-09   -5.483010244930D-10

        n-Index:     6                     7                     8
  T-Index:
        0    6.846238436472D-05   -2.242955329604D-06    2.890437688072D-08
        1   -3.130184159149D-06   -3.051994601527D-08    1.888148175469D-09
        2    4.646310029498D-06   -1.479612391848D-07    2.852251258320D-09
        3   -1.835196889733D-06    9.474014343303D-08   -2.342505583774D-09
        4   -1.001032516512D-08   -1.476839184318D-08    6.047700368169D-10
        5    5.176265845225D-08    1.291551676860D-09   -9.685157340473D-11
        6   -2.622921686955D-09   -2.259663431436D-10    1.161438990709D-11
        7   -4.197728680251D-10    3.032260338723D-11   -8.911076930014D-13
        8    3.064611702159D-11   -1.355903284487D-12    2.935080031599D-14
T1MIN =   0.10000D 00 EV
T1MAX =   2.00000D 04 EV
N2MIN =   1.00000D 08 1/CM3
N2MAX =   1.00000D 16 1/CM3
\end{verbatim}\end{small}

One alternative option is again provided in the form of the ionization rate as implemented in SD1D originating from the work by Havlickova \cite{havlickova2013}. It is slightly modified to remove the discontinuity at 20~eV and is given by \cite{SD1D}
\begin{equation}\label{ionization_SD1D}
    \langle\sigma_{\rm ion} v\rangle = \begin{cases} 7.638 \times 10^{-21} {\rm m}^3/{\rm s}             & \mbox{for } T \le 1 {\rm eV} \\
                                        10^{-6.0d+0} T^{-2.987} 10^{-15.72 \exp(-\log_{10}T) + 1.603*exp(-\log_{10}^2T)} {\rm m}^3/{\rm s} & \mbox{for } 1 {\rm eV} < T \le 20 {\rm eV} \\
                                        5.875 \times 10^{-12} T^{-0.5151} 10^{-2.563/\log_{10}T} {\rm m}^3/{\rm s}. \end{cases}
\end{equation}
This option is selected by setting {\tt case\_ion = "Havlickova"}

A further option again comes from the work of Nakazawa et al.~\cite{nakazawa2000} who use the ionization rate given in Table 3 of the report by Freeman and Jones~\cite{freeman1974}. They provide a fit to the electron impact ionization given by
\begin{equation}\label{ionization_FJ}
    \langle\sigma_{\rm ion} v\rangle = 10^{-6} \exp\left( \sum_{i=0}^6 b_i (\ln T)^i \right)  {\rm m}^3/{\rm s}
\end{equation}
with fitting coefficients $b_i$ given by
\begin{small}\begin{verbatim}
   b0 -0.3173850E+02  b1  0.1143818E+02  b2 -0.3833998E+01
   b3  0.7046692E+00  b4 -0.7431486E-01  b5  0.4153749E-02
   b6 -0.9486967E-04
\end{verbatim}\end{small}
This option is selected by setting {\tt case\_ion = "Freeman"}


\subsection{Excitation and ionization energy losses}

By default, the sum of energy losses from ionization and excitation is obtained from the AMJUEL data base, which provides an effective excitation rate in terms of an averaged effective energy loss per ionization. This is calculated from a fit function of the same form as defined above for the ionization rate (\ref{ionization_AMJUEL}), i.e.
\begin{equation}
  \langle E_{\rm ion} \sigma_{\rm ion} v\rangle + \langle E_{\rm exc} \sigma_{\rm exc} v\rangle =  10^{-6} \exp\left( \sum_{i=0}^8\sum_{j=0}^8 \alpha_{{\rm exc,} ij} (\ln \bar n)^j (\ln T)^i \right)  {\rm eV} {\rm m}^3/{\rm s}
\end{equation}
with the coefficients $\alpha_{{\rm exc,} ij}$ as tabulated in section 10.2 for reaction 2.5.1 of the AMJUEL document for the case of the total energy loss rate associated with Hydrogen ionization and excitation radiation:
\begin{small}\begin{verbatim}
        n-Index:     0                     1                     2
  T-Index:
        0   -2.497580168306D+01    1.081653961822D-03   -7.358936044605D-04
        1    1.004448839974D+01   -3.189474633369D-03    2.510128351932D-03
        2   -4.867952931298D+00   -5.852267850690D-03    2.867458651322D-03
        3    1.689422238067D+00    7.744372210287D-03   -3.087364236497D-03
        4   -4.103532320100D-01   -3.622291213236D-03    1.327415215304D-03
        5    6.469718387357D-02    8.268567898126D-04   -2.830939623802D-04
        6   -6.215861314764D-03   -9.836595524255D-05    3.017296919092D-05
        7    3.289809895460D-04    5.845697922558D-06   -1.479323780613D-06
        8   -7.335808238917D-06   -1.367574486885D-07    2.423236476442D-08

        n-Index:     3                     4                     5
  T-Index:
        0    4.122398646951D-04   -1.408153300988D-04    2.469730836220D-05
        1   -7.707040988954D-04    1.031309578578D-04   -3.716939423005D-06
        2   -8.328668093987D-04    2.056134355492D-04   -3.301570807523D-05
        3    4.707676288420D-04   -5.508611815406D-05    7.305867762241D-06
        4   -1.424078519508D-04    3.307339563081D-06    5.256679519499D-09
        5    2.411848024960D-05    5.707984861100D-07   -1.016945693300D-07
        6   -1.474253805845D-06   -2.397868837417D-07    1.518743025531D-08
        7   -4.633029022577D-08    3.337390374041D-08   -1.770252084837D-09
        8    5.733871119707D-09   -1.512777532459D-09    8.733801272834D-11

        n-Index:     6                     7                     8
  T-Index:
        0   -2.212823709798D-06    9.648139704737D-08   -1.611904413846D-09
        1   -4.249704742353D-07    4.164960852522D-08   -9.893423877739D-10
        2    2.831739755462D-06   -1.164969298033D-07    1.785440278790D-09
        3   -6.000115718138D-07    2.045211951761D-08   -1.790312871690D-10
        4    7.597020291557D-10    1.799505288362D-09   -9.280890205774D-11
        5    3.517154874443D-09   -4.453195673947D-10    2.002478264932D-11
        6    4.149084521319D-10   -6.803200444549D-12   -1.151855939531D-12
        7   -5.289806153651D-11    3.864394776250D-12   -8.694978774411D-15
        8    7.196798841269D-13   -1.441033650378D-13    1.734769090475D-15
T1MIN =   0.10000D 00 EV
T1MAX =   2.00000D 04 EV
N2MIN =   1.00000D 08 1/CM3
N2MAX =   1.00000D 16 1/CM3
\end{verbatim}\end{small}


Energy losses from ionization and due to excitation radiation are accounted for either as an effective energy loss constant per ionization (typically $E_{eff,ion} = 30$~eV as in Nakazawa et al. \cite{nakazawa2000}) which can be specified at input er by an effective excitation rate which accounts for the temperature and density dependent total energy loss from ionization and excitation radiation. In the SD1D code the latter is given by a simple fit function as
\begin{equation}\label{excitation_SD1D}
    \langle\sigma_{\rm exc} v\rangle = {4.90 \times 10^{-13} \over 0.28 + Y} \exp(-Y) \sqrt{Y (1.0 + Y)} {\rm eV} {\rm m}^3/{\rm s},
\end{equation}
where $Y = 10.2 / \max( 1, T)$. Added to this is the 13.6~eV of energy loss per ionization.



\subsection{Recombination}

In their 1D code Nakazawa et al.~\cite{nakazawa2000} use a combination of the radiative recombination rate from Gordeev et al. plus the three body recombination rate given by Hinnov et al., where the former is given by the expression~\cite{gordeev1977}
\begin{equation}\label{radiative_recombination}
    \langle\sigma_{\rm rad.rec} v\rangle = 1.27 \times 10^{-19} {(13.6./T)^{1.5} \over (13.6/T) + 0.59} {\rm m}^3 / s,
\end{equation}
and the latter is expressed as~\cite{hinnov1962}
\begin{equation}\label{three_body_recombination}
    \langle\sigma_{\rm 3bodyrec} v\rangle = 5.6 \times 10^{-39} \; T^{-4.5} n {\rm m}^3 / s.
\end{equation}
This model is selected by setting {\tt case\_rec = "Nakazawa"}.

Otherwise, the recombination rate is used from the AMJUEL data providing the total effective recombination rate including 3 body recombination using again a fit function (\ref{ionization_AMJUEL}) as defined above for the ionization rate, now with the coefficients $\alpha_{{\rm rec,} ij}$ as tabulated in section 4.6 reaction 2.1.8 of the AMJUEL document:
\begin{small}\begin{verbatim}
        n-Index:     0                     1                     2
  T-Index:
        0   -2.858858570847D+01    2.068671746773D-02   -7.868331504755D-03
        1   -7.676413320499D-01    1.278006032590D-02   -1.870326896978D-02
        2    2.823851790251D-03   -1.907812518731D-03    1.121251125171D-02
        3   -1.062884273731D-02   -1.010719783828D-02    4.208412930611D-03
        4    1.582701550903D-03    2.794099401979D-03   -2.024796037098D-03
        5   -1.938012790522D-04    2.148453735781D-04    3.393285358049D-05
        6    6.041794354114D-06   -1.421502819671D-04    6.143879076080D-05
        7    1.742316850715D-06    1.595051038326D-05   -7.858419208668D-06
        8   -1.384927774988D-07   -5.664673433879D-07    2.886857762387D-07

        n-Index:     3                     4                     5
  T-Index:
        0    3.843362133859D-03   -7.411492158905D-04    9.273687892997D-05
        1    3.828555048890D-03   -3.627770385335D-04    4.401007253801D-07
        2   -3.711328186517D-03    6.617485083301D-04   -6.860774445002D-05
        3   -1.005744410540D-03    1.013652422369D-04   -2.044691594727D-06
        4    6.250304936976D-04   -9.224891301052D-05    7.546853961575D-06
        5   -3.746423753955D-05    7.509176112468D-06   -8.688365258514D-07
        6   -1.232549226121D-05    1.394562183496D-06   -6.434833988001D-08
        7    1.774935420144D-06   -2.187584251561D-07    1.327090702659D-08
        8   -6.591743182569D-08    8.008790343319D-09   -4.805837071646D-10

        n-Index:     6                     7                     8
  T-Index:
        0   -7.063529824805D-06    3.026539277057D-07   -5.373940838104D-09
        1    1.932701779173D-06   -1.176872895577D-07    2.215851843121D-09
        2    4.508046989099D-06   -1.723423509284D-07    2.805361431741D-09
        3   -4.431181498017D-07    3.457903389784D-08   -7.374639775683D-10
        4   -3.682709551169D-07    1.035928615391D-08   -1.325312585168D-10
        5    7.144767938783D-08   -3.367897014044D-09    6.250111099227D-11
        6   -2.746804724917D-09    3.564291012995D-10   -8.551708197610D-12
        7   -1.386720240985D-10   -1.946206688519D-11    5.745422385081D-13
        8    6.459706573699D-12    5.510729582791D-13   -1.680871303639D-14
T1MIN =   0.10000D 00 EV
T1MAX =   2.00000D 04 EV
N2MIN =   1.00000D 08 1/CM3
N2MAX =   1.00000D 16 1/CM3
\end{verbatim}\end{small}

When the default option is selected to use the data from AMJUEL also the energy lost and gained by the electrons due to radiative and three-body recombination is taken into account. In this case the effective electron cooling rate from the associated processes is obtained from the fit specified in AMJUEL section 10.4 for the recombination reaction and taking into account the 13.6~eV potential enrgy gain per effective recombination event, i.e.
\begin{equation}
  \langle E_{\rm el} \sigma_{\rm rec} v\rangle =  10^{-6} \exp\left( \sum_{i=0}^8\sum_{j=0}^8 \alpha_{{\rm rec el,} ij} (\ln \bar n)^j (\ln T)^i \right) - 13.6 \langle\sigma_{\rm rec} v\rangle  {\rm eV} {\rm m}^3/{\rm s}
\end{equation}
where the effective recombination rate $\langle\sigma_{\rm rec} v\rangle$ is obtained from the fit specified before and the coefficients $\alpha_{{\rm rec el,} ij}$ of the current fit are given by
\begin{small}\begin{verbatim}

        E-Index:     0                     1                     2
  T-Index:
        0   -2.592450349909D+01    1.222097271874D-02    4.278499401907D-05
        1   -7.290670236493D-01   -1.540323930666D-02   -3.406093779190D-03
        2    2.363925869096D-02    1.164453346305D-02   -5.845209334594D-03
        3    3.645333930947D-03   -1.005820792983D-03    6.956352274249D-04
        4    1.594184648757D-03   -1.582238007548D-05    4.073695619272D-04
        5   -1.216668033378D-03   -3.503070140126D-04    1.043500296633D-04
        6    2.376115895241D-04    1.172709777146D-04   -6.695182045674D-05
        7   -1.930977636766D-05   -1.318401491304D-05    8.848025453481D-06
        8    5.599257775146D-07    4.977823319311D-07   -3.615013823092D-07

        E-Index:     3                     4                     5
  T-Index:
        0    1.943967743593D-03   -7.123474602102D-04    1.303523395892D-04
        1    1.532243431817D-03   -4.658423772784D-04    5.972448753445D-05
        2    2.854145868307D-03   -5.077485291132D-04    4.211106637742D-05
        3   -9.305056373739D-04    2.584896294384D-04   -3.294643898894D-05
        4   -9.379169243859D-05    1.490890502214D-06    2.245292872209D-06
        5    9.536162767321D-06   -6.908681884097D-06    8.232019008169D-07
        6    1.188184006210D-05   -4.381514364966D-07   -6.936267173079D-08
        7   -2.072370711390D-06    2.055919993599D-07   -7.489632654212D-09
        8    9.466989306497D-08   -1.146485227699D-08    6.772338917155D-10

        E-Index:     6                     7                     8
  T-Index:
        0   -1.186560752561D-05    5.334455630031D-07   -9.349857887253D-09
        1   -4.070843294052D-06    1.378709880644D-07   -1.818079729166D-09
        2   -1.251436618314D-06   -1.626555745259D-08    1.073458810743D-09
        3    2.112924018518D-06   -6.544682842175D-08    7.810293075700D-10
        4   -3.150901014513D-07    1.631965635818D-08   -2.984093025695D-10
        5   -2.905331051259D-08   -3.169038517749D-10    2.442765766167D-11
        6    6.592249255001D-09   -1.778887958831D-10    1.160762106747D-12
        7   -7.073797030749D-11    1.047087505147D-11   -1.877446271350D-13
        8   -1.776496344763D-11    7.199195061382D-14    3.929300283002D-15
T1MIN =   0.10000D 00 EV
T1MAX =   2.00000D 04 EV
N2MIN =   1.00000D 08 1/CM3
N2MAX =   1.00000D 16 1/CM3

  Max. rel. Error:   0.930E+01 %
  Mean rel. Error:   0.127E+01 %
\end{verbatim}\end{small}


\subsection{Impurity Radiation Losses}

The impurity radiation losses are typically given in the form of a radiative cooling function $L_Z(T)$ where $Z$ stands for the impurity species under consideration. For a given cooling rate function, the energy losses from impurity radiation are then given by
\begin{equation}\label{impurity_radiation}
    Q_{\rm imp} = - n^2 \xi_Z L_Z(T),
\end{equation}
where $\xi_Z$ is the concentration of the impurity. Post et al.\cite{post1977} tabulate fit functions for the most relevant impurities in fusion plasmas using the general functional form of \cite{post1977}
\begin{equation}\label{cooling_rate}
    \log_{10} L_Z = \sum_{i=0}^5 \; A(i) (\log_{10} T_{\rm keV})^i [{\rm cm}^3 {\rm erg/s}]
\end{equation}
where $T_{\rm keV}$ is the temperature in keV.

In the present code version the radiation cooling function of Carbon and Nitrogen are implemented, but any impurity tabulated by Post et al \cite{post1977} can readily be added. The coefficients of the Carbon fit function are given in the following table for various temperature ranges
%\begin{small}\begin{verbatim}
%  TMIN   TMAX       A(0)          A(1)          A(2)          A(3)          A(4)          A(5)
%     3     20   1.965300E+03  4.572039E+03  4.159590E+03  1.871560E+03  4.173889E+02  3.699382E+01
%    20    200   7.467599E+01  4.549038E+02  8.372937E+02  7.402515E+02  3.147607E+02  5.164578E+01
%   200   2000  -2.120151E+01 -3.668933E-01  7.295099E-01 -1.944827E-01 -1.263576E-01 -1.491027E-01
%\end{verbatim}\end{small}
%The Carbon table for energy losses from Post et al.
\begin{small}\begin{verbatim}
  TMAX         A(0)           A(1)            A(2)           A(3)           A(4)         A(5)
     20   1.965300E+03, 	4.572035E+03, 	4.159590E+03,   1.871560E+03, 	 4.173889E+02,  3.699382E+01, 
    200   7.467599E+01, 	4.549038E+02, 	8.372937E+02,   7.402515E+02, 	 3.147607E+02,  5.164578E+01, 
   2000  -2.120151E+01, -3.668933E-01, 	7.295099E-01,  -1.944827E-01,  -1.263576E-01, -1.491027E-01, 
  20000   2.121979E+01, -2.346986E-01,  4.093794E-01,   7.874548E-02,  -1.841379E-01,  5.590744E-02, 
 100000  -2.476796E+01,  9.408181E+00, -9.657446E+00,   4.999161E+00,  -1.237382E+00,  1.160610E-01, 
	  \end{verbatim}\end{small}
The Nitrogen table for energy losses from Post et al. is given in the following table
      \begin{small} \begin{verbatim}
      TMAX	 20,           200,          2000,         20000,        100000,
      A(0)	-1.967182E+02,-3.615155E+01,-2.093912E+01,-2.093039E+01,-9.452522E+00, 
      A(1)	-2.429049E+02,-3.943802E+01,-5.677397E-01,-6.617905E-01,-3.583144E+01,
      A(2)	-7.454123E+01,-5.564129E+00, 7.664689E-01, 1.146777E+00, 4.386446E+01,
      A(3)	 3.126366E+01, 5.140343E+01,-2.610450E-01,-7.390625E-01,-2.639331E+01,
      A(4)	 2.166881E+01, 4.369243E+01, 3.464473E-01, 3.042676E-01, 7.890268E+00,
      A(5)	 3.300054E+00, 1.027448E+01, 6.723385E-01,-6.024562E-02,-9.366682E-01,
      	\end{verbatim}  \end{small}
Note that in order to convert to units of [W m$^3$] requires multiplication of $L_{Z}$ with $10^{-13}/e_{charge}$ . %is to be subtracted from the coefficient A(0).

An alternative for the cooling rate of Carbon according to Post et al. is provided in an expression from Havlickova \cite{havlickova2013}
\begin{equation}
    L_C(T) = 2.0 \times 10^{-31} {(T/10)^3 \over 1 + (T/10)^{4.5}} [{\rm W m}^3]
\end{equation}
for $T$ in eV.


\section{Edge-localised modes}
\label{section:elms}

Edge-localised modes (ELMs) are sudden disruptions in the core plasma, which release large amounts of heat and particles. To simulate ELMs in DIV1D, the subroutine simulate\_elm is included. This subroutine adds an additional term to the boundary conditions for the parallel heat flux and the plasma density at the X-point, which simulates the short-term  heat flux and particle density surges that are characteristic of ELMs. Two models have been included. The DIV1D ELMs are characterised by a steep initial increase, and a longer decrease. The first is a model by Eich et al~\cite{eich2017}, which is a simple triangular waveform written as

\begin{equation}
\label{triangular_elm}
q_{\mathrm{ELM}}(t)=\frac{2}{3}\frac{Q_{\mathrm{ELM}}}{\tau_{\mathrm{r}}}\times
\left\{
\begin{matrix}
\frac{t}{\tau_{\mathrm{\mathrm{r}}}} &
0 \leq t \leq \tau_{\mathrm{r}} \\
1-\frac{t-\tau_{\mathrm{r}}}{2\tau_{\mathrm{r}}} &
\tau_{\mathrm{r}} < t \leq 3\tau_{\mathrm{r}}
\end{matrix}
\right.
\end{equation}

that has the useful property that

\begin{equation}
\int_{0}^{3\tau_{\mathrm{r}}} q_{\mathrm{ELM}} \mathrm{d}t=Q_{\mathrm{ELM}}.
\end{equation}

Here $\tau_{\mathrm{r}}$ is the time during which the ELM-induced heat flux grows, and $Q_{\mathrm{ELM}}$ the total fluence (joules per square meter) that is added on the upstream side of the domain over the course of the ELM. Equation~\ref{triangular_elm} is added to the steady-state heat flux ($q_{\parallel}$) as

\begin{equation}
\label{elm_boundary_condition}
q_{\parallel,\mathrm{tot}}=q_{\parallel}+q_{\mathrm{ELM}}(t).
\end{equation}

An alternative formulation to equation~\ref{triangular_elm}, to remove discontinuities in the first derivative of $q_{\mathrm{ELM}}$ with respect to time, is formulated as
\begin{equation}
\label{maxwell_boltzmann_elm}
q_{\mathrm{ELM}} = Q_{\mathrm{ELM}}\times
\frac{1}{N}\sqrt{\frac{2}{\pi}}
\sum_{n=1}^{N} \frac{t^2}{a_{n}^3}\exp{\frac{-t^2}{2a_n^2}}
\end{equation}
which is recognizable as a sum over Maxwell-Boltzmann distributions in time $t$ where the factors $a_i$ are given by
\begin{equation}
\label{maxwell_boltzmann_elm_a}
a_n = \Delta^{n-1}\frac{\tau_{\mathrm{r}}}{\sqrt{2}}.
\end{equation}

This follows from the fact that the maximum of equation~(\ref{maxwell_boltzmann_elm}) for $n=1$ is set to coincide with some ELM ramp time $\tau_{\mathrm{r}}$, which is an input parameter. Note that the ELM ramp time does not exactly reflect the maximum value for $n=1$, which is shifted to a slightly later time. Consecutive-order terms in $a_n$ differ by a multiplication factor of $\Delta$. Values of $\Delta=1.4$ and $N=2$ were chosen to match the triangular shape proposed by~\cite{eich2017}. It was initially hoped that this type of ELM would give smoother results than its triangular counterpart, but it would appear that this is not the case. This is why the triangular ELM is the default.\\

Single-pulse ELMs are not always realistic, particularly for the smaller (but higher-frequency) Type-II and Type-III ELMs. To this end, DIV1D includes the possibility of repeating ELMs on a fixed time interval.\\

Finally, the ELM functionality was written such that it can also include a density effect. Using the same equation \ref{triangular_elm} (or \ref{maxwell_boltzmann_elm}), but now considering $Q_{\mathrm{ELM}}$ to be the total amount of particles that is dispelled during an ELM we can take the time derivative of equation \ref{triangular_elm}/\ref{maxwell_boltzmann_elm}) and thus obtain the rate of change of the density due to this ELM. This effect can be summed with the density ramp feature in DIV1D to simulate the rise in the upstream density.\\

The full list of input parameters related to ELMs in DIV1D is
elm\_start\_time,
elm\_ramp\_time,
elm\_time\_between,
elm\_expelled\_heat,
elm\_expelled\_particles,
switch\_elm\_density,
switch\_elm\_heat\_flux,
switch\_elm\_series,
gaussian\_elm
explained more fully in table~2.

\section{Discretization and numerical implementation}\label{numerics}

The equations are discretized on a nonequidistant grid as employed also in the SD1D code \cite{SD1D}. For  $N$ grid cells the boundaries $x_{{\rm cb},i}$ counting from $i = 0$ at the X-point on the left to $i =N$ at the target on the right are given by \cite{SD1D}
\begin{equation}\label{cell_boundaries}
   x_{{\rm cb},i} = L \left( {(2 - \delta_{x,{\rm min}}) i \over N} - {(1 - \delta_{x,{\rm min}}) i^2 \over N^2} \right)
\end{equation}
where $L$ is the total distance from X-point to the target and $\delta_{x,{\rm min}}$ is a parameter that sets the ration between the smallest grid cell at the target to the average grid cell size. The cell centres are defined as
\begin{equation}\label{cell_centres}
    x_i = {x_{{\rm cb},i} + x_{{\rm cb},i-1} \over 2}, \quad\quad\hbox{for $i = 1 \cdots N$}.
\end{equation}
The widths of the grid cells are given by
\begin{equation}\label{cell_width}
    \Delta x_{{\rm cb},i} = x_{{\rm cb},i} - x_{{\rm cb},i-1}, \quad\quad\hbox{for $i = 1 \cdots N$}.
\end{equation}
Similarly, the distance between cell centres defines
\begin{equation}\label{deltax}
    \Delta x_i = x_{i+1} - x_{i}, \quad\quad\hbox{for $i = 1 \cdots N-1$},
\end{equation}
while for the final cell centre the distance to a virtual mirror point beyond the target is used to obtain
\begin{equation}
    \Delta x_N = 2 ( L - x_N ).
\end{equation}
All variables are calculated on the cell centres, while fluxes are calculated on the cell boundaries.

The primary variables that are evolved in the code are the plasma density $n$, the plasma momentum $P \equiv n m v_\parallel$, the total internal energy $E \equiv 3 n k T$, and the neutral density $n_n$. For evaluation of the ODE solver the normalized variables are stacked in a vector $Y$ of length $4N$ as
\begin{eqnarray}
    Y_i        &\equiv {\displaystyle n_i \over n_{\rm norm} } &\quad\quad\hbox{for $i = 1 \cdots N$}, \nonumber \\ \nonumber \\
    Y_{N + i}  &\equiv {\displaystyle P_i \over n_{\rm norm} m c_{\rm norm} } &\quad\quad\hbox{for $i = 1 \cdots N$}, \nonumber \\ \label{solution_vector} \\
    Y_{2N + i} &\equiv {\displaystyle E_i \over 3 n_{\rm norm} k T_{\rm norm} } &\quad\quad\hbox{for $i = 1 \cdots N$},  \nonumber\\ \nonumber \\
    Y_{3N + i} &\equiv {\displaystyle n_n \over n_{\rm norm} } &\quad\quad\hbox{for $i = 1 \cdots N$}, \nonumber
\end{eqnarray}
where $c_{\rm norm} = \sqrt{2 k T_{\rm norm}/m}$ is the sound speed at the normalizing temperature. Typically the normalizing density will be set equal to the initial density at the X-point $n_{\rm norm} = n_1$, while the normalizing temperature has the default value $T_{\rm norm} = 1$~eV.

The advected part of the fluxes are calculated according to the numerical scheme as used in SD1D \cite{SD1D}. First the parallel velocities on the cell boundaries are obtained from a simple average of the parallel velocities at the adjacent cell centres. The advected quantity is the obtained using a slope limiter and a upwind scheme scheme in case of supersonic flow. For subsonic flow a Lax flux is used to damp discontinuities. For details see Dudson et al. \cite{dudson2019} and references therein. Pressure gradients are discretized using a downwind scheme.


\section{Input and Output}\label{IO}

Here we provide a complete list of the input and output parameters of the code. The input is read from a file named {\tt input.txt} in the directory in which the code is executed. The input is read in the form of two FORTRAN namelists. The first namelist contains the settings of the numerical parameters for the code as listed in table~\ref{tab:input_numerics}. The second namelist contains the settings for the physics parameters as listed in table~~\ref{tab:input_physics}. Some inputs can be switched to be dynamic, which requires additional files in the execution director. These files should only contain columns with numerical values and sizes are indicated. \\

The output is written in the execution directory in {\tt div1doutput.txt}. Firstly, most inputs that were provided via {\tt input.txt} and the {\tt *.dat} files have been appended to {\tt div1doutput.txt} and represent the actual values used by DIV1D for input verification. Secondly, the outputs in terms of states: {\tt density, velocity, temperature, neutral }, fluxes: {\tt Gamma\_n, Gamma\_mom, q\_parallel, neutral\_flux }, and sources: {\tt Source\_n, Source\_v, Source\_Q, Source\_neutral} are written.

\begin{table}[h]
\begin{center}
  \caption{{\bf Namelist {\tt div1d\_numerics} setting parameters controlling the numerical implementation}.}
  \label{tab:input_numerics}
  \begin{tabular}{|| l  | l ||}
    \hline\hline
    variable name                & description \\ \hline\hline
    {\tt Nx}                     & integer: number of points in the grid along the field line (default: 1000) \\ \hline
    {\tt ntime}                  & integer: number of times steps (default: 1000)  \\ \hline
    {\tt delta\_t}               & real: time step size (default: 1.0D-06)  \\ \hline
    {\tt nout}                   & integer: output period in time steps (default: 100) \\ \hline
    {\tt method}                 & integer: set method flag for {\tt dvode} integrator (default: 227)   \\ \hline
	{\tt istate\_mod}            & integer: number of time steps between calls to {\tt dvode} with {\tt istate = 1} (default: 0) \\ \hline
	{\tt max\_step}              & integer: maximum number of internal steps in {\tt dvode} (default: 100000) \\ \hline
	{\tt max\_attempts}          & integer: maximum number of restarts of {\tt dvode} after failed integration (default: 1000) \\ \hline
    {\tt abstol}                 & real: absolute value of tolerance in numerical integration (default: 1.0D-4)  \\ \hline
    {\tt reltol}                 & real: relative value of tolerance in numerical integration (default: 1.0D-4)  \\ \hline
	{\tt evolve\_density}        & integer: multiplier for density evolution (default: 1) \\ \hline
	{\tt evolve\_momentum}       & integer: multiplier for momentum evolution (default: 1) \\ \hline
	{\tt evolve\_energy}         & integer: multiplier for energy evolution (default: 1) \\ \hline
	{\tt evolve\_neutral}        & integer: multiplier for neutral density evolution (default: 1) \\ \hline
	{\tt density\_norm}          & real: normalization of densities in solution vector {\tt Y} (default: {\tt initial\_n}) \\ \hline
	{\tt temperature\_norm}      & real: temperature used to normalize solution vector {\tt Y} (default: 1.0D0 eV) \\ \hline
	{\tt velocity\_norm}         & real: velocity used to normalize solution vector {\tt Y} (default: sound speed at {\tt temperature\_norm}) \\ \hline
	{\tt switch\_density\_source} & real: multiplier of source term in particle balance (default: 1.0D0) \\ \hline
	{\tt switch\_momentum\_source}& real: multiplier of source term in momentum balance (default: 1.0D0) \\ \hline
	{\tt switch\_energy\_source}  & real: multiplier of source term in energy balance (default: 1.0D0) \\ \hline
	{\tt switch\_neutral\_source} & real: multiplier of source term in neutral density equation (default: 1.0D0) \\ \hline
	{\tt switch\_charge\_exchange}& real: multiplier of the charge exchange rate (default: 1.0D0) \\ \hline
	{\tt switch\_recombination}  & real: multiplier of the recombination rate (default: 1.0D0) \\ \hline
	{\tt switch\_ionization}     & real: multiplier of the ionization rate (default: 1.0D0) \\ \hline
	{\tt switch\_excitation}     & real: multiplier of the excitation rate (default: 1.0D0) \\ \hline
	{\tt switch\_impurity\_radiation}& real: multiplier of the impurity radiation losses (default: 1.0) \\ \hline
	{\tt viscosity}              & real: numerical viscosity (default: 0.0D0)  \\ \hline
	{\tt central\_differencing}   & real: fraction of central differencing in pressure gradient term (default: 0.5D0) \\ \hline
	{\tt restart}                & logical: when {\tt .true.} initial conditions are read from restart file (default: {\tt .false.}) \\ \hline
    \hline
  \end{tabular}
\end{center}
\end{table}
\begin{table}[h]
\begin{center}
  \caption{{\bf Namelist {\tt div1d\_physics} setting parameters controlling the physics of the problem solved}. }
  \label{tab:input_physics}
  \begin{tabular}{|| l  | l ||}
    \hline\hline
    variable name                & description \\ \hline\hline
    {\tt L}                      & real: length of the field line between X-point and target (default: 5.0D1~m) \\ \hline
    {\tt mass}                   & real: mass of the main plasma ion (default: 3.3436D-17~kg for D) \\ \hline
    {\tt gamma}                  & real: sheath heat transmission factor (default: 6.5D0) \\ \hline
    {\tt sintheta}               & real: sinus of angle theta between B-field and target plate (default: 0.1D0) \\ \hline
    {\tt q\_parX}                & real: parallel heat flux at the X-point (default: 1.0D8~W/m$^2$)	\\
    				& (dynamic: set to -1 and provide {\tt dyn\_qpar.dat(ntime) }) \\ \hline
    {\tt flux\_expansion}        & real: expansion of total flux from X-point to target (default: 1.0D0) \\ \hline
    {\tt initial\_n}             & real: initial, homogeneous particle density cq. density at X-point (default: 1.0D20~m$^{-3}$) \\ 
    				& (dynamic: set to -1 and provide {\tt dyn\_nu.dat(ntime) }) \\ \hline
    {\tt density\_ramp\_rate}    & real: rate of change of X-point density (default: 0.0D0~/m$^3$/s) \\ \hline
    {\tt initial\_v}             & real: initial, homogeneous plasma velocity (default: 0.0D0~m/s) \\ \hline
    {\tt initial\_T}             & real: initial, homogeneous plasma temperature (default: 1.0D2~eV) \\ \hline
    {\tt initial\_a}             & real: (initial, homogeneous neutral density) and (homogeneous neutral background density) \\ 
    &  (default: 0.0D0~m$^{-3}$) (dynamic: set to -1 and provide {\tt dyn\_nb.dat(ntime) }) \\ \hline
    {\tt impurity\_Z}  & real: atomic number impurity ions (default: [6 0 0 0 0]) \\ \hline
    {\tt impurity\_concentration}  & real: concentration of impurity ions (default: [1 1 1 1 1]0.0D+1) \\ 
     & (dynamic: set to [-1 0 -1 0 0] and provide {\tt dyn\_imp\_con.dat(5,ntime) }) \\ \hline
    {\tt recycling}              & real: fraction of ion flux on target recycled as neutral (default: 1.0D0) \\ 
     & (dynamic: set to -1 and provide {\tt dyn\_rec.dat(ntime) }) \\ \hline
    {\tt redistributed\_fraction}& real: fraction of recycling neutral that is evenly distributed along the divertor leg (default: 0.0D0) \\
     & (dynamic: set to -1 and provide {\tt dyn\_red\_frc.dat(ntime) }) \\ \hline
    {\tt neutral\_residence\_time}& real: time scale in which neutral are lost from the SOL (default: 1.0D20~s) \\ \hline
    {\tt gas\_puff\_source}      & real: total particle source from gass puff along flux tube (default: 0.0D0~m$^{-2}$) \\ 
     & (dynamic: set to -1 and provide {\tt dyn\_gas.dat(ntime) }) \\ \hline
    {\tt gas\_puff\_location}    & real: location of gass puff along flux tube (default: 0.0D0~m) \\ \hline
    {\tt gas\_puff\_width}       & real: Gaussian width of gass puff region along flux tube (default: 1.0D20~m) \\ \hline
    {\tt charge\_exchange\_model}& character string selecting the charge exchange model (default: {\tt AMJUEL}) \\ \hline
    {\tt ionization\_model}      & character string selecting the ionization model (default: {\tt AMJUEL}) \\ \hline
    {\tt recombination\_model}   & character string selecting the recombination model (default: {\tt AMJUEL}) \\ \hline
    {\tt energy\_loss\_ion}      & real: effective energy loss per ionization (only used when {\tt switch\_excitation .eq. 0.0d+0})) \\ \hline
    {\tt minimum\_density}       & real: minimum value allowed for densities (default: 1.0D4~m$^{-3}$) \\ \hline
    {\tt minimum\_temperature}   & real: minimum value allowed for temperatures (default: 0.1D0~eV) \\ \hline
    {\tt elm\_start\_time}   	 & integer: outer time step at which the first ELM starts, in units of delta\_t (default: 0) \\ \hline
    {\tt elm\_ramp\_time}   	 & integer: ramp time of ELM (not exact for gaussian, see section~\ref{section:elms}), in units of delta\_t (default: 0) \\ \hline
    {\tt elm\_time\_between}    & integer: outer time steps between two the start of two ELMs, in units of delta\_t (default: 2d+08) \\ \hline
    {\tt elm\_expelled\_heat}    & real: total amount of expelled heat per unit area due to an ELM (default: 0 J m$^{-2}$) \\ \hline
    {\tt elm\_expelled\_particles}& real: total number of expelled particles due to an ELM (default: 0 m$^{-3}$) \\ \hline
    {\tt switch\_elm\_heat\_flux}& integer: turns off (0) or on (1) the ELM contribution to the heat flux (default: 0) \\ \hline
    {\tt switch\_elm\_density}   & integer: turns off (0) or on (1) the ELM contribution to the particle flux (default: 0) \\ \hline
    {\tt switch\_elm\_series}    & integer: turns off (0) or on (1) the multi-ELM sequence (default: 0) \\ \hline
    {\tt gaussian\_elm}      & integer: chooses triangular (0) or gaussian (1) ELM (default: 0) \\ \hline
    {\tt radial\_loss\_factor} & real: set fraction of qparX which is lost (except for radial\_loss\_gaussian=-1) (default: 0) \\ 
     & (dynamic: set to -1 and provide {\tt dyn\_rad\_los.dat(ntime) }) \\ \hline
    {\tt radial\_loss\_gaussian} & integer: chooses constant (0), gaussian (1) or locally dependent (-1) loss profile (default: 0) \\ \hline
    {\tt radial\_loss\_width}  & real: sets width of gaussian, only necessary if radial\_loss\_gaussian=1 (default: 1d+20) \\ \hline
    {\tt radial\_loss\_location} & real: sets peak position of gaussian, only necessary if radial\_loss\_gaussian=1 (default: 0) \\ \hline
   % {\tt switch\_dyn\_nu} & integer: switch dynamic X-point density (default: 0)   \\ \hline
   % {\tt switch\_dyn\_gas} & integer: switch dynamic gas puff source  (default: 0) (requires: {\tt dyn\_nu.dat(ntime) }) \\ \hline
   % {\tt switch\_dyn\_rec} & integer: switch dynamic recycling coefficient (default: 0) (requires: {\tt dyn\_rec.dat(ntime) })   \\ \hline
   % {\tt switch\_dyn\_rad\_los} & integer:  switch dynamic radial los factor (default: 0) (requires: {\tt dyn\_rad\_los.dat(ntime) })  \\ \hline
   % {\tt switch\_car\_con\_prf} & integer: switch carbon concentration profile (default: 0) (requires: {\tt car\_con\_prf.dat(Nx) }) \\ \hline
   % {\tt switch\_dyn\_qpar} & integer: switch dynamic X-point heat flux (default: 0) (requires: {\tt dyn\_qpar.dat(ntime) })  \\ \hline
   % {\tt switch\_dyn\_red\_frc} & integer: switch dynamic redistribution fraction (default: 0)  (requires: {\tt dyn\_red\_frc.dat(ntime) })  \\ \hline
    \hline
  \end{tabular}
\end{center}
\end{table}

\cleardoublepage
\section{Benchmark with 2 Point Model}
The 2 Point Model for an attached divertor plasma assuming purely conductive heat transport and no losses is discussed in Chapter 5.2 of~\cite{stangeby}). For a given parallel heat flux $q_{\parallel,\rm X}$ entering the divertor leg upstream, a given upstream plasma density $n_{\rm X}$, the length of the divertor leg $L$, and neglecting possible effects of flux expansion a set of three equations is obtained from the energy balance, the sheath heat transmission, and pressure balance respectively~\cite{stangeby}:
\begin{equation}\label{2PM_energy_balance}
    T_{\rm X} = \left(T_{\rm L}^{7/2} + {7 \over 2} {q_{\parallel,\rm X} L \over \kappa_0}\right)^{2/7},
\end{equation}
\begin{equation}\label{2PM_sheath_heat_transmission}
    q_{\parallel,\rm X} = \gamma n_{\rm L} e T_{\rm L} c_{\rm s},
\end{equation}
\begin{equation}\label{2PM_pressure_balance}
    2 n_{\rm L} T_{\rm L} = n_{\rm X} T_{\rm X}.
\end{equation}
with $\kappa_0 = 2000$~W/eV$^{7/2}$m. These equations can be used to solve for the three unknowns, upstream temperature $T_{\rm X}$, target temperature $T_{\rm L}$ and target density $n_{\rm L}$. The solution can be obtained iteratively starting from $T_{\rm L} = 0$. In the first iteration the solution for $T_{\rm X}$ then is
\begin{equation}\label{2PM_upstream_temperature}
    T_{\rm X} = \left({7 \over 2} {q_{\parallel,\rm X} L \over \kappa_0}\right)^{2/7}.
\end{equation}
Next an equation for $T_{\rm L}$ is obtained by substituting equation (\ref{2PM_pressure_balance}) in the sheath heat transmission (\ref{2PM_sheath_heat_transmission}):
\begin{equation}\label{2PM_target_temperature}
    T_{\rm L} = {m \over e} {2 q_{\parallel,\rm X}^2 \over \gamma^2 e^2 n_{\rm X}^2 T_{\rm X}^2}.
\end{equation}
Finally the pressure balance (\ref{2PM_pressure_balance}) then allows to obtain the target density as
\begin{equation}\label{2PM_target_density}
    n_{\rm L} = {\gamma^2 e^3 n_{\rm X}^3 T_{\rm X}^3 \over 4 m q_{\parallel,\rm X}^2}.
\end{equation}
In fact, these last three equations represent the simple 2 Point Model as obtained under the condition that $T_{\rm X} \gg T_{\rm L}$ \cite{stangeby}.

An extension of the simple 2 Point Model is obtained introducing the effects of a finite convective contribution to the parallel heat transport, as well as power and momentum losses (see Chapter 5.4 of~\cite{stangeby}). Assuming that the convective contribution is evenly distributed over the distance between the X-point and target, the solution to the energy balance equation \ref{2PM_energy_balance} is simply modified as
\begin{equation}\label{extended_2PM_energy_balance}
    T_{\rm X} = \left(T_{\rm L}^{7/2} + {7 \over 2} {q_{\parallel,\rm X} (1 - f_{\rm conv}) L \over \kappa_0}\right)^{2/7},
\end{equation}
introducing the factor $f_{\rm conv}$ representing the relative contribution from convection to the heat transport. When the power and momentum losses are not negligible yet sufficiently localized near the target, the energy balance solution for the upstream temperature (\ref{extended_2PM_energy_balance}) is not altered further, while the sheath heat transmission is simply modified by introducing the factor $(1-f_{\rm pwr})$ on the left hand side of equation (\ref{2PM_sheath_heat_transmission}) and the pressure balance is modified similarly by introducing the factor $(1-f_{\rm mom})$ on the right hand side of equation (\ref{2PM_pressure_balance}) to account for the power an momentum losses, respectively. As a result the equation for the target temperature (\ref{2PMF_target temperature}) is modified to
\begin{equation}\label{extended_2PM_target_temperature}
    T_{\rm L} = {m \over e} {2 q_{\parallel,\rm X}^2 (1-f_{\rm pwr})^2 \over \gamma^2 e^2 n_{\rm X}^2 T_{\rm X}^2 (1-f_{\rm mom})^2}.
\end{equation}
while the target density is modified to
\begin{equation}\label{extended_2PM_target_density}
    n_{\rm L} = {\gamma^2 e^3 n_{\rm X}^3 T_{\rm X}^3 (1-f_{\rm mom})^3 \over 4 m q_{\parallel,\rm X}^2 (1-f_{\rm pwr})^2}.
\end{equation}
The power and momentum loss fractions $f_{\rm pwr}$ and $f_{\rm mom}$ are defined as
%\begin{equation}\label{power_loss_fraction}
%    f_{\rm pwr} \equiv \int_{x=L}^0 Q {\rm d}x \Big/ q_{\parallel,\rm X},
%\end{equation}
\begin{equation} \label{power_loss_fraction}
f_{\rm pwr} \equiv \left(1-\frac{q_{\mathrm{\|,L}} }{q_{\mathrm{\|,X}}} \frac{B_{\mathrm{X}}}{B_{\mathrm{L}}}\right)=\int_{x=L}^{0}\left[\frac{v_{\|}}{B} \frac{\partial p}{\partial x}+\frac{Q}{B}\right] d x \frac{B_{\rm X}}{q_{\mathrm{\|,X}}}
\end{equation}
and
%\begin{equation}\label{momentum_loss_fraction}
%    f_{\rm mom} \equiv \int_{x=L}^0 S_{\rm mom} {\rm d}x \Big/ p_{\rm tot,X},
%\end{equation}
\begin{equation}\label{momentum_loss_fraction}
f_{\rm mom} \equiv 1-\frac{p_{\rm tot,L}}{p_{\rm tot,X}}=\int_{x=L}^{0}\left(-n m v_{\|}^{2} \frac{\partial B}{\partial x}+S_{\rm mom}\right) d x / p_{\rm tot,X}
\end{equation}

respectively \cite{stangeby2018}. Also the equations for the extended 2 Point Model (\ref{extended_2PM_energy_balance}-\ref{extended_2PM_target_density}) are easily solved by iteration when the convective fraction, and the power and momentum loss fractions are known.

Also flux expansion can be accounted for in the 2 Point model. In case of finite flux expansion, neglecting as usual the convective transport and assuming that all sources and sinks are well localized near the target the stationary energy equation becomes
\begin{equation}
    {q_\parallel \over B} = {q_{\parallel, \rm X} \over B_{\rm X}}.
\end{equation}
The solution of this equation now results in
\begin{equation}\label{extended_2PM_energy_balance_FE}
    T_{\rm X} = \left(T_{\rm L}^{7/2} + {7 \over 2} {q_{\parallel,\rm X} (1 - f_{\rm conv}) L  \Phi_{\rm f} \over \kappa_0}\right)^{2/7},
\end{equation}
introducing the flux expansion integral $\Phi_{\rm f}$ defined as
\begin{equation}\label{FE_integral}
    \Phi_{\rm f} \equiv \int_0^1 {B\over B_{\rm X}} {\rm d}{\bar x}.
\end{equation}
Here, ${\bar x} \equiv x /L$ is the normalized $x$ coordinate. Using that $B \propto 1/R$ and assuming that the major radius varies linearly along the divertor leg this integral evaluates to
\begin{equation}\label{FE_integral_}
    \Phi_{\rm f} = \int_0^1 {1 \over 1 + \varepsilon_{\rm f} {\bar x}} {\rm d}{\bar x} = {\ln( 1 + \varepsilon_{\rm f} ) \over \varepsilon_{\rm f}},
\end{equation}
with the flux expansion parameter $\varepsilon_{\rm f}$ defined by $\varepsilon_{\rm f} \equiv (R_{\rm L}/R_{\rm X}) - 1$. A second modification to the extended 2 Point Model equations comes from the sheath boundary condition, which accounting to the flux expansion near the target is modified to
\begin{equation}\label{2PM_sheath_heat_transmission_FE}
    q_{\parallel,\rm X} {(1-f_{\rm pwr}) \over (1 + \varepsilon_{\rm f})} = \gamma n_{\rm L} e T_{\rm L} c_{\rm s}.
\end{equation}
This results in a modification of the equation for the target temperature to
\begin{equation}\label{extended_2PM_target_temperature_FE}
    T_{\rm L} = {m \over e} {2 q_{\parallel,\rm X}^2 (1-f_{\rm pwr})^2 \over (1 + \varepsilon_{\rm f})^2 \gamma^2 e^2 n_{\rm X}^2 T_{\rm X}^2 (1-f_{\rm mom})^2},
\end{equation}
and of the target density is to
\begin{equation}\label{extended_2PM_target_density_FE}
    n_{\rm L} = {\gamma^2 e^3 n_{\rm X}^3 T_{\rm X}^3 (1-f_{\rm mom})^3 (1 + \varepsilon_{\rm f})^2 \over 4 m q_{\parallel,\rm X}^2 (1-f_{\rm pwr})^2}.
\end{equation}
Again these equations (\ref{extended_2PM_energy_balance_FE}), (\ref{extended_2PM_target_temperature_FE}), and (\ref{extended_2PM_target_density_FE}) are easily solved iteratively.

For the case that the losses along the divertor leg and/or the convective heat flux are not negligible, Kotov and Reiter \cite{kotov2009} derived an exact 2 Point Model Formulation for the target temperature and density, which under the conditions of the present model (equal electron and ion temperature, sound speed at the sheath boundary, and constant magnetic field) is simplified to
\begin{equation}\label{2PMF_target temperature}
    T_{\rm L} = {8m \over e \gamma^2} {q_{\parallel,\rm X}^2 \left(1-f_{\rm pwr}\right)^2 \over p_{\rm tot,X}^2 \left(1-f_{\rm mom}\right)^2}
\end{equation}
\begin{equation}\label{2PMF_target density}
    n_{\rm L} = {\gamma^2 \over 32 m} {p_{\rm tot,X}^3 \left(1-f_{\rm mom}\right)^3 \over q_{\parallel,\rm X} \left(1-f_{\rm pwr}\right)^2}
\end{equation}
where $p_{\rm tot,X} = 2 n_{\rm X} e T_{\rm X} + n_{\rm X} m v_{\parallel\rm X}^2$ is the total upstream pressure. When the upstream parallel momentum is negligible, these equations are identical to the extended 2 Point Model equations (\ref{extended_2PM_target_temperature}) and (\ref{extended_2PM_target_density}) for the target temperature and density. These equations can serve as an accuracy check for the numerical code.


\subsection{Analytical Extended 2 Point Model}

When one would have simple closed analytical expressions for the power and momentum loss fraction, an instructive version of the extended 2 Point Model would be obtained. Such a model might provide further insight into conditions for divertor detachment and can be used for extremely fast calculations as needed in detachment control. Noting that both power losses and momentum losses typically increase with decreasing target temperature or increasing target density, we suggest the simple scaling
\begin{equation}\label{scaling_of_losses}
    f_{\rm pwr, mom} = exp( - \alpha_{\rm pwr, mom} T_{\rm L} / n_{\rm L} ),
\end{equation}
where the coefficients $\alpha_{\rm pwr, mom}$ must be chosen appropriately. Note that detachment (a decreasing target temperature) requires that $f_{\rm mom} < f_{\rm pwr}$, i.e. $\alpha_{\rm mom} < \alpha_{\rm pwr}$.

Some other remarks on the solutions of the 1D problem. Momentum dissipation requires a finite velocity, which only appears near the target where a finite particle source from ionization must be compensated by a convective particle flux towards the target. When recycling is incomplete or in case of efficient pumping in the divertor region, a significant particle flux must be present at the X-point already and there will be much more possibilities for momentum dissipation along the entire divertor leg. Impurity puffing will mostly effect power dissipation without significant momentum dissipation. In contrast, increasing the neutral density in the divertor by normal gass puff will increase the charge exchange rate and thereby increases momentum dissipation. Note that increasing the length of the divertor leg increases the X-point temperature. This must have a limit somewhere which sets a maximum distance from the X-point to the recombination front, when the divertor leg length becomes larger that this the detachment front must separate from the target. The question then arises about the stability of the detachment front position. How is this determined??


\section{Further analysis of the 1D divertor model: stationary solutions}

Important insight into the existence and stability of stationary solutions can be obtained by some further assumptions concerning the 1D model. In particular, ionization and recombination are often located in a relatively narrow zone. For stationary solutions of the equations, the pressure can then be assumed to be constant outside this narrow ionization and recombination zone and the parallel velocity and the convective heat flux can be neglected. The stationary solutions are then entirely specified by the stationary energy balance
\begin{equation}\label{stationary_energy_balance}
    {\partial \over \partial x} q_\parallel = Q \quad\hbox{with}\quad q_\parallel = - \kappa_\parallel {\partial \over \partial x} T.
\end{equation}
Using a transformation first introduced by Wim Goedheer [check citation] in the context of radiation losses from the edge of the confined plasma, allows to recast this equation into an integral equation for the change in the heat flux over the divertor leg or, in other words, an integral equation for the radiative energy losses in the divertor region. To this end Eq.~({\ref{stationary_energy_balance}) is multiplied by the parallel heat flux and integrated to obtain
\begin{equation}\label{flux_change}
    q_X^2 - q_L^2 = - 2 \int_{T_L}^{T_X} \kappa_\parallel Q {\rm d}T.
\end{equation}
Next, the parallel heat conductivity is substituted as $\kappa_\parallel = \kappa_0 T^{5/2}$, the energy loss term is taken in the form of the radiative impurity losses (\ref{impurity_radiation}), and the constant pressure, $ p = 2 neT$, approximation is applied to obtain~\cite{lengyel1981,capes1992,hutchinson1994}
\begin{equation}\label{flux_change_2}
    q_X^2 - q_L^2 = {1\over2 e^2} \kappa_0 p^2 \xi_Z \int_{T_L}^{T_X} \sqrt{T} L_Z(T) {\rm d}T,
\end{equation}
where the impurity concentration $xi_Z$ has also been assumed to be constant over the radiative layer.

A number of conclusions can be drawn from this equation. Several authors apply it to quantify the maximum flux that can be radiated away in the divertor leg \cite{lengyel1981,lackner1993,kallenbach2013,siccinio2016}
\begin{equation}\label{maximum_radiated_flux}
    q_{\rm max} = \sqrt{{1\over2 e^2} \kappa_0 p^2 \xi_Z \int_0^\infty \sqrt{T} L_Z(T) {\rm d}T}.
\end{equation}
Hutchinson and Lipschultz use it to formulate the stability condition for a radiative front~\cite{hutchinson1994,lipschultz2016}:
\begin{equation}
    {{\rm d} \over {\rm d}x_f} (q_i - q_{\rm max}) \le 0,
\end{equation}
where $x_f$ is the front position and $q_i$ is the parallel energy flux entering the front. We can apply this to analyse the stability of a detachment front. In our model, $q_i$ is fixed by the boundary condition at the X-point, while the pressure is determined by the equilibrium condition that $q_X = q_{\rm max}$. For a given X-point density then also the X-point temperature is determined, which through the heat conductivity equation sets the position of the radiation front (i.e. the position where the temperature is reached at which $L_Z(T)$ is maximum). Perturbing the front position in the direction of the X-point would lower the X-point temperature and, consequently the X-pint temperature and the pressure in the divertor. This will reduce the losses in the radiation front and push it back to its original position. Similarly a movement in the direction of the target will result in an increase of the losses again pushing the radiative layer back to its original position. As a result, given the boundary conditions in our model, the position of the radiation front for a fully detached case is always stable. We should add here that the radiation losses in our model also include significant contributions from ionization and recombination. However, because for a given impurity concentration the detachment front in our model becomes effectively independent of the position of the front so that the analysis given above also applies when considering the total losses which effectively have become a function of the temperature only.

Capes et al. use this integral relation to study the stability of attached or semi-detached divertor solutions for which $q_X$ is larger than $q_{\rm max}$~\cite{capes1992}. The residual energy flux to the target then determines the properties at the target according to the sheath boundary condition. Substituting the target heat flux (\ref{sheath_heat_transmission}) in Eq.~(\ref{flux_change_2}), while taking into account the condition (\ref{2PM_pressure_balance}) that follows from the acceleration of the plasma in the presheath to the sound velocity, then results after some algebraic manipulation in a nonlinear relation between the target temperature $T_L$ to the conditions at the X-point:
\begin{equation}
    X(T_L) \equiv {8 m q_X^2 \over e \gamma^2 p^2} = T_L + {4 \kappa_0 m \xi_Z \over e^3 \gamma^2} \int_{T_L}^\infty \sqrt{T} L_Z(T) {\rm d}T,
\end{equation}
where for convenience the upper limit of the integration has been extended to infinity. Generalizing this to include all radiative losses this relation could be written as
\begin{equation}
    X(T_L) \equiv {8 m q_X^2 \over e \gamma^2 p^2} = T_L - {16 m \over e \gamma^2 p^2}  \int_{T_L}^\infty \kappa_\parallel Q {\rm d}T.
\end{equation}
Note that $X(T_L)$ is identical to the two-point-model target temperature (\ref{2PM_target_temperature}) in case there are no radiative losses. In case radiative losses are not negligible the integral on the right hand side is related to the dissipated power fraction in the two point model. For a given condition $X$ at the X-point, there can be one, multiple, or no solution for $T_L$. Bifurcation points can be found by differentiating $X(T_L)$ with respect to $T_L$ and finding the zeros of the resulting relation, i.e.
\begin{equation}
    {{\rm d} X(T_L) \over {\rm d}T_L} = 1 - {4 \kappa_0 m \xi_Z \over e^3 \gamma^2} \sqrt{T_L} L_Z(T_L) = 0,
\end{equation}
or
\begin{equation}
    {{\rm d} X(T_L) \over {\rm d}T_L} = 1 + {16 m \over e \gamma^2 p^2} [\kappa_\parallel Q](T_L) = 0.
\end{equation}
The impurity concentration $\xi_Z$ above which the solutions for the stationary divertor solution are bifurcated is then determined by the maximum of $\sqrt{T} L_Z(T)$
\begin{equation}
    \xi_\star = {e^3 \gamma^2 \over 4 \kappa_0 m [\sqrt{T} L_Z(T)]_{\rm max}}.
\end{equation}
Substituting the numbers from our simulations would predict a critical Carbon concentration of about 2.5\% [MUST BE CHECKED]. In case of a 5\% Carbon concentration our simulations indeed confirm the existence of a bifurcation. However, impurity radiation losses are not the only losses that must be accounted for. Also ionization and charge exchange losses are significant and should be accounted for. Apart from leading to energy losses these are also responsible for momentum losses that are not accounted for in the analysis given above. The simulations show tat the momentum losses are non-negligible in the parameter range where the bifurcation is observed.

Note that in their paper Capes et al. considered only the electron fluid and neglected the effect of the acceleration of the plasma in the presheath resulting in slightly different expressions. Also the values used for the parallel heat conductivity $\kappa_0$ and the heat transmission coefficient $\gamma$ differ from the present work.

\section*{Acknowledgments}
\noindent This project was carried out with financial support from NWO. The work has been carried out within the framework of the EUROfusion Consortium and has received funding from the Euratom research and training programme 2014-2018 under grant agreement No 633053. The views and opinions expressed herein do not necessarily reflect those of the European Commission.

%\section*{References}
\begin{thebibliography}{99}
\bibitem{nakazawa2000}
Shinji Nakazawa, Noriyoshi Nakajima, Masao Okamoto and Nobuyoshi Ohyabu, Plasma Phys. Control. Fusion {\bf 42} (2000) 401.

\bibitem{dudson2019}
B.D. Dudson, J. Allen, T. Body, B. Chapman, C. Lau, L. Townley, D. Moulton, J. Harrison and B. Lipschultz, Plasma Phys. Control. Fusion {\bf 61} (2019) 065008.

\bibitem{SD1D}
Ben Dudson, {\it SD1D: 1D divertor model for detachment studies}, 19 December 2016.

\bibitem{stangeby}
Peter C. Stangeby, 2000, {\it The Plasma Boundary of Magnetic Fusion Devices}, Institute of Physics Publishing, Dirac House, Temple Back, Bristol BS1 6BE,
UK.

\bibitem{stangeby2018}
P.C. Stangeby, Plasma Phys. Control. Fusion {\bf 60} (2018) 44022. 


\bibitem{EIRENE}
{\tt http://www.eirene.de/}.

\bibitem{havlickova2013}
E. Havlickova, et al., Plasma Phys. Control. Fusion {\bf 55} (2013) 065004.

\bibitem{freeman1974}
R.L. Freeman, and E.M. Jones, {\it 'Atomic collision processes in plasma physics experiments'}, Culham Report CLM-R137

\bibitem{gordeev1977}
Yu.S. Gordeev, A.N. Zinov'ev, and M.P. Petrov, J. Exp. Theor. Physics {\bf 25} (1977) 204.

\bibitem{post1977}
D.E. Post, et al., Atomic Data and Nuclear Data Tables {\bf 20} (1977) 397.

\bibitem{veres2009}
G. Veres, et al., J. Nucl. Mater. {\bf 390–391} (2009) 835.

\bibitem{eich2017}
T. Eich, et al., Nuclear Materials and Energy {\bf 12} (2017) 84-90l

\bibitem{hinnov1962}
E. Hinnov,. and J.G. Hirschberg, Phys. Review {\bf 125} (1962) 795.

\bibitem{kotov2009}
V. Kotov and D. reiter, Plasma Phys. Control. Fusion {\bf 51} (2009) 115002.

\bibitem{lengyel1981}
L.L. Lengyel, {\it 'Analysis of radiation Plasma Boundary Layers'}, Report IPP 1/191 (1981).

\bibitem{capes1992}
H. Capes, Ph. Ghendrih, and A. Samain, Phys. Fluids B {\bf 4} (1992) 1287.

\bibitem{hutchinson1994}
I.H. Hutchinson, Nucl. Fusion {bf 34} (1994) 1337.

\bibitem{lackner1993}
K. Lackner, and R. Schneider, Fusion Eng. Design {\bf 22} (1993) 107.

\bibitem{kallenbach2013}
A. Kallenbach, et al., Plasma Phys. Control. Fusion {\bf 55} (2013) 124041.

\bibitem{siccinio2016}
M. Siccinio, et al., Plasma Phys. Control. Fusion {\bf 58} (2016) 125011.

\bibitem{lipschultz2016}
B. Lipschultz, F.I. Parra, and I.H. Hutchinson, Nucl. Fusion {\bf 56} (2016) 056007.

\end{thebibliography}

\end{document}

\begin{eqnarray}
    &&2n_t T_t = n_u T_u ( 1 - f_{\rm mom}) \nonumber \\
    &&T_u^{7/2} = T_t^{7/2} + {\displaystyle {7 \over 2} {q_\parallel (1 - f_{\rm conv}) L \over \kappa_{0e}} } \nonumber \\
    &&q_\parallel (1 - f_{\rm pwr}) = \gamma n_t k_B T_t c_{st} \nonumber
\end{eqnarray}
  This explains that we find a bifurcation in the stationary divertor solutions already for Carbon impurity concentrations of 5\%.
