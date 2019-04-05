module constants
implicit none
integer, parameter, private :: wp = KIND(1.0D0)
real( wp ), parameter :: e_charge = 1.6022d-19 ! elementary charge [Coulomb]
real( wp ), parameter :: c        = 2.9979d+8  ! velocity of light [m/s]
real( wp ), parameter :: K_B      = 1.3807d-23 ! Boltzmann constant [Coulomb]
real( wp ), parameter :: amu      = 1.6605d-27 ! atomic mass unit [kg]
real( wp ), parameter :: me       = 9.1094d-31 ! electron mass [kg]

end module constants
