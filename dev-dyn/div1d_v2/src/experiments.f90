module experiments

    use constants, only : pi
    use physics_parameters, only : elm_start_time, elm_ramp_time, elm_time_between ,elm_expelled_heat, elm_expelled_particles, &
                                   switch_elm_density, switch_elm_heat_flux, switch_elm_series, gaussian_elm, &
                                   q_parX, L, radial_loss_factor, radial_loss_gaussian, radial_loss_width, radial_loss_location
    use numerics_parameters, only: delta_t
    use grid_data, only: x, delta_xcb

    implicit none
    integer, parameter, private :: wp = KIND(1.0D0)

contains

    subroutine simulate_elm(elm_heat_load, elm_density_change, time)

        ! This subroutine translates the ELM that comes out of the maxwell_boltzmann_elm() subroutine into
        ! parameters that can be added to the boundary conditions for density and heat flux at the X-point
        ! in the right_hand_side routine. Density and heat fluxes can be turned on independently through a
        ! pair of switches, which should be set to either 0 or 1. Two if statements test against this prior to execution
        ! to save computational time. The elm_count and time_mod combined with the while loop act as a counter to keep
        ! track of how many ELMs are currently in the system. Setting switch_elm_series (accepts 0 and 1) to 0 zero
        ! effectively kills the while loop at the end of the first iteration, resulting in only a single ELM in that case.

        implicit none
        ! Input parameters
        real(wp)                    :: time
        ! Internal parameters
        integer                     :: elm_count
        real(wp)                    :: time_relative, time_mod
        ! Parameters required by maxwell_boltzmann_elm()
        real(wp)                    :: time_integrated_quantity, quantity, dquantity
        ! Ouput parameters
        real(wp)                    :: elm_heat_load, elm_density_change

        ! Check against illegal input (which only happens if a user actively alters either of the switches,
        ! since they are set to zero by default in physics_parameters).

        if(switch_elm_density.ne.0.and.switch_elm_density.ne.1) then
            print *,"Error: switch_elm_density only accepts 0 or 1 as inputs. Terminating div1d.exe."
            call exit
        elseif(switch_elm_heat_flux.ne.0.and.switch_elm_heat_flux.ne.1) then
            print *,"Error: switch_elm_heat_flux only accepts 0 or 1 as inputs. Terminating div1d.exe."
            call exit
        elseif(switch_elm_series.ne.0.and.switch_elm_series.ne.1) then
            print *,"Error: switch_elm_series only accepts 0 or 1 as inputs. Terminating div1d.exe."
            call exit
        endif

        time_relative = time - elm_start_time*delta_t
        time_mod = time_relative
        elm_heat_load = 0d0
        elm_density_change = 0d0
        elm_count = 0

        do while (time_mod.gt.0d0)

            if (switch_elm_heat_flux.eq.1) then
                time_integrated_quantity = elm_expelled_heat
                if (gaussian_elm.eq.1) then
                    call maxwell_boltzmann_elm(time_mod,time_integrated_quantity,quantity,dquantity)
                else
                    call triangular_elm(time_mod,time_integrated_quantity,quantity,dquantity)
                endif
                elm_heat_load = elm_heat_load + quantity
                if (isnan(elm_heat_load)) then
                    elm_heat_load = 0d0
                endif
            endif

            if (switch_elm_density.eq.1) then
                time_integrated_quantity = elm_expelled_particles
                if (gaussian_elm.eq.1) then
                    call maxwell_boltzmann_elm(time_mod,time_integrated_quantity,quantity,dquantity)
                else
                    call triangular_elm(time_mod,time_integrated_quantity,quantity,dquantity)
                endif
                elm_density_change = elm_density_change + dquantity
                if (isnan(elm_density_change)) then
                    elm_density_change = 0d0
                endif
            endif

            elm_count = elm_count + 1
            time_mod = (time_relative - elm_count*elm_time_between*delta_t)*switch_elm_series

        end do

    end subroutine simulate_elm

    subroutine triangular_elm(time_relative,time_integrated_quantity,quantity,dquantity)

        ! Time-dependent ELM model from Eich et al (2017) [https://doi.org/10.1016/j.nme.2017.04.014]
        ! that simulates ELMs as triangular waveforms.

        implicit none

        real(wp)        :: t, tr, norm, quantity, dquantity, time_relative, time_integrated_quantity

        t = time_relative
        tr = elm_ramp_time*delta_t
        norm = time_integrated_quantity*(2d0/3d0)*(1/tr)

        if (t.lt.tr) then
            quantity = norm * t/tr
            dquantity= norm * 1/tr
        elseif (t.gt.tr.and.t.le.3d0*tr) then
            quantity = norm * (1d0-(t-tr)/(2d0*tr))
            dquantity= norm * -1d0/(2d0*tr)
        endif
        
    end subroutine triangular_elm

    subroutine maxwell_boltzmann_elm(time_relative,time_integrated_quantity,quantity,dquantity)

        ! This subroutine calculates the temporal evolution of an ELM at a fixed point in space,
        ! typically at the X-point boundary, as the sum of four Maxwell-Boltzmann distributions,
        ! MB0(t) through MB3(t), which are normalised and then multiplied by the input parameter
        ! time_integrated_quantity. It takes and elm_ramp_time and delta_t to calculate the steepness  
        ! of the initial growth. The mb2mb_distance is a fixed scaling factor that helps with keeping
        ! the ELM relatively smooth, minimising bumps in the tail of the profile.
        !
        ! The MBs are used in calculating, for example, the heat flux boundary q_parX in the 
        ! right_hand_side subroutine when the time_integrated_quantity is chosen to reflect to total
        ! energy that is dumped into the flux tube. The analytical derivates of the MBs, c.q. the dMBs, 
        ! are used in that same subrountine to calculate the density boundary, which requires the
        ! change in density rather than the absolute density, when time_integrated_quantity is set to
        ! the total number of particles dumped into the flux tube.

        implicit none

        real(wp)        :: t, tr, norm, mb2mb_distance
        real(wp)        :: a0, a1, a2, a3, MB0, MB1, MB2, MB3, dMB0, dMB1, dMB2, dMB3
        real(wp)        :: time_relative, time_integrated_quantity, quantity, dquantity
        integer         :: istep
        
        ! t = time_relative
        ! tr = elm_ramp_time*delta_t
        ! norm    = 1d0/4d0*sqrt(2d0/pi)
        ! mb2mb_distance = 1.9d0

        ! a0      = mb2mb_distance**0d0*tr/sqrt(2d0)
        ! a1      = mb2mb_distance**1d0*tr/sqrt(2d0)
        ! a2      = mb2mb_distance**2d0*tr/sqrt(2d0)
        ! a3      = mb2mb_distance**3d0*tr/sqrt(2d0)
        ! MB0     = norm/a0**3d0*t**2d0*exp(-t**2d0/(2d0*a0**2d0))
        ! MB1     = norm/a1**3d0*t**2d0*exp(-t**2d0/(2d0*a1**2d0))
        ! MB2     = norm/a2**3d0*t**2d0*exp(-t**2d0/(2d0*a2**2d0))
        ! MB3     = norm/a3**3d0*t**2d0*exp(-t**2d0/(2d0*a3**2d0))
        ! quantity    = (MB0 + MB1 + MB2 + MB3)*time_integrated_quantity

        ! dMB0    = (2d0/t-t/a0**2d0)*MB0
        ! dMB1    = (2d0/t-t/a1**2d0)*MB1
        ! dMB2    = (2d0/t-t/a2**2d0)*MB2
        ! dMB3    = (2d0/t-t/a3**2d0)*MB3
        ! dquantity   = (dMB0 + dMB1 + dMB2 + dMB3)*time_integrated_quantity

        t = time_relative
        tr = elm_ramp_time*delta_t
        norm    = 1d0/2d0*sqrt(2d0/pi)
        mb2mb_distance = 1.4d0

        a0      = mb2mb_distance**0d0*tr/sqrt(2d0)
        a1      = mb2mb_distance**1d0*tr/sqrt(2d0)
        MB0     = norm/a0**3d0*t**2d0*exp(-t**2d0/(2d0*a0**2d0))
        MB1     = norm/a1**3d0*t**2d0*exp(-t**2d0/(2d0*a1**2d0))
        quantity    = (MB0 + MB1)*time_integrated_quantity

        dMB0    = (2d0/t-t/a0**2d0)*MB0
        dMB1    = (2d0/t-t/a1**2d0)*MB1
        dquantity   = (dMB0 + dMB1)*time_integrated_quantity


    end subroutine maxwell_boltzmann_elm

!    subroutine calculate_radial_losses(Nx,radial_sink,q_parallel)
        
        ! This subroutine captures the radial losses as a volumetric energy sink with a gaussian
        ! profile. Inputs are the gaussian width and peak location, given by radial_loss_width and 
        ! radial_loss_location respectively. Depending on whether radial_loss_gaussian is
        ! positive, zero or negative, the radial loss profile is a bell curve, constant or dependent
        ! on the local heat flux, respectively. The normalisation of the gaussian is calculated numerically,
        ! so that the total lost heat flux is always fixed by q_parX.

!        implicit none
!        integer         :: Nx
!        real(wp)        :: radial_sink(Nx), a0, x0, norm, gaussian(Nx), normalisation,q_parallel(Nx)

!        if (radial_loss_gaussian.gt.0) then
!            a0 = radial_loss_width
!            x0 = radial_loss_location
!            gaussian = exp(-(x-x0)**2/(2*a0**2))
!            normalisation = sum(gaussian * delta_xcb)
!            radial_sink = radial_loss_factor * q_parX * gaussian / normalisation
!        elseif (radial_loss_gaussian.lt.0) then
!            radial_sink = radial_loss_factor *q_parallel / L
!        else
!            radial_sink = radial_loss_factor * q_parX / L
!        endif

!    end subroutine calculate_radial_losses
        

end module experiments
