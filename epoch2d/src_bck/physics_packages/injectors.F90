! Copyright (C) 2009-2019 University of Warwick
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

MODULE injectors

  USE shared_data
  USE partlist
  USE particle_temperature
  USE evaluator
  USE random_generator
  USE utilities

  IMPLICIT NONE

  REAL(num) :: flow_limit_val = 10.0_num

CONTAINS

  SUBROUTINE init_injector(boundary, injector)

    INTEGER, INTENT(IN) :: boundary
    TYPE(injector_block), INTENT(INOUT) :: injector

    injector%npart_per_cell = -1.0_num
    injector%species = -1
    injector%boundary = boundary
    injector%t_start = 0.0_num
    injector%t_end = t_end
    injector%has_t_end = .FALSE.
    injector%density_min = 0.0_num
    injector%density_max = HUGE(1.0_num)
    injector%use_flux_injector = .TRUE.
    injector%v_inject_z = 0.0_num
    injector%power_law = -1.0_num
    injector%g1 = 10.0_num*1e3*q0 !10*kev - 1Mev
    injector%g2 = 1.0e3_num*1e3*q0
    injector%wein_temperature = -1.0_num
    injector%monoenergy = -1.0_num
    injector%angle_theta = -100.0_num
    injector%angle_phi = 0.0_num
    injector%use_time_function = .FALSE.
    NULLIFY(injector%depth)
    NULLIFY(injector%depth_xy)
    NULLIFY(injector%next)

    IF (boundary == c_bd_x_min .OR. boundary == c_bd_x_max) THEN
      ALLOCATE(injector%depth(1-ng:ny+ng))
      injector%depth = 1.0_num
    END IF

    IF (boundary == c_bd_y_min .OR. boundary == c_bd_y_max) THEN
      ALLOCATE(injector%depth(1-ng:nx+ng))
      injector%depth = 1.0_num
    END IF

    IF (boundary == c_bd_z_min .OR. boundary == c_bd_z_max) THEN
      ALLOCATE(injector%depth_xy(1-ng:nx+ng,1-ng:ny+ng))
      injector%depth_xy = 1.0_num
    END IF

    need_random_state = .TRUE.

  END SUBROUTINE init_injector



  SUBROUTINE attach_injector(injector)

    TYPE(injector_block), POINTER :: injector
    INTEGER :: boundary

    boundary = injector%boundary

    IF (boundary == c_bd_x_min) THEN
      CALL attach_injector_to_list(injector_x_min, injector)
    ELSE IF (boundary == c_bd_x_max) THEN
      CALL attach_injector_to_list(injector_x_max, injector)
    ELSE IF (boundary == c_bd_y_min) THEN
      CALL attach_injector_to_list(injector_y_min, injector)
    ELSE IF (boundary == c_bd_y_max) THEN
      CALL attach_injector_to_list(injector_y_max, injector)
    ELSE IF (boundary == c_bd_z_min) THEN
      CALL attach_injector_to_list(injector_z_min, injector)
    ELSE IF (boundary == c_bd_z_max) THEN
      CALL attach_injector_to_list(injector_z_max, injector)
    END IF

  END SUBROUTINE attach_injector



  ! Actually does the attaching of the injector to the correct list
  SUBROUTINE attach_injector_to_list(list, injector)

    TYPE(injector_block), POINTER :: list
    TYPE(injector_block), POINTER :: injector
    TYPE(injector_block), POINTER :: current

    NULLIFY(injector%next)

    IF (ASSOCIATED(list)) THEN
      current => list
      DO WHILE(ASSOCIATED(current%next))
        current => current%next
      END DO
      current%next => injector
    ELSE
      list => injector
    END IF

  END SUBROUTINE attach_injector_to_list



  SUBROUTINE deallocate_injectors

    CALL deallocate_injector_list(injector_x_min)
    CALL deallocate_injector_list(injector_x_max)
    CALL deallocate_injector_list(injector_y_min)
    CALL deallocate_injector_list(injector_y_max)
    CALL deallocate_injector_list(injector_z_min)
    CALL deallocate_injector_list(injector_z_max)

  END SUBROUTINE deallocate_injectors



  SUBROUTINE deallocate_injector_list(list)

    TYPE(injector_block), POINTER :: list
    TYPE(injector_block), POINTER :: current, next
    INTEGER :: i

    current => list
    DO WHILE(ASSOCIATED(current))
      next => current%next
      IF (current%density_function%init) &
          CALL deallocate_stack(current%density_function)
      IF (ASSOCIATED(current%depth)) DEALLOCATE(current%depth)
      IF (ASSOCIATED(current%depth_xy)) DEALLOCATE(current%depth_xy)
      DO i = 1, 3
        IF (current%temperature_function(i)%init) &
            CALL deallocate_stack(current%temperature_function(i))
        IF (current%drift_function(i)%init) &
            CALL deallocate_stack(current%drift_function(i))
      END DO
      IF (current%use_time_function) &
            CALL deallocate_stack(current%time_function)
      DEALLOCATE(current)
      current => next
    END DO

  END SUBROUTINE deallocate_injector_list



  SUBROUTINE run_injectors

    TYPE(injector_block), POINTER :: current

    !WRITE(*,*) 'Begin run injectors z',z_min_boundary, ASSOCIATED(injector_z_min)
    !WRITE(*,*) 'Begin run injectors x',x_min_boundary, ASSOCIATED(injector_x_min)
    IF (x_min_boundary) THEN
      current => injector_x_min
      DO WHILE(ASSOCIATED(current))
        CALL run_single_injector(current, c_bd_x_min)
        current => current%next
      END DO
    END IF

    IF (x_max_boundary) THEN
      current => injector_x_max
      DO WHILE(ASSOCIATED(current))
        CALL run_single_injector(current, c_bd_x_max)
        current => current%next
      END DO
    END IF

    IF (y_min_boundary) THEN
      current => injector_y_min
      DO WHILE(ASSOCIATED(current))
        CALL run_single_injector(current, c_bd_y_min)
        current => current%next
      END DO
    END IF

    IF (y_max_boundary) THEN
      current => injector_y_max
      DO WHILE(ASSOCIATED(current))
        CALL run_single_injector(current, c_bd_y_max)
        current => current%next
      END DO
    END IF

    IF (z_min_boundary) THEN
      current => injector_z_min
      DO WHILE(ASSOCIATED(current))
        CALL run_single_injector(current, c_bd_z_min)
        current => current%next
      END DO
    END IF

    IF (z_max_boundary) THEN
      current => injector_z_max
      DO WHILE(ASSOCIATED(current))
        CALL run_single_injector(current, c_bd_z_max)
        current => current%next
      END DO
    END IF

  END SUBROUTINE run_injectors



  SUBROUTINE run_single_injector(injector, direction)

    TYPE(injector_block), POINTER :: injector
    INTEGER, INTENT(IN) :: direction
    REAL(num) :: bdy_pos, cell_size
    TYPE(particle), POINTER :: new
    TYPE(particle_list) :: plist
    REAL(num) :: mass, typical_mc2, p_therm, p_inject_drift
    REAL(num) :: gamma_mass, v_inject, density, vol, p_drift, p_ratio
    REAL(num) :: npart_ideal, itemp, v_inject_s, density_correction, dir_mult
    REAL(num) :: v_inject_dt
#ifndef PER_SPECIES_WEIGHT
    REAL(num) :: weight_fac
#endif
    REAL(num), DIMENSION(3) :: temperature, drift
    INTEGER :: parts_this_time, ipart, idir, dir_index, flux_dir, flux_dir_cell
    INTEGER :: ii
    INTEGER :: perp_dir_index, nperp
    REAL(num) :: perp_cell_size, cur_cell
    TYPE(parameter_pack) :: parameters
    REAL(num), PARAMETER :: sqrt2 = SQRT(2.0_num)
    REAL(num), PARAMETER :: sqrt2_inv = 1.0_num / sqrt2
    REAL(num), PARAMETER :: sqrt2pi_inv = 1.0_num / SQRT(2.0_num * pi)

    IF (time < injector%t_start .OR. time > injector%t_end) RETURN

    ! If you have a moving window that has started moving then unless you
    ! EXPLICITLY give a t_end value to the injector stop the injector
    IF (move_window .AND. window_started .AND. .NOT. injector%has_t_end) &
        RETURN

    IF (direction == c_bd_x_min) THEN
      bdy_pos = x_min
      parameters%pack_ix = 0
      dir_mult = 1.0_num
      ! x-direction
      dir_index = 1
      cell_size = dx
      perp_dir_index = 2
      perp_cell_size = dy
      nperp = ny
    ELSE IF (direction == c_bd_x_max) THEN
      bdy_pos = x_max
      parameters%pack_ix = nx
      dir_mult = -1.0_num
      ! x-direction
      dir_index = 1
      cell_size = dx
      perp_dir_index = 2
      perp_cell_size = dy
      nperp = ny
    ELSE IF (direction == c_bd_y_min) THEN
      bdy_pos = y_min
      parameters%pack_iy = 0
      dir_mult = 1.0_num
      ! y-direction
      dir_index = 2
      cell_size = dy
      perp_dir_index = 1
      perp_cell_size = dx
      nperp = nx
    ELSE IF (direction == c_bd_y_max) THEN
      bdy_pos = y_max
      parameters%pack_iy = ny
      dir_mult = -1.0_num
      ! y-direction
      dir_index = 2
      cell_size = dy
      perp_dir_index = 1
      perp_cell_size = dx
      nperp = nx
      ! XY for injector in z direction
    ELSE IF (direction == c_bd_z_min) THEN
        CALL injector_in_z_direction(injector=injector)
       RETURN
    END IF

    IF (injector%use_flux_injector) THEN
      flux_dir = dir_index
    ELSE
      flux_dir = -1
    END IF

    vol = dx * dy
    bdy_pos = bdy_pos - 0.5_num * dir_mult * cell_size * png

    mass = species_list(injector%species)%mass
    typical_mc2 = (mass * c)**2
    cur_cell = 0.0_num
#ifndef PER_SPECIES_WEIGHT
    weight_fac = vol / injector%npart_per_cell
#endif

    CALL create_empty_partlist(plist)

    DO ii = 1, nperp
      IF (perp_dir_index == 1) THEN
        cur_cell = x(ii)
        parameters%pack_ix = ii
      ELSE
        cur_cell = y(ii)
        parameters%pack_iy = ii
      END IF

      parameters%use_grid_position = .TRUE.

      CALL populate_injector_properties(injector, parameters, density=density)

      IF (density < injector%density_min) CYCLE

      CALL populate_injector_properties(injector, parameters, &
          temperature=temperature, drift=drift)

      ! Assume agressive maximum thermal momentum, all components
      ! like hottest component
      p_therm = SQRT(mass * kb * MAXVAL(temperature))
      p_inject_drift = drift(dir_index)
      flux_dir_cell = flux_dir

      IF (flux_dir_cell /= -1) THEN
        ! Drift adjusted so that +ve is 'inwards' through boundary
        p_drift = p_inject_drift * dir_mult

        ! Average momentum of inflowing part
        ! For large inwards drift, is asymptotic to drift
        ! Otherwise it is a complicated expression
        ! Inwards drift - lhs terms are same sign -> +ve
        IF (p_drift > flow_limit_val * p_therm) THEN
          ! For sufficiently large drifts, net inflow -> p_drift
          gamma_mass = SQRT(p_inject_drift**2 + typical_mc2) / c
          v_inject_s = p_inject_drift / gamma_mass
          density_correction = 1.0_num
          ! Large drift flux Maxwellian can be approximated by a
          ! non-flux Maxwellian
          flux_dir_cell = -1
        ELSE IF (p_drift < -flow_limit_val * p_therm) THEN
          ! Net is outflow - inflow velocity is zero
          CYCLE
        ELSE IF (ABS(p_therm) < c_tiny) THEN
          CYCLE
        ELSE IF (ABS(p_drift) < p_therm * 1.0e-9_num) THEN
          v_inject_s = 2.0_num * sqrt2pi_inv * p_therm &
              + (1.0_num - 2.0_num * sqrt2 / pi) * p_drift
          gamma_mass = SQRT(v_inject_s**2 + typical_mc2) / c
          v_inject_s = v_inject_s / gamma_mass
          density_correction = 0.5_num
        ELSE
          p_ratio = sqrt2_inv * p_drift / p_therm

          ! Fraction of the drifting Maxwellian distribution inflowing
          density_correction = 0.5_num * (1.0_num + erf_func(p_ratio))
          IF (density_correction < c_tiny) CYCLE

          ! Below is actually MOMENTUM, will correct on next line
          v_inject_s = dir_mult * (p_drift &
              + sqrt2pi_inv * p_therm * EXP(-p_ratio**2) / density_correction)

          gamma_mass = SQRT(v_inject_s**2 + typical_mc2) / c
          v_inject_s = v_inject_s / gamma_mass
        END IF
      ELSE
        ! User asked for Maxwellian only - no correction to apply
        gamma_mass = SQRT(p_inject_drift**2 + typical_mc2) / c
        v_inject_s = p_inject_drift / gamma_mass
        density_correction = 1.0_num
      END IF

      v_inject = ABS(v_inject_s)
      v_inject_dt = dt * v_inject_s

      npart_ideal = injector%npart_per_cell * v_inject * density_correction &
          * dt / cell_size
      itemp = random_box_muller(0.5_num * SQRT(npart_ideal &
          * (1.0_num - npart_ideal / injector%npart_per_cell))) + npart_ideal
      injector%depth(ii) = injector%depth(ii) - itemp

      IF (injector%depth(ii) >= 0.0_num) CYCLE

      parts_this_time = FLOOR(ABS(injector%depth(ii) - 1.0_num))
      injector%depth(ii) = injector%depth(ii) + REAL(parts_this_time, num)

      DO ipart = 1, parts_this_time
        CALL create_particle(new)
        new%part_pos(perp_dir_index) = &
            (random() - 0.5_num) * perp_cell_size + cur_cell

        new%part_pos(dir_index) = bdy_pos - random() * v_inject_dt
        parameters%pack_pos = new%part_pos
        parameters%use_grid_position = .FALSE.

#ifdef PER_SPECIES_WEIGHT
        CALL populate_injector_properties(injector, parameters, &
            temperature=temperature, drift=drift)
#else
        CALL populate_injector_properties(injector, parameters, density, &
            temperature, drift)
#endif

        DO idir = 1, 3
          IF (idir == flux_dir_cell) THEN
            ! Drift is signed - dir mult is the direciton we want to get
            new%part_p(idir) = flux_momentum_from_temperature(&
                mass, temperature(idir), drift(idir), dir_mult)
          ELSE
            new%part_p(idir) = momentum_from_temperature(mass, &
                temperature(idir), drift(idir))
          END IF
        END DO
#ifdef PER_PARTICLE_CHARGE_MASS
        new%charge = species_list(injector%species)%charge
        new%mass = mass
#endif
#ifndef PER_SPECIES_WEIGHT
        density = MIN(density, injector%density_max)
        new%weight = weight_fac * density
#endif
        CALL add_particle_to_partlist(plist, new)
      END DO
      END DO ! Do ipart

    CALL append_partlist(species_list(injector%species)%attached_list, plist)

  END SUBROUTINE run_single_injector

  SUBROUTINE injector_in_z_direction(injector)
      TYPE(parameter_pack) :: parameters
      TYPE(injector_block), POINTER :: injector
      TYPE(particle_list) :: plist
      INTEGER :: ix,iy,ipart,idir
      INTEGER :: npart_per_cell,parts_this_time
      REAL(num) :: density_max, density_min,mass
      REAL(num) :: density_local, weight_fac, vol
      LOGICAL :: load
      TYPE(particle), POINTER :: new
      REAL(num), DIMENSION(3) :: temperature, drift
      REAL(num) :: v_inject, density_correction,npart_ideal,itemp,cell_size,n,g1,g2,P,en,theta,phi
      REAL(num) :: inject_theta,inject_phi,px,py,pz,gam
      REAL(num) :: energy_pho,en_kbT,p_mag,p_tau

      weight_fac = 1.0_num
      mass = species_list(injector%species)%mass
      ! npart_per_cell 
      npart_per_cell = injector%npart_per_cell
      density_max = injector%density_max
      density_min = injector%density_min
      density_local = 0.0_num

      CALL create_empty_partlist(plist)
      DO iy = 1,ny
        parameters%pack_iy = iy
        DO ix = 1,nx
            load = .FALSE.
        parameters%pack_ix = ix
      ! Calculate density,temeperature and drift of inject
        CALL populate_injector_properties(injector, parameters, density=density_local)
        !IF (density < injector%density_min) CYCLE
        CALL populate_injector_properties(injector, parameters, &
           temperature=temperature, drift=drift)

        IF (density_local > density_max) density_local = density_max
        IF (density_local >= density_min) THEN
            load = .TRUE.
        END IF
        IF (.NOT. load) CYCLE
        !Calculate parts_this_time, some assumption see REAMME.txt
        cell_size = 1.0_num ! 1m in z direction
        v_inject = injector_time_profile(injector,parameters) * injector%v_inject_z ! must > 0.0

        density_correction = 1.0_num
        npart_ideal = injector%npart_per_cell * v_inject * density_correction * dt / cell_size
        !  * dt / cell_size
        !delete dt /cell_size

        itemp = random_box_muller(0.5_num * SQRT(npart_ideal &
          * (1.0_num - npart_ideal / injector%npart_per_cell))) + npart_ideal
        injector%depth_xy(ix,iy) = injector%depth_xy(ix,iy) - itemp


      IF (injector%depth_xy(ix,iy) >= 0.0_num) CYCLE

      parts_this_time = MIN(100,FLOOR(ABS(injector%depth_xy(ix,iy) - 1.0_num)))
      !IF (debug_mod) THEN
      !  WRITE(*,*) 'part_this_time,depth_xy',parts_this_time,injector%depth_xy(ix,iy),npart_ideal,itemp
      !END IF

      injector%depth_xy(ix,iy) = injector%depth_xy(ix,iy) + REAL(parts_this_time, num)
!
!        parts_this_time = npart_per_cell
        ipart = 0
        vol = dx * dy
#ifndef PER_SPECIES_WEIGHT
        weight_fac = vol / injector%npart_per_cell
#endif

        DO WHILE(ipart < parts_this_time)
            CALL create_particle(new)
            new%part_pos(1) = x(ix) + (random() - 0.5_num) * dx
            new%part_pos(2) = y(iy) + (random() - 0.5_num) * dy

#ifdef PER_PARTICLE_CHARGE_MASS
            new%charge = species_list(injector%species)%charge
            new%mass = mass
#endif
#ifndef PER_SPECIES_WEIGHT
            new%weight = weight_fac * density_local
#endif
        ! momentum
        IF (species_list(injector%species)%species_type == c_species_id_photon) THEN
            IF (injector%monoenergy > 0.0_num) THEN
                !theta is (p,x), phi 是p在yz平面投影与y夹角 (px,y)
                inject_theta = injector%angle_theta
                inject_phi = injector%angle_phi
                new%part_p = injector%monoenergy/c*(/cos(inject_theta),&
                                                    sin(inject_theta)*cos(inject_phi),&
                                                    sin(inject_theta)*sin(inject_phi)/)
                new%particle_energy = injector%monoenergy
                p_tau = random()
                new%optical_depth = -LOG(1.0_num - p_tau)
            ELSE IF (injector%wein_temperature > 0.0_num) THEN
                en = injector%wein_temperature*q0
                P = random()
                en_kbT = find_value_from_table_1d_inject(P, 1000, wein_table, log_wein_xi)
                energy_pho = en_kbT*en;
                p_mag = energy_pho/c

                inject_theta = injector%angle_theta
                inject_phi = injector%angle_phi
                IF (inject_theta < -2.0*pi) THEN
                  !for uniform solid angle
                  theta = 2*pi*random() ! dOmega is uniform 
                  phi = ACOS(2*random() - 1)
                  px = p_mag*SIN(theta)*SIN(phi)
                  py = p_mag*COS(theta)*SIN(phi)
                  pz = p_mag*COS(phi)
                  new%part_p = (/px,py,pz/)
                ELSE
                  !one direction
                  new%part_p = p_mag* (/cos(inject_theta),sin(inject_theta)*cos(inject_phi),sin(inject_theta)*sin(inject_phi)/)
                END IF
                new%particle_energy = energy_pho
                !Optical depth
                p_tau = random()
                new%optical_depth = -LOG(1.0_num - p_tau)
            END IF !monoenergy
        ELSE ! not photon
!           IF (species_list(injector%species)%species_type == c_species_id_electron .OR. &
!               species_list(injector%species)%species_type == c_species_id_positron) THEN
            IF (injector%power_law > 0.0) THEN
                g1 = injector%g1/1e3/q0
                g2 = injector%g2/1e3/q0
                n = injector%power_law
                P = random()
!               en = 10**(sqrt(P*(log10_g2**2-log10_g1**2) + log10_g1**2))*1e3*q0
                en = g1/(((g1/g2)**(n-1)-1)*P + 1)**(1/(n-1))*1e3*q0
                gam = en/mass/c**2 + 1
! angle
                p_mag = sqrt(gam**2 - 1)*mass*c
                inject_theta = injector%angle_theta
                inject_phi = injector%angle_phi
                IF (inject_theta < -2.0*pi) THEN
                    !for uniform solid angle
                    theta = 2*pi*random() ! dOmega is uniform 
                    phi = ACOS(2*random() - 1)
                    px = p_mag*SIN(theta)*SIN(phi)
                    py = p_mag*COS(theta)*SIN(phi)
                    pz = p_mag*COS(phi)
                    new%part_p = (/px,py,pz/)
                ELSE
                  !one direction
                  new%part_p = p_mag* (/cos(inject_theta),sin(inject_theta)*cos(inject_phi),sin(inject_theta)*sin(inject_phi)/)
                END IF
                new%particle_energy = en
            ELSE ! not power_law
              DO idir = 1, 3
                new%part_p(idir) = momentum_from_temperature(mass, &
                    temperature(idir), drift(idir))
                p_mag = SQRT(DOT_PRODUCT(new%part_p/m0/c,new%part_p/m0/c))
                new%particle_energy = (SQRT(p_mag**2 + 1) - 1)*m0*c**2
              END DO
            END IF ! for injector%power_law
        END IF!photon

           CALL add_particle_to_partlist(plist, new)
           ipart = ipart + 1
          ! One particle sucessfully placed
           !npart_left = npart_left - 1
      END DO ! Do While for ipart
    END DO ! ix
    END DO ! iy

    CALL append_partlist(species_list(injector%species)%attached_list, plist)

   END SUBROUTINE injector_in_z_direction

  FUNCTION injector_time_profile(injector,parameters)

    TYPE(injector_block), POINTER :: injector
    REAL(num) :: injector_time_profile
    INTEGER :: err
    TYPE(parameter_pack),INTENT(IN) :: parameters

    err = 0
    IF (injector%use_time_function) THEN
      injector_time_profile = evaluate_with_parameters(injector%time_function, &
          parameters, err)
      RETURN
    END IF

    injector_time_profile = 1.0_num

  END FUNCTION injector_time_profile


  SUBROUTINE populate_injector_properties(injector, parameters, density, &
      temperature, drift)

    TYPE(injector_block), POINTER :: injector
    TYPE(parameter_pack), INTENT(IN) :: parameters
    REAL(num), INTENT(OUT), OPTIONAL :: density
    REAL(num), DIMENSION(3), INTENT(OUT), OPTIONAL :: temperature, drift
    INTEGER :: errcode, i

    errcode = 0
    IF (PRESENT(density)) THEN
      density = 0.0_num
      IF (injector%density_function%init) THEN
        density = MAX(evaluate_with_parameters(injector%density_function, &
            parameters, errcode), 0.0_num)
      END IF
    END IF

    ! Stack can only be time varying if valid. Change if this isn't true
    IF (PRESENT(temperature)) THEN
      temperature(:) = 0.0_num
      DO i = 1, 3
        IF (injector%temperature_function(i)%init) THEN
          temperature(i) = &
              MAX(evaluate_with_parameters(injector%temperature_function(i), &
                  parameters, errcode), 0.0_num)
        END IF
      END DO
    END IF

    IF (PRESENT(drift)) THEN
      drift(:) = 0.0_num
      DO i = 1, 3
        IF (injector%drift_function(i)%init) THEN
          drift(i) = &
              evaluate_with_parameters(injector%drift_function(i), &
                                       parameters, errcode)
        END IF
      END DO
    END IF

    IF (errcode /= c_err_none) CALL abort_code(errcode)

  END SUBROUTINE populate_injector_properties



  SUBROUTINE finish_injector_setup

    TYPE(injector_block), POINTER :: current

    IF (x_min_boundary) THEN
      current => injector_x_min
      DO WHILE(ASSOCIATED(current))
        CALL finish_single_injector_setup(current, c_bd_x_min)
        current => current%next
      END DO
    END IF

    IF (x_max_boundary) THEN
      current => injector_x_max
      DO WHILE(ASSOCIATED(current))
        CALL finish_single_injector_setup(current, c_bd_x_max)
        current => current%next
      END DO
    END IF

    IF (y_min_boundary) THEN
      current => injector_y_min
      DO WHILE(ASSOCIATED(current))
        CALL finish_single_injector_setup(current, c_bd_y_min)
        current => current%next
      END DO
    END IF

    IF (y_max_boundary) THEN
      current => injector_y_max
      DO WHILE(ASSOCIATED(current))
        CALL finish_single_injector_setup(current, c_bd_y_max)
        current => current%next
      END DO
    END IF

    IF (z_min_boundary) THEN
      current => injector_z_min
      DO WHILE(ASSOCIATED(current))
        CALL finish_single_injector_setup(current, c_bd_z_min)
        current => current%next
      END DO
    END IF


  END SUBROUTINE finish_injector_setup



  SUBROUTINE finish_single_injector_setup(injector, boundary)

    TYPE(injector_block), POINTER :: injector
    INTEGER, INTENT(IN) :: boundary
    TYPE(particle_species), POINTER :: species
    INTEGER :: i

    species => species_list(injector%species)
    IF (injector%npart_per_cell < 0.0_num) THEN
      injector%npart_per_cell = species%npart_per_cell
    END IF

    IF (.NOT.injector%density_function%init) THEN
      CALL copy_stack(species%density_function, injector%density_function)
    END IF

    DO i = 1, 3
      IF (.NOT.injector%drift_function(i)%init) THEN
        CALL copy_stack(species%drift_function(i), injector%drift_function(i))
      END IF
      IF (.NOT.injector%temperature_function(i)%init) THEN
        CALL copy_stack(species%temperature_function(i), &
            injector%temperature_function(i))
      END IF
    END DO

  END SUBROUTINE finish_single_injector_setup



  SUBROUTINE create_boundary_injector(ispecies, bnd)

    INTEGER, INTENT(IN) :: ispecies, bnd
    TYPE(injector_block), POINTER :: working_injector

    IF (bnd <= 2*c_ndims) THEN 
    species_list(ispecies)%bc_particle(bnd) = c_bc_open
    END IF
    use_injectors = .TRUE.

    ALLOCATE(working_injector)

    CALL init_injector(bnd, working_injector)
    working_injector%species = ispecies

    CALL attach_injector(working_injector)

    WRITE(*,*) 'create_boundary:', ispecies,bnd, ASSOCIATED(working_injector)

  END SUBROUTINE create_boundary_injector



  SUBROUTINE setup_injector_depths(inj_init, depths, inj_count)

    TYPE(injector_block), POINTER :: inj_init
    REAL(num), DIMENSION(:,:), INTENT(IN) :: depths
    INTEGER, INTENT(OUT) :: inj_count
    TYPE(injector_block), POINTER :: inj
    INTEGER :: iinj

    iinj = 1
    inj => inj_init

    DO WHILE(ASSOCIATED(inj))
      ! Exclude ghost cells
      IF (inj%boundary == c_bd_x_min .OR. inj%boundary == c_bd_x_max) THEN
        inj%depth(1:ny) = depths(:,iinj)
      ELSE IF (inj%boundary == c_bd_y_min .OR. inj%boundary == c_bd_y_max) THEN
        inj%depth(1:nx) = depths(:,iinj)
      ELSE
        inj%depth_xy(:,:) = depths(:,:) ! Some problems for restart, to revise one day
      END IF
      iinj = iinj + 1
      inj => inj%next
    END DO

    inj_count = iinj - 1

  END SUBROUTINE setup_injector_depths

  FUNCTION find_value_from_table_1d_inject(x_in, nx, x, values)

    REAL(num) :: find_value_from_table_1d_inject
    REAL(num), INTENT(IN) :: x_in
    INTEGER, INTENT(IN) :: nx
    REAL(num), INTENT(IN) :: x(nx), values(nx)
    REAL(num) :: fx, x_value, value_interp, xdif1, xdif2, xdifm
    INTEGER :: i1, i2, im
    LOGICAL, SAVE :: warning = .TRUE.

    !No log10
    x_value = x_in

    i1 = 1
    i2 = nx
    xdif1 = x(i1) - x_value
    xdif2 = x(i2) - x_value
    IF (xdif1 * xdif2 < 0) THEN
      ! Use bisection to find the nearest cell
      DO
        im = (i1 + i2) / 2
        xdifm = x(im) - x_value
        IF (xdif1 * xdifm < 0) THEN
          i2 = im
        ELSE
          i1 = im
          xdif1 = xdifm
        END IF
        IF (i2 - i1 == 1) EXIT
      END DO
      ! Interpolate in x
      fx = (x_value - x(i1)) / (x(i2) - x(i1))
    ELSE
      IF (warning .AND. rank == 0) THEN
        PRINT*,'*** WARNING ***'
        PRINT*,'Argument to "find_value_from_table_1d" outside the range ', &
            'of the table.'
        PRINT*,'Using truncated value. No more warnings will be issued.'
        warning = .FALSE.
      END IF
      IF (xdif1 >= 0) THEN
        fx = 0.0_num
      ELSE
        fx = 1.0_num
      END IF
    END IF

    value_interp = (1.0_num - fx) * values(i1) + fx * values(i2)

    find_value_from_table_1d_inject = 10.0_num**value_interp
  END FUNCTION find_value_from_table_1d_inject

END MODULE injectors
