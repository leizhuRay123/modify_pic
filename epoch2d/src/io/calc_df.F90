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

MODULE calc_df

  USE boundary
#ifdef PHOTONS
  USE photons
#endif


  IMPLICIT NONE

CONTAINS

  SUBROUTINE calc_boundary(data_array, species)

    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN), OPTIONAL :: species

    CALL processor_summation_bcs(data_array, ng, species=species)

  END SUBROUTINE calc_boundary

  SUBROUTINE calc_mass_density(data_array, current_species, direction)

    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species
    INTEGER, INTENT(IN), OPTIONAL :: direction
    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_m
    ! The data to be weighted onto the grid
    REAL(num) :: wdata
    REAL(num) :: fac, idx
    INTEGER :: ispecies, ix, iy, spec_start, spec_end
    TYPE(particle), POINTER :: current
    LOGICAL :: spec_sum
#include "particle_head.inc"

    data_array = 0.0_num
    part_m = 0.0_num
    fac = 0.0_num

    idx = 1.0_num / dx / dy

    spec_start = current_species
    spec_end = current_species
    spec_sum = .FALSE.

    IF (current_species <= 0) THEN
      spec_start = 1
      spec_end = n_species
      spec_sum = .TRUE.
    END IF

    DO ispecies = spec_start, spec_end
      IF (io_list(ispecies)%species_type == c_species_id_photon) CYCLE
#ifndef NO_TRACER_PARTICLES
      IF (spec_sum .AND. io_list(ispecies)%zero_current) CYCLE
#endif
      current => io_list(ispecies)%attached_list%head
      part_m  = io_list(ispecies)%mass
      fac = io_list(ispecies)%weight
      wdata = part_m * fac

      DO WHILE (ASSOCIATED(current))
        ! Copy the particle properties out for speed
#ifdef PER_PARTICLE_CHARGE_MASS
        part_m  = current%mass
#ifndef PER_SPECIES_WEIGHT
        fac = current%weight
#endif
        wdata = part_m * fac
#else
#ifndef PER_SPECIES_WEIGHT
        fac = current%weight
        wdata = part_m * fac
#endif
#endif

#include "particle_to_grid.inc"

        DO iy = sf_min, sf_max
        DO ix = sf_min, sf_max
          data_array(cell_x+ix, cell_y+iy) = &
              data_array(cell_x+ix, cell_y+iy) + gx(ix) * gy(iy) * wdata
        END DO
        END DO

        current => current%next
      END DO
      CALL calc_boundary(data_array, ispecies)
    END DO

    CALL calc_boundary(data_array)

    data_array = data_array * idx
    DO ix = 1, 2*c_ndims
      CALL field_zero_gradient(data_array, c_stagger_centre, ix)
    END DO

  END SUBROUTINE calc_mass_density

  SUBROUTINE calc_ekbar(data_array, current_species, direction)

    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species
    INTEGER, INTENT(IN), OPTIONAL :: direction
    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_ux, part_uy, part_uz, part_mc, part_u2
    ! The weight of a particle
    REAL(num) :: part_w
    ! The data to be weighted onto the grid
    REAL(num) :: wdata
    REAL(num) :: fac, gamma_rel, gamma_rel_m1
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: wt
    INTEGER :: ispecies, ix, iy, spec_start, spec_end
    TYPE(particle), POINTER :: current
    LOGICAL :: spec_sum
#include "particle_head.inc"

    ALLOCATE(wt(1-ng:nx+ng,1-ng:ny+ng))
    data_array = 0.0_num
    wt = 0.0_num
    part_mc  = 1.0_num
    part_w = 1.0_num

    spec_start = current_species
    spec_end = current_species
    spec_sum = .FALSE.

    IF (current_species <= 0) THEN
      spec_start = 1
      spec_end = n_species
      spec_sum = .TRUE.
    END IF

    DO ispecies = spec_start, spec_end
      IF (spec_sum) CYCLE
#ifndef NO_TRACER_PARTICLES
      IF (spec_sum .AND. io_list(ispecies)%zero_current) CYCLE
#endif
      current => io_list(ispecies)%attached_list%head
      part_mc = c * io_list(ispecies)%mass
      part_w = io_list(ispecies)%weight
      fac = part_mc * part_w * c

      DO WHILE (ASSOCIATED(current))
        ! Copy the particle properties out for speed
#ifdef PER_PARTICLE_CHARGE_MASS
        part_mc = c * current%mass
#ifndef PER_SPECIES_WEIGHT
        part_w = current%weight
#endif
        fac = part_mc * part_w * c
#else
#ifndef PER_SPECIES_WEIGHT
        part_w = current%weight
        fac = part_mc * part_w * c
#endif
#endif

        IF (io_list(ispecies)%species_type /= c_species_id_photon) THEN
          part_ux = current%part_p(1) / part_mc
          part_uy = current%part_p(2) / part_mc
          part_uz = current%part_p(3) / part_mc

          part_u2 = part_ux**2 + part_uy**2 + part_uz**2
          gamma_rel = SQRT(part_u2 + 1.0_num)
          gamma_rel_m1 = part_u2 / (gamma_rel + 1.0_num)

          wdata = gamma_rel_m1 * fac
        ELSE
#if defined(PHOTONS) || defined(BREMSSTRAHLUNG)
          wdata = SQRT(SUM(current%part_p**2))*c*part_w
          !wdata = current%particle_energy * part_w

          !IF (current%particle_energy > 10000*1e6*q0) THEN
          !    WRITE(*,*) 'ispecies,particle_energy',ispecies,current%particle_energy/1e6/q0
          !    WRITE(*,*) 'id',current%id
          !END IF

#else
          wdata = 0.0_num
#endif
        END IF

#include "particle_to_grid.inc"

        DO iy = sf_min, sf_max
        DO ix = sf_min, sf_max
          data_array(cell_x+ix, cell_y+iy) = &
              data_array(cell_x+ix, cell_y+iy) + gx(ix) * gy(iy) * wdata
          wt(cell_x+ix, cell_y+iy) = &
              wt(cell_x+ix, cell_y+iy) + gx(ix) * gy(iy) * part_w
        END DO
        END DO

        current => current%next
      END DO
      CALL calc_boundary(data_array, ispecies)
      CALL calc_boundary(wt, ispecies)
    END DO

    CALL calc_boundary(data_array)
    CALL calc_boundary(wt)

    data_array = data_array / MAX(wt, c_tiny)
    DO ix = 1, 2*c_ndims
      CALL field_zero_gradient(data_array, c_stagger_centre, ix)
    END DO

    DEALLOCATE(wt)

  END SUBROUTINE calc_ekbar



  SUBROUTINE calc_ekflux(data_array, current_species, direction)

    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species
    INTEGER, INTENT(IN), OPTIONAL :: direction
    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_ux, part_uy, part_uz, part_mc, part_u2
    ! The weight of a particle
    REAL(num) :: part_w
    ! The data to be weighted onto the grid
    REAL(num) :: wdata
    REAL(num) :: fac, gamma_rel, gamma_rel_m1, part_flux, xfac, yfac, zfac
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: wt
    INTEGER :: ispecies, ix, iy, spec_start, spec_end
    TYPE(particle), POINTER :: current
    LOGICAL :: spec_sum
#include "particle_head.inc"

    ALLOCATE(wt(1-ng:nx+ng,1-ng:ny+ng))
    data_array = 0.0_num
    wt = 0.0_num
    part_mc  = 1.0_num
    part_w = 1.0_num

    xfac = c * dy
    yfac = c * dx
    zfac = c * dx * dy

    spec_start = current_species
    spec_end = current_species
    spec_sum = .FALSE.

    IF (current_species <= 0) THEN
      spec_start = 1
      spec_end = n_species
      spec_sum = .TRUE.
    END IF

    DO ispecies = spec_start, spec_end
      IF (spec_sum .AND. &
          io_list(ispecies)%species_type == c_species_id_photon) CYCLE
#ifndef NO_TRACER_PARTICLES
      IF (spec_sum .AND. io_list(ispecies)%zero_current) CYCLE
#endif
      current => io_list(ispecies)%attached_list%head
      part_mc = c * io_list(ispecies)%mass
      part_w = io_list(ispecies)%weight
      fac = part_mc * part_w * c

      DO WHILE (ASSOCIATED(current))
        ! Copy the particle properties out for speed
#ifdef PER_PARTICLE_CHARGE_MASS
        part_mc = c * current%mass
#ifndef PER_SPECIES_WEIGHT
        part_w = current%weight
#endif
        fac = part_mc * part_w * c
#else
#ifndef PER_SPECIES_WEIGHT
        part_w = current%weight
        fac = part_mc * part_w * c
#endif
#endif

        IF (io_list(ispecies)%species_type /= c_species_id_photon) THEN
          part_ux = current%part_p(1) / part_mc
          part_uy = current%part_p(2) / part_mc
          part_uz = current%part_p(3) / part_mc

          part_u2 = part_ux**2 + part_uy**2 + part_uz**2
          gamma_rel = SQRT(part_u2 + 1.0_num)
          gamma_rel_m1 = part_u2 / (gamma_rel + 1.0_num)

          wdata = gamma_rel_m1 * fac
        ELSE
#if defined(PHOTONS) || defined(BREMSSTRAHLUNG)
          fac = c / current%particle_energy
          part_ux = current%part_p(1) * fac
          part_uy = current%part_p(2) * fac
          part_uz = current%part_p(3) * fac
          gamma_rel = 1.0_num
          wdata = current%particle_energy * part_w
#else
          gamma_rel = 1.0_num
          wdata = 0.0_num
#endif
        END IF

        SELECT CASE(direction)
        CASE(-c_dir_x)
          ! negative flux in x
          part_flux = xfac * part_ux / gamma_rel
          wdata = -wdata * MIN(part_flux, 0.0_num)
        CASE( c_dir_x)
          ! positive flux in x
          part_flux = xfac * part_ux / gamma_rel
          wdata =  wdata * MAX(part_flux, 0.0_num)
        CASE(-c_dir_y)
          ! negative flux in y
          part_flux = yfac * part_uy / gamma_rel
          wdata = -wdata * MIN(part_flux, 0.0_num)
        CASE( c_dir_y)
          ! positive flux in y
          part_flux = yfac * part_uy / gamma_rel
          wdata =  wdata * MAX(part_flux, 0.0_num)
        CASE(-c_dir_z)
          ! negative flux in z
          part_flux = zfac * part_uz / gamma_rel
          wdata = -wdata * MIN(part_flux, 0.0_num)
        CASE( c_dir_z)
          ! positive flux in z
          part_flux = zfac * part_uz / gamma_rel
          wdata =  wdata * MAX(part_flux, 0.0_num)
        END SELECT

#include "particle_to_grid.inc"

        DO iy = sf_min, sf_max
        DO ix = sf_min, sf_max
          data_array(cell_x+ix, cell_y+iy) = &
              data_array(cell_x+ix, cell_y+iy) + gx(ix) * gy(iy) * wdata
          wt(cell_x+ix, cell_y+iy) = &
              wt(cell_x+ix, cell_y+iy) + gx(ix) * gy(iy) * part_w
        END DO
        END DO

        current => current%next
      END DO
      CALL calc_boundary(data_array, ispecies)
      CALL calc_boundary(wt, ispecies)
    END DO

    CALL calc_boundary(data_array)
    CALL calc_boundary(wt)

    data_array = data_array / MAX(wt, c_tiny)
    DO ix = 1, 2*c_ndims
      CALL field_zero_gradient(data_array, c_stagger_centre, ix)
    END DO

    DEALLOCATE(wt)

  END SUBROUTINE calc_ekflux



  SUBROUTINE calc_poynt_flux(data_array, current_species, direction)

    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species
    INTEGER, INTENT(IN), OPTIONAL :: direction
    INTEGER :: ix, iy
    REAL(num) :: ex_cc, ey_cc, ez_cc, bx_cc, by_cc, bz_cc

    SELECT CASE(direction)
    CASE(c_dir_x)
      DO iy = 1, ny
      DO ix = 1, nx
        ey_cc = 0.5_num  * (ey(ix  , iy-1) + ey(ix, iy  ))
        ez_cc = ez(ix, iy)
        by_cc = 0.5_num  * (by(ix-1, iy  ) + by(ix, iy  ))
        bz_cc = 0.25_num * (bz(ix-1, iy-1) + bz(ix, iy-1) &
                         +  bz(ix-1, iy  ) + bz(ix, iy  ))
        data_array(ix,iy) = (ey_cc * bz_cc - ez_cc * by_cc) / mu0
      END DO
      END DO
    CASE(c_dir_y)
      DO iy = 1, ny
      DO ix = 1, nx
        ex_cc = 0.5_num  * (ex(ix-1, iy  ) + ex(ix, iy  ))
        ez_cc = ez(ix, iy)
        bx_cc = 0.5_num  * (bx(ix  , iy-1) + bx(ix, iy  ))
        bz_cc = 0.25_num * (bz(ix-1, iy-1) + bz(ix, iy-1) &
                         +  bz(ix-1, iy  ) + bz(ix, iy  ))
        data_array(ix,iy) = (ez_cc * bx_cc - ex_cc * bz_cc) / mu0
      END DO
      END DO
    CASE(c_dir_z)
      DO iy = 1, ny
      DO ix = 1, nx
        ex_cc = 0.5_num  * (ex(ix-1, iy  ) + ex(ix, iy))
        ey_cc = 0.5_num  * (ey(ix  , iy-1) + ey(ix, iy))
        bx_cc = 0.5_num  * (bx(ix  , iy-1) + bx(ix, iy))
        by_cc = 0.5_num  * (by(ix-1, iy  ) + by(ix, iy))
        data_array(ix,iy) = (ex_cc * by_cc - ey_cc * bx_cc) / mu0
      END DO
      END DO
    END SELECT

  END SUBROUTINE calc_poynt_flux



  SUBROUTINE calc_charge_density(data_array, current_species, direction)

    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species
    INTEGER, INTENT(IN), OPTIONAL :: direction
    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_q
    ! The data to be weighted onto the grid
    REAL(num) :: wdata
    REAL(num) :: fac, idx
    INTEGER :: ispecies, ix, iy, spec_start, spec_end
    TYPE(particle), POINTER :: current
    LOGICAL :: spec_sum
#include "particle_head.inc"

    data_array = 0.0_num
    part_q = 0.0_num
    fac = 0.0_num

    idx = 1.0_num / dx / dy

    spec_start = current_species
    spec_end = current_species
    spec_sum = .FALSE.

    IF (current_species <= 0) THEN
      spec_start = 1
      spec_end = n_species
      spec_sum = .TRUE.
    END IF

    DO ispecies = spec_start, spec_end
      IF (io_list(ispecies)%species_type == c_species_id_photon) CYCLE
#ifndef NO_TRACER_PARTICLES
      IF (spec_sum .AND. io_list(ispecies)%zero_current) CYCLE
#endif
      current => io_list(ispecies)%attached_list%head
      part_q  = io_list(ispecies)%charge
      fac = io_list(ispecies)%weight
      wdata = part_q * fac

      DO WHILE (ASSOCIATED(current))
        ! Copy the particle properties out for speed
#ifdef PER_PARTICLE_CHARGE_MASS
        part_q  = current%charge
#ifndef PER_SPECIES_WEIGHT
        fac = current%weight
#endif
        wdata = part_q * fac
#else
#ifndef PER_SPECIES_WEIGHT
        fac = current%weight
        wdata = part_q * fac
#endif
#endif

#include "particle_to_grid.inc"

        DO iy = sf_min, sf_max
        DO ix = sf_min, sf_max
          data_array(cell_x+ix, cell_y+iy) = &
              data_array(cell_x+ix, cell_y+iy) + gx(ix) * gy(iy) * wdata
        END DO
        END DO

        current => current%next
      END DO
      CALL calc_boundary(data_array, ispecies)
    END DO

    CALL calc_boundary(data_array)

    data_array = data_array * idx
    DO ix = 1, 2*c_ndims
      CALL field_zero_gradient(data_array, c_stagger_centre, ix)
    END DO

  END SUBROUTINE calc_charge_density



  SUBROUTINE calc_number_density(data_array, current_species, direction)

    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species
    INTEGER, INTENT(IN), OPTIONAL :: direction
    ! The data to be weighted onto the grid
    REAL(num) :: wdata
    REAL(num) :: idx, vol
    INTEGER :: ispecies, ix, iy, spec_start, spec_end
    TYPE(particle), POINTER :: current
    LOGICAL :: spec_sum
#include "particle_head.inc"

    data_array = 0.0_num

    vol = dx * dy
    idx = 1.0_num / vol

    spec_start = current_species
    spec_end = current_species
    spec_sum = .FALSE.

    IF (current_species <= 0) THEN
      spec_start = 1
      spec_end = n_species
      spec_sum = .TRUE.
    END IF

    DO ispecies = spec_start, spec_end
      IF (spec_sum &
          .AND. io_list(ispecies)%species_type == c_species_id_photon) CYCLE
#ifndef NO_TRACER_PARTICLES
      IF (spec_sum .AND. io_list(ispecies)%zero_current) CYCLE
#endif
      IF (io_list(ispecies)%background_species) THEN
        data_array(1:nx, 1:ny) = data_array(1:nx, 1:ny) &
            + io_list(ispecies)%background_density(1:nx, 1:ny) * vol
      ELSE
        current => io_list(ispecies)%attached_list%head
        wdata = io_list(ispecies)%weight

        DO WHILE (ASSOCIATED(current))
#ifndef PER_SPECIES_WEIGHT
          wdata = current%weight
#endif

#include "particle_to_grid.inc"

          DO iy = sf_min, sf_max
          DO ix = sf_min, sf_max
            data_array(cell_x+ix, cell_y+iy) = &
                data_array(cell_x+ix, cell_y+iy) + gx(ix) * gy(iy) * wdata
          END DO
          END DO

          current => current%next
        END DO
      END IF
      CALL calc_boundary(data_array, ispecies)
    END DO

    CALL calc_boundary(data_array)

    data_array = data_array * idx
    DO ix = 1, 2*c_ndims
      CALL field_zero_gradient(data_array, c_stagger_centre, ix)
    END DO

  END SUBROUTINE calc_number_density



  SUBROUTINE calc_ppc(data_array, current_species, direction)

    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species
    INTEGER, INTENT(IN), OPTIONAL :: direction
    INTEGER :: ispecies, spec_start, spec_end
    TYPE(particle), POINTER :: current
    LOGICAL :: spec_sum
    REAL(num) :: cell_x_r, cell_y_r
    INTEGER :: cell_x, cell_y

    data_array = 0.0_num

    spec_start = current_species
    spec_end = current_species
    spec_sum = .FALSE.

    IF (current_species <= 0) THEN
      spec_start = 1
      spec_end = n_species
      spec_sum = .TRUE.
    END IF

    DO ispecies = spec_start, spec_end
#ifndef NO_TRACER_PARTICLES
      IF (spec_sum .AND. io_list(ispecies)%zero_current) CYCLE
#endif
      current => io_list(ispecies)%attached_list%head

      DO WHILE (ASSOCIATED(current))
#ifdef PARTICLE_SHAPE_TOPHAT
        cell_x_r = (current%part_pos(1) - x_grid_min_local) / dx
        cell_y_r = (current%part_pos(2) - y_grid_min_local) / dy
#else
        cell_x_r = (current%part_pos(1) - x_grid_min_local) / dx + 0.5_num
        cell_y_r = (current%part_pos(2) - y_grid_min_local) / dy + 0.5_num
#endif
        cell_x = FLOOR(cell_x_r) + 1
        cell_y = FLOOR(cell_y_r) + 1

        data_array(cell_x,cell_y) = data_array(cell_x,cell_y) + 1.0_num

        current => current%next
      END DO
    END DO

  END SUBROUTINE calc_ppc



  SUBROUTINE calc_average_weight(data_array, current_species, direction)

    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species
    INTEGER, INTENT(IN), OPTIONAL :: direction
    ! The data to be weighted onto the grid
    REAL(num) :: wdata
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: part_count
    INTEGER :: ispecies, spec_start, spec_end
    TYPE(particle), POINTER :: current
    LOGICAL :: spec_sum
    REAL(num) :: cell_x_r, cell_y_r
    INTEGER :: cell_x, cell_y

    data_array = 0.0_num

    spec_start = current_species
    spec_end = current_species
    spec_sum = .FALSE.

    IF (current_species <= 0) THEN
      spec_start = 1
      spec_end = n_species
      spec_sum = .TRUE.
    END IF

    ALLOCATE(part_count(1-ng:nx+ng,1-ng:ny+ng))
    part_count = 0.0_num

    DO ispecies = spec_start, spec_end
#ifndef NO_TRACER_PARTICLES
      IF (spec_sum .AND. io_list(ispecies)%zero_current) CYCLE
#endif
      current => io_list(ispecies)%attached_list%head
      wdata = io_list(ispecies)%weight

      DO WHILE (ASSOCIATED(current))
#ifndef PER_SPECIES_WEIGHT
        wdata = current%weight
#endif

#ifdef PARTICLE_SHAPE_TOPHAT
        cell_x_r = (current%part_pos(1) - x_grid_min_local) / dx
        cell_y_r = (current%part_pos(2) - y_grid_min_local) / dy
#else
        cell_x_r = (current%part_pos(1) - x_grid_min_local) / dx + 0.5_num
        cell_y_r = (current%part_pos(2) - y_grid_min_local) / dy + 0.5_num
#endif
        cell_x = FLOOR(cell_x_r) + 1
        cell_y = FLOOR(cell_y_r) + 1

        data_array(cell_x,cell_y) = data_array(cell_x,cell_y) + wdata
        part_count(cell_x,cell_y) = part_count(cell_x,cell_y) + 1.0_num

        current => current%next
      END DO
    END DO

    data_array = data_array / MAX(part_count, c_tiny)

    DEALLOCATE(part_count)

  END SUBROUTINE calc_average_weight



  SUBROUTINE calc_temperature(sigma, current_species, direction)

    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: sigma
    INTEGER, INTENT(IN) :: current_species
    INTEGER, INTENT(IN), OPTIONAL :: direction
    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_pmx, part_pmy, part_pmz, sqrt_part_m
    ! The weight of a particle
    REAL(num) :: part_w
    REAL(num) :: gf
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: part_count, meanx, meany, meanz
    INTEGER :: ispecies, ix, iy, spec_start, spec_end
    TYPE(particle), POINTER :: current
    LOGICAL :: spec_sum
    REAL(num) :: dof, wdata
    INTEGER :: dir

#include "particle_head.inc"

    ALLOCATE(meanx(1-ng:nx+ng,1-ng:ny+ng))
    ALLOCATE(meany(1-ng:nx+ng,1-ng:ny+ng))
    ALLOCATE(meanz(1-ng:nx+ng,1-ng:ny+ng))
    ALLOCATE(part_count(1-ng:nx+ng,1-ng:ny+ng))
    meanx = 0.0_num
    meany = 0.0_num
    meanz = 0.0_num
    part_count = 0.0_num
    sigma = 0.0_num

    IF (PRESENT(direction)) THEN
      dir = direction
      dof = 1.0_num
    ELSE
      dir = -1
      dof = 3.0_num
    END IF

    spec_start = current_species
    spec_end = current_species
    spec_sum = .FALSE.

    IF (current_species <= 0) THEN
      spec_start = 1
      spec_end = n_species
      spec_sum = .TRUE.
    END IF

    DO ispecies = spec_start, spec_end
      IF (io_list(ispecies)%species_type == c_species_id_photon) CYCLE
#ifndef NO_TRACER_PARTICLES
      IF (spec_sum .AND. io_list(ispecies)%zero_current) CYCLE
#endif
      current => io_list(ispecies)%attached_list%head
      sqrt_part_m  = SQRT(io_list(ispecies)%mass)
      part_w = io_list(ispecies)%weight

      DO WHILE(ASSOCIATED(current))
#ifdef PER_PARTICLE_CHARGE_MASS
        sqrt_part_m  = SQRT(current%mass)
#endif
#ifndef PER_SPECIES_WEIGHT
        part_w = current%weight
#endif
        ! Copy the particle properties out for speed
        part_pmx = current%part_p(1) / sqrt_part_m
        part_pmy = current%part_p(2) / sqrt_part_m
        part_pmz = current%part_p(3) / sqrt_part_m

#include "particle_to_grid.inc"

        SELECT CASE(dir)
          CASE(c_dir_x)
            DO iy = sf_min, sf_max
            DO ix = sf_min, sf_max
              gf = gx(ix) * gy(iy) * part_w
              meanx(cell_x+ix, cell_y+iy) = &
                  meanx(cell_x+ix, cell_y+iy) + gf * part_pmx
              part_count(cell_x+ix, cell_y+iy) = &
                  part_count(cell_x+ix, cell_y+iy) + gf
            END DO
            END DO
          CASE(c_dir_y)
            DO iy = sf_min, sf_max
            DO ix = sf_min, sf_max
              gf = gx(ix) * gy(iy) * part_w
              meany(cell_x+ix, cell_y+iy) = &
                  meany(cell_x+ix, cell_y+iy) + gf * part_pmy
              part_count(cell_x+ix, cell_y+iy) = &
                  part_count(cell_x+ix, cell_y+iy) + gf
            END DO
            END DO
          CASE(c_dir_z)
            DO iy = sf_min, sf_max
            DO ix = sf_min, sf_max
              gf = gx(ix) * gy(iy) * part_w
              meanz(cell_x+ix, cell_y+iy) = &
                  meanz(cell_x+ix, cell_y+iy) + gf * part_pmz
              part_count(cell_x+ix, cell_y+iy) = &
                  part_count(cell_x+ix, cell_y+iy) + gf
            END DO
            END DO
          CASE DEFAULT
            DO iy = sf_min, sf_max
            DO ix = sf_min, sf_max
              gf = gx(ix) * gy(iy) * part_w
              meanx(cell_x+ix, cell_y+iy) = &
                  meanx(cell_x+ix, cell_y+iy) + gf * part_pmx
              meany(cell_x+ix, cell_y+iy) = &
                  meany(cell_x+ix, cell_y+iy) + gf * part_pmy
              meanz(cell_x+ix, cell_y+iy) = &
                  meanz(cell_x+ix, cell_y+iy) + gf * part_pmz
              part_count(cell_x+ix, cell_y+iy) = &
                  part_count(cell_x+ix, cell_y+iy) + gf
            END DO
            END DO
        END SELECT
        current => current%next
      END DO

      SELECT CASE(dir)
        CASE(c_dir_x)
          CALL calc_boundary(meanx, ispecies)
        CASE(c_dir_y)
          CALL calc_boundary(meany, ispecies)
        CASE(c_dir_z)
          CALL calc_boundary(meanz, ispecies)
        CASE DEFAULT
          CALL calc_boundary(meanx, ispecies)
          CALL calc_boundary(meany, ispecies)
          CALL calc_boundary(meanz, ispecies)
      END SELECT

      CALL calc_boundary(part_count, ispecies)
    END DO

    SELECT CASE(dir)
      CASE(c_dir_x)
        CALL calc_boundary(meanx)
      CASE(c_dir_y)
        CALL calc_boundary(meany)
      CASE(c_dir_z)
        CALL calc_boundary(meanz)
      CASE DEFAULT
        CALL calc_boundary(meanx)
        CALL calc_boundary(meany)
        CALL calc_boundary(meanz)
    END SELECT
    CALL calc_boundary(part_count)

    part_count = MAX(part_count, 1.e-6_num)

    meanx = meanx / part_count
    meany = meany / part_count
    meanz = meanz / part_count

    ! Restore ghost cell values for means
    SELECT CASE(dir)
      CASE(c_dir_x)
        CALL field_bc(meanx, ng)
      CASE(c_dir_y)
        CALL field_bc(meany, ng)
      CASE(c_dir_z)
        CALL field_bc(meanz, ng)
      CASE DEFAULT
        CALL field_bc(meanx, ng)
        CALL field_bc(meany, ng)
        CALL field_bc(meanz, ng)
    END SELECT

    part_count = 0.0_num
    DO ispecies = spec_start, spec_end
      IF (io_list(ispecies)%species_type == c_species_id_photon) CYCLE
#ifndef NO_TRACER_PARTICLES
      IF (spec_sum .AND. io_list(ispecies)%zero_current) CYCLE
#endif
      current => io_list(ispecies)%attached_list%head
      sqrt_part_m  = SQRT(io_list(ispecies)%mass)

      DO WHILE(ASSOCIATED(current))
#ifdef PER_PARTICLE_CHARGE_MASS
        sqrt_part_m  = SQRT(current%mass)
#endif

#ifndef PER_SPECIES_WEIGHT
        part_w = current%weight
#endif
        ! Copy the particle properties out for speed
        part_pmx = current%part_p(1) / sqrt_part_m
        part_pmy = current%part_p(2) / sqrt_part_m
        part_pmz = current%part_p(3) / sqrt_part_m

#include "particle_to_grid.inc"

        SELECT CASE(dir)
          CASE(c_dir_x)
            DO iy = sf_min, sf_max
            DO ix = sf_min, sf_max
              !gf = gx(ix) * gy(iy) 
              gf = gx(ix) * gy(iy) * part_w
              wdata = (part_pmx - meanx(cell_x+ix, cell_y+iy))**2
              sigma(cell_x+ix, cell_y+iy) = &
                  sigma(cell_x+ix, cell_y+iy) + gf * wdata
              part_count(cell_x+ix, cell_y+iy) = &
                  part_count(cell_x+ix, cell_y+iy) + gf
            END DO
            END DO
          CASE(c_dir_y)
            DO iy = sf_min, sf_max
            DO ix = sf_min, sf_max
              gf = gx(ix) * gy(iy) * part_w
              wdata = (part_pmy - meany(cell_x+ix, cell_y+iy))**2
              sigma(cell_x+ix, cell_y+iy) = &
                  sigma(cell_x+ix, cell_y+iy) + gf * wdata
              part_count(cell_x+ix, cell_y+iy) = &
                  part_count(cell_x+ix, cell_y+iy) + gf
            END DO
            END DO
          CASE(c_dir_z)
            DO iy = sf_min, sf_max
            DO ix = sf_min, sf_max
              gf = gx(ix) * gy(iy) * part_w
              wdata = (part_pmz - meanz(cell_x+ix, cell_y+iy))**2
              sigma(cell_x+ix, cell_y+iy) = &
                  sigma(cell_x+ix, cell_y+iy) + gf * wdata
              part_count(cell_x+ix, cell_y+iy) = &
                  part_count(cell_x+ix, cell_y+iy) + gf
            END DO
            END DO
          CASE DEFAULT
            DO iy = sf_min, sf_max
            DO ix = sf_min, sf_max
              gf = gx(ix) * gy(iy) * part_w
              wdata = (part_pmx - meanx(cell_x+ix, cell_y+iy))**2 &
                    + (part_pmy - meany(cell_x+ix, cell_y+iy))**2 &
                    + (part_pmz - meanz(cell_x+ix, cell_y+iy))**2
              sigma(cell_x+ix, cell_y+iy) = &
                  sigma(cell_x+ix, cell_y+iy) + gf * wdata
              part_count(cell_x+ix, cell_y+iy) = &
                  part_count(cell_x+ix, cell_y+iy) + gf
            END DO
            END DO
        END SELECT
        current => current%next
      END DO
      CALL calc_boundary(sigma, ispecies)
      CALL calc_boundary(part_count, ispecies)
    END DO

    CALL calc_boundary(sigma)
    CALL calc_boundary(part_count)

    ! 3/2 kT = <p^2>/(2m)
    sigma = sigma / MAX(part_count, 1.e-6_num) / kb / dof

    DEALLOCATE(part_count, meanx, meany, meanz)

  END SUBROUTINE calc_temperature



  SUBROUTINE calc_per_species_current(data_array, current_species, direction)

    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species
    INTEGER, INTENT(IN), OPTIONAL :: direction
    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_q, part_mc
    REAL(num) :: part_px, part_py, part_pz
    ! The data to be weighted onto the grid
    REAL(num) :: wdata
    REAL(num) :: fac, idx, root
    INTEGER :: ispecies, ix, iy, spec_start, spec_end
    TYPE(particle), POINTER :: current
    LOGICAL :: spec_sum
#include "particle_head.inc"

    IF (.NOT. PRESENT(direction)) THEN
      IF (rank == 0) THEN
        PRINT*, 'Error: No direction argument supplied to ', &
            'calc_per_species_current'
        CALL abort_code(c_err_bad_value)
      END IF
    END IF

    data_array = 0.0_num
    part_q = 0.0_num
    fac = 0.0_num

    idx = 1.0_num / dx / dy

    spec_start = current_species
    spec_end = current_species
    spec_sum = .FALSE.

    IF (current_species <= 0) THEN
      spec_start = 1
      spec_end = n_species
      spec_sum = .TRUE.
    END IF

    DO ispecies = spec_start, spec_end
      IF (io_list(ispecies)%species_type == c_species_id_photon) CYCLE
#ifndef NO_TRACER_PARTICLES
      IF (spec_sum .AND. io_list(ispecies)%zero_current) CYCLE
#endif
      current => io_list(ispecies)%attached_list%head
      part_mc = c * io_list(ispecies)%mass
      part_q  = io_list(ispecies)%charge
      fac = io_list(ispecies)%weight
      wdata = part_q * fac

      DO WHILE (ASSOCIATED(current))
        ! Copy the particle properties out for speed
#ifdef PER_PARTICLE_CHARGE_MASS
        part_mc = c * current%mass
        part_q  = current%charge
#ifndef PER_SPECIES_WEIGHT
        fac = current%weight
#endif
        wdata = part_q * fac
#else
#ifndef PER_SPECIES_WEIGHT
        fac = current%weight
        wdata = part_q * fac
#endif
#endif

        ! Copy the particle properties out for speed
        part_px = current%part_p(1)
        part_py = current%part_p(2)
        part_pz = current%part_p(3)
        ! XY
        root = 1.0_num / SQRT(part_mc**2 + part_px**2 + part_py**2 + part_pz**2)
        SELECT CASE (direction)
          CASE(c_dir_x)
            wdata = wdata * part_px * root
          CASE(c_dir_y)
            wdata = wdata * part_py * root
          CASE(c_dir_z)
            wdata = wdata * part_pz * root
        END SELECT

#include "particle_to_grid.inc"

        DO iy = sf_min, sf_max
        DO ix = sf_min, sf_max
          data_array(cell_x+ix, cell_y+iy) = &
              data_array(cell_x+ix, cell_y+iy) + gx(ix) * gy(iy) * wdata
        END DO
        END DO

        current => current%next
      END DO
      CALL calc_boundary(data_array, ispecies)
    END DO

    CALL calc_boundary(data_array)

    idx = c * idx
    data_array = data_array * idx
    DO ix = 1, 2*c_ndims
      CALL field_zero_gradient(data_array, c_stagger_centre, ix)
    END DO

  END SUBROUTINE calc_per_species_current



  SUBROUTINE calc_average_momentum(data_array, current_species, direction)

    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species
    INTEGER, INTENT(IN), OPTIONAL :: direction
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: part_count
    ! The data to be weighted onto the grid
    REAL(num) :: wdata, weight
    INTEGER :: ispecies, ix, iy, spec_start, spec_end
    TYPE(particle), POINTER :: current
    LOGICAL :: spec_sum
#include "particle_head.inc"

    IF (.NOT. PRESENT(direction)) THEN
      IF (rank == 0) THEN
        PRINT*, 'Error: No direction argument supplied to calc_average_momentum'
        CALL abort_code(c_err_bad_value)
      END IF
    END IF

    ALLOCATE(part_count(1-ng:nx+ng,1-ng:ny+ng))

    part_count = 0.0_num
    data_array = 0.0_num

    spec_start = current_species
    spec_end = current_species
    spec_sum = .FALSE.

    IF (current_species <= 0) THEN
      spec_start = 1
      spec_end = n_species
      spec_sum = .TRUE.
    END IF

    DO ispecies = spec_start, spec_end
      IF (io_list(ispecies)%species_type == c_species_id_photon) CYCLE
#ifndef NO_TRACER_PARTICLES
      IF (spec_sum .AND. io_list(ispecies)%zero_current) CYCLE
#endif
      current => io_list(ispecies)%attached_list%head

      weight = io_list(ispecies)%weight

      DO WHILE (ASSOCIATED(current))
        ! Copy the particle properties out for speed
#ifndef PER_SPECIES_WEIGHT
        weight = current%weight
#endif
        wdata = weight * current%part_p(direction)

#include "particle_to_grid.inc"

        DO iy = sf_min, sf_max
        DO ix = sf_min, sf_max
          data_array(cell_x+ix, cell_y+iy) = &
              data_array(cell_x+ix, cell_y+iy) + gx(ix) * gy(iy) * wdata
          part_count(cell_x+ix, cell_y+iy) = &
              part_count(cell_x+ix, cell_y+iy) + gx(ix) * gy(iy) * weight
        END DO
        END DO

        current => current%next
      END DO
      CALL calc_boundary(data_array, ispecies)
      CALL calc_boundary(part_count, ispecies)
    END DO

    CALL calc_boundary(data_array)
    CALL calc_boundary(part_count)

    data_array = data_array / MAX(part_count, c_tiny)
    DO ix = 1, 2*c_ndims
      CALL field_zero_gradient(data_array, c_stagger_centre, ix)
    END DO

    DEALLOCATE(part_count)

  END SUBROUTINE calc_average_momentum



  SUBROUTINE calc_total_energy_sum(per_species)

    LOGICAL, INTENT(IN) :: per_species
    REAL(num) :: particle_energy, field_energy
    REAL(num) :: part_ux, part_uy, part_uz, part_u2
    REAL(num) :: part_mc, part_w, fac, gamma_rel, gamma_rel_m1, part_energy
    REAL(num), ALLOCATABLE :: sum_out(:), sum_in(:)
    REAL(num), ALLOCATABLE :: species_energy(:)
    REAL(num), PARAMETER :: c2 = c**2
    INTEGER :: ispecies, i, j, nsum
    TYPE(particle), POINTER :: current

    particle_energy = 0.0_num
    IF (per_species) THEN
      ALLOCATE(species_energy(n_species))
      nsum = 1 + n_species
    ELSE
      nsum = 2
    END IF

    ! Sum over all particles to calculate total kinetic energy
    DO ispecies = 1, n_species
#ifndef NO_TRACER_PARTICLES
      !IF (species_list(ispecies)%zero_current) CYCLE
#endif
      current => species_list(ispecies)%attached_list%head
      part_mc = c * species_list(ispecies)%mass
      part_w = species_list(ispecies)%weight
      fac = part_mc * part_w * c
      part_energy = 0.0_num

      DO WHILE (ASSOCIATED(current))
        ! Copy the particle properties out for speed
#ifdef PER_PARTICLE_CHARGE_MASS
        part_mc = c * current%mass
#ifndef PER_SPECIES_WEIGHT
        part_w = current%weight
#endif
        fac = part_mc * part_w * c
#else
#ifndef PER_SPECIES_WEIGHT
        part_w = current%weight
        fac = part_mc * part_w * c
#endif
#endif

        IF (species_list(ispecies)%species_type /= c_species_id_photon) THEN
          part_ux = current%part_p(1) / part_mc
          part_uy = current%part_p(2) / part_mc
          part_uz = current%part_p(3) / part_mc

          part_u2 = part_ux**2 + part_uy**2 + part_uz**2
          gamma_rel = SQRT(part_u2 + 1.0_num)
          gamma_rel_m1 = part_u2 / (gamma_rel + 1.0_num)

          part_energy = part_energy + gamma_rel_m1 * fac
#if defined(PHOTONS) || defined(BREMSSTRAHLUNG)
        ELSE
          part_energy = part_energy + current%particle_energy * part_w
#endif
        END IF

        current => current%next
      END DO
        !XY
      IF (per_species) THEN 
          !这个用于计算光子能量,而不用输出光子,但是一旦生成正负电子对的话，会高估
          !解决办法是在正负电子对那里把remove光子的时候同时减掉这个光子能量
          !但是太麻烦了，因为正负电子对本来就需要光子生成.
        IF (species_list(ispecies)%species_type == c_species_id_photon &
            .AND. (.NOT. produce_pairs)) THEN
        part_energy = photon_energy_record
        END IF
        species_energy(ispecies) = part_energy
      END IF
      particle_energy = particle_energy + part_energy
    END DO

    ! EM field energy
    field_energy = 0.0_num
    DO j = 1, ny
    DO i = 1, nx
      field_energy = field_energy + ex(i,j)**2 + ey(i,j)**2 &
          + ez(i,j)**2 + c2 * (bx(i,j)**2 + by(i,j)**2 + bz(i,j)**2)
    END DO
    END DO
    field_energy = 0.5_num * epsilon0 * field_energy * dx * dy

    ALLOCATE(sum_out(nsum))
    ALLOCATE(sum_in(nsum))
    sum_out(1) = field_energy
    IF (per_species) THEN
      sum_out(2:1+n_species) = species_energy(:)
    ELSE
      sum_out(2) = particle_energy
    END IF

    CALL MPI_REDUCE(sum_out, sum_in, nsum, mpireal, MPI_SUM, 0, comm, errcode)
    total_field_energy = sum_in(1)
    IF (per_species) THEN
      total_particle_energy_species(:) = sum_in(2:1+n_species)
      total_particle_energy = SUM(sum_in(2:1+n_species))
    ELSE
      total_particle_energy = sum_in(2)
    END IF

  END SUBROUTINE calc_total_energy_sum



  SUBROUTINE calc_initial_current

    REAL(num), DIMENSION(:,:), ALLOCATABLE :: jx, jy, jz
    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_q, part_mc
    REAL(num) :: part_px, part_py, part_pz
    ! The data to be weighted onto the grid
    REAL(num) :: wdata
    REAL(num) :: part_jx, part_jy, part_jz
    REAL(num) :: fac, idx, root
    REAL(num) :: sum_in(3), sum_out(3)
    INTEGER :: ispecies, ix, iy
    TYPE(particle), POINTER :: current
#include "particle_head.inc"

    ALLOCATE(jx(1-jng:nx+jng, 1-jng:ny+jng))
    ALLOCATE(jy(1-jng:nx+jng, 1-jng:ny+jng))
    ALLOCATE(jz(1-jng:nx+jng, 1-jng:ny+jng))

    jx = 0.0_num
    jy = 0.0_num
    jz = 0.0_num

    idx = 1.0_num / dx / dy

    DO ispecies = 1, n_species
      IF (species_list(ispecies)%species_type == c_species_id_photon) CYCLE
#ifndef NO_TRACER_PARTICLES
      IF (species_list(ispecies)%zero_current) CYCLE
#endif
      current => species_list(ispecies)%attached_list%head
      part_mc = c * species_list(ispecies)%mass
      part_q  = species_list(ispecies)%charge
      fac = species_list(ispecies)%weight
      wdata = part_q * fac

      DO WHILE (ASSOCIATED(current))
        ! Copy the particle properties out for speed
#ifdef PER_PARTICLE_CHARGE_MASS
        part_mc = c * current%mass
        part_q  = current%charge
#ifndef PER_SPECIES_WEIGHT
        fac = current%weight
#endif
        wdata = part_q * fac
#else
#ifndef PER_SPECIES_WEIGHT
        fac = current%weight
        wdata = part_q * fac
#endif
#endif

        ! Copy the particle properties out for speed
        part_px = current%part_p(1)
        part_py = current%part_p(2)
        part_pz = current%part_p(3)
        !XY  four current vector \rho_0 *U = (c\rho,jx,jy,jz)
        !假设密度是实验室系下看到的密度已经乘上了Gamma_rel
        root = 1.0_num / SQRT(part_mc**2 + part_px**2 + part_py**2 + part_pz**2)

        part_jx = wdata * part_px * root
        part_jy = wdata * part_py * root
        part_jz = wdata * part_pz * root
!        part_jx = wdata * part_px 
!        part_jy = wdata * part_py 
!        part_jz = wdata * part_pz

#include "particle_to_grid.inc"

        DO iy = sf_min, sf_max
        DO ix = sf_min, sf_max
          jx(cell_x+ix, cell_y+iy) = &
              jx(cell_x+ix, cell_y+iy) + gx(ix) * gy(iy) * part_jx
          jy(cell_x+ix, cell_y+iy) = &
              jy(cell_x+ix, cell_y+iy) + gx(ix) * gy(iy) * part_jy
          jz(cell_x+ix, cell_y+iy) = &
              jz(cell_x+ix, cell_y+iy) + gx(ix) * gy(iy) * part_jz
        END DO
        END DO

        current => current%next
      END DO
    END DO

    sum_out(1) = SUM(jx)
    sum_out(2) = SUM(jy)
    sum_out(3) = SUM(jz)
    DEALLOCATE(jx, jy, jz)

    CALL MPI_ALLREDUCE(sum_out, sum_in, 3, mpireal, MPI_SUM, comm, errcode)

    fac = c * idx / nx_global / ny_global

    initial_jx = sum_in(1) * fac
    initial_jy = sum_in(2) * fac
    initial_jz = sum_in(3) * fac

  END SUBROUTINE calc_initial_current


  !For optical intensity and polarization
    SUBROUTINE calc_optical_polarization(data_array, current_species, direction)
    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN), OPTIONAL :: direction
    INTEGER, INTENT(IN) :: current_species
    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_m, part_w,part_mc
    ! The data to be weighted onto the grid
    REAL(num) :: wdata
    REAL(num) :: fac, idx
    INTEGER :: ispecies, ix, iy, spec_start, spec_end
    TYPE(particle), POINTER :: current
    LOGICAL :: spec_sum
    !For polarization
    REAL(num), DIMENSION(3) :: dir1,dir2,e_at_part,b_at_part
    REAL(num) :: part_ux, part_uy, part_uz, part_u2,part_x,part_y
    REAL(num) :: I1,I2,Uo,Vo,I1n,I2n,Un,Vn, I1dirx,I1diry,I1dirz,px,py,pz
    REAL(num) :: cos_phi,Omega,cos_theta,theta
    
#include "particle_head.inc"

    data_array = 0.0_num
    part_w = 1.0_num
    part_m = 0.0_num
    fac = 0.0_num

    !idx = 1.0_num / dx / dy

    spec_start = current_species
    spec_end = current_species
    spec_sum = .FALSE.

    !Calculate_detection_direction for I(phi) and I(phi + pi/2)
    !fix a direction for (1,0,0) in the x direction
    !de_dir1 = (/ 1.0_num,0.0_num,0.0_num /)

    IF (current_species <= 0) THEN
      spec_start = 1
      spec_end = n_species
      spec_sum = .TRUE.
    END IF

    DO ispecies = spec_start, spec_end
      IF ((.NOT. spec_sum) .AND. (io_list(ispecies)%species_type == c_species_id_electron &
            .OR. io_list(ispecies)%species_type == c_species_id_positron)) THEN ! optical radiation from electron and positron
      current => io_list(ispecies)%attached_list%head
      part_mc = c * species_list(ispecies)%mass
      part_w = species_list(ispecies)%weight
      fac = io_list(ispecies)%weight

      DO WHILE (ASSOCIATED(current))
        ! Copy the particle properties out for speed
        part_ux = current%part_p(1) / part_mc
        part_uy = current%part_p(2) / part_mc
        part_uz = current%part_p(3) / part_mc
        part_u2 = part_ux**2 + part_uy**2 + part_uz**2
        !judgement the direction
        cos_phi = (part_ux*dir_ob(1) + part_uy*dir_ob(2) + part_uz*dir_ob(3))/SQRT(part_u2)
        Omega = 2*pi*(1 - cos_phi)
        IF  (Omega < dOmega) THEN 
            !calculate optical_part_energy
            part_x  = current%part_pos(1) - x_grid_min_local
            part_y  = current%part_pos(2) - y_grid_min_local
            CALL field_at_particle(part_x,part_y,e_at_part, b_at_part)
            CALL calculate_optical_stokes2(p = current%part_p,b = b_at_part,wei = part_w, &
                                           I1 = I1,I2 = I2,U = Uo,V = Vo,dir1 = dir1,dir2 = dir2)
            cos_theta = dir1(1)*de_dir1(1) + dir1(2)*de_dir1(2) + dir1(3)*de_dir1(3)
            theta = ACOS(cos_theta)
            ! Rotate the proper axes to the detection direction
            CALL calculate_rotate_I(theta,I1,I2,Uo,Vo,I1n,I2n,Un,Vn)
            IF (debug_mod) THEN
              WRITE(*,*) 'I1,I1n:',I1,'->',I1n
              WRITE(*,*) 'I2,I2n:',I2,'->',I2n
              WRITE(*,*) 'Uo,Un:',Uo,'->',Un
              WRITE(*,*) 'Vo,Vn:',Vo,'->',Vn
            END IF

        SELECT Case(direction)
            Case(c_dir_stokes_I)
                wdata = (I1n + I2n)
            Case(c_dir_stokes_Q)
                wdata = (I1n - I2n)
            Case(c_dir_stokes_U)
                wdata = Un
            Case(c_dir_stokes_V)
                wdata = Vn
            Case(c_dir_stokes_N)
                wdata = 1 
            !for calculate critical value
        END SELECT
#include "particle_to_grid.inc"
        DO iy = sf_min, sf_max
        DO ix = sf_min, sf_max
          data_array(cell_x+ix, cell_y+iy) = &
              data_array(cell_x+ix, cell_y+iy) + gx(ix) * gy(iy) * wdata
        END DO
        END DO
       END IF! For Omega
        current => current%next
      END DO ! for current 

      !need thinking maybe clamp_zero
      CALL calc_boundary(data_array, ispecies)
    END IF  !if electron and positron
  END DO ! for ispecies

    CALL calc_boundary(data_array)

    DO ix = 1, 2*c_ndims
      CALL field_zero_gradient(data_array, c_stagger_centre, ix)
    END DO
  END SUBROUTINE calc_optical_polarization


  SUBROUTINE calc_total_optical_intensity()

    REAL(num) :: optical_I1,optical_I2,optical_U,optical_V
    REAL(num) :: I1,I2,Uo,Vo,I1n,I2n,Un,Vn
    REAL(num), DIMENSION(3) :: dir1,dir2,e_at_part,b_at_part
    REAL(num) :: part_ux, part_uy, part_uz, part_u2,part_x,part_y
    REAL(num) :: part_mc, part_w, fac, gamma_e
    REAL(num) :: cos_phi,Omega   ! angle between radiation and observation
    REAL(num) :: cos_theta,theta ! for proper rotation
    REAL(num), ALLOCATABLE :: sum_out(:), sum_in(:)
    REAL(num), PARAMETER :: c2 = c**2
    INTEGER :: ispecies, i, j, nsum
    TYPE(particle), POINTER :: current

    optical_I1 = 0.0_num
    optical_I2 = 0.0_num
    optical_U = 0.0_num
    optical_V = 0.0_num
    nsum = 4
    !total_optical_I1
    !total_optical_I2
    !total_optical_U
    !total_optical_V
    ! Sum over all particles to calculate total kinetic energy
    DO ispecies = 1, n_species
#ifndef NO_TRACER_PARTICLES
      IF (species_list(ispecies)%zero_current) CYCLE
#endif
      IF (species_list(ispecies)%species_type == c_species_id_electron &
          .OR. species_list(ispecies)%species_type == c_species_id_positron) THEN
      current => species_list(ispecies)%attached_list%head
      part_mc = c * species_list(ispecies)%mass
      part_w = species_list(ispecies)%weight
      fac = part_mc * part_w * c

      DO WHILE (ASSOCIATED(current))
        ! Copy the particle properties out for speed
#ifdef PER_PARTICLE_CHARGE_MASS
        part_mc = c * current%mass
#ifndef PER_SPECIES_WEIGHT
        part_w = current%weight
#endif
        fac = part_mc * part_w * c
#else
#ifndef PER_SPECIES_WEIGHT
        part_w = current%weight
        fac = part_mc * part_w * c
#endif
#endif
        part_ux = current%part_p(1) / part_mc
        part_uy = current%part_p(2) / part_mc
        part_uz = current%part_p(3) / part_mc
        part_u2 = part_ux**2 + part_uy**2 + part_uz**2
        !judgement the direction
        cos_phi = (part_ux*dir_ob(1) + part_uy*dir_ob(2) + part_uz*dir_ob(3))/SQRT(part_u2)
        Omega = 2*pi*(1 - cos_phi)
        IF  (Omega < dOmega) THEN 
            !calculate optical_part_energy
            part_x  = current%part_pos(1) - x_grid_min_local
            part_y  = current%part_pos(2) - y_grid_min_local
            CALL field_at_particle(part_x,part_y,e_at_part, b_at_part)

            CALL calculate_optical_stokes2(current%part_p,b_at_part,part_w, &
                                            I1, I2,Uo,Vo,dir1,dir2)
!            I1 = 0.0_num
!            I2 = 0.0_num
!            dir1 = (/1.0,0.0,0.0/)
!            dir2 = (/0.0,1.0,0.0/)
!
            cos_theta = dir1(1)*de_dir1(1) + dir1(2)*de_dir1(2) + dir1(3)*de_dir1(3)
            theta = ACOS(cos_theta)
            CALL calculate_rotate_I(theta,I1,I2,Uo,Vo,I1n,I2n,Un,Vn)
            !test for no rotation
            !CALL calculate_rotate_I(0.0_num,I1,I2,Uo,Vo,I1n,I2n,Un,Vn)
            !WRITE(*,*) I1,I2,I1n,I2n,Un

            optical_I1 = optical_I1 + I1n
            optical_I2 = optical_I2 + I2n
            optical_U  = optical_U  + Un
            optical_V  = optical_V  + Vn
        END IF
            current => current%next
      END DO ! Current
    END IF !electron and positron
    END DO !do species

    ALLOCATE(sum_out(nsum))
    ALLOCATE(sum_in(nsum))
    sum_out(1) = optical_I1
    sum_out(2) = optical_I2
    sum_out(3) = optical_U
    sum_out(4) = optical_V


    CALL MPI_REDUCE(sum_out, sum_in, nsum, mpireal, MPI_SUM, 0, comm, errcode)
    total_optical_I1 = sum_in(1)
    total_optical_I2 = sum_in(2)
    total_optical_U  = sum_in(3)
    total_optical_V  = sum_in(4)

  END SUBROUTINE calc_total_optical_intensity


    SUBROUTINE calc_polarization(data_array, current_species, direction)
    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN), OPTIONAL :: direction
    INTEGER, INTENT(IN) :: current_species
    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_m, part_w
    ! The data to be weighted onto the grid
    REAL(num) :: wdata
    REAL(num) :: fac, idx
    INTEGER :: ispecies, ix, iy, spec_start, spec_end
    TYPE(particle), POINTER :: current
    LOGICAL :: spec_sum
    !For polarization
    REAL(num) :: I1o,I2o,Uo,Vo,I1n,I2n,Un,Vn, I1dirx,I1diry,I1dirz,px,py,pz
    REAL(num) :: cos_phi,Omega,cos_theta,theta
    
#include "particle_head.inc"

    data_array = 0.0_num
    part_w = 1.0_num
    part_m = 0.0_num
    fac = 0.0_num

    !idx = 1.0_num / dx / dy

    spec_start = current_species
    spec_end = current_species
    spec_sum = .FALSE.

    !Calculate_detection_direction for I(phi) and I(phi + pi/2)
    !fix a direction for (1,0,0) in the x direction
    !de_dir1 = (/ 1.0_num,0.0_num,0.0_num /)

    IF (current_species <= 0) THEN
      spec_start = 1
      spec_end = n_species
      spec_sum = .TRUE.
    END IF

    DO ispecies = spec_start, spec_end
      IF ((.NOT. spec_sum) .AND. io_list(ispecies)%species_type == c_species_id_photon) THEN
      current => io_list(ispecies)%attached_list%head
      fac = io_list(ispecies)%weight

      DO WHILE (ASSOCIATED(current))
        ! Copy the particle properties out for speed
        px = current%part_p(1)/m0/c
        py = current%part_p(2)/m0/c
        pz = current%part_p(3)/m0/c
        I1o = current%stokes(1)
        I2o = current%stokes(2)
        part_w = current%weight
        wdata = 0.0_num

        !judge the photon moving direction
        cos_phi = (px*dir_ob(1) + py*dir_ob(2) + pz*dir_ob(3))/SQRT(px**2 + py**2 + pz**2)
        Omega = 2*pi*(1 - cos_phi)
        IF  (Omega < dOmega) THEN 
            I1o = current%stokes(1)
            I2o = current%stokes(2)
            Uo = 0.0_num
            Vo = 0.0_num
            I1dirx = current%stokes_dir(1)
            I1diry = current%stokes_dir(2)
            I1dirz = current%stokes_dir(3)
            cos_theta = I1dirx*de_dir1(1) + I1diry*de_dir1(2) + I1dirz*de_dir1(3)
            theta = ACOS(cos_theta)
            CALL calculate_rotate_I(theta,I1o,I2o,Uo,Vo,I1n,I2n,Un,Vn)
        SELECT Case(direction)
            Case(c_dir_stokes_I)
                wdata = (I1n + I2n)
            Case(c_dir_stokes_Q)
                wdata = (I1n - I2n)
            Case(c_dir_stokes_U)
                wdata = Un
            Case(c_dir_stokes_N)
                wdata = 1 
            !for calculate critical value
        END SELECT
#include "particle_to_grid.inc"
        DO iy = sf_min, sf_max
        DO ix = sf_min, sf_max
          data_array(cell_x+ix, cell_y+iy) = &
              data_array(cell_x+ix, cell_y+iy) + gx(ix) * gy(iy) * wdata
        END DO
        END DO
       END IF! For Omega
        current => current%next
      END DO ! for current 

      !need thinking maybe clamp_zero
      CALL calc_boundary(data_array, ispecies)
  END IF  !if photon
  END DO ! for ispecies

    CALL calc_boundary(data_array)

    DO ix = 1, 2*c_ndims
      CALL field_zero_gradient(data_array, c_stagger_centre, ix)
    END DO
  END SUBROUTINE calc_polarization

  SUBROUTINE calculate_rotate_I(theta, I1o, I2o,Uo,Vo,I1n,I2n,Un,Vn)
    REAL(num), INTENT(IN) :: I1o, I2o, Uo, Vo, theta
    REAL(num), INTENT(OUT) :: I1n,I2n,Un,Vn
    REAL(num):: Io, Qo,eta 
    I1n = cos(theta)**2 * I1o + sin(theta)**2 * I2o + 0.5*sin(2*theta)*Uo
    I2n = cos(theta)**2 * I2o + sin(theta)**2 * I1o - 0.5*sin(2*theta)*Uo
    Un =  sin(2*theta)*(-I1o) + sin(2*theta)*(I2o) + cos(2*theta)*Uo
    Vn = Vo
  END SUBROUTINE calculate_rotate_I


  !XY  for polarization
  !calculate polarizition degree in each cell by considering the magnetic field is uniform in each cell and need decide the
  !direction of observation
    SUBROUTINE calc_polarization_degree(data_array, current_species, direction)
    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: wt,I1,I2
    INTEGER, INTENT(IN), OPTIONAL :: direction
    INTEGER, INTENT(IN) :: current_species
    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_m, part_w
    ! The data to be weighted onto the grid
    REAL(num) :: wdata
    REAL(num) :: fac, idx
    INTEGER :: ispecies, ix, iy, spec_start, spec_end
    TYPE(particle), POINTER :: current
    LOGICAL :: spec_sum
    !For polarization
    REAL(num) :: I1o,I2o,Uo,Vo, I1_de1, I1_de2, I1dirx,I1diry,I1dirz,px,py,pz, I,Q
    REAL(num) :: cos_phi,Omega,cos_theta,theta
    
#include "particle_head.inc"

    ALLOCATE(I1(1-ng:nx+ng,1-ng:ny+ng))
    ALLOCATE(I2(1-ng:nx+ng,1-ng:ny+ng))
    data_array = 0.0_num
    I1 = 0.0_num
    I2 = 0.0_num

    part_w = 1.0_num
    part_m = 0.0_num
    fac = 0.0_num

    !idx = 1.0_num / dx / dy

    spec_start = current_species
    spec_end = current_species
    spec_sum = .FALSE.

    !Calculate_detection_direction for I(phi) and I(phi + pi/2)
    !fix a direction for (1,0,0) in the x direction
    !de_dir1 = (/ 1.0_num,0.0_num,0.0_num /)

    IF (current_species <= 0) THEN
      spec_start = 1
      spec_end = n_species
      spec_sum = .TRUE.
    END IF

    DO ispecies = spec_start, spec_end
      IF ((.NOT. spec_sum) .AND. io_list(ispecies)%species_type == c_species_id_photon) THEN
      current => io_list(ispecies)%attached_list%head
      fac = io_list(ispecies)%weight

      DO WHILE (ASSOCIATED(current))
        ! Copy the particle properties out for speed
        px = current%part_p(1)/m0/c
        py = current%part_p(2)/m0/c
        pz = current%part_p(3)/m0/c
        I1o = current%stokes(1)
        I2o = current%stokes(2)
        part_w = current%weight
        wdata = current%particle_energy*part_w

        !judge the photon moving direction
        cos_phi = (px*dir_ob(1) + py*dir_ob(2) + pz*dir_ob(3))/SQRT(px**2 + py**2 + pz**2)
        Omega = 2*pi*(1 - cos_phi)
        IF  (Omega < dOmega) THEN 
            I1o = current%stokes(1)
            I2o = current%stokes(2)
            Uo = 0.0_num
            Vo = 0.0_num
            I1dirx = current%stokes_dir(1)
            I1diry = current%stokes_dir(2)
            I1dirz = current%stokes_dir(3)
            cos_theta = I1dirx*de_dir1(1) + I1diry*de_dir1(2) + I1dirz*de_dir1(3)
            theta = ACOS(cos_theta)
            CALL calculate_detect_I(theta,I1o,I2o,Uo,Vo,I1_de1)
            cos_theta = I1dirx*de_dir2(1) + I1diry*de_dir2(2) + I1dirz*de_dir2(3)
            theta = ACOS(cos_theta)
            CALL calculate_detect_I(theta,I1o,I2o,Uo,Vo,I1_de2)
#include "particle_to_grid.inc"
        DO iy = sf_min, sf_max
        DO ix = sf_min, sf_max
          I1(cell_x+ix, cell_y+iy) = &
              I1(cell_x+ix, cell_y+iy) + gx(ix) * gy(iy) * I1_de1
          I2(cell_x+ix, cell_y+iy) = &
              I2(cell_x+ix, cell_y+iy) + gx(ix) * gy(iy) * I1_de2
        END DO
        END DO
       END IF! For Omega
        current => current%next
      END DO ! for current 

      !need thinking maybe clamp_zero
      CALL calc_boundary(data_array, ispecies)
  END IF  !if photon
  END DO ! for current
    Do ix = 1, nx
    Do iy = 1, ny
        I = I1(ix,iy) + I2(ix,iy)
        Q = I1(ix,iy) - I2(ix,iy)
        I = MAX(c_tiny,I)
        data_array(ix,iy) =  Q/I
    ENDDO
    ENDDO

    CALL calc_boundary(data_array)

    DO ix = 1, 2*c_ndims
      CALL field_zero_gradient(data_array, c_stagger_centre, ix)
    END DO

    DEALLOCATE(I1)
    DEALLOCATE(I2)
  END SUBROUTINE calc_polarization_degree

    SUBROUTINE calc_polar_intensity(data_array, current_species, direction)
    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: wt,I1,I2
    INTEGER, INTENT(IN), OPTIONAL :: direction
    INTEGER, INTENT(IN) :: current_species
    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_m, part_w
    ! The data to be weighted onto the grid
    REAL(num) :: wdata
    REAL(num) :: fac, idx
    INTEGER :: ispecies, ix, iy, spec_start, spec_end
    TYPE(particle), POINTER :: current
    LOGICAL :: spec_sum
    !For polarization
    REAL(num) :: I1o,I2o,Uo,Vo, I1_de1, I1_de2, I1dirx,I1diry,I1dirz,px,py,pz, I
    REAL(num) :: cos_phi,Omega,cos_theta,theta
    
#include "particle_head.inc"

    data_array = 0.0_num

    part_w = 1.0_num
    part_m = 0.0_num
    fac = 0.0_num

    !idx = 1.0_num / dx / dy

    spec_start = current_species
    spec_end = current_species
    spec_sum = .FALSE.

    !Calculate_detection_direction for I(phi) and I(phi + pi/2)
    !fix a direction for (1,0,0) in the x direction
    !de_dir1 = (/ 1.0_num,0.0_num,0.0_num /)

    IF (current_species <= 0) THEN
      spec_start = 1
      spec_end = n_species
      spec_sum = .TRUE.
    END IF

    DO ispecies = spec_start, spec_end
      IF ((.NOT. spec_sum) .AND. io_list(ispecies)%species_type == c_species_id_photon) THEN
      current => io_list(ispecies)%attached_list%head
      fac = io_list(ispecies)%weight

      DO WHILE (ASSOCIATED(current))
        ! Copy the particle properties out for speed
        px = current%part_p(1)/m0/c
        py = current%part_p(2)/m0/c
        pz = current%part_p(3)/m0/c
        I1o = current%stokes(1)
        I2o = current%stokes(2)
        part_w = current%weight
        wdata = current%particle_energy*part_w

        !judge the photon moving direction
        cos_phi = (px*dir_ob(1) + py*dir_ob(2) + pz*dir_ob(3))/SQRT(px**2 + py**2 + pz**2)
        Omega = 2*pi*(1 - cos_phi)
        IF  (Omega < dOmega) THEN 
            I1o = current%stokes(1)
            I2o = current%stokes(2)
            Uo = 0.0_num
            Vo = 0.0_num
            I1dirx = current%stokes_dir(1)
            I1diry = current%stokes_dir(2)
            I1dirz = current%stokes_dir(3)

            SELECT CASE(direction)
                CASE(c_dir_polar1)
                cos_theta = I1dirx*de_dir1(1) + I1diry*de_dir1(2) + I1dirz*de_dir1(3)
                theta = ACOS(cos_theta)
                CALL calculate_detect_I(theta,I1o,I2o,Uo,Vo,I1_de1)
                CASE(c_dir_polar2)
                cos_theta = I1dirx*de_dir2(1) + I1diry*de_dir2(2) + I1dirz*de_dir2(3)
                theta = ACOS(cos_theta)
                CALL calculate_detect_I(theta,I1o,I2o,Uo,Vo,I1_de1)
            END SELECT
            I1_de1 = MAX(c_tiny, I1_de1)
!
#include "particle_to_grid.inc"
        DO iy = sf_min, sf_max
        DO ix = sf_min, sf_max
          data_array(cell_x+ix, cell_y+iy) = &
              data_array(cell_x+ix, cell_y+iy) + gx(ix) * gy(iy) * I1_de1
        END DO
        END DO
       END IF! For Omega
        current => current%next
      END DO ! for current 

      !need thinking maybe clamp_zero
      CALL calc_boundary(data_array, ispecies)
  END IF  !if photon
    END DO ! for current

    CALL calc_boundary(data_array)

    !data_array = data_array / Max(wt,c_tiny)

    DO ix = 1, 2*c_ndims
      CALL field_zero_gradient(data_array, c_stagger_centre, ix)
    END DO

  END SUBROUTINE calc_polar_intensity


  SUBROUTINE calculate_detect_I(theta, I1o, I2o,Uo,Vo,I1_de)
    REAL(num), INTENT(IN) :: I1o, I2o, Uo, Vo, theta
    REAL(num), INTENT(OUT) :: I1_de
    REAL(num):: Io, Qo,eta 
    eta = 0.0_num
    Io = I1o+I2o
    Qo = I1o-I2o
    I1_de = 0.5*(Io + Qo*cos(2*theta))
    !+ (Uo*cos(eta) - Vo*sin(eta))*sin(2*theta))
  END SUBROUTINE calculate_detect_I
!Fluid statistics

  SUBROUTINE calc_per_species_rel_temperature(data_array,current_species,direction)
    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: u1,u2,u3
    !Boost 
    REAL(num) :: A1,A2,A3,A4,A5,A6,A7,A8,A9,A10
    REAL(num) :: Gamma_bulk, uu1, uu2, uu3, u_2
    !Particle viriable in labatory frame
    REAL(num) :: part_px, part_py, part_pz, part_m, gamma_rel_c,idx
    !Particle viriable in rest frame
    REAL(num) :: part_px_rest, part_py_rest, part_pz_rest, gamma_rest
    ! The weight of a particle
    REAL(num) :: part_w
    REAL(num) :: gf
    INTEGER, INTENT(IN) :: current_species
    INTEGER, INTENT(IN),OPTIONAL :: direction 
    INTEGER :: ispecies, ix, iy, spec_start, spec_end
    TYPE(particle), POINTER :: current
    LOGICAL :: spec_sum
#include "particle_head.inc"
    

    data_array = 0.0_num

    idx = 1.0_num / dx / dy

    spec_start = current_species
    spec_end = current_species
    spec_sum = .FALSE.

    IF (current_species <= 0) THEN
      spec_start = 1
      spec_end = n_species
      spec_sum = .TRUE.
    ENDIF
    ! if current_species <=0 means calculate the whole velocity
    ALLOCATE(u1(1-ng:nx+ng, 1-ng:ny+ng))
    u1 = 0.0_num
    CALL calc_per_species_velocity(u1, current_species, c_dir_x)
    ALLOCATE(u2(1-ng:nx+ng, 1-ng:ny+ng))
    u2 = 0.0_num
    CALL calc_per_species_velocity(u2, current_species, c_dir_y)
    ALLOCATE(u3(1-ng:nx+ng, 1-ng:ny+ng))
    u3 = 0.0_num
    CALL calc_per_species_velocity(u3, current_species, c_dir_z)
    u1 = u1/c 
    u2 = u2/c 
    u3 = u3/c 

    DO ispecies = spec_start, spec_end
#ifndef NO_TRACER_PARTICLES
      IF (spec_sum .AND. io_list(ispecies)%zero_current) CYCLE
#endif
      current => io_list(ispecies)%attached_list%head
      part_m  = io_list(ispecies)%mass
      part_w  = io_list(ispecies)%weight

      DO WHILE(ASSOCIATED(current))
#ifdef PER_PARTICLE_CHARGE_MASS
        part_m  = current%mass
#endif

#ifndef PER_SPECIES_WEIGHT
        part_w = current%weight
#endif
        ! Copy the particle properties out for speed
        part_px = current%part_p(1) / part_m / c
        part_py = current%part_p(2) / part_m / c
        part_pz = current%part_p(3) / part_m / c 
        gamma_rel_c = SQRT(1.0_num + part_px**2 + part_py**2 + part_pz**2) 
        ! Copy gamma_rel vx, gamma_rel vy, gamma_rel vz

#include "particle_to_grid.inc"

        DO iy = sf_min, sf_max
        DO ix = sf_min, sf_max
            ! Boost see jackson's book P547
            uu1 = u1(cell_x+ix,cell_y+iy)
            uu2 = u2(cell_x+ix,cell_y+iy)
            uu3 = u3(cell_x+ix,cell_y+iy)
            Gamma_bulk = SQRT(1.0/(1.0 - uu1*uu1 - uu2*uu2 - uu3*uu3))
            !IF (rank == 0 ) WRITE(*,*) Gamma_bulk
            !1.0/(1.0 - (uu1*uu1 + uu2*uu2 + uu3*uu3))
            A1 = Gamma_bulk
            A2 = -A1*uu1
            A3 = -A1*uu2
            A4 = -A1*uu3
            u_2 = uu1**2 + uu2**2+ uu3**2
            A5 = 1 + (A1-1) * uu1*uu1/u_2
            A6 =     (A1-1) * uu1*uu2/u_2
            A7 =     (A1-1) * uu1*uu3/u_2
            A8 = 1 + (A1-1) * uu2*uu2/u_2
            A9 =     (A1-1) * uu2*uu3/u_2
            A10 =1 + (A1-1) * uu3*uu3/u_2
            ! Boost rule matrix = (1,-1,-1,-1) 
            gamma_rest =   A1 * gamma_rel_c + A2* part_px + A3 * part_py + A4 * part_pz
            part_px_rest = A2 * gamma_rel_c + A5* part_px + A6 * part_py + A7 * part_pz
            part_py_rest = A3 * gamma_rel_c + A6* part_px + A8 * part_py + A9 * part_pz
            part_pz_rest = A4 * gamma_rel_c + A7* part_px + A9 * part_py + A10 * part_pz
            !IF (rank == 0) WRITE(*,*) Gamma_bulk, gamma_rel_c, gamma_rest,
            !IF (rank == 0) WRITE(*,*) gamma_rest**2 - part_px_rest**2 - part_py_rest**2 - part_pz_rest**2
            ! Statistics T_mu_nu = c*int f dp3/p0 p_mu*p_nu 
            gf = gx(ix) * gy(iy) *part_w
            SELECT CASE (direction)
                CASE (c_dir_t)
                data_array(cell_x+ix, cell_y+iy) = &
                    data_array(cell_x+ix, cell_y+iy) + gf * &
               gamma_rest*part_m*c**2 
                CASE (c_dir_x)
                data_array(cell_x+ix, cell_y+iy) = data_array(cell_x+ix, cell_y+iy) + gf &
              * (part_px_rest)**2/gamma_rest*part_m*c**2
                CASE (c_dir_y)
                data_array(cell_x+ix, cell_y+iy) = data_array(cell_x+ix, cell_y+iy) + gf &
              * (part_py_rest)**2/gamma_rest*part_m*c**2 
                CASE (c_dir_z)
                data_array(cell_x+ix, cell_y+iy) = data_array(cell_x+ix, cell_y+iy) + gf &
              * (part_pz_rest)**2/gamma_rest*part_m*c**2 
            END SELECT
        ENDDO
        ENDDO
        current => current%next
      ENDDO
      CALL calc_boundary(data_array, ispecies)
    ENDDO
    CALL calc_boundary(data_array)
    data_array = data_array * idx

    DO ix = 1, 2*c_ndims
      CALL field_zero_gradient(data_array, c_stagger_centre, ix)
    END DO

  END SUBROUTINE calc_per_species_rel_temperature

   SUBROUTINE calc_per_species_rel_velocity(data_array, current_species, direction)
    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species
    INTEGER, INTENT(IN), OPTIONAL :: direction
    ! Properties of the current particles. Copy out of particle arrays for speed
    REAL(num) :: part_mc, part_w
    REAL(num) :: part_ux, part_uy, part_uz, part_u2
    REAL(num) :: gamma_rel, part_p_gamma
    ! The data to be weighted onto the grid
    REAL(num) :: wdata
    REAL(num) :: gf
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: part_count
    INTEGER :: ispecies, spec_start, spec_end
    INTEGER :: ix, iy
    TYPE(particle), POINTER :: current
    LOGICAL :: spec_sum
#include "particle_head.inc"


    data_array = 0.0_num
    ALLOCATE (part_count(1-ng:nx+ng, 1-ng:ny+ng))
    part_count = 0.0_num

    spec_start = current_species
    spec_end = current_species
    spec_sum = .FALSE.

    IF (current_species <= 0) THEN
        spec_start = 1
        spec_end = n_species
        spec_sum = .TRUE.
    END IF

    DO ispecies = spec_start, spec_end
      IF (spec_sum .AND. &
          io_list(ispecies)%species_type == c_species_id_photon) CYCLE
#ifndef NO_TRACER_PARTICLES
      IF (spec_sum .AND. io_list(ispecies)%zero_current) CYCLE
#endif
      current => io_list(ispecies)%attached_list%head
      part_mc = c * io_list(ispecies)%mass
      part_w = io_list(ispecies)%weight
      wdata = part_mc * part_w
      
      DO WHILE (ASSOCIATED(current))
#ifdef PER_PARTICLE_CHARGE_MASS
        part_mc = c * current%mass
#ifndef PER_SPECIES_WEIGHT
        part_w = current%weight
#endif
        wdata = part_mc * part_w
#else
#ifndef PER_SPECIES_WEIGHT
        part_w = current%weight
        wdata = part_mc * part_w
#endif
#endif
        ! Copy the particle properties out for speed
        part_ux = current%part_p(1) / part_mc
        part_uy = current%part_p(2) / part_mc
        part_uz = current%part_p(3) / part_mc

        part_u2 = part_ux**2 + part_uy**2 + part_uz**2
        gamma_rel = SQRT(part_u2 + 1.0_num)

        SELECT CASE (direction)
          CASE (c_dir_x)
            part_p_gamma = part_ux 
          CASE (c_dir_y)
            part_p_gamma = part_uy
          CASE (c_dir_z)
            part_p_gamma = part_uz
        END SELECT

#include "particle_to_grid.inc"

        DO iy = sf_min, sf_max
        DO ix = sf_min, sf_max
          gf = gx(ix) * gy(iy) * wdata
          data_array(cell_x+ix, cell_y+iy) = &
              data_array(cell_x+ix, cell_y+iy) + gf * part_p_gamma
          part_count(cell_x+ix, cell_y+iy) = &
              part_count(cell_x+ix, cell_y+iy) + gf
        END DO
        END DO
        current => current%next
      END DO
      CALL calc_boundary(data_array, ispecies)
      CALL calc_boundary(part_count, ispecies)
    END DO
    CALL calc_boundary(data_array)
    CALL calc_boundary(part_count)

    part_count = MAX(part_count, 1.e-31_num)
    data_array = c * data_array / part_count

    DO ix = 1, 2*c_ndims
      CALL field_zero_gradient(data_array, c_stagger_centre, ix)
    END DO

    DEALLOCATE(part_count)

  END SUBROUTINE calc_per_species_rel_velocity

  SUBROUTINE calc_per_species_velocity(data_array, current_species, direction)
    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species, direction
    ! Properties of the current particles. Copy out of particle arrays for speed
    REAL(num) :: part_mc, part_w
    REAL(num) :: part_ux, part_uy, part_uz, part_u2
    REAL(num) :: gamma_rel, part_p_gamma
    ! The data to be weighted onto the grid
    REAL(num) :: wdata
    REAL(num) :: gf
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: part_count
    INTEGER :: ispecies, spec_start, spec_end
    INTEGER :: ix, iy
    TYPE(particle), POINTER :: current
    LOGICAL :: spec_sum
#include "particle_head.inc"


    data_array = 0.0_num
    ALLOCATE (part_count(1-ng:nx+ng, 1-ng:ny+ng))
    part_count = 0.0_num

    spec_start = current_species
    spec_end = current_species
    spec_sum = .FALSE.

    IF (current_species <= 0) THEN
        spec_start = 1
        spec_end = n_species
        spec_sum = .TRUE.
    END IF

    DO ispecies = spec_start, spec_end
      IF (spec_sum .AND. &
          io_list(ispecies)%species_type == c_species_id_photon) CYCLE
#ifndef NO_TRACER_PARTICLES
      IF (spec_sum .AND. io_list(ispecies)%zero_current) CYCLE
#endif
      current => io_list(ispecies)%attached_list%head
      part_mc = c * io_list(ispecies)%mass
      part_w = io_list(ispecies)%weight
      wdata = part_mc * part_w
      
      DO WHILE (ASSOCIATED(current))
#ifdef PER_PARTICLE_CHARGE_MASS
        part_mc = c * current%mass
#ifndef PER_SPECIES_WEIGHT
        part_w = current%weight
#endif
        wdata = part_mc * part_w
#else
#ifndef PER_SPECIES_WEIGHT
        part_w = current%weight
        wdata = part_mc * part_w
#endif
#endif
        ! Copy the particle properties out for speed
        part_ux = current%part_p(1) / part_mc
        part_uy = current%part_p(2) / part_mc
        part_uz = current%part_p(3) / part_mc

        part_u2 = part_ux**2 + part_uy**2 + part_uz**2
        gamma_rel = SQRT(part_u2 + 1.0_num)

        SELECT CASE (direction)
          CASE (c_dir_x)
            part_p_gamma = part_ux / gamma_rel
          CASE (c_dir_y)
            part_p_gamma = part_uy / gamma_rel
          CASE (c_dir_z)
            part_p_gamma = part_uz / gamma_rel
        END SELECT

#include "particle_to_grid.inc"

        DO iy = sf_min, sf_max
        DO ix = sf_min, sf_max
          gf = gx(ix) * gy(iy) * wdata
          data_array(cell_x+ix, cell_y+iy) = &
              data_array(cell_x+ix, cell_y+iy) + gf * part_p_gamma
          part_count(cell_x+ix, cell_y+iy) = &
              part_count(cell_x+ix, cell_y+iy) + gf
        END DO
        END DO
        current => current%next
      END DO
      CALL calc_boundary(data_array, ispecies)
      CALL calc_boundary(part_count, ispecies)
    END DO
    CALL calc_boundary(data_array)
    CALL calc_boundary(part_count)

    part_count = MAX(part_count, 1.e-31_num)
    data_array = c * data_array / part_count

    DO ix = 1, 2*c_ndims
      CALL field_zero_gradient(data_array, c_stagger_centre, ix)
    END DO

    DEALLOCATE(part_count)

  END SUBROUTINE calc_per_species_velocity


  SUBROUTINE calc_per_species_ux(data_array, current_species)

    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species
    
        CALL calc_per_species_velocity(data_array, current_species, c_dir_x)

  END SUBROUTINE calc_per_species_ux



  SUBROUTINE calc_per_species_uy(data_array, current_species)

    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species

        CALL calc_per_species_velocity(data_array, current_species, c_dir_y)

  END SUBROUTINE calc_per_species_uy



  SUBROUTINE calc_per_species_uz(data_array, current_species)

    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species
        CALL calc_per_species_velocity(data_array, current_species, c_dir_z)

  END SUBROUTINE calc_per_species_uz


  SUBROUTINE calc_per_species_pressure_tensor(data_array, current_species, u1, u2, directions)
    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species
    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(IN) :: u1, u2
    INTEGER, DIMENSION(2), INTENT(IN) :: directions
    ! Properties of the current particles. Copy out of particle arrays for speed
    REAL(num) :: part_mc, part_w
    REAL(num) :: part_ux, part_uy, part_uz, part_u2
    REAL(num) :: gamma_rel_c
    ! The data to be weighted onto the grid
    REAL(num), DIMENSION(2):: wdata
    REAL(num) :: gf
    !REAL(num), DIMENSION(:,:), ALLOCATABLE :: part_count
    INTEGER :: ispecies, spec_start, spec_end
    INTEGER :: ix, iy
    REAL(num) :: idx
    TYPE(particle), POINTER :: current
    LOGICAL :: spec_sum
#include "particle_head.inc"

    data_array = 0.0_num
    !ALLOCATE (part_count(1-ng:nx+ng, 1-ng:ny+ng))
    !part_count = 0.0_num
    idx = 1.0_num / dx / dy

    spec_start = current_species
    spec_end = current_species
    spec_sum = .FALSE.

    IF (current_species <= 0) THEN
        spec_start = 1
        spec_end = n_species
        spec_sum = .TRUE.
    END IF

    DO ispecies = spec_start, spec_end
      IF (spec_sum .AND. &
          io_list(ispecies)%species_type == c_species_id_photon) CYCLE
#ifndef NO_TRACER_PARTICLES
      IF (spec_sum .AND. io_list(ispecies)%zero_current) CYCLE
#endif
      current => io_list(ispecies)%attached_list%head
      part_mc = c * io_list(ispecies)%mass
      part_w = io_list(ispecies)%weight
      
      DO WHILE (ASSOCIATED(current))
#ifdef PER_PARTICLE_CHARGE_MASS
        part_mc = c * current%mass
#ifndef PER_SPECIES_WEIGHT
        part_w = current%weight
#endif
#else
#ifndef PER_SPECIES_WEIGHT
        part_w = current%weight
#endif
#endif
        ! Copy the particle properties out for speed
        part_ux = current%part_p(1) / part_mc
        part_uy = current%part_p(2) / part_mc
        part_uz = current%part_p(3) / part_mc

        part_u2 = part_ux**2 + part_uy**2 + part_uz**2
        gamma_rel_c = SQRT(part_u2 + 1.0_num) / c

        SELECT CASE (directions(1))
          CASE (c_dir_x)
            wdata(1) = part_ux
          CASE (c_dir_y)
            wdata(1) = part_uy
          CASE (c_dir_z)
            wdata(1) = part_uz
        END SELECT
        SELECT CASE (directions(2))
          CASE (c_dir_x)
            wdata(2) = part_ux
          CASE (c_dir_y)
            wdata(2) = part_uy
          CASE (c_dir_z)
            wdata(2) = part_uz
        END SELECT


#include "particle_to_grid.inc"
        

        DO iy = sf_min, sf_max
        DO ix = sf_min, sf_max
          !Gamma1 = 1.0/(1.0- (u1(cell_x+ix,cell_y+iy)**2 + u2(cell_x+ix,cell_y+iy)**2 + u3(cell_x+ix,cell_y+iy)**2)/c**2)
          gf = gx(ix) * gy(iy) * part_w
          data_array(cell_x+ix, cell_y+iy) = &
              data_array(cell_x+ix, cell_y+iy) + gf * (part_mc) * &
              !Xiey 
              (wdata(1) - gamma_rel_c * u1(cell_x+ix, cell_y+iy)) * (wdata(2) - gamma_rel_c * u2(cell_x+ix, cell_y+iy)) / & 
              gamma_rel_c
        END DO
        END DO
        current => current%next
      END DO
      CALL calc_boundary(data_array, ispecies)
    END DO
    CALL calc_boundary(data_array)

    data_array = data_array * idx

    DO ix = 1, 2*c_ndims
      CALL field_zero_gradient(data_array, c_stagger_centre, ix)
    END DO

  END SUBROUTINE calc_per_species_pressure_tensor

  SUBROUTINE calc_per_species_Pxx(data_array, current_species,direction)

    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species

    REAL(num), DIMENSION(:,:), ALLOCATABLE :: velocity_1
    INTEGER, INTENT(IN), OPTIONAL :: direction
    INTEGER, DIMENSION(2) :: directions

    ALLOCATE(velocity_1(1-ng:nx+ng, 1-ng:ny+ng))


    velocity_1 = 0.0_num
    CALL calc_per_species_velocity(velocity_1, current_species, c_dir_x)

    directions = (/ c_dir_x, c_dir_x /)
    CALL calc_per_species_pressure_tensor(data_array, current_species, velocity_1, velocity_1, directions)

    DEALLOCATE(velocity_1)
  END SUBROUTINE calc_per_species_Pxx



  SUBROUTINE calc_per_species_Pxy(data_array, current_species,direction)

    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species

    REAL(num), DIMENSION(:,:), ALLOCATABLE :: velocity_1
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: velocity_2
    INTEGER, INTENT(IN), OPTIONAL :: direction
    INTEGER, DIMENSION(2) :: directions

    ALLOCATE(velocity_1(1-ng:nx+ng, 1-ng:ny+ng))
    ALLOCATE(velocity_2(1-ng:nx+ng, 1-ng:ny+ng))


    velocity_1 = 0.0_num
    CALL calc_per_species_velocity(velocity_1, current_species, c_dir_x)
    velocity_2 = 0.0_num
    CALL calc_per_species_velocity(velocity_2, current_species, c_dir_y)

    directions = (/ c_dir_x, c_dir_y /)
    CALL calc_per_species_pressure_tensor(data_array, current_species, velocity_1, velocity_2, directions)

    DEALLOCATE(velocity_1)
    DEALLOCATE(velocity_2)
  END SUBROUTINE calc_per_species_Pxy



  SUBROUTINE calc_per_species_Pxz(data_array, current_species,direction)

    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species

    REAL(num), DIMENSION(:,:), ALLOCATABLE :: velocity_1
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: velocity_2
    INTEGER, INTENT(IN), OPTIONAL :: direction
    INTEGER, DIMENSION(2) :: directions

    ALLOCATE(velocity_1(1-ng:nx+ng, 1-ng:ny+ng))
    ALLOCATE(velocity_2(1-ng:nx+ng, 1-ng:ny+ng))


    velocity_1 = 0.0_num
    CALL calc_per_species_velocity(velocity_1, current_species, c_dir_x)
    velocity_2 = 0.0_num
    CALL calc_per_species_velocity(velocity_2, current_species, c_dir_z)

    directions = (/ c_dir_x, c_dir_z /)
    CALL calc_per_species_pressure_tensor(data_array, current_species, velocity_1, velocity_2, directions)

    DEALLOCATE(velocity_1)
    DEALLOCATE(velocity_2)
  END SUBROUTINE calc_per_species_Pxz


  SUBROUTINE calc_per_species_Pyy(data_array, current_species,direction)

    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species

    REAL(num), DIMENSION(:,:), ALLOCATABLE :: velocity_1
    INTEGER, INTENT(IN), OPTIONAL :: direction
    INTEGER, DIMENSION(2) :: directions

    ALLOCATE(velocity_1(1-ng:nx+ng, 1-ng:ny+ng))


    velocity_1 = 0.0_num
    CALL calc_per_species_velocity(velocity_1, current_species, c_dir_y)

    directions = (/ c_dir_y, c_dir_y /)
    CALL calc_per_species_pressure_tensor(data_array, current_species, velocity_1, velocity_1, directions)

    DEALLOCATE(velocity_1)
  END SUBROUTINE calc_per_species_Pyy


  SUBROUTINE calc_per_species_Pyz(data_array, current_species,direction)

    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species

    REAL(num), DIMENSION(:,:), ALLOCATABLE :: velocity_1
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: velocity_2
    INTEGER, INTENT(IN), OPTIONAL :: direction
    INTEGER, DIMENSION(2) :: directions

    ALLOCATE(velocity_1(1-ng:nx+ng, 1-ng:ny+ng))
    ALLOCATE(velocity_2(1-ng:nx+ng, 1-ng:ny+ng))


    velocity_1 = 0.0_num
    CALL calc_per_species_velocity(velocity_1, current_species, c_dir_y)
    velocity_2 = 0.0_num
    CALL calc_per_species_velocity(velocity_2, current_species, c_dir_z)

    directions = (/ c_dir_y, c_dir_z /)
    CALL calc_per_species_pressure_tensor(data_array, current_species, velocity_1, velocity_2, directions)

    DEALLOCATE(velocity_1)
    DEALLOCATE(velocity_2)
  END SUBROUTINE calc_per_species_Pyz


  SUBROUTINE calc_per_species_Pzz(data_array, current_species,direction)

    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species

    REAL(num), DIMENSION(:,:), ALLOCATABLE :: velocity_1
    INTEGER, INTENT(IN), OPTIONAL :: direction
    INTEGER, DIMENSION(2) :: directions

    ALLOCATE(velocity_1(1-ng:nx+ng, 1-ng:ny+ng))


    velocity_1 = 0.0_num
    CALL calc_per_species_velocity(velocity_1, current_species, c_dir_z)

    directions = (/ c_dir_z, c_dir_z /)
    CALL calc_per_species_pressure_tensor(data_array, current_species, velocity_1, velocity_1, directions)

    DEALLOCATE(velocity_1)
  END SUBROUTINE calc_per_species_Pzz

  SUBROUTINE calc_per_species_thermal_flux(data_array, current_species, u1, u2, u3, direction)
    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species
    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(IN) :: u1, u2, u3
    INTEGER, INTENT(IN),OPTIONAL :: direction
    ! Properties of the current particles. Copy out of particle arrays for speed
    REAL(num) :: part_mc, part_w
    REAL(num) :: part_ux, part_uy, part_uz, part_u2
    REAL(num) :: gamma_rel_c, part_w2
    ! The data to be weighted onto the grid
    REAL(num), DIMENSION(3):: part_w1
    REAL(num) :: gf
    !REAL(num), DIMENSION(:,:), ALLOCATABLE :: part_count
    INTEGER :: ispecies, spec_start, spec_end
    INTEGER :: ix, iy, w_direction
    REAL(num) :: idx
    TYPE(particle), POINTER :: current
    LOGICAL :: spec_sum
#include "particle_head.inc"

    data_array = 0.0_num
    !ALLOCATE (part_count(1-ng:nx+ng, 1-ng:ny+ng))
    !part_count = 0.0_num
    idx = 1.0_num / dx / dy

    spec_start = current_species
    spec_end = current_species
    spec_sum = .FALSE.

    IF (current_species <= 0) THEN
        spec_start = 1
        spec_end = n_species
        spec_sum = .TRUE.
    END IF

    DO ispecies = spec_start, spec_end
      IF (spec_sum .AND. &
          io_list(ispecies)%species_type == c_species_id_photon) CYCLE
#ifndef NO_TRACER_PARTICLES
      IF (spec_sum .AND. io_list(ispecies)%zero_current) CYCLE
#endif
      current => io_list(ispecies)%attached_list%head
      part_mc = c * io_list(ispecies)%mass
      part_w = io_list(ispecies)%weight
      
      DO WHILE (ASSOCIATED(current))
#ifdef PER_PARTICLE_CHARGE_MASS
        part_mc = c * current%mass
#ifndef PER_SPECIES_WEIGHT
        part_w = current%weight
#endif
#else
#ifndef PER_SPECIES_WEIGHT
        part_w = current%weight
#endif
#endif
        ! Copy the particle properties out for speed
        part_ux = current%part_p(1) / part_mc
        part_uy = current%part_p(2) / part_mc
        part_uz = current%part_p(3) / part_mc

        part_u2 = part_ux**2 + part_uy**2 + part_uz**2
        gamma_rel_c = SQRT(part_u2 + 1.0_num) / c


#include "particle_to_grid.inc"

        part_w1 = 0.0_num
        DO iy = sf_min, sf_max
        DO ix = sf_min, sf_max
          gf = gx(ix) * gy(iy) * part_w
          part_w1(1) = part_ux - gamma_rel_c * u1(cell_x+ix, cell_y+iy)
          part_w1(2) = part_uy - gamma_rel_c * u2(cell_x+ix, cell_y+iy)
          part_w1(3) = part_uz - gamma_rel_c * u3(cell_x+ix, cell_y+iy)
          part_w2 = DOT_PRODUCT(part_w1, part_w1)
          data_array(cell_x+ix, cell_y+iy) = &
              data_array(cell_x+ix, cell_y+iy) + gf * (part_mc) * &
              part_w2 * part_w1(direction) / (gamma_rel_c**2)
        END DO
        END DO
        current => current%next
      END DO
      CALL calc_boundary(data_array, ispecies)
    END DO
    CALL calc_boundary(data_array)

    data_array = 0.5 * data_array * idx

    DO ix = 1, 2*c_ndims
      CALL field_zero_gradient(data_array, c_stagger_centre, ix)
    END DO

  END SUBROUTINE calc_per_species_thermal_flux



  SUBROUTINE calc_thermal_flux(data_array, current_species, direction)
    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species
    INTEGER, INTENT(IN),OPTIONAL :: direction
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: velocity_1, velocity_2, velocity_3

    ALLOCATE(velocity_1(1-ng:nx+ng, 1-ng:ny+ng))
    velocity_1 = 0.0_num
    CALL calc_per_species_velocity(velocity_1, current_species, c_dir_x)
    ALLOCATE(velocity_2(1-ng:nx+ng, 1-ng:ny+ng))
    velocity_2 = 0.0_num
    CALL calc_per_species_velocity(velocity_2, current_species, c_dir_y)
    ALLOCATE(velocity_3(1-ng:nx+ng, 1-ng:ny+ng))
    velocity_3 = 0.0_num
    CALL calc_per_species_velocity(velocity_3, current_species, c_dir_z)

    CALL calc_per_species_thermal_flux(data_array, current_species, velocity_1, velocity_2, velocity_3, direction)

    DEALLOCATE(velocity_1)
    DEALLOCATE(velocity_2)
    DEALLOCATE(velocity_3)
  END SUBROUTINE calc_thermal_flux



  SUBROUTINE calc_per_species_gyrotropy(data_array, current_species)
    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species
    ! Define array to evaluate bulk velocity.
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: u1, u2, u3
    ! Define array to evaluate pressure tensor components.
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: pxx, pyy, pzz
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: pxy, pxz, pyz
    ! Define array to hold B unit vector.
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: Babs
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: Btheta
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: Pparallel
    ! Invariant number
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: I1, I2
    INTEGER, DIMENSION(2) :: directions
    REAL(num) :: denominator_local
    INTEGER :: ix, iy

    ! Compute bulk velocity
    ALLOCATE(u1(1-ng:nx+ng, 1-ng:ny+ng))
    CALL calc_per_species_ux(u1, current_species)
    ALLOCATE(u2(1-ng:nx+ng, 1-ng:ny+ng))
    CALL calc_per_species_ux(u2, current_species)
    ALLOCATE(u3(1-ng:nx+ng, 1-ng:ny+ng))
    CALL calc_per_species_ux(u3, current_species)

    ! Compute pressure tensor
    ALLOCATE(pxx(1-ng:nx+ng, 1-ng:ny+ng))
    directions = (/ c_dir_x, c_dir_x /)
    CALL calc_per_species_pressure_tensor(pxx, current_species, u1, u1, directions)
    ALLOCATE(pyy(1-ng:nx+ng, 1-ng:ny+ng))
    directions = (/ c_dir_y, c_dir_y /)
    CALL calc_per_species_pressure_tensor(pyy, current_species, u2, u2, directions)
    ALLOCATE(pzz(1-ng:nx+ng, 1-ng:ny+ng))
    directions = (/ c_dir_z, c_dir_z /)
    CALL calc_per_species_pressure_tensor(pzz, current_species, u3, u3, directions)
    ALLOCATE(pxy(1-ng:nx+ng, 1-ng:ny+ng))
    directions = (/ c_dir_x, c_dir_y /)
    CALL calc_per_species_pressure_tensor(pxy, current_species, u1, u2, directions)
    ALLOCATE(pxz(1-ng:nx+ng, 1-ng:ny+ng))
    directions = (/ c_dir_x, c_dir_z /)
    CALL calc_per_species_pressure_tensor(pxz, current_species, u1, u3, directions)
    ALLOCATE(pyz(1-ng:nx+ng, 1-ng:ny+ng))
    directions = (/ c_dir_y, c_dir_z /)
    CALL calc_per_species_pressure_tensor(pyz, current_species, u2, u3, directions)

    ! Compute magnetic field
    ALLOCATE(Babs(1-ng:nx+ng, 1-ng:ny+ng))
    Babs = SQRT(bx**2 + by**2 + bz**2)
    ALLOCATE(Btheta(3, 1-ng:nx+ng, 1-ng:ny+ng))
    DO iy = 1-ng, ny+ng
    DO ix = 1-ng, nx+ng
        IF (Babs(ix,iy) .NE. 0.0_num) THEN
            Btheta(1, ix, iy) = bx(ix, iy) / Babs(ix, iy)
            Btheta(2, ix, iy) = by(ix, iy) / Babs(ix, iy)
            Btheta(3, ix, iy) = bz(ix, iy) / Babs(ix, iy)
        ELSE
            Btheta(1, ix, iy) = 0.0_num
            Btheta(2, ix, iy) = 0.0_num
            Btheta(3, ix, iy) = 0.0_num
        END IF
    END DO
    END DO

    ! Compute P parallel
    ALLOCATE(Pparallel(1-ng:nx+ng, 1-ng:ny+ng))
    DO iy = 1-ng, ny+ng
    DO ix = 1-ng, nx+ng
        Pparallel(ix, iy) = &
            (Btheta(1,ix,iy)**2)*pxx(ix,iy) + (Btheta(2,ix,iy)**2)*pyy(ix,iy) + (Btheta(3,ix,iy)**2)*pzz(ix,iy) + &
            2.0*(Btheta(1,ix,iy)*Btheta(2,ix,iy)*pxy(ix,iy) + Btheta(1,ix,iy)*Btheta(3,ix,iy)*pxz(ix,iy) + &
                 Btheta(2,ix,iy)*Btheta(3,ix,iy)*pyz(ix,iy))
    END DO
    END DO

    ! Compute I1 and I2, which are invariant under coordinate ratations.
    ! I1, the trace of pressure tensor
    ALLOCATE(I1(1-ng:nx+ng, 1-ng:ny+ng))
    I1 = pxx + pyy + pzz
    ! I2, sums of principal minors
    ALLOCATE(I2(1-ng:nx+ng, 1-ng:ny+ng))
    I2 = (pxx*pyy + pxx*pzz + pyy*pzz) - (pxy**2 + pxz**2 + pyz**2)

    ! Compute gyrotropy 
    DO iy = 1-ng, ny+ng
    DO ix = 1-ng, nx+ng
        denominator_local = (I1(ix,iy) - Pparallel(ix,iy)) * (I1(ix,iy) + 3.0*Pparallel(ix,iy))
        data_array(ix,iy) = 1.0_num - 4.0_num * I2(ix,iy) / denominator_local
        !IF (denominator_local .NE. 0.0_num) THEN
        !    data_array(ix,iy) = 1.0_num - 4.0_num * I2(ix,iy) / denominator_local
        !ELSE
        !    data_array(ix,iy) = 0.0_num
        !END IF
    END DO
    END DO

    ! Deallocate array
    DEALLOCATE(u1)
    DEALLOCATE(u2)
    DEALLOCATE(u3)
    DEALLOCATE(pxx)
    DEALLOCATE(pxy)
    DEALLOCATE(pxz)
    DEALLOCATE(pyy)
    DEALLOCATE(pyz)
    DEALLOCATE(pzz)
    DEALLOCATE(Babs)
    DEALLOCATE(Btheta)
    DEALLOCATE(Pparallel)
    DEALLOCATE(I1)
    DEALLOCATE(I2)
  END SUBROUTINE calc_per_species_gyrotropy


  SUBROUTINE calc_magnetic_vector_potential(data_array,current_species,direction)
    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species
    INTEGER, INTENT(IN), OPTIONAL :: direction
    REAL(num) :: integral_x, integral_y
    INTEGER :: ix, iy

    integral_y = 0.0_num
    DO iy = 1-ng, ny+ng
        integral_y = integral_y + bx(1-ng, iy)
        integral_x = 0.0_num
        DO ix = 1-ng, nx+ng
            integral_x = integral_x - by(ix, iy)
            data_array(ix, iy) = integral_x * dx + integral_y * dy
        END DO
    END DO
        
  END SUBROUTINE calc_magnetic_vector_potential


  SUBROUTINE calc_per_species_Hall_electric_field(data_array, current_species, direction)
    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species
    INTEGER, INTENT(IN),OPTIONAL :: direction
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: charge_density

    ! Compute charge density
    ALLOCATE(charge_density(1-ng:nx+ng, 1-ng:ny+ng))
    CALL calc_charge_density(charge_density, current_species)

    ! Select direction
    SELECT CASE(direction)
    CASE(c_dir_x)
        data_array = -(jy*bz - jz*by) / charge_density
    CASE(c_dir_y)
        data_array = -(jz*bx - jx*bz) / charge_density
    CASE(c_dir_z)
        data_array = -(jx*by - jy*bx) / charge_density
    END SELECT
  END SUBROUTINE calc_per_species_Hall_electric_field


  SUBROUTINE calc_per_species_convection_electric_field(data_array, current_species, direction)
    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species
    INTEGER, INTENT(IN),OPTIONAL :: direction
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: velocity_1
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: velocity_2

    ALLOCATE(velocity_1(1-ng:nx+ng, 1-ng:ny+ng))
    ALLOCATE(velocity_2(1-ng:nx+ng, 1-ng:ny+ng))

    ! Select direction
    SELECT CASE(direction)
    CASE(c_dir_x)
        CALL calc_per_species_velocity(velocity_1, current_species, c_dir_y)
        CALL calc_per_species_velocity(velocity_2, current_species, c_dir_z)
        data_array = -(velocity_1*bz - velocity_2*by)
    CASE(c_dir_y)
        CALL calc_per_species_velocity(velocity_1, current_species, c_dir_z)
        CALL calc_per_species_velocity(velocity_2, current_species, c_dir_x)
        data_array = -(velocity_1*bx - velocity_2*bz)
    CASE(c_dir_z)
        CALL calc_per_species_velocity(velocity_1, current_species, c_dir_x)
        CALL calc_per_species_velocity(velocity_2, current_species, c_dir_y)
        data_array = -(velocity_1*by - velocity_2*bx)
    END SELECT
  END SUBROUTINE calc_per_species_convection_electric_field


  SUBROUTINE calc_per_species_Nernst_electric_field(data_array, current_species, direction)
    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species
    INTEGER, INTENT(IN),OPTIONAL :: direction
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: thermal_qx
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: thermal_qy
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: thermal_qz
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: pressure_p

    ALLOCATE(pressure_p(1-ng:nx+ng, 1-ng:ny+ng))

    SELECT CASE(direction)
    CASE(c_dir_x)
        CALL calc_per_species_Pxx(pressure_p, current_species)

        ALLOCATE(thermal_qy(1-ng:nx+ng, 1-ng:ny+ng))
        CALL calc_thermal_flux(thermal_qy, current_species, c_dir_y)
        ALLOCATE(thermal_qz(1-ng:nx+ng, 1-ng:ny+ng))
        CALL calc_thermal_flux(thermal_qz, current_species, c_dir_z)

        data_array = -0.4_num * (thermal_qy*bz - thermal_qz*by) / pressure_p
    CASE(c_dir_y)
        CALL calc_per_species_Pyy(pressure_p, current_species)

        ALLOCATE(thermal_qx(1-ng:nx+ng, 1-ng:ny+ng))
        CALL calc_thermal_flux(thermal_qx, current_species, c_dir_x)
        ALLOCATE(thermal_qz(1-ng:nx+ng, 1-ng:ny+ng))
        CALL calc_thermal_flux(thermal_qz, current_species, c_dir_z)

        data_array = -0.4_num * (thermal_qz*bx - thermal_qx*bz) / pressure_p
    CASE(c_dir_z)
        CALL calc_per_species_Pzz(pressure_p, current_species)

        ALLOCATE(thermal_qx(1-ng:nx+ng, 1-ng:ny+ng))
        CALL calc_thermal_flux(thermal_qx, current_species, c_dir_x)
        ALLOCATE(thermal_qy(1-ng:nx+ng, 1-ng:ny+ng))
        CALL calc_thermal_flux(thermal_qy, current_species, c_dir_y)

        data_array = -0.4_num * (thermal_qx*by - thermal_qy*bx) / pressure_p
    END SELECT

  END SUBROUTINE calc_per_species_Nernst_electric_field

  SUBROUTINE calc_per_species_Biermann_electric_field(data_array, current_species, direction)
    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species
    INTEGER, INTENT(IN),OPTIONAL :: direction
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: pressure_p
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: charge_density

    ALLOCATE(pressure_p(1-ng:nx+ng, 1-ng:ny+ng))

    SELECT CASE(direction)
    CASE(c_dir_x)
        CALL calc_per_species_Pxx(pressure_p, current_species)
    CASE(c_dir_y)
        CALL calc_per_species_Pyy(pressure_p, current_species)
    CASE(c_dir_z)
        CALL calc_per_species_Pzz(pressure_p, current_species)
    END SELECT

    ! Compute derivation
    CALL calc_derivate(pressure_p, direction)

    ! Compute charge density
    ALLOCATE(charge_density(1-ng:nx+ng, 1-ng:ny+ng))
    CALL calc_charge_density(charge_density, current_species)

    data_array = pressure_p / charge_density

  END SUBROUTINE calc_per_species_Biermann_electric_field

  SUBROUTINE calc_per_species_viscous_electric_field(data_array, current_species, directions)
    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species
    INTEGER, DIMENSION(3), INTENT(IN) :: directions
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: pressure_p
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: charge_density

    ! Compute charge density
    ALLOCATE(charge_density(1-ng:nx+ng, 1-ng:ny+ng))
    CALL calc_charge_density(charge_density, current_species)

    ! Compute off diagnal pressure tensor
    ALLOCATE(pressure_p(1-ng:nx+ng, 1-ng:ny+ng))
    IF ((directions(1) == c_dir_x) .AND. (directions(2) == c_dir_y)) THEN
        CALL calc_per_species_Pxy(pressure_p, current_species)
    ELSE IF ((directions(1) == c_dir_x) .AND. (directions(2) == c_dir_z)) THEN
        CALL calc_per_species_Pxz(pressure_p, current_species)
    ELSE IF ((directions(1) == c_dir_y) .AND. (directions(2) == c_dir_z)) THEN
        CALL calc_per_species_Pyz(pressure_p, current_species)
    END IF

    ! Compute derivate
    CALL calc_derivate(pressure_p, directions(3))

    data_array = pressure_p / charge_density
  END SUBROUTINE calc_per_species_viscous_electric_field


  SUBROUTINE calc_viscous_electric_field_xyx(data_array, current_species)
    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species
    INTEGER, DIMENSION(3) :: directions

    directions = (/ c_dir_x, c_dir_y, c_dir_x /)
    CALL calc_per_species_viscous_electric_field(data_array, current_species, directions)
  END SUBROUTINE calc_viscous_electric_field_xyx


  SUBROUTINE calc_viscous_electric_field_xyy(data_array, current_species)
    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species
    INTEGER, DIMENSION(3) :: directions

    directions = (/ c_dir_x, c_dir_y, c_dir_y /)
    CALL calc_per_species_viscous_electric_field(data_array, current_species, directions)
  END SUBROUTINE calc_viscous_electric_field_xyy


  SUBROUTINE calc_viscous_electric_field_xzx(data_array, current_species)
    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species
    INTEGER, DIMENSION(3) :: directions

    directions = (/ c_dir_x, c_dir_z, c_dir_x /)
    CALL calc_per_species_viscous_electric_field(data_array, current_species, directions)
  END SUBROUTINE calc_viscous_electric_field_xzx


  SUBROUTINE calc_viscous_electric_field_xzz(data_array, current_species)
    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species
    INTEGER, DIMENSION(3) :: directions

    directions = (/ c_dir_x, c_dir_z, c_dir_z /)
    CALL calc_per_species_viscous_electric_field(data_array, current_species, directions)
  END SUBROUTINE calc_viscous_electric_field_xzz


  SUBROUTINE calc_viscous_electric_field_yzy(data_array, current_species)
    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species
    INTEGER, DIMENSION(3) :: directions

    directions = (/ c_dir_y, c_dir_z, c_dir_y /)
    CALL calc_per_species_viscous_electric_field(data_array, current_species, directions)
  END SUBROUTINE calc_viscous_electric_field_yzy


  SUBROUTINE calc_viscous_electric_field_yzz(data_array, current_species)
    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species
    INTEGER, DIMENSION(3) :: directions

    directions = (/ c_dir_y, c_dir_z, c_dir_z /)
    CALL calc_per_species_viscous_electric_field(data_array, current_species, directions)
  END SUBROUTINE calc_viscous_electric_field_yzz

  SUBROUTINE calc_derivate(data_array, direction)
    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: direction
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: data_array_derivate
    REAL(num) :: idx_interval
    INTEGER :: ix, iy

    ALLOCATE(data_array_derivate(1-ng:nx+ng, 1-ng:ny+ng))
    data_array_derivate = 0.0_num
    SELECT CASE(direction)
    CASE(c_dir_x)
        ! Centre data
        DO iy = 1, ny
        DO ix = 2, nx-1
            data_array_derivate(ix, iy) = 0.5_num * (data_array(ix+1, iy) - data_array(ix-1, iy))
        END DO
        END DO
        ! Edge data
        DO iy = 1, ny
            data_array_derivate(1, iy) = 0.5_num*(data_array_derivate(2, iy) + data_array_derivate(3, iy))
            data_array_derivate(nx, iy) = 0.5_num*(data_array_derivate(nx-1, iy) + data_array_derivate(nx-2, iy))
        !    data_array_derivate(1-ng, iy) = -0.5_num * &
        !        (3.0_num*data_array(1-ng, iy) - 4.0_num*data_array(2-ng, iy) + data_array(3-ng, iy))
        !    data_array_derivate(nx+ng, iy) = 0.5_num * &
        !        (3.0_num*data_array(nx+ng, iy) - 4.0_num*data_array(nx+ng-1, iy) + data_array(nx+ng-2, iy))
        END DO

        idx_interval = 1.0_num / dx
        data_array = data_array_derivate * idx_interval
    CASE(c_dir_y)
        ! Centre data
        DO iy = 2, ny-1
        DO ix = 1, nx
            data_array_derivate(ix, iy) = 0.5_num * (data_array(ix, iy+1) - data_array(ix, iy-1))
        END DO
        END DO
        ! Edge data
        DO ix = 1, nx
            data_array_derivate(ix, 1) = 0.5_num*(data_array_derivate(ix, 2) + data_array_derivate(ix, 3))
            data_array_derivate(ix, ny) = 0.5_num*(data_array_derivate(ix, ny-1) + data_array_derivate(ix, ny-2))
            !data_array_derivate(ix, 1) = -0.5_num * &
            !    (3.0_num*data_array(ix, 1) - 4.0_num*data_array(ix, 2) + data_array(ix, 3))
            !data_array_derivate(ix, ny) = 0.5_num * &
            !    (3.0_num*data_array(ix, ny) - 4.0_num*data_array(ix, ny-1) + data_array(ix, ny-2))
        END DO

        idx_interval = 1.0_num / dy
        data_array = data_array_derivate * idx_interval
    CASE(c_dir_z)
        data_array = 0.0_num

    END SELECT

  END SUBROUTINE calc_derivate



END MODULE calc_df
