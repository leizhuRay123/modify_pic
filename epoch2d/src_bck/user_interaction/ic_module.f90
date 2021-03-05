! Copyright (C) 2010-2014 Keith Bennett <K.Bennett@warwick.ac.uk>
! Copyright (C) 2009      Chris Brady <C.S.Brady@warwick.ac.uk>
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

MODULE ic_module

  USE shared_data
  USE helper
  USE random_generator

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: manual_load_delta_f,manual_load_wein

CONTAINS

  SUBROUTINE manual_load_delta_f
      TYPE(particle_species), POINTER :: species
      TYPE(particle), POINTER :: current
      REAL(num) :: rpos, dx,g1,g2,gam, P,theta,phi,px,py,pz,mass,n, energy_pho,p_mag,stdev,mu,en_kbT,energy0
      REAL(num), PARAMETER :: eps = EPSILON(1.0_num)
      INTEGER :: ispecies
!      energy_pho = ml_energy_pho

      energy0 = delta_f_temp*q0
      stdev = DBLE(delta_f_stdev)
      mu = delta_f_mu

      IF (rank == 0) THEN
        WRITE(*,*) '####-----manual_load_delta_f for photon---#####'
        WRITE(*,*) 'energy0=',delta_f_temp,'*q0'
        WRITE(*,*) 'stdev=',delta_f_stdev
        WRITE(*,*) 'mu=',delta_f_mu
      END IF

      DO ispecies = 1, n_species
        IF (species_list(ispecies)%species_type == c_species_id_photon) THEN
        current => species_list(ispecies)%attached_list%head
!        mass = species_list(ispecies)%mass
            DO WHILE(ASSOCIATED(current))
                en_kbT = random_box_muller(stdev,mu)
                IF (en_kbT < 1e-3_num) THEN 
                    en_kbT = 1e-3_num
                END IF
                energy_pho = en_kbT*energy0;

                p_mag = energy_pho/c
            !P = random()
!            gam = g1/(((g1/g2)**(n-1)-1)*P + 1)**(1/(n-1))
!!            !for uniform solid angle
            theta = 2*pi*random() ! dOmega is uniform 
            phi = ACOS(2*random() - 1)
            px = p_mag*SIN(theta)*SIN(phi)
            py = p_mag*COS(theta)*SIN(phi)
            pz = p_mag*COS(phi)
!
            ! for unifrom projection angle
!            theta = 2*pi*random() ! <p,B> is uniform
!            phi = 2*pi*random()
!            px = sqrt(gam**2 - 1)*SIN(theta)*COS(phi)
!            py = sqrt(gam**2 - 1)*COS(theta)
!            pz = sqrt(gam**2 - 1)*SIN(theta)*SIN(phi)
!
            current%part_p = (/px,py,pz/)
            current%particle_energy = energy_pho
!            current%part_p(2) = py*mass*c
!            current%part_p(3) = pz*mass*c
            current => current%next
            ENDDO ! do while current
        ENDIF !if photon
     ENDDO !ispecies
     !WRITE(*,*) 'ipart',ipart
  END SUBROUTINE manual_load_delta_f

  SUBROUTINE manual_load_wein
      TYPE(particle_species), POINTER :: species
      TYPE(particle), POINTER :: current
      REAL(num) :: rpos, dx,g1,g2,gam, P,theta,phi,px,py,pz,mass,n, energy_pho,p_mag,stdev,mu,en_kbT,energy0
      REAL(num), PARAMETER :: eps = EPSILON(1.0_num)
      INTEGER :: ispecies,ipart
!      energy_pho = ml_energy_pho
      energy0 = wein_temperature*q0

      DO ispecies = 1, n_species
        IF (species_list(ispecies)%species_type == c_species_id_photon) THEN
        IF (rank == 0) THEN
            WRITE(*,*) 'manual_load_wein: for photon'
        END IF
        current => species_list(ispecies)%attached_list%head
!        mass = species_list(ispecies)%mass
            DO WHILE(ASSOCIATED(current))
                P = random()
                en_kbT = find_value_from_table_1d_ic(P, 1000, wein_table, log_wein_xi)
                energy_pho = en_kbT*energy0;
                p_mag = energy_pho/c
!!            !for uniform solid angle
                theta = 2*pi*random() ! dOmega is uniform 
                phi = ACOS(2*random() - 1)
                px = p_mag*SIN(theta)*SIN(phi)
                py = p_mag*COS(theta)*SIN(phi)
                pz = p_mag*COS(phi)
                current%part_p = (/px,py,pz/)
                current%particle_energy = energy_pho
                current => current%next
            ENDDO ! do while current
        ENDIF !if photon
     ENDDO !ispecies
     !WRITE(*,*) 'ipart',ipart

  END SUBROUTINE manual_load_wein

  FUNCTION find_value_from_table_1d_ic(x_in, nx, x, values)

    REAL(num) :: find_value_from_table_1d_ic
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

    find_value_from_table_1d_ic = 10.0_num**value_interp

  END FUNCTION find_value_from_table_1d_ic


END MODULE ic_module
