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

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: manual_load

CONTAINS

  SUBROUTINE manual_load
      TYPE(particle_species), POINTER :: species
      TYPE(particle), POINTER :: current
      REAL(num) :: rpos, dx,g1,g2,gam, P,theta,phi,px,py,pz,mass,n, energy_pho
      INTEGER :: ispecies
      energy_pho = ml_energy_pho

      IF (rank == 0) THEN
        WRITE(*,*) 'manual_load: for photon'
      END IF
      DO ispecies = 1, n_species
        IF (species_list(ispecies)%species_type == c_species_id_photon) THEN
        current => species_list(ispecies)%attached_list%head
!        mass = species_list(ispecies)%mass
            DO WHILE(ASSOCIATED(current))
!            P = random()
!            gam = g1/(((g1/g2)**(n-1)-1)*P + 1)**(1/(n-1))
!!            !for uniform solid angle
!            theta = 2*pi*random() ! dOmega is uniform 
!            phi = ACOS(2*random() - 1)
!            px = sqrt(gam**2 - 1)*SIN(theta)*SIN(phi)
!            py = sqrt(gam**2 - 1)*COS(theta)*SIN(phi)
!            pz = sqrt(gam**2 - 1)*COS(phi)
!
            ! for unifrom projection angle
!            theta = 2*pi*random() ! <p,B> is uniform
!            phi = 2*pi*random()
!            px = sqrt(gam**2 - 1)*SIN(theta)*COS(phi)
!            py = sqrt(gam**2 - 1)*COS(theta)
!            pz = sqrt(gam**2 - 1)*SIN(theta)*SIN(phi)
!
            current%part_p(1) = energy_pho/c
            current%particle_energy = energy_pho
!            current%part_p(2) = py*mass*c
!            current%part_p(3) = pz*mass*c
            current => current%next
            ENDDO ! do while current
        ENDIF !if photon
     ENDDO !ispecies

  END SUBROUTINE manual_load

END MODULE ic_module
