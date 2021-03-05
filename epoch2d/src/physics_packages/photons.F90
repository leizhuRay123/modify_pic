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

MODULE photons
#ifdef PHOTONS

  USE partlist

  IMPLICIT NONE

  SAVE
  REAL(num) :: part_pos_global, gamma_global, eta_global
CONTAINS

  SUBROUTINE setup_qed_module

    INTEGER :: ispecies, iu, io
    LOGICAL :: found

    ! Sanity check
    !IF (produce_photons .AND. .NOT.ic_from_restart) THEN
    IF (.NOT.ic_from_restart) THEN
      ! Identify if there exists any *populated* electron/positron species
      found = .FALSE.
      DO ispecies = 1, n_species
        IF (species_list(ispecies)%count <= 0) CYCLE

        IF (species_list(ispecies)%species_type == c_species_id_electron &
            .OR. species_list(ispecies)%species_type == c_species_id_positron &
            .OR. ispecies == photon_species ) THEN
          found = .TRUE.
          EXIT
        END IF
      END DO

      IF (.NOT.found) THEN
        IF (rank == 0) THEN
          DO iu = 1, nio_units ! Print to stdout and to file
            io = io_units(iu)
            WRITE(io,*)
            WRITE(io,*) '*** WARNING ***'
            WRITE(io,*) 'Electron, positron and photon species are either ', &
                'unspecified or contain no'
            WRITE(io,*) 'particles. QED routines will do nothing.'
          END DO
        END IF
      END IF
    END IF

    ! Load the tables for the QED routines
    CALL setup_tables_qed
    IF (.NOT.ic_from_restart) THEN
      DO ispecies = 1, n_species
        CALL initialise_optical_depth(species_list(ispecies))
      END DO
    END IF
  END SUBROUTINE setup_qed_module



  SUBROUTINE shutdown_qed_module

    CALL deallocate_tables_qed

  END SUBROUTINE shutdown_qed_module



  FUNCTION check_qed_variables()

    INTEGER :: check_qed_variables
    INTEGER :: io, iu, ispecies
    INTEGER :: first_electron = -1, first_positron = -1

    check_qed_variables = c_err_none

    IF (.NOT.use_qed) RETURN

    ! If you're only doing radiation reaction force then don't need any special
    ! species, so don't do any checking here
    !XY close produce_photon limit
    !IF (.NOT.produce_photons) RETURN

    ! Identify if there exists any electron species
    DO ispecies = 1, n_species
      IF (species_list(ispecies)%species_type == c_species_id_electron &
          .AND. first_electron == -1) THEN
        first_electron = ispecies
      ELSE IF (species_list(ispecies)%species_type == c_species_id_positron &
          .AND. first_positron == -1) THEN
        first_positron = ispecies
      END IF
    END DO

    IF (first_electron < 0 .AND. first_positron < 0) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'No electron or positron species specified.'
          WRITE(io,*) 'Specify using "identify:electron" or "identify:positron"'
          WRITE(io,*) 'QED routines require at least one species of ', &
              'electrons or positrons.'
        END DO
      END IF
      check_qed_variables = c_err_missing_elements
      RETURN
    END IF

    IF (photon_species < 0) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'No photon species specified. Specify using ', &
              '"identify:photon"'
        END DO
      END IF
      check_qed_variables = c_err_missing_elements
      RETURN
    END IF

    ! If you're not producing pairs then you don't have to designate special
    ! electron or positron species so just return
    IF (.NOT.produce_pairs) RETURN

    IF (first_electron < 0 .OR. first_positron < 0) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'To use pair production routines need at least one ', &
              'positron species and one ', 'electron species. Specify ', &
              'using "identify:electron" or "identify:positron"'
        END DO
      END IF
      check_qed_variables = c_err_missing_elements
      RETURN
    END IF

    IF (breit_wheeler_positron_species < 0) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** WARNING ***'
          WRITE(io,*) 'No Breit-Wheeler positron species specified.'
          WRITE(io,*) 'Specify using "identify:breit_wheeler_positron".'
          WRITE(io,*) 'Using species ', &
              TRIM(species_list(first_positron)%name), ' instead.'
        END DO
      END IF
      breit_wheeler_positron_species = first_positron
    END IF

#ifdef TRIDENT_PHOTONS
    IF (trident_positron_species < 0) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** WARNING ***'
          WRITE(io,*) 'No trident positron species specified.'
          WRITE(io,*) 'Specify using "identify:trident_positron".'
          WRITE(io,*) 'Using species ', &
              TRIM(species_list(first_positron)%name), ' instead.'
        END DO
      END IF
      trident_positron_species = first_positron
    END IF
#endif

    IF (breit_wheeler_electron_species < 0) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** WARNING ***'
          WRITE(io,*) 'No Breit-Wheeler electron species specified.'
          WRITE(io,*) 'Specify using "identify:breit_wheeler_electron".'
          WRITE(io,*) 'Using species ', &
              TRIM(species_list(first_electron)%name), ' instead.'
        END DO
      END IF
      breit_wheeler_electron_species = first_electron
    END IF

#ifdef TRIDENT_PHOTONS
    IF (trident_electron_species < 0) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** WARNING ***'
          WRITE(io,*) 'No trident electron species specified.'
          WRITE(io,*) 'Specify using "identify:trident_electron".'
          WRITE(io,*) 'Using species ', &
              TRIM(species_list(first_electron)%name), ' instead.'
        END DO
      END IF
      trident_electron_species = first_electron
    END IF
#endif
  END FUNCTION check_qed_variables



  SUBROUTINE setup_tables_qed

    ! Reads files epsilon.table, log_chi.table, energy_split.table
    ! and sets up appropriate tables

    REAL(num) :: etalog_min = 0.0_num, etalog_max = 0.0_num
    REAL(num) :: etalog_dx, chi_min, chi_dx
    REAL(num), ALLOCATABLE :: realbuf(:)
    INTEGER :: i, n, ichi2, iepsilon, ieta, ichi, bufsize
    !XY 
    INTEGER :: ixi, intbuf(14),iy,igp,ix

    IF (rank == 0) THEN
      OPEN(unit=lu, file=TRIM(qed_table_location)//'/hsokolov.table', &
          status='OLD')
      READ(lu,*) n_sample_h
      ALLOCATE(log_hsokolov(n_sample_h,2))
      DO ieta = 1, n_sample_h
        READ(lu,*) log_hsokolov(ieta,1), log_hsokolov(ieta,2)
      END DO
      CLOSE(unit=lu)

      OPEN(unit=lu, file=TRIM(qed_table_location)//'/pairprod.table', &
          status='OLD')
      READ(lu,*) n_sample_t
      ALLOCATE(log_tpair(n_sample_t,2))
      ALLOCATE(log_omegahat(n_sample_t,2))
      DO ichi = 1, n_sample_t
        READ(lu,*) log_tpair(ichi,1), log_omegahat(ichi,2), log_tpair(ichi,2)
      END DO
      CLOSE(unit=lu)

      OPEN(unit=lu, file=TRIM(qed_table_location)//'/ksi_sokolov.table', &
          status='OLD')
      READ(lu,*) n_sample_eta, n_sample_chi, etalog_min, etalog_max
      ALLOCATE(p_photon_energy(n_sample_eta,n_sample_chi))
      DO ieta = 1, n_sample_eta
        READ(lu,*) (p_photon_energy(ieta,ichi), ichi=1,n_sample_chi)
      END DO
      CLOSE(unit=lu)

      OPEN(unit=lu, file=TRIM(qed_table_location)//'/chimin.table', &
          status='OLD')
      ALLOCATE(chimin_table(n_sample_eta))
      READ(lu,*) (chimin_table(ieta), ieta=1,n_sample_eta)
      CLOSE(unit=lu)

      OPEN(unit=lu, file=TRIM(qed_table_location)//'/log_chi2.table', &
          status='OLD')
      READ(lu,*) n_sample_chi2
      ALLOCATE(log_chi2(n_sample_chi2))
      DO ichi2 = 1, n_sample_chi2
        READ(lu,*) log_chi2(ichi2)
      END DO
      CLOSE(unit=lu)

      OPEN(unit=lu, file=TRIM(qed_table_location)//'/epsilon.table', &
          status='OLD')
      READ(lu,*) n_sample_epsilon
      ALLOCATE(epsilon_split(n_sample_epsilon))
      DO iepsilon = 1, n_sample_epsilon
        READ(lu,*) epsilon_split(iepsilon)
      END DO
      CLOSE(unit=lu)

      OPEN(unit=lu, file=TRIM(qed_table_location)//'/energy_split.table', &
          status='OLD')
      ALLOCATE(p_energy(n_sample_chi2,n_sample_epsilon))
      DO ichi2 = 1, n_sample_chi2
        DO iepsilon = 1, n_sample_epsilon
          READ(lu,*) p_energy(ichi2,iepsilon)
        END DO
      END DO
      CLOSE(unit=lu)
      !XY for stokes
      ! for kv13/kv23
   OPEN(unit=lu, file=TRIM(qed_table_location)//'/kv13_div_kv23.table', &
       status='OLD')
   READ(lu,*) n_sample_xi
   ALLOCATE(kv13_div_kv23_table(n_sample_xi))
   DO ixi = 1, n_sample_xi
       READ(lu,*) kv13_div_kv23_table(ixi)
     ENDDO
   CLOSE(unit=lu)


     ! For kv(2/3,x)/Int(kv(5/3,x,infty))
     OPEN(unit=lu, file=TRIM(qed_table_location)//'/kv23_Intkv53.table', &
         status='OLD')
     READ(lu,*) n_sample_kv23divIntkv53
     ALLOCATE(kv23_div_Intkv53_table(n_sample_kv23divIntkv53))
     DO ixi = 1, n_sample_kv23divIntkv53
         READ(lu,*) kv23_div_Intkv53_table(ixi)
       ENDDO
     CLOSE(unit=lu)


     OPEN(unit=lu, file=TRIM(qed_table_location)//'/Psi23.table', &
         status='OLD')
     READ(lu,*) n_sample_y, n_sample_gp
     ALLOCATE(Psi23_table(n_sample_y,n_sample_gp))
     DO iy = 1, n_sample_y
        DO igp = 1, n_sample_gp
           READ(lu,*) Psi23_table(iy,igp)
       ENDDO
     ENDDO
     CLOSE(unit=lu)

     OPEN(unit=lu, file=TRIM(qed_table_location)//'/Psi13.table', &
         status='OLD')
     ALLOCATE(Psi13_table(n_sample_y,n_sample_gp))
     DO iy = 1, n_sample_y
       Do igp = 1, n_sample_gp
           READ(lu,*) Psi13_table(iy,igp)
       ENDDO
     ENDDO
     CLOSE(unit=lu)

     ! for xIntkv53
   OPEN(unit=lu, file=TRIM(qed_table_location)//'/xIntkv53.table', &
       status='OLD')
   READ(lu,*) n_sample_xIntkv53
   ALLOCATE(xIntkv53_table(n_sample_xIntkv53))
   DO ixi = 1, n_sample_xIntkv53
       READ(lu,*) xIntkv53_table(ixi)
   ENDDO
   CLOSE(unit=lu)


     ! for ComptonScatter
     OPEN(unit=lu, file=TRIM(qed_table_location)//'/table_for_compton.table', &
         status='OLD')
     READ(lu,*) n_sample_epo,n_sample_mu
     ALLOCATE(Compton_table(n_sample_epo,n_sample_mu))
     DO ix = 1, n_sample_epo
       Do iy = 1, n_sample_mu
           READ(lu,*) Compton_table(ix,iy)
       ENDDO
     ENDDO
   CLOSE(unit=lu)
     ! for photon wein distribution
     OPEN(unit=lu, file=TRIM(qed_table_location)//'/wein.table', &
         status='OLD')
     ALLOCATE(wein_table(1000))
     DO ix = 1, 1000
           READ(lu,*) wein_table(ix)
     ENDDO
   CLOSE(unit=lu)

     ! for xkv13
!   OPEN(unit=lu, file=TRIM(qed_table_location)//'/x_kv13.table', &
!       status='OLD')
!   READ(lu,*) n_sample_xkv13
!   ALLOCATE(xkv13_table(n_sample_xkv13))
!   DO ixi = 1, n_sample_xkv13
!       READ(lu,*) xkv13_table(ixi)
!   ENDDO
!   CLOSE(unit=lu)

      !XY 
       WRITE(*,*) 'n_sample_xi='  ,n_sample_xi
       WRITE(*,*) 'n_sample_y='   ,n_sample_y
       WRITE(*,*) 'n_sample_gp='  ,n_sample_gp
       WRITE(*,*) 'n_sample_xIntkv53=',n_sample_xIntkv53
       WRITE(*,*) 'n_sample_kv23_div_Intkv53=',n_sample_kv23divIntkv53
!       WRITE(*,*) 'n_sample_xkv13=',n_sample_xkv13
       WRITE(*,*) 'n_sample_epo=%d, n_sample_mu=%d',n_sample_epo,n_sample_mu

      ! Pack data for broadcasting to other processes

      intbuf(1) = n_sample_h
      intbuf(2) = n_sample_t
      intbuf(3) = n_sample_eta
      intbuf(4) = n_sample_chi
      intbuf(5) = n_sample_chi2
      intbuf(6) = n_sample_epsilon
        !XY
     intbuf(7) = n_sample_xi
     intbuf(8) = n_sample_y
     intbuf(9) = n_sample_gp
     intbuf(10) = n_sample_xIntkv53
     intbuf(11) = n_sample_kv23divIntkv53
     intbuf(12) = n_sample_epo
     intbuf(13) = n_sample_mu
!     intbuf(14) = n_sample_xkv13

      CALL MPI_BCAST(intbuf, 13, MPI_INTEGER, 0, comm, errcode)

      bufsize = n_sample_h * 2
      bufsize = bufsize + n_sample_t * 3
      bufsize = bufsize + 2 + n_sample_eta * n_sample_chi
      bufsize = bufsize + n_sample_eta
      bufsize = bufsize + n_sample_chi2
      bufsize = bufsize + n_sample_epsilon
      bufsize = bufsize + n_sample_chi2 * n_sample_epsilon
       !For polarization
      bufsize = bufsize + n_sample_xi
      bufsize = bufsize + n_sample_kv23divIntkv53 !for kv(2/3,x)/Int(kv(5/3,x,infty))
      bufsize = bufsize + n_sample_y*n_sample_gp
      bufsize = bufsize + n_sample_y*n_sample_gp
      bufsize = bufsize + n_sample_xIntkv53
      ! For compton scatter
      bufsize = bufsize + n_sample_epo*n_sample_mu
      ! For wein table
      bufsize = bufsize + 1000
      ! For xkv13 for synchrotron circular polarization
 !     bufsize = bufsize + n_sample_xkv13

      ALLOCATE(realbuf(bufsize))
      n = 1

      DO i = 1, 2
        DO ieta = 1, n_sample_h
          realbuf(n) = log_hsokolov(ieta,i)
          n = n + 1
        END DO
      END DO

      DO ichi = 1, n_sample_t
        realbuf(n) = log_tpair(ichi,1)
        n = n + 1
        realbuf(n) = log_tpair(ichi,2)
        n = n + 1
        realbuf(n) = log_omegahat(ichi,2)
        n = n + 1
      END DO

      realbuf(n) = etalog_min
      n = n + 1
      realbuf(n) = etalog_max
      n = n + 1

      DO ichi = 1, n_sample_chi
        DO ieta = 1, n_sample_eta
          realbuf(n) = p_photon_energy(ieta,ichi)
          n = n + 1
        END DO
      END DO

      DO ieta = 1, n_sample_eta
        realbuf(n) = chimin_table(ieta)
        n = n + 1
      END DO

      DO ichi2 = 1, n_sample_chi2
        realbuf(n) = log_chi2(ichi2)
        n = n + 1
      END DO

      DO iepsilon = 1, n_sample_epsilon
        realbuf(n) = epsilon_split(iepsilon)
        n = n + 1
      END DO

      DO iepsilon = 1, n_sample_epsilon
        DO ichi2 = 1, n_sample_chi2
          realbuf(n) = p_energy(ichi2,iepsilon)
          n = n + 1
        END DO
      END DO

      ! XY for stokes
      DO ixi = 1,n_sample_xi
            realbuf(n) = kv13_div_kv23_table(ixi)
            n = n+1
        ENDDO


      DO ixi = 1,n_sample_kv23divIntkv53
            realbuf(n) = kv23_div_Intkv53_table(ixi)
            n = n+1
        ENDDO

   DO iy = 1, n_sample_y
    DO igp = 1, n_sample_gp
      realbuf(n) = Psi23_table(iy,igp)
      n = n + 1
    ENDDO
  ENDDO

   DO iy = 1, n_sample_y
    DO igp = 1, n_sample_gp
      realbuf(n) = Psi13_table(iy,igp)
      n = n + 1
    ENDDO
  ENDDO

      DO ixi = 1,n_sample_xIntkv53
            realbuf(n) = xIntkv53_table(ixi)
            n = n+1
        ENDDO

        !compton scatter
   DO ix = 1, n_sample_epo
    DO iy = 1, n_sample_mu
      realbuf(n) = Compton_table(ix,iy)
      n = n + 1
    ENDDO
  ENDDO
    ! wein table
   DO ix = 1, 1000
      realbuf(n) = wein_table(ix)
      n = n + 1
  ENDDO

!   DO ix = 1, n_sample_xkv13
!      realbuf(n) = xkv13_table(ix)
!      n = n + 1
!  ENDDO
!
      CALL MPI_BCAST(realbuf, bufsize, mpireal, 0, comm, errcode)

      DEALLOCATE(realbuf)
    ELSE
      ! Unpack data from rank zero

      CALL MPI_BCAST(intbuf, 13, MPI_INTEGER, 0, comm, errcode)

      n_sample_h       = intbuf(1)
      n_sample_t       = intbuf(2)
      n_sample_eta     = intbuf(3)
      n_sample_chi     = intbuf(4)
      n_sample_chi2    = intbuf(5)
      n_sample_epsilon = intbuf(6)
       ! XY
      n_sample_xi      = intbuf(7)
      n_sample_y      = intbuf(8)
      n_sample_gp      = intbuf(9)
      n_sample_xIntkv53      = intbuf(10)
      n_sample_kv23divIntkv53 = intbuf(11)
      n_sample_epo     = intbuf(12)
      n_sample_mu      = intbuf(13)
 !     n_sample_xkv13   = intbuf(14)

      bufsize = n_sample_h * 2
      bufsize = bufsize + n_sample_t * 3
      bufsize = bufsize + 2 + n_sample_eta * n_sample_chi
      bufsize = bufsize + n_sample_eta
      bufsize = bufsize + n_sample_chi2
      bufsize = bufsize + n_sample_epsilon
      bufsize = bufsize + n_sample_chi2 * n_sample_epsilon

      bufsize = bufsize + n_sample_xi
      bufsize = bufsize + n_sample_kv23divIntkv53
      bufsize = bufsize + n_sample_y*n_sample_gp ! For Psi23
      bufsize = bufsize + n_sample_y*n_sample_gp ! For Psi13
      bufsize = bufsize + n_sample_xIntkv53 ! For optical detection
      !compton scatter
      bufsize = bufsize + n_sample_epo*n_sample_mu
      !wein table
      bufsize = bufsize + 1000
      !read bufsize of n_sample_xkv13 
 !     bufsize = bufsize + n_sample_xkv13

      ALLOCATE(realbuf(bufsize))
      n = 1

      CALL MPI_BCAST(realbuf, bufsize, mpireal, 0, comm, errcode)

      ALLOCATE(log_hsokolov(n_sample_h,2))
      DO i = 1, 2
        DO ieta = 1, n_sample_h
          log_hsokolov(ieta,i) = realbuf(n)
          n = n + 1
        END DO
      END DO

      ALLOCATE(log_tpair(n_sample_t,2))
      ALLOCATE(log_omegahat(n_sample_t,2))
      DO ichi = 1, n_sample_t
        log_tpair(ichi,1)    = realbuf(n)
        n = n + 1
        log_tpair(ichi,2)    = realbuf(n)
        n = n + 1
        log_omegahat(ichi,2) = realbuf(n)
        n = n + 1
      END DO

      etalog_min = realbuf(n)
      n = n + 1
      etalog_max = realbuf(n)
      n = n + 1

      ALLOCATE(p_photon_energy(n_sample_eta,n_sample_chi))
      DO ichi = 1, n_sample_chi
        DO ieta = 1, n_sample_eta
          p_photon_energy(ieta,ichi) = realbuf(n)
          n = n + 1
        END DO
      END DO

      ALLOCATE(chimin_table(n_sample_eta))
      DO ieta = 1, n_sample_eta
        chimin_table(ieta) = realbuf(n)
        n = n + 1
      END DO

      ALLOCATE(log_chi2(n_sample_chi2))
      DO ichi2 = 1, n_sample_chi2
        log_chi2(ichi2) = realbuf(n)
        n = n + 1
      END DO

      ALLOCATE(epsilon_split(n_sample_epsilon))
      DO iepsilon = 1, n_sample_epsilon
        epsilon_split(iepsilon) = realbuf(n)
        n = n + 1
      END DO

      ALLOCATE(p_energy(n_sample_chi2,n_sample_epsilon))
      DO iepsilon = 1, n_sample_epsilon
        DO ichi2 = 1, n_sample_chi2
          p_energy(ichi2,iepsilon) = realbuf(n)
          n = n + 1
        END DO
      END DO
      !XY for stokes
      ALLOCATE(kv13_div_kv23_table(n_sample_xi))
       DO ixi = 1, n_sample_xi
           kv13_div_kv23_table(ixi) = realbuf(n)
         n = n + 1
       ENDDO

     ALLOCATE(kv23_div_Intkv53_table(n_sample_kv23divIntkv53))
       DO ixi = 1, n_sample_kv23divIntkv53
           kv23_div_Intkv53_table(ixi) = realbuf(n)
         n = n + 1
       ENDDO

     ALLOCATE(Psi23_table(n_sample_y,n_sample_gp))
      DO iy = 1, n_sample_y
       DO igp = 1, n_sample_gp
         Psi23_table(iy,igp) = realbuf(n)
         n = n + 1
       ENDDO
     ENDDO

     ALLOCATE(Psi13_table(n_sample_y,n_sample_gp))
      DO iy = 1, n_sample_y
       DO igp = 1, n_sample_gp
         Psi13_table(iy,igp) = realbuf(n)
         n = n + 1
       ENDDO
     ENDDO

      ALLOCATE(xIntkv53_table(n_sample_xIntkv53))
       DO ixi = 1, n_sample_xIntkv53
           xIntkv53_table(ixi) = realbuf(n)
         n = n + 1
       ENDDO

!Compton scatter
     ALLOCATE(Compton_table(n_sample_epo,n_sample_mu))
      DO ix = 1, n_sample_epo
       DO iy = 1, n_sample_mu
         Compton_table(ix,iy) = realbuf(n)
         n = n + 1
       ENDDO
     ENDDO
     !wein table
     ALLOCATE(wein_table(1000))
      DO ix = 1, 1000
         wein_table(ix) = realbuf(n)
         n = n + 1
     ENDDO

!     ALLOCATE(xkv13_table(n_sample_xkv13))
!      DO ix = 1, n_sample_xkv13
!         xkv13_table(ix) = realbuf(n)
!         n = n + 1
!     ENDDO
!
      DEALLOCATE(realbuf)
    END IF

    log_omegahat(:,1) = log_tpair(:,1)

    ALLOCATE(log_eta(n_sample_eta))
    ALLOCATE(log_chi(n_sample_eta,n_sample_chi))
    !XY
    ALLOCATE(log_xi(n_sample_xi))
    ALLOCATE(log_y(n_sample_y))
    ALLOCATE(gp(n_sample_gp))
    ALLOCATE(log_xIntkv53(n_sample_xIntkv53))
    !For Compton scatter
    ALLOCATE(log_epo(n_sample_epo))
    ALLOCATE(com_mu(n_sample_epo,n_sample_mu))
    !For wein 
    ALLOCATE(log_wein_xi(1000))
    

    etalog_dx = (etalog_max - etalog_min) / REAL(n_sample_eta-1,num)
    DO ieta = 1, n_sample_eta
      log_eta(ieta) = etalog_min + REAL(ieta-1,num) * etalog_dx
      chi_min = LOG10(chimin_table(ieta))
      chi_dx  = (log_eta(ieta) - LOG10(2.0_num) - chi_min) &
          / REAL(n_sample_chi-1,num)
      DO ichi = 1, n_sample_chi
        log_chi(ieta,ichi) = chi_min + REAL(ichi-1,num) * chi_dx
      END DO
    END DO

      !XY
      log_xi_min = -5
      log_xi_max = 2.5
      log_xi_dx = (log_xi_max - log_xi_min)/(REAL(n_sample_xi,num) - 1)
      DO ixi = 1, n_sample_xi
          log_xi(ixi) = log_xi_min + REAL(ixi-1,num) * log_xi_dx
      ENDDO
      IF(rank == 0) THEN
          WRITE(*,*) 'xi:', n_sample_xi, log_xi(1), log_xi(n_sample_xi)
      END IF
 
      !For gp
      gp_min = -10
      gp_max = 10
      dgp = (gp_max - gp_min)/ REAL(n_sample_gp - 1,num)
      Do igp =1, n_sample_gp
          gp(igp) =  gp_min + REAL(igp - 1, num) * dgp
      ENDDO
      !For y = o_oc
      log_y_min = -4
      log_y_max = 2.0
      log_dy = (log_y_max - log_y_min)/ REAL(n_sample_y - 1, num)
 
      Do iy =1, n_sample_y
          log_y(iy) =  log_y_min + REAL(iy - 1, num) * log_dy
      ENDDO
      !For xIntkv53
      log_xInt_min = -9.0
      log_xInt_max = 1.0
      log_dxInt = (log_xInt_max - log_xInt_min)/ REAL(n_sample_xIntkv53 - 1, num)
      Do iy = 1, n_sample_xIntkv53
          log_xIntkv53(iy) =  log_xInt_min + REAL(iy - 1, num) * log_dxInt
      ENDDO
      !For kv23_div_intkv53
      !log_xIntkv53 is the same tas the log_kv23_div_intkv53
      !For log_epo
      log_epo_min = -3.0
      log_epo_max = 3.0 !1e3*Kev/me/c**2
      log_depo = (log_epo_max - log_epo_min)/ REAL(n_sample_epo - 1, num)
      Do ix = 1, n_sample_epo
          log_epo(ix) =  log_epo_min + REAL(ix - 1, num) * log_depo
      ENDDO

      !For compton theta
      com_theta_min = 0 + 1e-9
      com_theta_max = pi - 1e-9
      dcom_theta = (com_theta_max - com_theta_min)/ REAL(n_sample_mu - 1, num)
      Do ix = 1, n_sample_epo
      Do iy = 1, n_sample_mu
          com_theta =  com_theta_min + REAL(iy - 1, num) * dcom_theta
          com_mu(ix,iy) = COS(com_theta)
      ENDDO
      ENDDO
      !WRITE(*,*) com_mu(401,68)
      !wein 
      log_wein_xi_min = -1.5
      log_wein_xi_max = 1.5
      dlog_wein_xi = (log_wein_xi_max - log_wein_xi_min)/ REAL(1000 - 1, num)
      Do ix = 1, 1000
          log_wein_xi(ix) =  log_wein_xi_min + REAL(ix - 1, num) * dlog_wein_xi
      ENDDO
     
  END SUBROUTINE setup_tables_qed



  SUBROUTINE deallocate_tables_qed

    DEALLOCATE(log_chi2)
    DEALLOCATE(epsilon_split)
    DEALLOCATE(p_energy)
    DEALLOCATE(log_hsokolov)
    DEALLOCATE(log_eta)
    DEALLOCATE(log_chi)
    DEALLOCATE(p_photon_energy)
    DEALLOCATE(log_tpair)
    DEALLOCATE(log_omegahat)

          !XY
    DEALLOCATE(kv13_div_kv23_table)
    DEALLOCATE(kv23_div_Intkv53_table)
    DEALLOCATE(Psi23_table)
    DEALLOCATE(Psi13_table)
    DEALLOCATE(xIntkv53_table)
    DEALLOCATE(log_y)
    DEALLOCATE(gp)
    DEALLOCATE(log_xi)
    DEALLOCATE(log_xIntkv53)
    DEALLOCATE(chimin_table)
    !compton and wein
    DEALLOCATE(Compton_table)
    DEALLOCATE(log_epo)
    DEALLOCATE(com_mu)
    DEALLOCATE(wein_table)
    DEALLOCATE(log_wein_xi)

  END SUBROUTINE deallocate_tables_qed



  SUBROUTINE initialise_optical_depth(current_species)

    ! Resets optical depth (to random number) of all particles
    TYPE(particle_species) :: current_species
    TYPE(particle), POINTER :: current
    REAL(num) :: p_tau

    current => current_species%attached_list%head
    DO WHILE(ASSOCIATED(current))
      p_tau = random()
      current%optical_depth = -LOG(1.0_num - p_tau)

#ifdef TRIDENT_PHOTONS
      p_tau = random()
      current%optical_depth_tri = -LOG(1.0_num - p_tau)
#endif
      current => current%next
    END DO

  END SUBROUTINE initialise_optical_depth



  FUNCTION reset_optical_depth()

    ! Resets optical depth of particle
    REAL(num) :: reset_optical_depth
    REAL(num) :: p_tau

    p_tau = random()
    reset_optical_depth = -LOG(1.0_num - p_tau)

  END FUNCTION reset_optical_depth



  SUBROUTINE qed_update_optical_depth

    ! Updates the optical depth for electrons and photons
    INTEGER :: ispecies
    TYPE(particle), POINTER :: current, next_pt

    REAL(num) :: part_x, part_y
    REAL(num) :: part_ux, part_uy, part_uz
    REAL(num) :: dir_x, dir_y, dir_z
    REAL(num) :: eta, chi_val, part_e, gamma_rel, norm

    DO ispecies = 1, n_species

      ! First consider electrons and positrons
      IF (species_list(ispecies)%species_type == c_species_id_electron &
          .OR. species_list(ispecies)%species_type == c_species_id_positron) &
          THEN
        current => species_list(ispecies)%attached_list%head
        DO WHILE(ASSOCIATED(current))
          ! Find eta at particle position
          part_x  = current%part_pos(1) - x_grid_min_local
          part_y  = current%part_pos(2) - y_grid_min_local
          part_ux = current%part_p(1) / mc0
          part_uy = current%part_p(2) / mc0
          part_uz = current%part_p(3) / mc0
          gamma_rel = SQRT(part_ux**2 + part_uy**2 + part_uz**2 + 1.0_num)

          eta = calculate_eta(part_x, part_y, part_ux, part_uy, &
              part_uz, gamma_rel)

          current%optical_depth = &
              current%optical_depth - gfactor*delta_optical_depth(eta, gamma_rel)
#ifdef TRIDENT_PHOTONS
          current%optical_depth_tri = current%optical_depth_tri &
              - gfactor*delta_optical_depth_tri(eta, gamma_rel)
#endif
          ! If optical depth dropped below zero generate photon...
          IF (current%optical_depth <= 0.0_num) THEN
            CALL generate_photon(current, photon_species, eta)
            ! ... and reset optical depth
            current%optical_depth = reset_optical_depth()
          END IF

#ifdef TRIDENT_PHOTONS
          IF (current%optical_depth_tri <= 0.0_num) THEN
            CALL generate_pair_tri(current, trident_electron_species, &
                trident_positron_species)
            ! ... and reset optical depth
            current%optical_depth_tri = reset_optical_depth()
          END IF
#endif
          current => current%next
        END DO

      ! and finally photons
      ELSE IF (species_list(ispecies)%species_type == c_species_id_photon &
          .AND. produce_pairs) THEN
        current => species_list(ispecies)%attached_list%head
        DO WHILE(ASSOCIATED(current))
          ! Current may be deleted
          next_pt => current%next
          part_x  = current%part_pos(1) - x_grid_min_local
          part_y  = current%part_pos(2) - y_grid_min_local
          norm  = c / current%particle_energy
          dir_x = current%part_p(1) * norm
          dir_y = current%part_p(2) * norm
          dir_z = current%part_p(3) * norm
          part_e  = current%particle_energy / m0 / c**2

          chi_val = calculate_chi(part_x, part_y, dir_x, dir_y, &
              dir_z, part_e)

          current%optical_depth = current%optical_depth &
              - gfactor*delta_optical_depth_photon(chi_val, part_e)
          ! If optical depth dropped below zero generate pair...
          IF (current%optical_depth <= 0.0_num) THEN
            CALL generate_pair(current, chi_val, photon_species, &
                breit_wheeler_electron_species, breit_wheeler_positron_species)
          END IF
          current => next_pt
        END DO
      END IF
    END DO

  END SUBROUTINE qed_update_optical_depth



  FUNCTION delta_optical_depth(eta, gamma_rel)

    ! Function that calcualtes the change to the optical depth
    REAL(num) :: delta_optical_depth
    REAL(num), INTENT(IN) :: eta, gamma_rel
    REAL(num) :: hsokolov


    hsokolov = find_value_from_table_1d(eta, n_sample_h, log_hsokolov(:,1), &
        log_hsokolov(:,2))

    delta_optical_depth = dt * eta * alpha_f * SQRT(3.0_num) * hsokolov &
        / (2.0_num * pi * tau_c * gamma_rel)

  END FUNCTION delta_optical_depth



  FUNCTION delta_optical_depth_tri(eta, gamma_rel)

    ! Function that calcualtes the change to the optical depth
    REAL(num) :: delta_optical_depth_tri
    REAL(num), INTENT(IN) :: eta, gamma_rel
    REAL(num) :: omegahat

    omegahat = find_value_from_table_1d(eta, n_sample_t, log_omegahat(:,1), &
        log_omegahat(:,2))

    delta_optical_depth_tri = dt * eta * alpha_f**2 * 0.64_num * omegahat &
        / (2.0_num * pi * tau_c * gamma_rel)

  END FUNCTION delta_optical_depth_tri



  FUNCTION delta_optical_depth_photon(chi_val, part_e)

    ! Function that calcualtes the change to the optical depth
    REAL(num) :: delta_optical_depth_photon
    REAL(num), INTENT(IN) :: chi_val, part_e
    REAL(num) :: tpair

    tpair = find_value_from_table_1d(chi_val, n_sample_t, log_tpair(:,1), &
        log_tpair(:,2))

    delta_optical_depth_photon = dt / tau_c * alpha_f / part_e * chi_val * tpair

  END FUNCTION delta_optical_depth_photon



  FUNCTION calculate_eta(part_x, part_y, part_ux, part_uy, part_uz, &
      gamma_rel)

    REAL(num) :: calculate_eta
    REAL(num), INTENT(IN) :: part_x, part_y
    REAL(num), INTENT(IN) :: part_ux, part_uy, part_uz, gamma_rel
    REAL(num) :: e_at_part(3), b_at_part(3)
    REAL(num) :: beta_x, beta_y, beta_z
    REAL(num) :: flperp(3), i_e, tau0, roland_eta
    REAL(num) :: lambdac, coeff_eta, moduclip, moduclip2, u_dot_e

    CALL field_at_particle(part_x, part_y, e_at_part, b_at_part)

    moduclip2 = MAX(part_ux**2 + part_uy**2 + part_uz**2, c_tiny)
    moduclip = SQRT(moduclip2)

    beta_x = part_ux / gamma_rel
    beta_y = part_uy / gamma_rel
    beta_z = part_uz / gamma_rel

    lambdac = h_bar / mc0

    coeff_eta = SQRT(3.0_num * lambdac / (2.0_num * alpha_f * m0 * c**3))

    u_dot_e = (part_ux * e_at_part(1) + part_uy * e_at_part(2) &
        + part_uz * e_at_part(3)) / moduclip2

    flperp(1) = q0 * (e_at_part(1) - u_dot_e * part_ux &
        + c * (beta_y * b_at_part(3) - beta_z * b_at_part(2)))

    flperp(2) = q0 * (e_at_part(2) - u_dot_e * part_uy &
        + c * (beta_z * b_at_part(1) - beta_x * b_at_part(3)))

    flperp(3) = q0 * (e_at_part(3) - u_dot_e * part_uz &
        + c * (beta_x * b_at_part(2) - beta_y * b_at_part(1)))

    ! Dipole emission intensity

    tau0 = q0**2 / (6.0_num * pi * epsilon0 * m0 * c**3)

    i_e = tau0 * gamma_rel**2 * (flperp(1)**2 + flperp(2)**2 + flperp(3)**2 &
        + (q0 * (beta_x * e_at_part(1) + beta_y * e_at_part(2) &
        + beta_z * e_at_part(3)) / moduclip)**2) / m0

    roland_eta = coeff_eta * SQRT(i_e)

    ! Determine eta from fields
    calculate_eta = roland_eta

  END FUNCTION calculate_eta



  FUNCTION calculate_chi(part_x, part_y, dir_x, dir_y, dir_z, part_e)

    REAL(num) :: calculate_chi
    REAL(num), INTENT(IN) :: part_x, part_y
    REAL(num), INTENT(IN) :: dir_x, dir_y, dir_z, part_e
    REAL(num) :: e_at_part(3), b_at_part(3), q(3)
    REAL(num) :: e_dot_dir, calculate_chi_roland

    CALL field_at_particle(part_x, part_y, e_at_part, b_at_part)

    e_dot_dir = e_at_part(1) * dir_x + e_at_part(2) * dir_y &
        + e_at_part(3) * dir_z

    q(1) = e_at_part(1) - e_dot_dir * dir_x &
        + c * (dir_y * b_at_part(3) - dir_z * b_at_part(2))
    q(2) = e_at_part(2) - e_dot_dir * dir_y &
        + c * (dir_z * b_at_part(1) - dir_x * b_at_part(3))
    q(3) = e_at_part(3) - e_dot_dir * dir_z &
        + c * (dir_x * b_at_part(2) - dir_y * b_at_part(1))

    calculate_chi_roland = 0.5_num * SQRT(q(1)**2 + q(2)**2 + q(3)**2) &
        * part_e / e_s

    ! Determine chi from fields
    calculate_chi = calculate_chi_roland

  END FUNCTION calculate_chi



  SUBROUTINE field_at_particle(part_x, part_y, e_at_part, b_at_part)

    REAL(num), INTENT(IN) :: part_x, part_y
    REAL(num), INTENT(OUT) :: e_at_part(3), b_at_part(3)

    ! Contains the integer cell position of the particle in x, y, z
    INTEGER :: cell_x1, cell_x2, cell_y1, cell_y2

    ! Contains the floating point version of the cell number (never actually
    ! used)
    REAL(num) :: cell_x_r, cell_y_r

    ! The fraction of a cell between the particle position and cell boundary
    REAL(num) :: cell_frac_x, cell_frac_y

    ! Weighting factors as Eqn 4.77 page 25 of manual
    ! Eqn 4.77 would be written as
    ! F(j-1) * gmx + F(j) * g0x + F(j+1) * gpx
    ! Defined at the particle position
    REAL(num), DIMENSION(sf_min:sf_max) :: gx, gy

    ! Defined at the particle position - 0.5 grid cell in each direction
    ! This is to deal with the grid stagger
    REAL(num), DIMENSION(sf_min:sf_max) :: hx, hy
    ! Temporary variables
    INTEGER :: dcellx, dcelly
    REAL(num) :: ex_part, ey_part, ez_part
    REAL(num) :: bx_part, by_part, bz_part

    ! Particle weighting multiplication factor
#ifdef PARTICLE_SHAPE_BSPLINE3
    REAL(num) :: cf2
    REAL(num), PARAMETER :: fac = (1.0_num / 24.0_num)**c_ndims
#elif  PARTICLE_SHAPE_TOPHAT
    REAL(num), PARAMETER :: fac = (1.0_num)**c_ndims
#else
    REAL(num) :: cf2
    REAL(num), PARAMETER :: fac = (0.5_num)**c_ndims
#endif

    ! Grid cell position as a fraction.
#ifdef PARTICLE_SHAPE_TOPHAT
    cell_x_r = part_x / dx - 0.5_num
    cell_y_r = part_y / dy - 0.5_num
#else
    cell_x_r = part_x / dx
    cell_y_r = part_y / dy
#endif
    ! Round cell position to nearest cell
    cell_x1 = FLOOR(cell_x_r + 0.5_num)
    ! Calculate fraction of cell between nearest cell boundary and particle
    cell_frac_x = REAL(cell_x1, num) - cell_x_r
    cell_x1 = cell_x1 + 1

    cell_y1 = FLOOR(cell_y_r + 0.5_num)
    cell_frac_y = REAL(cell_y1, num) - cell_y_r
    cell_y1 = cell_y1 + 1

    ! Particle weight factors as described in the manual, page25
    ! These weight grid properties onto particles
    ! Also used to weight particle properties onto grid, used later
    ! to calculate J
    ! NOTE: These weights require an additional multiplication factor!
#ifdef PARTICLE_SHAPE_BSPLINE3
#include "bspline3/gx.inc"
#elif  PARTICLE_SHAPE_TOPHAT
#include "tophat/gx.inc"
#else
#include "triangle/gx.inc"
#endif

    ! Now redo shifted by half a cell due to grid stagger.
    ! Use shifted version for ex in X, ey in Y, ez in Z
    ! And in Y&Z for bx, X&Z for by, X&Y for bz
    cell_x2 = FLOOR(cell_x_r)
    cell_frac_x = REAL(cell_x2, num) - cell_x_r + 0.5_num
    cell_x2 = cell_x2 + 1

    cell_y2 = FLOOR(cell_y_r)
    cell_frac_y = REAL(cell_y2, num) - cell_y_r + 0.5_num
    cell_y2 = cell_y2 + 1

    dcellx = 0
    dcelly = 0
    ! NOTE: These weights require an additional multiplication factor!
#ifdef PARTICLE_SHAPE_BSPLINE3
#include "bspline3/hx_dcell.inc"
#elif  PARTICLE_SHAPE_TOPHAT
#include "tophat/hx_dcell.inc"
#else
#include "triangle/hx_dcell.inc"
#endif

    ! These are the electric and magnetic fields interpolated to the
    ! particle position. They have been checked and are correct.
    ! Actually checking this is messy.
#ifdef PARTICLE_SHAPE_BSPLINE3
#include "bspline3/e_part.inc"
#include "bspline3/b_part.inc"
#elif  PARTICLE_SHAPE_TOPHAT
#include "tophat/e_part.inc"
#include "tophat/b_part.inc"
#else
#include "triangle/e_part.inc"
#include "triangle/b_part.inc"
#endif

    ! update particle momenta using weighted fields
    ! ex_part etc are NOT fields at particle, but fac times
    ! field

    e_at_part(1) = fac * ex_part
    e_at_part(2) = fac * ey_part
    e_at_part(3) = fac * ez_part

    b_at_part(1) = fac * bx_part
    b_at_part(2) = fac * by_part
    b_at_part(3) = fac * bz_part

  END SUBROUTINE field_at_particle



  SUBROUTINE generate_photon(generating_electron, iphoton, eta)

    ! Generates a photon moving in same direction as electron
    ! (generates entirely new photon)

    TYPE(particle), POINTER :: generating_electron
    INTEGER, INTENT(IN) :: iphoton
    REAL(num), INTENT(IN) :: eta
    REAL(num) :: dir_x, dir_y, dir_z, mag_p, generating_gamma
    REAL(num) :: rand_temp, photon_energy
    TYPE(particle), POINTER :: new_photon
    !XY for stokes
    REAL(num) :: part_x,part_y, I1,I2,V,cos_phi,Omega
    REAL(num), DIMENSION(3) :: dir_observation,dir1,dir2,e_at_part,b_at_part

    mag_p = MAX(SQRT(generating_electron%part_p(1)**2 &
        + generating_electron%part_p(2)**2 &
        + generating_electron%part_p(3)**2), c_tiny)

    dir_x = generating_electron%part_p(1) / mag_p
    dir_y = generating_electron%part_p(2) / mag_p
    dir_z = generating_electron%part_p(3) / mag_p

    generating_gamma = SQRT(1.0_num + (mag_p / m0 / c)**2)

    ! Determine photon energy

    rand_temp = random()
    photon_energy = calculate_photon_energy(rand_temp, eta, generating_gamma)
    photon_energy_record = photon_energy_record +photon_energy*generating_electron%weight

    IF (use_radiation_reaction) THEN
      ! Calculate electron recoil
      mag_p = mag_p - photon_energy / c

      generating_electron%part_p(1) = dir_x * mag_p
      generating_electron%part_p(2) = dir_y * mag_p
      generating_electron%part_p(3) = dir_z * mag_p
    END IF

    ! This will only create photons that have energies above a user specified
    ! cutoff and if photon generation is turned on. E+/- recoil is always
    ! considered
    IF (photon_energy > photon_energy_min .AND. produce_photons) THEN
      IF (photon_energy < c_tiny) photon_energy = c_tiny

      CALL create_particle(new_photon)
      new_photon%part_pos = generating_electron%part_pos

      new_photon%part_p(1) = dir_x * photon_energy / c
      new_photon%part_p(2) = dir_y * photon_energy / c
      new_photon%part_p(3) = dir_z * photon_energy / c

      new_photon%optical_depth = reset_optical_depth()
      new_photon%particle_energy = photon_energy
      new_photon%weight = generating_electron%weight
      !XY 
      new_photon%generate_time = time
      
            !XY
      part_x  = new_photon%part_pos(1) - x_grid_min_local
      part_y  = new_photon%part_pos(2) - y_grid_min_local
      CALL field_at_particle(part_x,part_y,e_at_part, b_at_part)

      !CALL calculate_stokes(new_photon,b_at_part,dir_observation,generating_gamma,eta,     I1,I2,V,o_oc,gam_phi)

      CALL calculate_stokes2(new_photon,b_at_part,generating_gamma,eta,I1,I2,V,dir1,dir2)

     new_photon%stokes = (/ I1, I2, V /)
     new_photon%stokes_dir = dir1

      CALL add_particle_to_partlist(species_list(iphoton)%attached_list, &
          new_photon)
    END IF

  END SUBROUTINE generate_photon

  SUBROUTINE calculate_optical_stokes2(p,b,wei,I1,I2,U,V,dir1,dir2)
      ! calculate the optical radiation of electron 
      ! so p is the electron's momentum
    REAL(num), DIMENSION(3), INTENT(IN) ::b,p
    REAL(num), INTENT(IN) :: wei

    REAL(num), DIMENSION(3), INTENT(OUT) ::dir1,dir2
    REAL(num), INTENT(OUT) :: I1,I2,U,V
    REAL(num), DIMENSION(3):: part_u,temp_dir
    REAL(num) :: abs_b, energy_c,x,gam_phi,gamma_e,xIntkv53,part_u2,energy0
    REAL(num) :: sin_alpha,energy,y,phi,Psi_23,Psi_13,P1_Pt,P2_Pt,xi,kv23_div_Intkv53,log10_x,log10_y,dphi,dtheta,abs_p,cos_alpha
    REAL(num) :: sin_theta, cos_theta
    REAL(num) :: kv13_div_kv23,coef,g_theta

      part_u = p/m0/c
      abs_b = MAX(c_tiny,SQRT(b(1)**2 + b(2)**2 + b(3)**2))
      abs_p = MAX(c_tiny,SQRT(part_u(1)**2 + part_u(2)**2 + part_u(3)**2))
      !cos_theta =  (b(1)*n(1) + b(2) * n(2) + b(3)*n(3))/ abs_b
      !sin_theta = SQRT(1 - cos_theta**2)

      part_u2 = part_u(1)**2+part_u(2)**2+part_u(3)**2
      gamma_e = SQRT(part_u2 + 1.0_num)

      ! energy for specific wavelength 500nm: 
      energy = detect_energy

      !cos_alpha = (b(1)*part_u(1) + b(2)*part_u(2)+b(3)*part_u(3))/abs_b/abs_p
      !sin_alpha = SQRT(1 - cos_alpha**2)
      !磁场与观测方向的夹角
      cos_theta = DOT_PRODUCT(b,dir_ob)/abs_b
      sin_theta = SQRT(1 - cos_theta**2)
      energy_c = hbar_qe_div_me*3/2*abs_b*sin_theta*gamma_e**2

      x = energy/energy_c

      !x = 1.767447518097506E-002
      xIntkv53 = find_value_from_table_1d(x,n_sample_xIntkv53,log_xIntkv53,xIntkv53_table)
!      !xIntkv53 = 1.0_num
      !energy0 = SQRT(3.0_num)*q0**2/c*q0*abs_b/m0*sin_alpha*xIntkv53

      energy0 = sin_theta*xIntkv53*abs_b*q0/m0
      !改成sin_theta
      kv23_div_Intkv53 = find_value_from_table_1d(x,n_sample_kv23divIntkv53,log_xIntkv53,kv23_div_Intkv53_table)
      kv13_div_kv23 = find_value_from_table_1d(x,n_sample_xi,log_xi,kv13_div_kv23_table)
!     
      I1 = MAX(c_tiny, 0.5*energy0*wei*(1 + kv23_div_Intkv53))
      I2 = MAX(c_tiny, 0.5*energy0*wei*(1 - kv23_div_Intkv53))
      !I1 = MAX(c_tiny, 0.5*(1 + kv23_div_Intkv53))
      !I2 = MAX(c_tiny, 0.5*(1 - kv23_div_Intkv53))
      !这里I1(+)是垂直于磁场投影方向的强度,I2(-)是平行于磁场投影方向的强度
      U = 0
      !coef = cot(theta)
      !theta 是磁场同观测方向的夹角
      CALL cross(b,dir_ob,temp_dir);
      !coef = MIN(100.0_num,DOT_PRODUCT(b,dir_ob)/SQRT(SUM(temp_dir**2)))
      !coef = MAX(EPSILON(1.0_num),DOT_PRODUCT(b,dir_ob)/SQRT(SUM(temp_dir**2)))
      coef = 1.0_num
      g_theta = 1.0_num
      IF (sin_theta > 1e-6_num) THEN
      V = MAX(c_tiny,energy0*wei*4.0_num/3.0_num*coef*(kv23_div_Intkv53*kv13_div_kv23 + &
        (1 + g_theta)/x*(kv23_div_Intkv53 - 0.5_num))*sqrt(3.0_num*x/2.0_num)*x**(-1/2)*gamma_e*cos_theta/sin_theta)
      ELSE
      V = 0.0_num
      END IF
      !WRITE(*,*) x, I1 + I2

      IF (debug_mod) THEN
        WRITE(*,*) 'energy,energy_c,b,sin_alpha,xIntkv53',energy/q0,energy_c/q0,abs_b,sin_theta,xIntkv53
        WRITE(*,*) 'x,gamma_e,kv23/Intkv53,kv13/kv23',x,gamma_e,kv23_div_Intkv53,kv13_div_kv23
        WRITE(*,*) 'I1',I1
        WRITE(*,*) 'I2',I2
        WRITE(*,*) 'U',U
        WRITE(*,*) 'V',V
      END IF
      CALL calculate_stokes3_dir(dir_ob,b,dir1,dir2)
      !dir1 = (/1.0,0.0,0.0/)
      !dir2 = (/0.0,1.0,0.0/)
 END SUBROUTINE calculate_optical_stokes2

  SUBROUTINE calculate_stokes_angle(n,b,p,phi,sin_alpha)
    REAL(num), DIMENSION(3), INTENT(IN) :: n,b,p
    REAL(num), INTENT(OUT) :: phi,sin_alpha
    REAL(num), DIMENSION(3) :: l3, p2, p1,p1xn
    REAL(num) :: abs_p1,abs_b,cos_alpha,cos_phi,sgn
    
      abs_b = MAX(c_tiny, SQRT(b(1)**2 + b(2)**2 + b(3)**2))
      l3 = (/b(2)*n(3) - b(3)*n(2),&
             b(3)*n(1) - b(1)*n(3),&
             b(1)*n(2) - b(2)*n(1)/)
      l3 = l3/MAX(c_tiny,SQRT(l3(1)**2 + l3(2)**2 + l3(3)**2))
      p2 = (p(1)*l3(1) + p(2)*l3(2) + p(3)*l3(3))*l3
      p1 = (/p(1) - p2(1),p(2) - p2(2),p(3) - p2(3)/)
      abs_p1 = MAX(c_tiny, SQRT(p1(1)**2 + p1(2)**2 + p1(3)**2))
      

      cos_alpha = (p1(1)*b(1) + p1(2)*b(2) + p1(3)*b(3))/abs_p1/abs_b
      sin_alpha = SQRT(1 - cos_alpha**2)
      cos_phi = (p1(1)*n(1) + p1(2)*n(2)+p1(3)*n(3))/abs_p1

      !sgn is dependent (n cross p1) * l3
      p1xn = (/p1(2)*n(3) - p1(3)*n(2),&
             p1(3)*n(1) - p1(1)*n(3),&
             p1(1)*n(2) - p1(2)*n(1)/)
      
      sgn = p1xn(1)*l3(1) + p1xn(2)*l3(2) + p1xn(3)*l3(3)
      phi = SIGN(ACOS(cos_phi),sgn)

  END SUBROUTINE

  SUBROUTINE calculate_stokes2_dir(p1,b,dir1,dir2)
    REAL(num), DIMENSION(3), INTENT(IN) ::b,p1
    REAL(num), DIMENSION(3), INTENT(OUT) ::dir1,dir2
    REAL(num), DIMENSION(3) :: p
    !可能有问题，按照 M.P.C.Legg & K.C.Westfold的文章[1959,1968].
    !单电子辐射的主轴是在平行或垂直于（磁场k在垂直于观测方向n的平面上的投影,i1是平行,i2是垂直)
    !这里写的是b \times p 是以电子运动方向作为观测方向，需要修正[2021.01.26 XY],见calculate_stokes3_dir()

    !dir1 = dir_parallel = b \times p
    p = p1/MAX(c_tiny,SQRT(p1(1)**2 + p1(2)**2 + p1(3)**2))
    dir1  = (/b(2)*p(3) - b(3)*p(2),&
             b(3)*p(1) - b(1)*p(3),&
             b(1)*p(2) - b(2)*p(1)/)
    dir1 = dir1/MAX(c_tiny,SQRT(dir1(1)**2 + dir1(2)**2 + dir1(3)**2))
    !dir2 = p \times dir1
    dir2  =(/p(2)*dir1(3) - p(3)*dir1(2),&
             p(3)*dir1(1) - p(1)*dir1(3),&
             p(1)*dir1(2) - p(2)*dir1(1)/)
    dir2 = dir2/MAX(c_tiny,SQRT(dir2(1)**2 + dir2(2)**2 + dir2(3)**2))
  END SUBROUTINE calculate_stokes2_dir

  SUBROUTINE calculate_stokes3_dir(n,b,dir1,dir2)
    !这里是计算完主轴上的能量后，确定主轴方向的函数
    !按照Westfold(1959a,b),M.P.C.Legg(1968),其主轴方向应该这样来确定
    !假设观测方向是n,磁场方向是k,由n \times k 可以确定一个方向是i2,
    !再用n \times i2 = i1,即i2是垂直于磁场在n确定平面上的投影的方向
    !而i1是平行与磁场投影的方向.
    !但是疑惑的是，可能和电子运动方向无关。可能的解释是需要在回旋时间内进行积分,还有就是如何排除n // b的情况?
    !dir1 是垂直于磁场的投影
    !dir2 是平行于磁场的投影
    REAL(num), DIMENSION(3), INTENT(IN) ::b,n
    REAL(num), DIMENSION(3), INTENT(OUT) ::dir1,dir2
    REAL(num), DIMENSION(3) :: p

    !dir1 = dir_parallel = b \times p
    dir1  = (/b(2)*n(3) - b(3)*n(2),&
             b(3)*n(1) - b(1)*n(3),&
             b(1)*n(2) - b(2)*n(1)/)

    dir1 = dir1/MAX(c_tiny,SQRT(dir1(1)**2 + dir1(2)**2 + dir1(3)**2))
    !dir2 = n \times dir1
    dir2  =(/n(2)*dir1(3) - n(3)*dir1(2),&
             n(3)*dir1(1) - n(1)*dir1(3),&
             n(1)*dir1(2) - n(2)*dir1(1)/)

    dir2 = dir2/MAX(c_tiny,SQRT(dir2(1)**2 + dir2(2)**2 + dir2(3)**2))
  END SUBROUTINE calculate_stokes3_dir


  SUBROUTINE calculate_stokes2(photon,b,gamma_e,eta,I1,I2,V,dir1,dir2)
    TYPE(particle), POINTER, INTENT(IN) :: photon
    REAL(num), DIMENSION(3), INTENT(IN) ::b
    REAL(num), DIMENSION(3), INTENT(OUT) ::dir1,dir2

    REAL(num), INTENT(IN) :: gamma_e, eta
    REAL(num), INTENT(OUT) :: I1,I2,V
    REAL(num) :: abs_b, energy_c,x,gam_phi
    REAL(num) :: sin_alpha,energy,y,phi,Psi_23,Psi_13,P1_Pt,P2_Pt,xi,kv23_div_Intkv53,log10_x,log10_y,dphi,dtheta,abs_p,cos_alpha
    REAL(num), DIMENSION(3) :: p
      !photon_motion is the electron motion
      !n is the observation direction
      !cos(alpha) = b \cdot n/abs_b
      p = (/photon%part_p(1)/m0/c,&
            photon%part_p(2)/m0/c,&
            photon%part_p(3)/m0/c/)

      abs_b = MAX(c_tiny,SQRT(b(1)**2 + b(2)**2 + b(3)**2))
      abs_p = MAX(c_tiny,SQRT(p(1)**2 + p(2)**2 + p(3)**2))
      !cos_theta =  (b(1)*n(1) + b(2) * n(2) + b(3)*n(3))/ abs_b
      !sin_theta = SQRT(1 - cos_theta**2)

      energy = photon%particle_energy

      cos_alpha = (b(1)*p(1) + b(2)*p(2)+b(3)*p(3))/abs_b/abs_p
      sin_alpha = SQRT(1 - cos_alpha**2)
      energy_c = hbar_qe_div_me*3/2*abs_b*sin_alpha*gamma_e**2
      x = energy/energy_c

      gam_phi = 0.0_num
      kv23_div_Intkv53 = find_value_from_table_1d(x,n_sample_kv23divIntkv53,log_xIntkv53,kv23_div_Intkv53_table)
      !kv23_div_Intkv53 = find_value_from_table_1d(x,n_sample_xi,log_xi,kv23_div_Intkv53_table)
     
     I1 = 0.5*energy*photon%weight*(1 + kv23_div_Intkv53)
     I2 = 0.5*energy*photon%weight*(1 - kv23_div_Intkv53)
     !U = 0
     V =  gamma_e
     CALL calculate_stokes2_dir(p,b,dir1,dir2)
 END SUBROUTINE calculate_stokes2

  SUBROUTINE calculate_stokes(photon,b,n,gamma_e,eta,I1,I2,V)
    TYPE(particle), POINTER, INTENT(IN) :: photon
    REAL(num), DIMENSION(3), INTENT(IN) ::b,n
    REAL(num), INTENT(IN) :: gamma_e, eta
    REAL(num), INTENT(OUT) :: I1,I2,V
    REAL(num) :: abs_b,x,gam_phi
    REAL(num) :: sin_alpha,energy,y,phi,Psi_23,Psi_13,P1_Pt,P2_Pt,xi,kv13_div_kv23,log10_x,log10_y,dphi,dtheta
    REAL(num), DIMENSION(3) :: p
      !photon_motion is the electron motion
      !n is the observation direction
      !cos(alpha) = b \cdot n/abs_b
      p = (/photon%part_p(1)/m0/c,&
            photon%part_p(2)/m0/c,&
            photon%part_p(3)/m0/c/)

      abs_b = MAX(c_tiny,SQRT(b(1)**2 + b(2)**2 + b(3)**2))
      !cos_theta =  (b(1)*n(1) + b(2) * n(2) + b(3)*n(3))/ abs_b
      !sin_theta = SQRT(1 - cos_theta**2)

      energy = photon%particle_energy

      !First step calculate the projection of momentum in k, B plane
      !p1 for the momentum in the plane, p2 for the perpendicular vector b \cross n
      !
      CALL calculate_stokes_angle(n,b,p,phi,sin_alpha)
      !sin_theta = SQRT(1 - cos_theta**2)
      !x = omega/omega_c
      x = energy/(h_bar*(3/2*abs_b*q0/m0*sin_alpha*gamma_e**2))
      !x = energy/(gamma_e*mc0*c)*(2/3+eta)/eta

      !phi is the angle between n and k, need +/- phi = \theta - \alpha
      !sin_phi = sin_theta*cos_alpha-cos_theta*sin_alpha
      !WRITE(*,*) 'cos_theta,cos_alpha,sin_phi',cos_theta, cos_alpha, sin_phi

     !Psi_{2/3}(x,\gamma/phi)
     log10_x = LOG10(x) 
     gam_phi = gamma_e*phi
     CALL interpolate_2d(log10_x, gam_phi, &
        n_sample_y, n_sample_gp, &
        log_y, gp, &
        Psi23_table, Psi_23)
     CALL interpolate_2d(log10_x, gam_phi, &
        n_sample_y, n_sample_gp, log_y, gp, Psi13_table, Psi_13)

    !    WRITE(*,*) 'log10_x,gamma_e*phi',log10_x,gamma_e*phi
    !    WRITE(*,*) 'Psi23,Psi13',Psi_23,Psi_13
     P1_Pt = MAX(SQRT(3.0)/pi**2/sin_alpha**3*gamma_e**2*Psi_23,c_tiny)
     P2_Pt = MAX(SQRT(3.0)/pi**2/sin_alpha**3*gamma_e**2*Psi_13,c_tiny)
     !WRITE(*,*) x,gam_phi,Psi_23

     !xi = ellipse
     y = x/2.0*(1.0 + gam_phi**2)**(3.0/2.0)

     kv13_div_kv23 = find_value_from_table_1d(y,n_sample_xi,log_xi,kv13_div_kv23_table)
     xi = gam_phi/SQRT(1 + gam_phi**2)*kv13_div_kv23

     !Then P1*cos(theta)*dphi*dtheta/Pt
     dphi = pi/100
     dtheta = pi/100
     I1 = Psi_23*energy*COS(phi)*dphi*dtheta
     I2 = Psi_13*energy*COS(phi)*dphi*dtheta
     !U = 0
     V = 2*xi/(1+xi**2)*(I1+I2)
     !IF (I1+I2 > c_tiny .AND. V/(I1+I2) > 0.9) THEN
     !WRITE(*,*) 'x,gp,y:',log10(x),gam_phi,y
     !WRITE(*,*) 'Psi23,Psi13,kv13_div_kv23',Psi_23,Psi_13,kv13_div_kv23
     !WRITE(*,*) 'I1,I2,V, Sum',I1,I2,V,(V**2 + (I2-I1)**2)
     !ENDIF

 END SUBROUTINE calculate_stokes

SUBROUTINE interpolate_2d(x_in, y_in, nx, ny, x, y, f_table,f)

    REAL(num), INTENT(IN) :: x_in, y_in
    INTEGER, INTENT(IN) :: nx, ny
    REAL(num), INTENT(IN) :: x(nx), y(ny), f_table(nx,ny)
    REAL(num), INTENT(OUT) :: f
    INTEGER :: ix, iy,ixp,iyp
    REAL(num) :: xr,yr,dx,dy
    LOGICAL, SAVE :: warning = .TRUE.
    !
    dx = x(2) - x(1);
    dy = y(2) - y(1);
    ix = FLOOR((x_in-x(1))/dx) + 1
    iy = FLOOR((y_in-y(1))/dy) + 1
    IF (ix < nx .AND. ix > 0 .AND. iy < ny .AND. iy > 0) THEN
    xr = (x_in - x(ix))/dx
    yr = (y_in - y(iy))/dy
    ixp = ix + 1
    iyp = iy + 1
!   IF (ix < 1) THEN 
!        ix = 1
!        xr = 0.0
!        ixp = ix
!        WRITE(*,*), 'ix < 1'
!   ENDIF
!   IF (ix >= nx) THEN
!        ix = nx
!        ixp = ix
!        xr = 1.0
!        WRITE(*,*), 'ix >=nx'
!   ENDIF
!   IF (iy < 1) THEN 
!        iy = 1
!        yr = 0.0
!        iyp = iy
!        WRITE(*,*), 'iy <1'
!   ENDIF
!   IF (iy >= ny) THEN
!        iy = ny
!        yr = 1.0
!        iyp = iy
!        WRITE(*,*), 'iy >=nx'
!   ENDIF
!
    f =   f_table(ix,iy)*(1-xr)*(1-yr) &
        + f_table(ixp,iy)*xr*(1-yr) &
        + f_table(ix,iyp)*(1-xr)*yr &
        + f_table(ixp,iyp)*xr*yr
    ELSE
    f = 0.0_num
ENDIF


    !f(x,y) = f(0,0)(1-x)*(1-y)+f(1.0)*x*(1-y)+f(0,1)*(1-x)*y+f(1,1)*x*y
  END SUBROUTINE interpolate_2d



  FUNCTION calculate_photon_energy(rand_seed, eta, generating_gamma)

    REAL(num) :: calculate_photon_energy
    REAL(num), INTENT(IN) :: rand_seed, eta, generating_gamma
    REAL(num) :: chi_final

    chi_final = find_value_from_table_alt(eta, rand_seed, &
        n_sample_eta, n_sample_chi, log_eta, log_chi, p_photon_energy)

    calculate_photon_energy = (2.0_num * chi_final / eta) * generating_gamma &
        * m0 * c**2

  END FUNCTION calculate_photon_energy



  SUBROUTINE generate_pair(generating_photon, chi_val, iphoton, ielectron, &
      ipositron)

    ! Generates a pair moving in same direction as photon
    TYPE(particle), POINTER :: generating_photon
    REAL(num), INTENT(IN) :: chi_val
    INTEGER, INTENT(IN) :: iphoton, ielectron, ipositron
    REAL(num) :: dir_x, dir_y, dir_z, mag_p
    REAL(num) :: probability_split, epsilon_frac, norm
    TYPE(particle), POINTER :: new_electron, new_positron

    CALL create_particle(new_electron)
    CALL create_particle(new_positron)

    new_electron%part_pos = generating_photon%part_pos
    new_positron%part_pos = generating_photon%part_pos

    norm  = c / generating_photon%particle_energy
    dir_x = generating_photon%part_p(1) * norm
    dir_y = generating_photon%part_p(2) * norm
    dir_z = generating_photon%part_p(3) * norm

    ! Determine how to split the energy amoung e-/e+
    ! IS CHI HERE SAME AS ROLAND'S? DEFINED BSinT/B_s

    probability_split = random()

    epsilon_frac = find_value_from_table(chi_val, probability_split, &
        n_sample_chi2, n_sample_epsilon, log_chi2, epsilon_split, p_energy)

    mag_p = MAX(generating_photon%particle_energy / c, c_tiny)

    new_electron%part_p(1) = epsilon_frac * mag_p * dir_x
    new_electron%part_p(2) = epsilon_frac * mag_p * dir_y
    new_electron%part_p(3) = epsilon_frac * mag_p * dir_z

    new_positron%part_p(1) = (1.0_num - epsilon_frac) * mag_p * dir_x
    new_positron%part_p(2) = (1.0_num - epsilon_frac) * mag_p * dir_y
    new_positron%part_p(3) = (1.0_num - epsilon_frac) * mag_p * dir_z

    new_electron%optical_depth = reset_optical_depth()
    new_positron%optical_depth = reset_optical_depth()

#ifdef TRIDENT_PHOTONS
    new_electron%optical_depth_tri = reset_optical_depth()
    new_positron%optical_depth_tri = reset_optical_depth()
#endif

    new_electron%weight = generating_photon%weight
    new_positron%weight = generating_photon%weight
    CALL add_particle_to_partlist(species_list(ielectron)%attached_list, &
        new_electron)
    CALL add_particle_to_partlist(species_list(ipositron)%attached_list, &
        new_positron)

    ! Remove photon
    CALL remove_particle_from_partlist(species_list(iphoton)%attached_list, &
        generating_photon)

    DEALLOCATE(generating_photon)

  END SUBROUTINE generate_pair



  SUBROUTINE generate_pair_tri(generating_electron, ielectron, ipositron)

    ! Generates a pair moving in same direction as photon
    TYPE(particle), POINTER :: generating_electron
    INTEGER, INTENT(IN) :: ielectron, ipositron
    TYPE(particle), POINTER :: new_electron, new_positron

    CALL create_particle(new_electron)
    CALL create_particle(new_positron)

    new_electron%part_pos = generating_electron%part_pos
    new_positron%part_pos = generating_electron%part_pos

    new_electron%part_p = 0.0_num
    new_positron%part_p = 0.0_num

    new_electron%optical_depth = reset_optical_depth()
    new_positron%optical_depth = reset_optical_depth()

#ifdef TRIDENT_PHOTONS
    new_electron%optical_depth_tri = reset_optical_depth()
    new_positron%optical_depth_tri = reset_optical_depth()
#endif

    new_electron%weight = generating_electron%weight
    new_positron%weight = generating_electron%weight

    CALL add_particle_to_partlist(species_list(ielectron)%attached_list, &
        new_electron)
    CALL add_particle_to_partlist(species_list(ipositron)%attached_list, &
        new_positron)

  END SUBROUTINE generate_pair_tri



  FUNCTION find_value_from_table_1d(x_in, nx, x, values)

    REAL(num) :: find_value_from_table_1d
    REAL(num), INTENT(IN) :: x_in
    INTEGER, INTENT(IN) :: nx
    REAL(num), INTENT(IN) :: x(nx), values(nx)
    REAL(num) :: fx, x_value, value_interp, xdif1, xdif2, xdifm
    INTEGER :: i1, i2, im
    LOGICAL, SAVE :: warning = .TRUE.

    x_value = LOG10(MAX(x_in,c_tiny))

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

    find_value_from_table_1d = 10.0_num**value_interp

  END FUNCTION find_value_from_table_1d



  FUNCTION find_value_from_table_alt(x_in, p_value, nx, ny, x, y, p_table)

    REAL(num) :: find_value_from_table_alt
    REAL(num), INTENT(IN) :: x_in, p_value
    INTEGER, INTENT(IN) :: nx, ny
    REAL(num), INTENT(IN) :: x(nx), y(nx,ny), p_table(nx,ny)
    INTEGER :: ix, index_lt, index_gt, i1, i2, im
    REAL(num) :: fx, fp, y_lt, y_gt, x_value, y_interp, xdif1, xdif2, xdifm
    LOGICAL, SAVE :: warning = .TRUE.

    x_value = LOG10(x_in)

    ! Scan through x to find correct row of table
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
        PRINT*,'Argument to "find_value_from_table_alt" outside the range ', &
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

    index_lt = i1
    index_gt = i2

    ix = index_lt
    ! Scan through table row to find p_value
    i1 = 1
    i2 = ny
    xdif1 = p_table(ix,i1) - p_value
    xdif2 = p_table(ix,i2) - p_value
    IF (xdif1 * xdif2 < 0) THEN
      ! Use bisection to find the nearest cell
      DO
        im = (i1 + i2) / 2
        xdifm = p_table(ix,im) - p_value
        IF (xdif1 * xdifm < 0) THEN
          i2 = im
        ELSE
          i1 = im
          xdif1 = xdifm
        END IF
        IF (i2 - i1 == 1) EXIT
      END DO
      ! Interpolate in x
      fp = (p_value - p_table(ix,i1)) / (p_table(ix,i2) - p_table(ix,i1))
    ELSE
      IF (warning .AND. rank == 0) THEN
        PRINT*,'*** WARNING ***'
        PRINT*,'Argument to "find_value_from_table_alt" outside the range ', &
            'of the table.'
        PRINT*,'Using truncated value. No more warnings will be issued.'
        warning = .FALSE.
      END IF
      IF (xdif1 >= 0) THEN
        fp = 0.0_num
      ELSE
        fp = 1.0_num
      END IF
    END IF

    y_lt = (1.0_num - fp) * y(ix,i1) + fp * y(ix,i2)

    ix = index_gt
    ! Scan through table row to find p_value
    i1 = 1
    i2 = ny
    xdif1 = p_table(ix,i1) - p_value
    xdif2 = p_table(ix,i2) - p_value
    IF (xdif1 * xdif2 < 0) THEN
      ! Use bisection to find the nearest cell
      DO
        im = (i1 + i2) / 2
        xdifm = p_table(ix,im) - p_value
        IF (xdif1 * xdifm < 0) THEN
          i2 = im
        ELSE
          i1 = im
          xdif1 = xdifm
        END IF
        IF (i2 - i1 == 1) EXIT
      END DO
      ! Interpolate in x
      fp = (p_value - p_table(ix,i1)) / (p_table(ix,i2) - p_table(ix,i1))
    ELSE
      IF (warning .AND. rank == 0) THEN
        PRINT*,'*** WARNING ***'
        PRINT*,'Argument to "find_value_from_table_alt" outside the range ', &
            'of the table.'
        PRINT*,'Using truncated value. No more warnings will be issued.'
        warning = .FALSE.
      END IF
      IF (xdif1 >= 0) THEN
        fp = 0.0_num
      ELSE
        fp = 1.0_num
      END IF
    END IF

    !WRITE(*,*) 'epo, p, ix,i1,i2', x_in,p_value, ix, i1 ,i2 
    y_gt = (1.0_num - fp) * y(ix,i1) + fp * y(ix,i2)

    ! Interpolate in x

    y_interp = (1.0_num - fx) * y_lt + fx * y_gt

    find_value_from_table_alt = 10.0_num**y_interp


  END FUNCTION find_value_from_table_alt



  FUNCTION find_value_from_table(x_in, p_value, nx, ny, x, y, p_table)

    REAL(num) :: find_value_from_table
    REAL(num), INTENT(IN) :: x_in, p_value
    INTEGER, INTENT(IN) :: nx, ny
    REAL(num), INTENT(IN) :: x(nx), y(ny), p_table(nx,ny)
    INTEGER :: ix, index_lt, index_gt, i1, i2, im
    REAL(num) :: fx, fp, y_lt, y_gt, x_value, y_interp, xdif1, xdif2, xdifm
    LOGICAL, SAVE :: warning = .TRUE.

    x_value = LOG10(x_in)

    ! Scan through x to find correct row of table
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
        PRINT*,'Argument to "find_value_from_table" outside the range ', &
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

    index_lt = i1
    index_gt = i2

    ix = index_lt
    ! Scan through table row to find p_value
    i1 = 1
    i2 = ny
    xdif1 = p_table(ix,i1) - p_value
    xdif2 = p_table(ix,i2) - p_value
    IF (xdif1 * xdif2 < 0) THEN
      ! Use bisection to find the nearest cell
      DO
        im = (i1 + i2) / 2
        xdifm = p_table(ix,im) - p_value
        IF (xdif1 * xdifm < 0) THEN
          i2 = im
        ELSE
          i1 = im
          xdif1 = xdifm
        END IF
        IF (i2 - i1 == 1) EXIT
      END DO
      ! Interpolate in x
      fp = (p_value - p_table(ix,i1)) / (p_table(ix,i2) - p_table(ix,i1))
    ELSE
      IF (warning .AND. rank == 0) THEN
        PRINT*,'*** WARNING ***'
        PRINT*,'Argument to "find_value_from_table" outside the range ', &
            'of the table.'
        PRINT*,'Using truncated value. No more warnings will be issued.'
        warning = .FALSE.
      END IF
      IF (xdif1 >= 0) THEN
        fp = 0.0_num
      ELSE
        fp = 1.0_num
      END IF
    END IF

    y_lt = (1.0_num - fp) * y(i1) + fp * y(i2)

    ix = index_gt
    ! Scan through table row to find p_value
    i1 = 1
    i2 = ny
    xdif1 = p_table(ix,i1) - p_value
    xdif2 = p_table(ix,i2) - p_value
    IF (xdif1 * xdif2 < 0) THEN
      ! Use bisection to find the nearest cell
      DO
        im = (i1 + i2) / 2
        xdifm = p_table(ix,im) - p_value
        IF (xdif1 * xdifm < 0) THEN
          i2 = im
        ELSE
          i1 = im
          xdif1 = xdifm
        END IF
        IF (i2 - i1 == 1) EXIT
      END DO
      ! Interpolate in x
      fp = (p_value - p_table(ix,i1)) / (p_table(ix,i2) - p_table(ix,i1))
    ELSE
      IF (warning .AND. rank == 0) THEN
        PRINT*,'*** WARNING ***'
        PRINT*,'Argument to "find_value_from_table" outside the range ', &
            'of the table.'
        PRINT*,'Using truncated value. No more warnings will be issued.'
        warning = .FALSE.
      END IF
      IF (xdif1 >= 0) THEN
        fp = 0.0_num
      ELSE
        fp = 1.0_num
      END IF
    END IF

    y_gt = (1.0_num - fp) * y(i1) + fp * y(i2)

    ! Interpolate in x

    y_interp = (1.0_num - fx) * y_lt + fx * y_gt

    find_value_from_table = y_interp

  END FUNCTION find_value_from_table

  SUBROUTINE particle_scatter
    INTEGER :: ispecies, jspecies,photon_species,electron_species,positron_species
    INTEGER(i8) :: ix, iy
    TYPE(particle_list), POINTER :: p_list1,p_list2
    LOGICAL :: isscatter
    REAL(num) :: w1,w2,B_field

    photon_species = -1
    electron_species = -1
    positron_species = -1
    DO ispecies = 1, n_species
      IF (species_list(ispecies)%species_type == c_species_id_photon) &
          photon_species = ispecies
      IF (species_list(ispecies)%species_type == c_species_id_electron &
          .AND. ispecies /= breit_wheeler_electron_species) &
          electron_species = ispecies
      IF (species_list(ispecies)%species_type == c_species_id_positron &
          .AND. ispecies /= breit_wheeler_positron_species) &
          positron_species = ispecies
    END DO 

    isscatter = .False.
    IF (photon_species > 0 .AND. (electron_species > 0 .OR. positron_species > 0)) THEN
        isscatter = .True.
    END IF


    IF (isscatter) THEN
        IF (electron_species > 0) THEN
            DO iy = 1, ny
            DO ix = 1, nx
            p_list1 => species_list(electron_species)%secondary_list(ix,iy)
            p_list2 => species_list(photon_species)%secondary_list(ix,iy)
            CALL shuffle_particle_list_random_CS(p_list1)
            CALL shuffle_particle_list_random_CS(p_list2)
    END DO ! ix
    END DO ! iy
        
    w1 = species_list(electron_species)%weight
    w2 = species_list(photon_species)%weight

    !WRITE(*,*) photon_species,w2
    !WRITE(*,*) electron_species,w1
    DO iy = 1, ny
    DO ix = 1, nx
        B_field = sqrt(bx(ix,iy)**2 + by(ix,iy)**2 + bz(ix,iy)**2)
        CALL compton_scatter(species_list(electron_species)%secondary_list(ix,iy), &
         species_list(photon_species)%secondary_list(ix,iy), &
         w1,w2,electron_species,photon_species,B_field)
    END DO ! ix
    END DO ! iy
        END IF ! electron_species

        IF (positron_species > 0) THEN
            DO iy = 1, ny
            DO ix = 1, nx
            p_list1 => species_list(positron_species)%secondary_list(ix,iy)
            p_list2 => species_list(photon_species)%secondary_list(ix,iy)
            CALL shuffle_particle_list_random_CS(p_list1)
            CALL shuffle_particle_list_random_CS(p_list2)
    END DO ! ix
    END DO ! iy
        
    w1 = species_list(positron_species)%weight
    w2 = species_list(photon_species)%weight

    !WRITE(*,*) photon_species,w2
    !WRITE(*,*) electron_species,w1
    DO iy = 1, ny
    DO ix = 1, nx
        B_field = sqrt(bx(ix,iy)**2 + by(ix,iy)**2 + bz(ix,iy)**2)
        CALL compton_scatter(species_list(positron_species)%secondary_list(ix,iy), &
         species_list(photon_species)%secondary_list(ix,iy), &
         w1,w2,positron_species,photon_species,B_field)
    END DO ! ix
    END DO ! iy
        END IF ! positron_species

    END IF !isscatter

  END SUBROUTINE particle_scatter


  ! Binary collision compton scattering operator based jointly on:
  ! Perez et al. PHYSICS OF PLASMAS 19, 083104 (2012), and
  ! K. Nanbu and S. Yonemura, J. Comput. Phys. 145, 639 (1998)

  SUBROUTINE compton_scatter(p_list1, p_list2, weight1,weight2,ep_species,photon_species,B_field)
!      !p_list1 is elecron for compon scatter and p_list2 must be photon
!
    TYPE(particle_list), INTENT(INOUT) :: p_list1
    TYPE(particle_list), INTENT(INOUT) :: p_list2
    INTEGER, INTENT(IN),OPTIONAL :: ep_species, photon_species
    REAL(num), INTENT(IN), OPTIONAL :: weight1, weight2
    REAL(num), INTENT(IN), OPTIONAL :: B_field
    INTEGER :: icount, jcount, pcount, k
    INTEGER(num) :: N_max
    REAL(num) ::np, factor, weight_max, P_max, sigma_pk,sigma_o_pk,Pij,gamma_e
    TYPE(particle), POINTER :: current, impact,new_particle
    REAL(num) :: w1,w2,wr,energy_pho, rnd,rand_seed,energy_ele_sc,rand_seed_phi
    REAL(num):: en1,en2,ph1,ph2,m0c2,delta_energy_ele,delta_energy_pho
    REAL(num), PARAMETER :: eps = EPSILON(1.0_num)
    REAL(num), PARAMETER :: one_m_2eps = 1.0_num - 1e-7
    REAL(num), PARAMETER :: one_p_2eps = 1.0_num + 1e-7
    REAL(num), DIMENSION(3) :: p1, p2, p3, p4, vc, v1, v2, p5, p6,pe3,pp4
    REAL(num), DIMENSION(4) :: pk_in,pk_out,pk_o_sc,pk_o_out,delta_p,pk_sc
    REAL(num), DIMENSION(3) :: p1_norm, p2_norm
    TYPE(particle_list) :: append_ele_list, append_pho_list
    REAL(num) :: jitter_x, jitter_y
    !Resonance scatter
    REAL(num) :: sigma_rs_max,sigma_rs,sigma_o_rs,Prs
    LOGICAL(num) :: is_rs
!    REAL(num), PARAMETER :: pi4_eps2_c4 = 4.0_num * pi * epsilon0**2 * c**4
!    REAL(num), PARAMETER :: two_thirds = 2.0_num / 3.0_num
!    REAL(num), PARAMETER :: pi_fac = &
!                                (4.0_num * pi / 3.0_num)**(1.0_num / 3.0_num)
!
    factor = 0.0_num
    np = 0.0_num
    m0c2 = m0*c**2
!
!    ! Inter-species collisions
    icount = p_list1%count
    jcount = p_list2%count
     ! pcount should be modified

!     IF (debug_mod) THEN
!         IF (rank == 0) THEN
!            WRITE(*,*) 'icount,jcount',icount,jcount
!         END IF
!    END IF
    IF (icount > 0 .AND. jcount > 0) THEN
      ! temporarily join tail to the head of the lists to make them circular
      p_list1%tail%next => p_list1%head
      p_list2%tail%next => p_list2%head

#ifdef PER_SPECIES_WEIGHT
      weight_max =  MAX(weight1, weight2)
#else
      weight_max = 0.0_num
      current => p_list1%head
      impact => p_list2%head

      ! calculate P_max
      DO k = 1, icount
          weight_max = MAX(weight_max, current%weight)
          current => current%next
      END DO

      DO k = 1, jcount
          weight_max = MAX(weight_max, impact%weight)
          impact => impact%next
      END DO

        !get field after Lorentz_transform
        !judge whether satisfy the resonance requirement (6) and (7) in Paper PRD. J.K.Daugherty 1978
        
      !artificial limitation
      IF (B_field > 1e-3_num*b_s .AND. use_resonance_scatter) THEN

            sigma_rs_max = 9.0_num/4.0_num*(b_s/alpha_f/B_field)**2
      ELSE 
          sigma_rs_max = 0.0_num
      END IF

      sigma_rs_max = MIN(sigma_rs_max,10000.0_num)

      P_max = gfactor*Thomson_cross_section*(1 + sigma_rs_max)*c*dt*weight_max
      !P_max = gfactor*Thomson_cross_section*(sigma_rs_max)*c*dt*weight_max

      N_max = FLOOR(MIN(P_max*icount*jcount/dx/dy,REAL(HUGE(N_max),num)),num)

      IF (N_max > MIN(icount,jcount)) N_max = MIN(icount,jcount)

      IF (debug_mod) THEN
      IF( rank == 0) THEN
!      WRITE(*,*) 'weight_max',weight_max
!      WRITE(*,*) 'dt',dt
!      WRITE(*,*) 'dx,dy',dx,dy
!      WRITE(*,*) 'i,jcount',icount,jcount
      WRITE(*,*) 'gfactor, dt,weight_max,sigma_rs_max,B',gfactor,dt,weight_max,sigma_rs_max,B_field
      WRITE(*,*) 'P_max,N_max1,N_max2,icount,jcount',P_max,P_max*icount*jcount/dx/dy, N_max, icount, jcount
!      WRITE(*,*) 'i,j,P_max,N_max',icount,jcount,P_max,N_max
      END IF
      END IF !debug_mod

      current => p_list1%head
      impact => p_list2%head

#endif

!      ! If possible, use per-species properties
      w1 = weight1
      w2 = weight2
      wr = w1 / w2
!
      current => p_list1%head
      impact => p_list2%head
!
!      ! Per-cell constant factors
!
      delta_energy_ele = 0.0_num
      delta_energy_pho = 0.0_num

      ! create_empty_
      CALL create_empty_partlist(append_ele_list)
      CALL create_empty_partlist(append_pho_list)
      DO k = 1, N_max
      !cs_times for test:
        IF (impact%cs_times > 0.0_num .AND. use_cs_times_limit) CYCLE
#ifndef PER_SPECIES_WEIGHT
        w1 = current%weight
        w2 = impact%weight
        wr = w1 / w2
#endif
!       write(*,*) wr
!
        p1 = current%part_p / c
        p2 = impact%part_p / c
!
        p1_norm = p1 / m0
        p2_norm = p2 / m0
!
        ! Two stationary particles can't collide, so don't try
        IF (DOT_PRODUCT(p1_norm, p1_norm) < eps &
            .AND. DOT_PRODUCT(p2_norm, p2_norm) < eps) CYCLE


        ! Ditto for two particles with the same momentum
        !vc = (p1_norm - p2_norm)
        !IF (DOT_PRODUCT(vc, vc) < eps) THEN
        !    WRITE(*,*) 'same_momentum'
        !    CYCLE
        !ENDIF

        ! Update particle properties
        energy_pho = SQRT(DOT_PRODUCT(p2_norm,p2_norm)) ! E/mc**2
        pk_in =(/ energy_pho, p2_norm(1), p2_norm(2), p2_norm(3)/)
        gamma_e = SQRT(1 + DOT_PRODUCT(p1_norm,p1_norm))

        en1 = gamma_e*m0c2
        ph1 = energy_pho*m0c2

        !p1_4v = (/gamma_e, p1_norm(1), p1_norm(2), p1_norm(3)/)

        CALL Lorentz_transform(p1_norm,pk_in,pk_out)

        !initial cross section
        sigma_o_pk = 0.0_num
        sigma_o_rs = 0.0_num
        is_rs = .FALSE.

        CALL Compton_cross_section(pk_out,sigma_o_pk)


        IF (use_resonance_scatter) THEN

        CALL Resonance_scatter_section(current,pk_out, sigma_o_rs, is_rs,sigma_rs_max)

        END IF

        IF (pk_in(1) < eps) THEN ! For low energy photon
            sigma_pk = sigma_o_pk/gamma_e
            sigma_rs = sigma_o_rs/gamma_e
        ELSE
            sigma_pk = sigma_o_pk*pk_out(1)/pk_in(1)/gamma_e
            sigma_rs = sigma_o_rs*pk_out(1)/pk_in(1)/gamma_e
        END IF


        Pij = sigma_pk*c*dt*weight_max/P_max
        Prs = sigma_rs*c*dt*weight_max/P_max

!        IF (debug_mod) THEN
!            IF (rank == 0) THEN
!                WRITE(*,*) 'Pij,Prs,is_rs',Pij, Prs,is_rs
!                WRITE(*,*) 'sigma_pk,sigma_rs,sigma_o_rs',sigma_pk,sigma_rs,sigma_o_rs
!                WRITE(*,*) 'pk_in,pk_out,gamma_e', pk_in,pk_out,gamma_e
!            END IF
!        END IF
!
        !Mento_carlo simulation Rejection method 
        rnd = random() !test for scatter
        IF ((rnd < Pij .AND. .NOT.(is_RS) ).OR. (rnd > Pij .AND. rnd < Prs+Pij .AND. is_rs)) THEN  ! Scatter in rest frame.
            !WRITE(*,*) 'rnd, Pij', rnd, Pij
            rand_seed = random()
            rand_seed_phi = random()
            !test 
            CALL compton_scatter_photon(rand_seed,rand_seed_phi, pk_out,pk_o_sc)
            !pk_o_sc = pk_out
            !inverse Lorentz transform
            CALL Lorentz_transform(-p1_norm,pk_o_sc,pk_sc)
        !p' = p + delta P_pho
        !Update momentum
            delta_p = pk_sc - pk_in
        !current%part_p
            pe3 = current%part_p - delta_p(2:4)*m0c2/c
            pp4 = pk_sc(2:4)*m0c2/c
            

            en2 = c*SQRT(DOT_PRODUCT(pe3,pe3) + (m0*c)**2)
            ph2 = pk_sc(1)*m0c2
            !supplement energy

            IF (wr > 10) THEN
                print*, 'wr > 10',wr
!                WRITE(*,*) 'A', wr, one_p_2eps
                ALLOCATE(new_particle)
                jitter_x = (2 * random() - 1) * 0.25_num * dx
                jitter_y = (2 * random() - 1) * 0.25_num * dy
                new_particle = current
!                !modified weight
                new_particle%weight = w2
                IF (use_compton_recoil) THEN
                new_particle%part_p = pe3
                new_particle%particle_energy = en2
                END IF
                !w1 > w2 w1 - w2
                current%weight = current%weight - w2
                impact%part_p = pp4
                impact%particle_energy =  ph2
#if defined(PARTICLE_ID) || defined(PARTICLE_ID4)
              new_particle%id = 0
#endif
              new_particle%part_pos(1) = current%part_pos(1) + jitter_x
              new_particle%part_pos(2) = current%part_pos(2) + jitter_y
              new_particle%cs_times = new_particle%cs_times + 1.0_num
              new_particle%generate_time = time
              impact%cs_times = impact%cs_times + 1.0_num
#ifdef PARTICLE_DEBUG
              ! If running with particle debugging, specify that this
              ! particle has been split
              new_particle%processor_at_t0 = -1
#endif
              CALL add_particle_to_partlist(append_ele_list, new_particle)
              NULLIFY(new_particle)

              current%part_pos(1) = current%part_pos(1) - jitter_x
              current%part_pos(2) = current%part_pos(2) - jitter_y
!          !electron
!                CALL weighted_particles_correction_CS(w2 / w1, current%part_p, pe3, gamma_e*m0*c**2, energy_ele_sc, m0) 
            ELSE IF (wr < 0.1) THEN
                print*, 'wr <0.1',wr
                !WRITE(*,*) 'B',wr
                ALLOCATE(new_particle)
                jitter_x = (2 * random() - 1) * 0.25_num * dx
                jitter_y = (2 * random() - 1) * 0.25_num * dy
                new_particle = impact
                !modified weight
                new_particle%weight = w1
                impact%weight = impact%weight - w1
#if defined(PARTICLE_ID) || defined(PARTICLE_ID4)
              new_particle%id = 0
#endif
              new_particle%part_pos(1) = impact%part_pos(1) + jitter_x
              new_particle%part_pos(2) = impact%part_pos(2) + jitter_y
              new_particle%cs_times = impact%cs_times + 1.0_num
              current%cs_times = current%cs_times + 1.0_num
              new_particle%part_p = pp4
              new_particle%generate_time = time

              new_particle%particle_energy = ph2

              CALL add_particle_to_partlist(append_pho_list, new_particle)
#ifdef PARTICLE_DEBUG
              ! If running with particle debugging, specify that this
              ! particle has been split
              new_particle%processor_at_t0 = -1
#endif
              NULLIFY(new_particle)

              impact%part_pos(1) = impact%part_pos(1) - jitter_x
              impact%part_pos(2) = impact%part_pos(2) - jitter_y

              IF (use_compton_recoil) THEN
                current%part_p = pe3
                current%particle_energy = en2
              END IF
            !Merge method
            ELSE IF (wr > one_p_2eps) THEN
                print*, 'wr > one_p_2eps merge new photon',wr
                CALL weighted_particles_correction_CS(w2 / w1, current%part_p, pe3, gamma_e*m0*c**2, energy_ele_sc, m0) 
                IF (use_compton_recoil) THEN
                current%part_p = pe3
                current%particle_energy = en2
                END IF
                impact%part_p = pp4
                impact%particle_energy =  ph2
                impact%cs_times = impact%cs_times + 1.0_num
                current%cs_times = current%cs_times + w2/w1
            ELSE IF (wr < one_m_2eps) THEN
                print*, 'wr < one_m_2eps merge new photon',wr
                CALL weighted_particles_correction_CS(w1 / w2, pk_in(2:4)*m0*c, pp4, pk_in(1)*m0*c**2, pk_sc(1)*m0*c**2, 0.0_num)
                IF (use_compton_recoil) THEN
                current%part_p = pe3
                current%particle_energy = en2
                END IF
                impact%part_p = pp4
                impact%particle_energy =  ph2
                current%cs_times = current%cs_times + 1.0_num
                impact%cs_times = impact%cs_times + w1/w2
            ELSE
                IF (use_compton_recoil) THEN
                current%part_p = pe3
                current%particle_energy = en2
                END IF
            current%cs_times = current%cs_times + 1.0_num
            !for test
            impact%part_p = pp4
            impact%particle_energy =  ph2
            impact%cs_times = impact%cs_times + 1.0_num
            END IF !wr
            !Update momentum and energy

!        IF (ph2 > 1e6*q0) THEN
!            WRITE(*,*) 'wr',wr
!            WRITE(*,*) 'en1,en2,ph1,ph2',en1,en2,ph1,ph2
!            WRITE(*,*) 'p1_norm,pe3',p1_norm,'|',pe3
!            WRITE(*,*) 'pk_in,pk_sc',pk_in,'|',pk_sc
!        END IF
            !delta_energy_ele = delta_energy_ele + ph2 - ph1
            !delta_energy_pho = delta_energy_pho + en2 - en1
        !WRITE(*,*) wr
        END IF !rnd


        current => current%next
        impact => impact%next
#ifdef PREFETCH
        CALL prefetch_particle(current)
        CALL prefetch_particle(impact)
#endif
      END DO ! do 
      !WRITE(*,*) 'Delta Energy electron(photon)', delta_energy_ele, delta_energy_pho ,1 - delta_energy_ele/delta_energy_pho

      ! restore the tail of the lists
      
      CALL append_partlist(species_list(ep_species)%attached_list, append_ele_list)
      CALL append_partlist(species_list(photon_species)%attached_list, append_pho_list)
      NULLIFY(p_list1%tail%next)
      NULLIFY(p_list2%tail%next)
    END IF !icount .AND. jcount

  END SUBROUTINE compton_scatter

  SUBROUTINE compton_scatter_photon(rand_seed,rand_seed_phi, pk_out,pk_o_sc)
      REAL(num), DIMENSION(4), INTENT(IN) :: pk_out
      REAL(num), DIMENSION(4), INTENT(OUT) :: pk_o_sc
      REAL(num), INTENT(IN) :: rand_seed, rand_seed_phi
      REAL(num) :: epo, mu,phi,epo_sc
      REAL(num), DIMENSION(3) :: e1,e2,e3,pk_o_sc3,direction_z
      !Calculate differential scatter cross section
      !e1 = pk_out is the direction of photon
      !e2 is porportional to the e1, and e3 = e1 x e2
      !theta is the angle with respect to e1
      !phi is the angle on the plane e2,e3
      !e1 is the direction of pk_out 
      !e2 is the direction of e1 x (0,0,1)
      !e3 is the direction of e1 x e2
      !pk_out and pk_o_sc is dimensionless paramter = p/m0c
      !mento-corlo method for compton scatter
      epo = pk_out(1)
      mu = find_value_from_table_alt(epo, rand_seed, &
        n_sample_epo, n_sample_mu, log_epo, com_mu, Compton_table)
      mu = LOG10(mu)

      !WRITE(*,*) 'epsilon_o', epo
      !WRITE(*,*) 'rand_seed, mu', rand_seed, mu
      epo_sc = epo/(1.0_num + epo*(1-mu))

      phi = rand_seed_phi*2*pi
      e1 = pk_out(2:4)/pk_out(1)

      direction_z = (/0.0_num,0.0_num,1.0_num/)
      CALL cross(e1,direction_z,e2)
      e2 = e2/SQRT(sum(e2**2))
      CALL cross(e1,e2,e3)
      !WRITE(*,*) 'e1:', e1, SQRT(SUM(e1**2))
      !WRITE(*,*) 'e2:', e2, SQRT(SUM(e2**2))
      !WRITE(*,*) 'e3:', e3, SQRT(SUM(e3**2))
      pk_o_sc3 = epo_sc*(mu*e1 + SQRT(1-mu**2)*COS(phi)*e2 + SQRT(1-mu**2)*SIN(phi)*e3)
      pk_o_sc(2:4) = pk_o_sc3
      pk_o_sc(1) = epo_sc 

  END SUBROUTINE compton_scatter_photon

  SUBROUTINE cross(A,B,AxB)
      REAL(num), DIMENSION(3),INTENT(IN) :: A,B 
      REAL(num), DIMENSION(3),INTENT(OUT) :: AxB
      AxB = (/A(2)*B(3) - A(3)*B(2),&
             A(3)*B(1) - A(1)*B(3),&
             A(1)*B(2) - A(2)*B(1)/)
  END SUBROUTINE


  SUBROUTINE shuffle_particle_list_random_CS(p_list)

    TYPE(particle_list), INTENT(INOUT) :: p_list
    TYPE(particle), POINTER :: particle1, particle2

    INTEGER :: i, idx, swap_idx
    INTEGER :: p_num

    p_num = INT(p_list%count)

    ! Nothing to be done
    IF (p_num <= 2) RETURN

    ! First make sure that the sorting array is large enough
    ! This should happen infrequently
    IF (p_num > coll_sort_array_size) THEN
      DEALLOCATE(coll_sort_array)

      ! make the sort array somewhat larger to avoid frequent deallocation
      ! and reallocation
      coll_sort_array_size = (11 * p_num) / 10 + 10
      ALLOCATE(coll_sort_array(coll_sort_array_size))
    END IF

    ! Copy all the particle pointers into the array and create random
    ! sort indices
    particle1 => p_list%head
    DO i = 1,p_num
      coll_sort_array(i)%particle => particle1
      particle1 => particle1%next
    END DO

    ! Shuffle particles using Durstenfeld's algorithm
    DO idx = p_num,2,-1
      swap_idx = FLOOR(idx * random()) + 1
      particle1 => coll_sort_array(idx)%particle
      coll_sort_array(idx)%particle => coll_sort_array(swap_idx)%particle
      coll_sort_array(swap_idx)%particle => particle1
    END DO

    ! Finally we have to copy back to the list
    ! Do head first
    particle1 => coll_sort_array(1)%particle
    p_list%head => particle1
    NULLIFY(particle1%prev)

    ! Then do all the particle links between head and tail.
    DO i = 2, p_num
      particle2 => particle1
      particle1 => coll_sort_array(i)%particle

      particle1%prev => particle2
      particle2%next => particle1
    END DO

    ! Finally set the tail (at the end of the loop, particle is pointing to
    ! the tail)
    p_list%tail => particle1
    NULLIFY(particle1%next)

  END SUBROUTINE shuffle_particle_list_random_CS

  SUBROUTINE weighted_particles_correction_CS(wtr, p, p_scat, en, en_scat, mass)

    ! This is the correction to the particle according to
    ! Sentoku and Kemp (2008) formulas 21 to 26.
    REAL(num), INTENT(INOUT) :: p_scat(3)
    REAL(num), INTENT(IN) :: p(3)
    REAL(num), INTENT(IN) :: wtr
    REAL(num), INTENT(IN) :: en, en_scat
    REAL(num), INTENT(IN) :: mass

    REAL(num) :: p_after(3)
    REAL(num) :: en_after
    REAL(num) :: gamma_en, gamma_p
    REAL(num) :: delta_p, p_mag, p_trans_mag,p_after_m0c2, en_after_m0c2
    REAL(num) :: c1(3), c2(3), c3(3)
    REAL(num) :: phi

    en_after = (1.0_num - wtr) * en + wtr * en_scat
    p_after  = (1.0_num - wtr) * p  + wtr * p_scat
    p_mag = SQRT(DOT_PRODUCT(p_after, p_after))
    en_after_m0c2 = en_after/(m0*c**2)
    p_after_m0c2 = SQRT((p_mag/m0/c)**2 + (mass/m0)**2)
    !p_after_mev = m0*c**2*SQRT((mass/m0)**2 + p_mag**2*c**2)/1e6/q0
    !gamma_en_mass = en_after / (mass * cc)
    !gamma_p_mass = SQRT(mass + (p_mag / mass / c)**2)

    ! This if-statement is just to take care of possible rounding errors
    ! gamma_p should always be smaller than gamma_en
    !WRITE(*,*) 'p_after,en_after',p_after_m0c2,en_after_m0c2
    IF (p_after_m0c2 < en_after_m0c2) THEN
      ! magnitude of the momentum correction
      delta_p = c * SQRT((mass + en_after/c**2)**2 - (mass**2 + p_mag**2/c**2))
      p_trans_mag = SQRT(p_after(2)**2 + p_after(3)**2)

      CALL new_coords_CS(p_after, c1, c2, c3)

      phi = 2.0_num * pi * random()

      ! Correcting for the loss in energy by adding a perpendicular
      ! momentum correction
      p_scat = p_after + delta_p * (c2 * COS(phi) + c3 * SIN(phi))
    END IF
    !Test 

  END SUBROUTINE weighted_particles_correction_CS

  SUBROUTINE new_coords_CS(vector, c1, c2, c3)

    ! New orthonormal basis vectors for calculating scattering angles:
    ! c1 = v / |v|
    ! c2 = (v x e1) / |v x e1|
    ! c3 = ((v x e1) x v) / |(v x e1) x v|
    ! where e1 is a unit vector parallel to x-axis. I.e., e1 = (1,0,0)
    ! New x-axis, c1, is parallel to the input vector
    ! New y-axis, c2, is perpendicular to c1 and the x-axis
    ! New z-axis ,c3, is perpendicular to c1 and c2
    REAL(num), DIMENSION(3), INTENT(IN) :: vector
    REAL(num), DIMENSION(3), INTENT(OUT) :: c1, c2, c3
    REAL(num) :: vtrans, vmag

    vmag = SQRT(DOT_PRODUCT(vector, vector))
    vtrans = SQRT(vector(2)**2 + vector(3)**2)

    IF (vtrans > c_tiny) THEN
      c1 = vector / vmag
      c2 = (/ 0.0_num, vector(3), -vector(2) /)
      c2 = c2 / vtrans
      c3 = (/ vtrans**2, &
          -(vector(1) * vector(2)), &
          -(vector(1) * vector(3)) /)
      c3 = c3 / (vmag * vtrans)
    ELSE
      c1 = (/ 1.0_num, 0.0_num, 0.0_num /)
      c2 = (/ 0.0_num, 1.0_num, 0.0_num /)
      c3 = (/ 0.0_num, 0.0_num, 1.0_num /)
    END IF

  END SUBROUTINE new_coords_CS

  SUBROUTINE Lorentz_transform(p1,pk_in,pk_out)
      !pk_in,pk_out is the four wavevector/m0/c
    REAL(num), DIMENSION(3), INTENT(IN) :: p1
    REAL(num), DIMENSION(4), INTENT(IN) :: pk_in
    REAL(num), DIMENSION(4), INTENT(OUT) :: pk_out
    REAL(num), DIMENSION(3) :: p1_norm
    REAL(num) :: gam, L11,L12,L13,L14,L22,L23,L24,L33,L34,L44
    !L(p1_norm)
    p1_norm = p1
    gam = SQRT(1 + DOT_PRODUCT(p1_norm,p1_norm)) 
    L11 = gam 
    L12 = -p1_norm(1)
    L13 = -p1_norm(2)
    L14 = -p1_norm(3)
    L22 = 1.0 + p1_norm(1)*p1_norm(1)/(1.0_num+gam)
    L23 = p1_norm(1)*p1_norm(2)/(1.0_num + gam)
    L24 = p1_norm(1)*p1_norm(3)/(1.0_num + gam)
    L33 = 1 + p1_norm(2)*p1_norm(2)/(1.0_num+gam)
    L34 = p1_norm(2)*p1_norm(3)/(1.0_num+gam)
    L44 = 1 + p1_norm(3)*p1_norm(3)/(1.0_num+gam)

    pk_out(1) = L11*pk_in(1) + L12*pk_in(2) + L13*pk_in(3) + L14*pk_in(4)
    pk_out(2) = L12*pk_in(1) + L22*pk_in(2) + L23*pk_in(3) + L24*pk_in(4)
    pk_out(3) = L13*pk_in(1) + L23*pk_in(2) + L33*pk_in(3) + L34*pk_in(4)
    pk_out(4) = L14*pk_in(1) + L24*pk_in(2) + L34*pk_in(3) + L44*pk_in(4)
  END SUBROUTINE 

  SUBROUTINE Compton_cross_section(pk_out,sigma_ij)
      REAL(num), DIMENSION(4), INTENT(IN) :: pk_out
      REAL(num), INTENT(OUT) :: sigma_ij
      REAL(num) :: en
      REAL(num), PARAMETER :: eps = EPSILON(1.0_num)
      en = pk_out(1)
      IF (en < eps) THEN
          sigma_ij = Thomson_cross_section
      ELSE
          sigma_ij = pi*re**2/en*((1.0-2.0/en-2.0/en**2)*LOG(1.0 + 2.0*en) + 0.5 + 4.0/en - 1.0/(2.0*(1.0 + 2.0*en)**2))
      ENDIF
  END SUBROUTINE

  SUBROUTINE Resonance_scatter_section(current, pk, sigma_o_rs, is_rs,sigma_rs_max)
      TYPE(particle), POINTER :: current
      REAL(num), DIMENSION(4), INTENT(IN) :: pk
      REAL(num), INTENT(IN) :: sigma_rs_max
      LOGICAL(num), INTENT(OUT) :: is_rs
      REAL(num), INTENT(OUT) :: sigma_o_rs
      REAL(num) , DIMENSION(3) :: b_local
      REAL(num) ::  part_x,part_y,part_ux,part_uy,part_uz,gamma_rel
      REAL(num) ::  mag_b,en_rs,tau,Gam,mass,sigma_o_rs_tcs,sigma_o_rs_max
      REAL(num) :: sin_theta,cos_theta,en_photon
      REAL(num), PARAMETER :: eps = EPSILON(1.0_num)
      !1. 转到电子静止坐标系下，计算光子的动量能量和坐标系下的磁场大小
      !pk 是电子静止坐标系下的光子动量
      b_local = (/0.0_num,0.0_num,0.0_num/)
      part_x  = current%part_pos(1) - x_grid_min_local
      part_y  = current%part_pos(2) - y_grid_min_local
      part_ux = current%part_p(1) / mc0
      part_uy = current%part_p(2) / mc0
      part_uz = current%part_p(3) / mc0
      gamma_rel = SQRT(part_ux**2 + part_uy**2 + part_uz**2 + 1.0_num)

      !得到电子静止坐标系下的磁场
      b_local = calculate_bfield(part_x, part_y, part_ux, part_uy, part_uz, gamma_rel)
      mag_b = SQRT(DOT_PRODUCT(b_local,b_local))
      mass = m0

      en_photon = pk(1)*mass*c**2
      cos_theta = DOT_PRODUCT(pk(2:4),b_local)/mag_b/pk(1)
      sin_theta = sqrt(1.0_num - cos_theta**2)
      !2. 判断共振条件 en_photon > 1.0e-5_num* hbar*omega_H 若不是共振直接返回。若是，可以继续计算截面

      !IF (sin_theta < eps) THEN
        en_rs = mag_b*q0/mass*h_bar
      !ELSE
      !  en_rs = (1.0_num - SQRT(1.0_num - 2.0_num*mag_b/b_s*sin_theta**2))/sin_theta**2*mass*c**2
      !END IF

      is_rs = .False.
      sigma_o_rs = 0.0_num

      !3. Gam = hbar/tau;tau = 3/4*hbar/alpha_f*(b_s/B_local)**2/me/c**2
      tau = 0.75_num/alpha_f*(b_s/mag_b)**2 
      Gam = 1.0_num/tau
      !4. 给出截面sigma_o_rs = sigma_ts*(en_photon**2)/((en_photon - en_rs)**2 + Gam**2/4) 
      sigma_o_rs_tcs = (en_photon**2)/((en_photon - en_rs)**2 + (Gam*mass*c**2)**2/4.0_num)

      sigma_o_rs_max = en_rs**2*4.0_num/(Gam*mass*c**2)**2
      sigma_o_rs = Thomson_cross_section*MIN(sigma_o_rs_tcs,sigma_rs_max)

!      IF (debug_mod) THEN
!      IF( rank == 0) THEN
!          WRITE(*,*) 'b_local',b_local
!          WRITE(*,*) 'pk',pk
!          WRITE(*,*) 'sin_theta,cos_theta',sin_theta,cos_theta
!          WRITE(*,*) 'mag_b,Gam,sigma_o_rs_tcs,sigma_o_rs_max,tau', mag_b,Gam,sigma_o_rs_tcs,sigma_o_rs_max,tau
!          WRITE(*,*) 'E_pho,en_rs,',en_photon/1e6/q0,en_rs/1e6/q0,mag_b*q0/mass*h_bar/1e6/q0
!      END IF
!      END IF !debug_mod
!
      IF (sigma_o_rs_tcs > 1e-1) THEN
          is_rs = .TRUE.
      END IF

  END SUBROUTINE 

  FUNCTION calculate_bfield(part_x, part_y, part_ux, part_uy, part_uz, &
      gamma_rel)

    REAL(num),DIMENSION(3) :: calculate_bfield
    REAL(num), INTENT(IN) :: part_x, part_y
    REAL(num), INTENT(IN) :: part_ux, part_uy, part_uz, gamma_rel
    REAL(num) :: e_at_part(3), b_at_part(3)
    REAL(num) :: beta_x, beta_y, beta_z, beta_dot_b,moduclip2,moduclip
    REAL(num), DIMENSION(3) :: beta,b2,beta_cross_E


    CALL field_at_particle(part_x, part_y, e_at_part, b_at_part)

    moduclip2 = MAX(part_ux**2 + part_uy**2 + part_uz**2, c_tiny)
    moduclip = SQRT(moduclip2)

    beta_x = part_ux / gamma_rel
    beta_y = part_uy / gamma_rel
    beta_z = part_uz / gamma_rel

    beta_dot_b = (beta_x * b_at_part(1) + beta_y * b_at_part(2) &
        + beta_z * b_at_part(3))

    beta = (/ beta_x, beta_y, beta_z/)
    CALL cross(beta,e_at_part,beta_cross_E)
    !B2 = gamma_rel*(B1-beta \cross E) - gamma_rel**2/(gamma_rel + 1)*beta*(beta \cdot B);
    b2 = gamma_rel*(b_at_part - beta_cross_E/c) - gamma_rel**2/(gamma_rel + 1)*beta*(beta_dot_b)

    calculate_bfield = b2
!    IF (rank == 0) THEN
!        WRITE(*,*) 'Field Transformation Begin'
!        WRITE(*,*) 'b at part', b_at_part
!        WRITE(*,*) 'e at part', e_at_part
!        WRITE(*,*) 'gamma_rel,beta',gamma_rel,beta
!        WRITE(*,*) 'beta_cross_E',beta_cross_E
!        WRITE(*,*) 'beta_dot_b',beta_dot_b
!        WRITE(*,*) 'b2',b2
!        WRITE(*,*) 'Field Transformation End'
!    END IF
  END FUNCTION calculate_bfield

#endif
END MODULE photons
