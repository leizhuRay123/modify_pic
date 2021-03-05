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

MODULE split_particle

  USE boundary

  IMPLICIT NONE

  SAVE

  INTEGER(i8) :: npart_per_cell_min = 5
  LOGICAL :: use_split = .FALSE.

CONTAINS

  SUBROUTINE reorder_particles_to_grid

    INTEGER :: ispecies, ix, iy
    INTEGER :: cell_x, cell_y
    TYPE(particle), POINTER :: current, next
    INTEGER(i8) :: local_count
    INTEGER :: i0, i1, n

    i0 = 1 - ng
    IF (use_field_ionisation) i0 = -ng
    i1 = 1 - i0

    DO ispecies = 1, n_species
      local_count = species_list(ispecies)%attached_list%count
      CALL MPI_ALLREDUCE(local_count, species_list(ispecies)%global_count, &
          1, MPI_INTEGER8, MPI_SUM, comm, errcode)
      ALLOCATE(species_list(ispecies)%secondary_list(i0:nx+i1,i0:ny+i1))
      DO iy = i0, ny + i1
        DO ix = i0, nx + i1
          CALL create_empty_partlist(&
              species_list(ispecies)%secondary_list(ix,iy))
        END DO
      END DO
      current => species_list(ispecies)%attached_list%head
      n = 0
      DO WHILE(ASSOCIATED(current))
        next => current%next
#ifdef PARTICLE_SHAPE_TOPHAT
        cell_x = FLOOR((current%part_pos(1) - x_grid_min_local) / dx) + 1
        cell_y = FLOOR((current%part_pos(2) - y_grid_min_local) / dy) + 1
#else
        cell_x = FLOOR((current%part_pos(1) - x_grid_min_local) / dx + 1.5_num)
        cell_y = FLOOR((current%part_pos(2) - y_grid_min_local) / dy + 1.5_num)
#endif

        CALL remove_particle_from_partlist(&
            species_list(ispecies)%attached_list, current)
        CALL add_particle_to_partlist(&
            species_list(ispecies)%secondary_list(cell_x,cell_y), current)
        current => next
        n = n + 1
      END DO
      !WRITE(*,*) ispecies, n, species_list(ispecies)%global_count,species_list(ispecies)%attached_list%count
      !WRITE(*,*) ispecies,species_list(ispecies)%secondary_list(10,10)%count
    END DO

  END SUBROUTINE reorder_particles_to_grid



  SUBROUTINE reattach_particles_to_mainlist

    INTEGER :: ispecies, ix, iy
    INTEGER :: i0, i1

    i0 = 1 - ng
    IF (use_field_ionisation) i0 = -ng
    i1 = 1 - i0

    DO ispecies = 1, n_species
      DO iy = i0, ny + i1
        DO ix = i0, nx + i1
          CALL append_partlist(species_list(ispecies)%attached_list, &
              species_list(ispecies)%secondary_list(ix,iy))
          !WRITE(*,*) species_list(ispecies)%secondary_list(ix,iy)%head%next%part_pos
        END DO
      END DO
      DEALLOCATE(species_list(ispecies)%secondary_list)
    END DO

    CALL particle_bcs

  END SUBROUTINE reattach_particles_to_mainlist



  SUBROUTINE setup_split_particles

#ifndef PER_SPECIES_WEIGHT
    INTEGER :: ispecies

    use_split = .FALSE.
    DO ispecies = 1, n_species
      IF (species_list(ispecies)%split) THEN
        use_split = .TRUE.
        EXIT
      END IF
    END DO

    use_particle_lists = use_particle_lists .OR. use_split
    IF (use_split) need_random_state = .TRUE.
#endif

  END SUBROUTINE setup_split_particles



  SUBROUTINE split_particles

#ifndef PER_SPECIES_WEIGHT
    INTEGER :: ispecies, ix, iy
    INTEGER(i8) :: count
    TYPE(particle), POINTER :: current, new_particle
    TYPE(particle_list) :: append_list
    REAL(num) :: jitter_x, jitter_y
    INTEGER :: i0, i1

    i0 = 1 - ng
    IF (use_field_ionisation) i0 = -ng
    i1 = 1 - i0

    DO ispecies = 1, n_species
      IF (.NOT. species_list(ispecies)%split) CYCLE
      IF (species_list(ispecies)%npart_max > 0 &
          .AND. species_list(ispecies)%global_count &
          >= species_list(ispecies)%npart_max) CYCLE

      CALL create_empty_partlist(append_list)

      DO iy = i0, ny + i1
        DO ix = i0, nx + i1
          count = species_list(ispecies)%secondary_list(ix,iy)%count
          IF (count > 0 .AND. count <= npart_per_cell_min) THEN
            current => species_list(ispecies)%secondary_list(ix,iy)%head
            DO WHILE(ASSOCIATED(current) .AND. count <= npart_per_cell_min &
                .AND. current%weight >= 1.0_num)
              count = &
                  species_list(ispecies)%secondary_list(ix,iy)%count
              jitter_x = (2 * random() - 1) * 0.25_num * dx
              jitter_y = (2 * random() - 1) * 0.25_num * dy
              current%weight = 0.5_num * current%weight
              ALLOCATE(new_particle)
              new_particle = current
#if defined(PARTICLE_ID) || defined(PARTICLE_ID4)
              new_particle%id = 0
#endif
              new_particle%part_pos(1) = current%part_pos(1) + jitter_x
              new_particle%part_pos(2) = current%part_pos(2) + jitter_y
              CALL add_particle_to_partlist(append_list, new_particle)
#ifdef PARTICLE_DEBUG
              ! If running with particle debugging, specify that this
              ! particle has been split
              new_particle%processor_at_t0 = -1
#endif
              NULLIFY(new_particle)

              current%part_pos(1) = current%part_pos(1) - jitter_x
              current%part_pos(2) = current%part_pos(2) - jitter_y
              current => current%next
            END DO
          END IF
        END DO
      END DO

      CALL append_partlist(species_list(ispecies)%attached_list, append_list)
    END DO
#endif

  END SUBROUTINE split_particles

!Merge Particles
  SUBROUTINE merge_particles
    INTEGER :: ispecies, ix, iy,i0,i1
    !INTEGER :: num_spx, num_spy, num_spz !number of momentum space grid
    INTEGER(i8) :: count, ii ,merge_count,del_number1,del_number2,del_number3
    TYPE(particle), POINTER :: current, new_particle1, new_particle2,next
    TYPE(particle_list) :: append_list
    REAL(num) :: max_px, max_py, max_pz, dpx, dpy, dpz
    REAL(num) :: min_px, min_py, min_pz, time_line

    !momentum grid
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE:: c_weight, c_px, c_py, c_pz, c_energy
    INTEGER(num), DIMENSION(:,:,:),ALLOCATABLE:: c_count 
    !use to check
    !REAL(num), DIMENSION(:,:,:), ALLOCATABLE:: c2_weight, c2_px, c2_py, c2_pz, c2_energy
    !INTEGER(num), DIMENSION(:,:,:),ALLOCATABLE:: c2_count 

    !new structure of photon's infomation for saving
    REAL(num), DIMENSION(:), ALLOCATABLE :: part_px, part_py, part_pz
    REAL(num), DIMENSION(:), ALLOCATABLE :: part_weight, part_energy
    REAL(num), DIMENSION(-1:1) :: gx,gy,gz
    INTEGER(num) :: mgx, mgy, mgz,cell_z,cell_x,cell_y,ipx,ipy,ipz
    REAL(num) :: jitter_x,jitter_y

    
    !variables for calculating
    REAL(num):: var_A, var_B, var_C, var_D
    REAL(num):: xx1, xx2, yy1, yy2
    REAL(num):: diag_x, diag_y, diag_z, diag_mag
    REAL(num):: pvec_x, pvec_y, pvec_z
    INTEGER :: num_spx,num_spy,num_spz,step_merge,tolerent_number,min_per_momentum
    REAL(num) :: start_merge_time,end_merge_time
    REAL(num), DIMENSION(3) :: norm_diag,norm_pvec,norm_diag_cross_pvec,diag,pvec
    REAL(num), DIMENSION(3) :: delta_p
    REAL(num) :: p_mag, pvec_mag, cos_theta,delta_p_mag,mc2,particle_energy,m0c2
    LOGICAL :: isbeam,isdel

    !print the merge info at start time

    DO ispecies = 1, n_species 
        ! XY
      IF ( .NOT. species_list(ispecies)%merging) cycle

      step_merge = species_list(ispecies)%step_merge
      start_merge_time = species_list(ispecies)%start_merge_time
      end_merge_time = species_list(ispecies)%end_merge_time

      IF( MOD(step, step_merge) /= 0 .or. time < start_merge_time .or. time > end_merge_time) cycle
      !Merging paramaters
      num_spx = species_list(ispecies)%num_spx
      num_spy = species_list(ispecies)%num_spy
      num_spz = species_list(ispecies)%num_spz
      tolerent_number = species_list(ispecies)%tolerent_number
      min_per_momentum = species_list(ispecies)%min_per_momentum
      mc2 = species_list(ispecies)%mass*c*c
      m0c2 = m0*c*c


      ALLOCATE(c_count(1:num_spx,1:num_spy,1:num_spz))
      ALLOCATE(c_weight(1:num_spx,1:num_spy,1:num_spz))
      ALLOCATE(c_energy(1:num_spx,1:num_spy,1:num_spz))
      ALLOCATE(c_px(1:num_spx,1:num_spy,1:num_spz))
      ALLOCATE(c_py(1:num_spx,1:num_spy,1:num_spz))
      ALLOCATE(c_pz(1:num_spx,1:num_spy,1:num_spz))

      !ALLOCATE(c2_count(1:num_spx,1:num_spy,1:num_spz))
      !ALLOCATE(c2_weight(1:num_spx,1:num_spy,1:num_spz))
      !ALLOCATE(c2_energy(1:num_spx,1:num_spy,1:num_spz))
      !ALLOCATE(c2_px(1:num_spx,1:num_spy,1:num_spz))
      !ALLOCATE(c2_py(1:num_spx,1:num_spy,1:num_spz))
      !ALLOCATE(c2_pz(1:num_spx,1:num_spy,1:num_spz))

    IF ( species_list(ispecies)%first_call_merge .and. rank .eq. 0) THEN
      PRINT*,"--------------------------------------------------------------------------"
      PRINT*, "Merging info, num_spx,num_spy,num_spz for ispecies",ispecies
      PRINT*, num_spx,num_spy,num_spz 
      PRINT*, "step_merge:",step_merge, "min_per_momentum",min_per_momentum
      PRINT*, "tolerent_number:",tolerent_number
      PRINT*, "Total limitation nx*ny*nz*npx*npy*npz*n",min_per_momentum*num_spx*num_spy*num_spz*nx*ny
      PRINT*,"--------------------------------------------------------------------------"
      species_list(ispecies)%first_call_merge = .FAlSE.
    ENDIF

    If (rank==0) PRINT*,"Call Merge", time ,'For', ispecies

      i0 = 1 - ng
      i1 = 1 - i0
      del_number1 = 0
      del_number2 = 0
      del_number3 = 0
      !CALL create_empty_partlist(append_list)
      DO iy = i0, ny + i1
        DO ix = i0, nx + i1
        count = species_list(ispecies)%secondary_list(ix,iy)%count
        IF (count < species_list(ispecies)%tolerent_number) THEN
            CYCLE
        END IF
        
        ALLOCATE(part_px(1:count-2))
        ALLOCATE(part_py(1:count-2))
        ALLOCATE(part_pz(1:count-2))
        ALLOCATE(part_weight(1:count-2))
        ALLOCATE(part_energy(1:count-2))

        part_px = 0.0_num
        part_py = 0.0_num
        part_pz = 0.0_num
        part_weight = 0.0_num
        part_energy = 0.0_num
        merge_count = 0

        !Max Momentum Space 
        max_px = -1.0e30_num
        max_py = -1.0e30_num
        max_pz = -1.0e30_num
        !xy revised
        min_px = 1.0e30_num
        min_py = 1.0e30_num
        min_pz = 1.0e30_num

        !do not deal with tail or head
        current => species_list(ispecies)%secondary_list(ix,iy)%head
        current => current%next 
        ii = 1 !initial for particle counting
        !WRITE(*,*) 'Begin Find P_max'
        !log scale?
        DO WHILE(ASSOCIATED(current) .and. ASSOCIATED(current%next))
          !copy particle info to new structure for speed 
          ! sort particle which need to merge
          IF (current%particle_energy < species_list(ispecies)%merge_max_energy) THEN 
          part_px(ii) = current%part_p(1)
          part_py(ii) = current%part_p(2)
          part_pz(ii) = current%part_p(3)
          part_weight(ii) = current%weight
          part_energy(ii) = SQRT(DOT_PRODUCT(current%part_p,current%part_p)*c**2 + mc2**2) - mc2
          max_px = max(max_px,part_px(ii))
          max_py = max(max_py,part_py(ii))
          max_pz = max(max_pz,part_pz(ii))
          min_px = min(min_px,part_px(ii))
          min_py = min(min_py,part_py(ii))
          min_pz = min(min_pz,part_pz(ii))
          ii = ii + 1
          END IF !merge_max_energy 
          current => current%next
        ENDDO
        merge_count = ii - 1
        !PRINT*, merge_count, count

        !core merge algorithm 
        !initialize
        c_count = 0
        c_energy = 0.0_num
        c_weight = 0.0_num
        c_px = 0.0_num
        c_py = 0.0_num
        c_pz = 0.0_num
        !XY for check
        !c2_count = 0
        !c2_energy = 0.0_num
        !c2_weight = 0.0_num
        !c2_px = 0.0_num
        !c2_py = 0.0_num
        !c2_pz = 0.0_num

        !XY revise
        dpx = max_px - min_px
        dpy = max_py - min_py
        dpz = max_pz - min_pz
        IF (dpx < EPSILON(1.0_num)*m0*c .AND. dpy < EPSILON(1.0_num)*m0*c &
            .AND. dpz < EPSILON(1.0_num) .AND. mc2 > EPSILON(1.0_num)) THEN
            !This mean a extreme case that a relativstic beam can not be merge for mass > 0 particle
            DEALLOCATE(part_px)
            DEALLOCATE(part_py)
            DEALLOCATE(part_pz)
            DEALLOCATE(part_weight)
            DEALLOCATE(part_energy)
            CYCLE ! for ix,iy
        ELSE
            dpx = MAX(dpx,EPSILON(1.0_num))
            dpy = MAX(dpy,EPSILON(1.0_num))
            dpz = MAX(dpz,EPSILON(1.0_num))
        END IF

        !put particle into the momentum grid
        !WRITE(*,*) 'Begin merge photon in to num_px,y,z'
        Do ii = 1, merge_count !do not deal with head or tail 
          mgx=floor((part_px(ii)-min_px)/dpx*(num_spx-1))+1
          mgy=floor((part_py(ii)-min_py)/dpy*(num_spy-1))+1
          mgz=floor((part_pz(ii)-min_pz)/dpz*(num_spz-1))+1
          c_count(mgx,mgy,mgz) = c_count(mgx,mgy,mgz) + 1
          c_weight(mgx,mgy,mgz) = c_weight(mgx,mgy,mgz)+ part_weight(ii)
          c_energy(mgx,mgy,mgz) = c_energy(mgx,mgy,mgz)+ part_energy(ii)*part_weight(ii)
          c_px(mgx,mgy,mgz) = c_px(mgx,mgy,mgz)+ part_px(ii)*part_weight(ii)
          c_py(mgx,mgy,mgz) = c_py(mgx,mgy,mgz)+ part_py(ii)*part_weight(ii)
          c_pz(mgx,mgy,mgz) = c_pz(mgx,mgy,mgz)+ part_pz(ii)*part_weight(ii)
          !WRITE(*,*) 'part_energy,ii',ii,part_energy(ii)
        ENDDO
        !WRITE(*,*) 'END merge photon in to num_px,y,z'
        !norm

        current => species_list(ispecies)%secondary_list(ix,iy)%head
        current => current%next !do not deal with head or tail
        ii = 1
        DO WHILE(ASSOCIATED(current) .and. associated(current%next))
          !WRITE(*,*) ii
          next => current%next !current maybe deleted
          isdel = .TRUE.
          IF (current%particle_energy < species_list(ispecies)%merge_max_energy) THEN 
          !maybe not be merged 
          !find the coordinate in momentum grid
          mgx=floor((current%part_p(1)-min_px)/dpx*(num_spx-1))+1
          mgy=floor((current%part_p(2)-min_py)/dpy*(num_spy-1))+1
          mgz=floor((current%part_p(3)-min_pz)/dpz*(num_spz-1))+1
          
          IF(c_count(mgx,mgy,mgz) .ge. min_per_momentum .or. & 
             c_count(mgx,mgy,mgz) .eq. -1)  THEN !  ...-1 means have done generating
                                                 ! 第一次合并就是 c_count > min_per_momentum，之后的合并都是c_count = -1
                                                 ! 直接删除粒子
            IF(c_count(mgx,mgy,mgz) .ge. min_per_momentum) THEN !num in momentum grid
            del_number1 = del_number1 - c_count(mgx,mgy,mgz)
              ALLOCATE(new_particle1)
              ALLOCATE(new_particle2)

#if defined(PARTICLE_ID) || defined(PARTICLE_ID4)
               new_particle1%id = 0
               new_particle2%id = 0
#endif
              !calculate momentum of two particles 
              IF (num_spx > 1) THEN
               diag_x = dpx/(num_spx-1)
               diag_y = dpy/(num_spy-1)
               diag_z = dpz/(num_spz-1)
              ELSE
               diag_x = dpx
               diag_y = dpy
               diag_z = dpz
              END IF
               !XY 没有除weight?
               pvec_x = c_px(mgx,mgy,mgz)/c_weight(mgx,mgy,mgz)
               pvec_y = c_py(mgx,mgy,mgz)/c_weight(mgx,mgy,mgz)
               pvec_z = c_pz(mgx,mgy,mgz)/c_weight(mgx,mgy,mgz)
               diag = (/diag_x,diag_y,diag_z/)
               norm_diag = diag/SQRT(SUM(diag**2)) !norm
               pvec = (/pvec_x,pvec_y,pvec_z/)
               norm_pvec = pvec /SQRT(SUM(pvec**2)) !norm
               
               norm_diag_cross_pvec = 0.0_num
               CALL cross(norm_diag,norm_pvec,norm_diag_cross_pvec)

               !正方体四个对角线，总能找到与pvec不平行的
               IF (SQRT(SUM(norm_diag_cross_pvec**2)) < EPSILON(1.0_num)) THEN 
                 !PRINT*,"pvec_vec colinear with diag_vec"
                 IF(norm_diag(3) > EPSILON(1.0_num)) THEN
                    diag_z = -diag_z
                 ELSEIF(norm_diag(2) > EPSILON(1.0_num)) THEN
                    diag_y = -diag_y
                 ELSE
                    diag_x = -diag_x
                 ENDIF
               ENDIF

               norm_diag = diag/SQRT(SUM(diag**2)) !norm
               CALL cross(norm_diag,norm_pvec,norm_diag_cross_pvec)
               !还是平行，一定是一条线的beam case
               isbeam = .False.

               !WRITE(*,*) 'norm_diag_cross_pvec',SQRT(SUM(norm_diag_cross_pvec)**2)
               IF (SQRT(SUM(norm_diag_cross_pvec**2)) < 0.001) THEN
                    isbeam = .TRUE.
               END IF

               !1.pt//d 
               !1-1.如果p_perp_mag < EPSILON(1.0_num)  beam case,对于光子直接合成一个粒子
               !1-2否则随便选一个phi值
               !2.pt与d构成平面,要选平面上与pt成phi角的矢量作为两个粒子中的一个动量。
               IF (isbeam ) THEN ! beam case 合成一个粒子
                   IF (mc2/m0c2 < EPSILON(1.0_num)) THEN ! 除光子外，例如正负电子即使动量平行也不能合并。
                 new_particle1%part_pos = current%part_pos
                 new_particle1%weight= c_weight(mgx,mgy,mgz) !equal weight
                 new_particle1%part_p= pvec
                 new_particle1%particle_energy= &
                        c_energy(mgx,mgy,mgz)/c_weight(mgx,mgy,mgz)
                 new_particle1%optical_depth = reset_optical_depth()
                 !CALL add_particle_to_partlist(append_list, new_particle1)
                 del_number3 = del_number3 + 1
                 CALL add_particle_to_partlist(&
                      species_list(ispecies)%attached_list,new_particle1)
                 c_count(mgx,mgy,mgz) = -1 ! -1 means have done generating
                    
                   WRITE(*,*) 'beam case: merge'
                   ELSE 
                   WRITE(*,*) 'beam case: can not merge'
                        c_count(mgx,mgy,mgz) = -2 ! -1 means have done generating
                        isdel = .False.
                   END IF
               ELSE

               !XY
               !如果diag与pvec可以成一个平面，先求平面的法线n，pvec绕法线n旋转theta角,theta 由cos_theta = pvec/pa决定
               !pa和能量相关满足c*pa = c_energy/c_weight/2
               !位置错开
                jitter_x = 0.0_num
                !(2.0_num * random() - 1.0_num) * 0.25_num * dx
                jitter_y = 0.0_num
                !(2.0_num * random() - 1.0_num) * 0.25_num * dy
               particle_energy = c_energy(mgx,mgy,mgz)/c_weight(mgx,mgy,mgz) 
               !为了避免数值截断，这里需要先/m0c归一化一下

               p_mag = SQRT((particle_energy/m0c2 + mc2/m0c2)**2 - mc2/m0c2)
               pvec_mag = SQRT(DOT_PRODUCT(pvec/m0/c,pvec/m0/c))
               delta_p_mag = SQRT(p_mag**2 - pvec_mag**2)*m0*c
               !WRITE(*,*) 'c_energy',c_energy(mgx,mgy,mgz)

               CALL merge_part_momentum(norm_pvec,norm_diag,delta_p)

               !WRITE(*,*) 'particle_energy,m0c2,mc2',particle_energy,m0c2,mc2
               !WRITE(*,*) 'p_mag,pvec_mag,pvec,delta_p_mag',p_mag,pvec_mag,delta_p_mag
               !WRITE(*,*) 'pvec',pvec
               !WRITE(*,*) 'norm_diag',norm_diag
               !WRITE(*,*) 'delta_p',delta_p
               !WRITE(*,*) 'weight',c_weight(mgx,mgy,mgz)
               !new particle 1
               !WRITE(*,*) 'BBB, no beam case'
               new_particle1%part_pos(1) = current%part_pos(1) + jitter_x
               new_particle1%part_pos(2) = current%part_pos(2) + jitter_y

               new_particle1%weight= c_weight(mgx,mgy,mgz)/2.0_num !equal weight
               new_particle1%part_p = pvec + delta_p_mag*delta_p
               new_particle1%particle_energy= particle_energy
               new_particle1%optical_depth = reset_optical_depth()

               !new particle 2
               new_particle2%part_pos(1) = current%part_pos(1) - jitter_x
               new_particle2%part_pos(2) = current%part_pos(2) - jitter_y

               new_particle2%weight= c_weight(mgx,mgy,mgz)/2.0_num !equal weight
               new_particle2%part_p = pvec - delta_p_mag*delta_p
               new_particle2%particle_energy= particle_energy
               new_particle2%optical_depth = reset_optical_depth()
!------------------add particle
               !CALL add_particle_to_partlist(append_list, new_particle1)
               !CALL add_particle_to_partlist(append_list, new_particle2)
                 !WRITE(*,*) 'ispecies',ispecies
                 !WRITE(*,*) 'no beam case: part_pos1',new_particle1%part_pos,new_particle1%part_p
                 !WRITE(*,*) 'no beam case: part_pos2',new_particle2%part_pos,new_particle2%part_p
               del_number3 = del_number3 + 2
               CALL add_particle_to_partlist(&
                      species_list(ispecies)%attached_list,new_particle1)
               CALL add_particle_to_partlist(&
                      species_list(ispecies)%attached_list,new_particle2)
               !WRITE(*,*) 'non-beam case: part_p1',new_particle1%part_p
               !WRITE(*,*) 'part_p2', new_particle2%part_p
               ENDIF !a beam or not

               !check momentum conservation
               !WRITE(*,*) 'Pt:', c_px(mgx,mgy,mgz),c_py(mgx,mgy,mgz),c_pz(mgx,mgy,mgz)
               !WRITE(*,*) 'P1:', new_particle1%part_p*new_particle2%weight 
               !WRITE(*,*) 'P2:', new_particle2%part_p*new_particle2%weight 
               !WRITE(*,*) 'P1+P2:', new_particle1%part_p*new_particle1%weight+ new_particle2%part_p*new_particle2%weight 
               !WRITE(*,*) 'new_particle',new_particle1%part_p,new_particle2%part_p
               NULLIFY(new_particle1)
               NULLIFY(new_particle2)
               c_count(mgx,mgy,mgz) = -1 ! -1 means have done generating
            ENDIF! first merge,之后直接删除粒子

!------------------delete particle
            !c2_px(mgx,mgy,mgz) = c2_px(mgx,mgy,mgz) + current%part_p(1)*current%weight
            !c2_py(mgx,mgy,mgz) = c2_py(mgx,mgy,mgz) + current%part_p(2)*current%weight
            !c2_pz(mgx,mgy,mgz) = c2_pz(mgx,mgy,mgz) + current%part_p(3)*current%weight
            !c2_count(mgx,mgy,mgz) = c2_count(mgx,mgy,mgz) + 1

            IF (isdel) THEN !
            del_number2 = del_number2 - 1
            CALL remove_particle_from_partlist(&
                  species_list(ispecies)%attached_list, current)
            END IF
            !WRITE(*,*) 'END remove_particle_from_partlist'
          ENDIF ! c_count > min_per_momentum
        ENDIF !merge_max_energy
        current => next
        ii = ii + 1
        ENDDO
            !write(*,*) 'ix,iy',ix,iy
            !write(*,*) 'px:',c_px(1,1,1),c2_px(1,1,1)
            !write(*,*) 'py:',c_py(1,1,1),c2_py(1,1,1)
            !write(*,*) 'pz:',c_pz(1,1,1),c2_pz(1,1,1)

        DEALLOCATE(part_px)
        DEALLOCATE(part_py)
        DEALLOCATE(part_pz)
        DEALLOCATE(part_weight)
        DEALLOCATE(part_energy)
        ENDDO ! ix
      ENDDO ! iy
      !CALL append_partlist(species_list(ispecies)%attached_list,append_list)
      !WRITE(*,*) 'END append_partlist'
      DEALLOCATE(c_count)
      DEALLOCATE(c_weight)
      DEALLOCATE(c_energy)
      DEALLOCATE(c_px)
      DEALLOCATE(c_py)
      DEALLOCATE(c_pz)

      !WRITE(*,*) 'ispecies,Want to Delete, Real delete, Add particle Numbers:', &
      !           ispecies,del_number1,del_number2,del_number3
    ENDDO ! nspecies
        
  END SUBROUTINE merge_particles

  FUNCTION reset_optical_depth()

    ! Resets optical depth of particle
    REAL(num) :: reset_optical_depth
    REAL(num) :: p_tau

    p_tau = random()
    reset_optical_depth = -LOG(1.0_num - p_tau)

  END FUNCTION reset_optical_depth

  SUBROUTINE cross(A,B,AxB)
      REAL(num), DIMENSION(3),INTENT(IN) :: A,B 
      REAL(num), DIMENSION(3),INTENT(OUT) :: AxB
      AxB = (/A(2)*B(3) - A(3)*B(2),&
             A(3)*B(1) - A(1)*B(3),&
             A(1)*B(2) - A(2)*B(1)/)
  END SUBROUTINE

  SUBROUTINE merge_part_momentum(a,b,delta_p)
      REAL(num), DIMENSION(3),INTENT(IN) :: a,b
      REAL(num), DIMENSION(3),INTENT(OUT) :: delta_p
      REAL(num), DIMENSION(3) :: n
      REAL(num) :: beta,alpha,cos_phi,AA,BB,CC
      !先求法线方向n = a x b 
      !在求\delta p 方向\delta p  = n x a 
      CALL cross(a,b,n)
      n = n/SQRT(SUM(n**2))
      CALL cross(n,a,delta_p)
      delta_p = delta_p/SQRT(SUM(delta_p**2))
  END SUBROUTINE

END MODULE split_particle
