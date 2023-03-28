!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
MODULE qs_tddfpt2_soc
   USE admm_types,                      ONLY: admm_type
   USE bibliography,                    ONLY: Grimme2013,&
                                              Grimme2016,&
                                              Iannuzzi2005,&
                                              Bussy2021a,&
                                              cite_reference
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_p_type,&
                                              gto_basis_set_type
   USE cp_array_utils,                  ONLY: cp_1d_i_p_type,&
                                              cp_1d_r_p_type,&
                                              cp_2d_i_p_type,&
                                              cp_2d_r_p_type,&
                                              cp_3d_r_p_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_cfm_diag,                     ONLY: cp_cfm_heevd
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_get_info,&
                                              cp_cfm_get_submatrix,&
                                              cp_cfm_to_cfm,&!_columns,&
                                              cp_cfm_release,&
                                              cp_cfm_type,&
                                              cp_fm_to_cfm,&
                                              cp_cfm_to_fm
   USE cp_control_types,                ONLY: dft_control_type,&
                                              qs_control_type,&
                                              tddfpt2_control_type
   USE cp_dbcsr_operations,             ONLY: dbcsr_allocate_matrix_set,&
                                              cp_dbcsr_sm_fm_multiply, &
                                              dbcsr_deallocate_matrix_set, &
                                              copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr
   USE cp_fm_basic_linalg,              ONLY: cp_fm_scale, &
                                              cp_fm_upper_to_full, &
                                              cp_fm_transpose
   USE cp_fm_pool_types,                ONLY: fm_pool_create_fm
   USE cp_fm_struct,                    ONLY: cp_fm_struct_type, &
                                              cp_fm_struct_create, &
                                              cp_fm_struct_release
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_get_diag,&
                                              cp_fm_to_fm,&
                                              cp_fm_to_fm_submat,&
                                              cp_fm_type, &
                                              cp_fm_get_submatrix, &
                                              cp_fm_set_element,&
                                              cp_fm_write_formatted
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_add_iter_level,&
                                              cp_iterate,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_unit_nr,&
                                              cp_rm_iter_level
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: dbcsr_copy,&
                                              dbcsr_create,&
                                              dbcsr_multiply, &
                                              dbcsr_type,&
                                              dbcsr_p_type,&
                                              dbcsr_get_block_p,&
                                              dbcsr_set, &
                                              dbcsr_print, &
                                              dbcsr_release, &
                                              dbcsr_iterator_start, &
                                              dbcsr_iterator_stop, &
                                              dbcsr_iterator_type, &
                                              dbcsr_iterator_blocks_left,&
                                              dbcsr_iterator_next_block, &
                                              dbcsr_distribution_new,&
                                              dbcsr_distribution_get,&
                                              dbcsr_distribution_release,&
                                              dbcsr_distribution_type,&
                                              dbcsr_type_no_symmetry,&
                                              dbcsr_reserve_all_blocks,&
                                              dbcsr_get_info
   USE exstates_types,                  ONLY: excited_energy_type
   USE header,                          ONLY: tddfpt_header, &
                                              tddfpt_soc_header
   USE hfx_types,                       ONLY: compare_hfx_sections
   USE input_constants,                 ONLY: tddfpt_dipole_velocity,&
                                              tddfpt_kernel_full,&
                                              tddfpt_kernel_none,&
                                              tddfpt_kernel_stda,&
                                              tddfpt_soc_dipol_len,&
                                              tddfpt_soc_dipol_vel
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp,&
                                              default_string_length
   USE lri_environment_methods,         ONLY: lri_print_stat
   USE lri_environment_types,           ONLY: lri_density_release,&
                                              lri_env_release
   USE machine,                         ONLY: m_flush
   USE mathlib,                         ONLY: get_diag
   USE mulliken,                        ONLY: compute_charges
   USE orbital_pointers,                ONLY: indco,&
                                              indso,&
                                              nco,&
                                              ncoset,&
                                              nso,&
                                              nsoset
   USE physcon,                         ONLY: evolt
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE qs_2nd_kernel_ao,                ONLY: build_dm_response
   USE particle_methods,                ONLY: get_particle_set
   USE particle_types,                  ONLY: particle_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              get_qs_kind_set,&
                                              qs_kind_type

   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kernel_methods,               ONLY: create_kernel_env
   USE qs_kernel_types,                 ONLY: full_kernel_env_type,&
                                              kernel_env_type,&
                                              release_kernel_env
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_type
   USE qs_operators_ao,                 ONLY: p_xyz_ao,&
                                              rRc_xyz_ao
   USE qs_scf_methods,                  ONLY: eigensolver
   USE qs_scf_types,                    ONLY: qs_scf_env_type
   USE qs_tddfpt2_densities,            ONLY: tddfpt_construct_aux_fit_density,&
                                              tddfpt_construct_ground_state_orb_density
   USE qs_tddfpt2_eigensolver,          ONLY: tddfpt_davidson_solver,&
                                              tddfpt_orthogonalize_psi1_psi0,&
                                              tddfpt_orthonormalize_psi1_psi1
   USE qs_tddfpt2_lri_utils,            ONLY: tddfpt2_lri_init
   USE qs_tddfpt2_properties,           ONLY: tddfpt_dipole_operator,&
                                              tddfpt_print_excitation_analysis,&
                                              tddfpt_print_nto_analysis,&
                                              tddfpt_print_summary
   USE qs_tddfpt2_restart,              ONLY: tddfpt_read_restart,&
                                              tddfpt_write_newtonx_output,&
                                              tddfpt_write_restart
   USE qs_tddfpt2_stda_types,           ONLY: allocate_stda_env,&
                                              deallocate_stda_env,&
                                              stda_env_type,&
                                              stda_init_param
   USE qs_tddfpt2_stda_utils,           ONLY: stda_init_matrices
   USE qs_tddfpt2_subgroups,            ONLY: tddfpt_sub_env_init,&
                                              tddfpt_sub_env_release,&
                                              tddfpt_subgroup_env_type
   USE qs_tddfpt2_types,                ONLY: stda_create_work_matrices,&
                                              tddfpt_create_work_matrices,&
                                              tddfpt_ground_state_mos,&
                                              tddfpt_release_work_matrices,&
                                              tddfpt_work_matrices
   USE qs_tddfpt2_utils,                ONLY: tddfpt_guess_vectors,&
                                              tddfpt_init_mos,&
                                              tddfpt_oecorr,&
                                              tddfpt_release_ground_state_mos
   USE string_utilities,                ONLY: integer_to_string
   USE xc_write_output,                 ONLY: xc_write
!!Under construction!!
   USE qs_tddfpt2_soc_types,            ONLY: soc_atom_env_type, &
                                              soc_atom_create, &
                                              soc_atom_release, &
                                              soc_env_type, &
                                              soc_env_create, &
                                              soc_env_release
   USE xas_tdp_atom,                   ONLY: integrate_soc_atoms, &
                                             calculate_density_coeffs,&
                                             compute_sphi_so, &
                                             truncate_radial_grid
   !USE xas_tdp_utils,                  ONLY: rcs_amew_soc_elements
   USE qs_tddfpt2_soc_utils,            ONLY: rcs_amew_soc_elements
   USE lebedev,                        ONLY: deallocate_lebedev_grids,&
                                             get_number_of_lebedev_grid,&
                                             init_lebedev_grids,&
                                             lebedev_grid
   USE qs_grid_atom,                   ONLY: allocate_grid_atom,&
                                             create_grid_atom,&
                                             grid_atom_type
   USE qs_interactions,                ONLY: init_interaction_radii_orb_basis
   USE qs_harmonics_atom,              ONLY: allocate_harmonics_atom,&
                                             create_harmonics_atom,&
                                             get_maxl_CG,&
                                             get_none0_cg_list,&
                                             harmonics_atom_type
   USE spherical_harmonics,            ONLY: clebsch_gordon,&
                                             clebsch_gordon_deallocate,&
                                             clebsch_gordon_init
   USE memory_utilities,               ONLY: reallocate
   USE string_utilities,                ONLY: integer_to_string
   USE cp_fm_struct,                    ONLY: cp_fm_struct_get



#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_tddfpt2_soc'

   LOGICAL, PARAMETER, PRIVATE          :: debug_this_module = .FALSE.

   INTEGER, PARAMETER, PRIVATE          :: nderivs = 3
   INTEGER, PARAMETER, PRIVATE          :: maxspins = 2

   !A helper type for SOC
   TYPE dbcsr_soc_package_type
      TYPE(dbcsr_type), POINTER     :: dbcsr_sg
      TYPE(dbcsr_type), POINTER     :: dbcsr_tp
      TYPE(dbcsr_type), POINTER     :: dbcsr_sc
      TYPE(dbcsr_type), POINTER     :: dbcsr_sf
      TYPE(dbcsr_type), POINTER     :: dbcsr_prod
      TYPE(dbcsr_type), POINTER     :: dbcsr_ovlp
      TYPE(dbcsr_type), POINTER     :: dbcsr_tmp
      TYPE(dbcsr_type), POINTER     :: dbcsr_work
   END TYPE dbcsr_soc_package_type

   PUBLIC :: tddfpt_soc, inititialize_soc

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief Perform TDDFPT-SOC calculation.
!> \param qs_env  Quickstep environment
!> \param calc_forces ...
!> \par History
!>    * 02.2023 created [Jan-Robert Vogt]
!> \note Based on tddfpt2_methods and xas_tdp_utils.
! **************************************************************************************************

   SUBROUTINE tddfpt_soc(qs_env, evals_a,&
                         evals_b, evects_a, evects_b)!, gs)!, gs_mos,&
                         !tddfpt_control)

      TYPE(qs_environment_type), POINTER               :: qs_env
      CHARACTER(LEN=*), PARAMETER                      :: routineN="tddfpt_soc"
      INTEGER                                          :: handle!,&
                                                        ! log_unit,&
                                                        ! istate
      TYPE(dft_control_type), POINTER                  :: dft_control
      REAL(kind=dp), INTENT(in), DIMENSION(:), TARGET  :: evals_a, &
                                                          evals_b
      TYPE(cp_fm_type),DIMENSION(:,:)                  :: evects_a, &
                                                          evects_b
      !TYPE(cp_fm_type), POINTER                        :: gs
      LOGICAL                                          :: do_os
      TYPE(cp_logger_type), POINTER                    :: logger
   !   TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
   !      POINTER                                         :: gs_mos
   !   TYPE(tddfpt2_control_type), POINTER                :: tddfpt_control

      NULLIFY(dft_control)

      CALL timeset(routineN, handle)
      logger => cp_get_default_logger()

      CALL get_qs_env(qs_env, dft_control=dft_control)
      do_os = dft_control%uks .OR. dft_control%roks

      IF (do_os) THEN
         CALL tddfpt_soc_os(qs_env, evals_a,&
                             evals_b, evects_a, evects_b)
      ELSE
         CALL tddfpt_soc_rcs(qs_env, evals_a,&
                             evals_b, evects_a, evects_b)!,gs)!,gs_mos,&
                             !tddfpt_control)
      END IF

      CALL timestop(handle)
      
   END SUBROUTINE

   SUBROUTINE tddfpt_soc_rcs(qs_env, evals_sing,&
                         evals_trip, evects_sing, evects_trip)!, gs)!,gs_mos,&
                       !  tddfpt_control)

      TYPE(qs_environment_type), POINTER        :: qs_env
      CHARACTER(LEN=*), PARAMETER               :: routineN ="tddfpt_soc_rcs"
      INTEGER                                   :: handle, &
                                                   i,&
                                                   ntot,&
                                                   isg,& 
                                                   itp,&
                                                   nsg,&
                                                   group,&
                                                   ntp,&
                                                   nex,&
                                                   iex,&
                                                   npcols,&
                                                   nprows,&
                                                   dipol_form       
      TYPE(cp_blacs_env_type), POINTER          :: blacs_env
      TYPE(cp_para_env_type), POINTER           :: para_env
      TYPE(cp_fm_type)                          :: real_fm,& 
                                                   img_fm
      TYPE(cp_fm_struct_type), POINTER          :: full_struct 
      TYPE(cp_logger_type), POINTER             :: logger
      TYPE(section_vals_type), POINTER          :: tddfpt_print_section
      !!Under Contsruction!!
      INTEGER, POINTER                          :: nmo, &
                                                   nao
      INTEGER, DIMENSION(:), POINTER            :: row_dist,&
                                                   row_dist_new,&
                                                   col_dist,&
                                                   col_blk_size,&
                                                   row_blk_size
      INTEGER, DIMENSION(:, :), POINTER         :: pgrid
      TYPE(cp_cfm_type)                         :: evecs_cfm,&
                                                   hami_cfm
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER :: matrix_s
      TYPE(cp_fm_type), POINTER                 :: gs_coeffs,&
      !                                             gs,&
                                                   sg_coeffs, &
                                                   tp_coeffs
      TYPE(cp_fm_type)                          :: prod_fm, &
                                                   tmp_fm, &
                                                   work_fm, &
                                                   gstp_fm, &
                                                   vec_soc_x, &
                                                   vec_soc_y, &
                                                   vec_soc_z
      TYPE(cp_fm_type), DIMENSION(:,:)          :: evects_sing, &
                                                   evects_trip
      TYPE(cp_fm_struct_type), POINTER          :: vec_struct,& 
                                                   work_struct,&
                                                   prod_struct,&
                                                   gstp_struct
      TYPE(mo_set_type), DIMENSION(:), POINTER  :: mos
      TYPE(dbcsr_distribution_type), POINTER    :: coeffs_dist,&
                                                   dbcsr_dist,&
                                                   prod_dist
      TYPE(soc_atom_env_type), POINTER          :: soc_atom_env
      TYPE(soc_env_type), TARGET                :: soc_env
      TYPE(dbcsr_type), POINTER                 :: dbcsr_ovlp,& 
                                                   dbcsr_prod,& 
                                                   dbcsr_sg, &
                                                   dbcsr_tmp, & 
                                                   dbcsr_tp, &
                                                   dbcsr_work
      TYPE(dbcsr_type), POINTER                 :: orb_soc_x, & 
                                                   orb_soc_y, &
                                                   orb_soc_z
      REAL(dp)                                  :: sqrt2, &
                                                   soc_gst, &
                                                   eps_filter
      TYPE(section_vals_type), POINTER          :: soc_section
      REAL(dp), ALLOCATABLE, DIMENSION(:)       :: diag, &
                                                   tmp_evals
      REAL(kind=dp), INTENT(in),&
         DIMENSION(:), TARGET                   :: evals_sing, &
                                                   evals_trip
      REAL(kind=dp), ALLOCATABLE,&
                     DIMENSION(:,:)             :: mo_soc_x, &
                                                   mo_soc_y, &
                                                   mo_soc_z, &
                                                   gstp_block
      TYPE(dbcsr_soc_package_type)              :: dbcsr_soc_package
    !  !!DEBUG
    !   TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: dipole_op_mos_occ
    !   TYPE(cp_fm_type), DIMENSION(:, :), ALLOCATABLE     :: soc_evects
    !   INTEGER                                            :: log_unit
    !    TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
    !      POINTER                                         :: gs_mos
    !   TYPE(tddfpt2_control_type), POINTER                :: tddfpt_control

      CALL timeset(routineN, handle)

       NULLIFY (logger, full_struct, pgrid, row_dist, row_dist_new)
       NULLIFY (col_dist, col_blk_size, row_blk_size, matrix_s)
       NULLIFY (gstp_struct, tddfpt_print_section, soc_section)
      logger => cp_get_default_logger()

      tddfpt_print_section => section_vals_get_subs_vals(qs_env%input, &
                                       "PROPERTIES%TDDFPT%PRINT")
      soc_section => section_vals_get_subs_vals(qs_env%input, &
                                       "PROPERTIES%TDDFPT%SOC")
      CALL section_vals_val_get(soc_section, "EPS_FILTER", r_val=eps_filter)
      CALL section_vals_val_get(soc_section, "DIPOLE_FORM",i_val=dipol_form)
      nsg = SIZE(evals_sing)  !! Number of excited singlet states
      ntp = SIZE(evals_trip) !! Number of excited triplett states
      nex = nsg               !! Number of excited states of each multiplicity
      ntot = 1 + nsg + 3 * ntp !! number of (GS + S + T^-1 + T^0 + T^1)


!! Initzialize Working envirment
      CALL inititialize_soc(qs_env, soc_atom_env, soc_env, & 
                            evects_sing,&
                            evects_trip,&
                            dipol_form)

      blacs_env => soc_env%blacs_env
      para_env => soc_env%para_env
      mos => soc_env%mos
      nmo => soc_env%nmo
      nao => soc_env%nao

      CALL get_qs_env(qs_env, matrix_s= matrix_s)

      gs_coeffs => soc_env%gs_coeff

       !! This will create the H^SOC, but in an atomic basis !!
      CALL integrate_soc_atoms(soc_env%orb_soc, qs_env, soc_atom_env=soc_atom_env)
     
     !! DEBUG::
      DO i = 1,3
      CALL print_soc(tddfpt_print_section, "SOC_PRINT", logger,".aosoc", &
                     dbcsr_matrix=soc_env%orb_soc(i) ,&
                       !dbcsr2fm=.TRUE., blacs_env=blacs_env, &
                        header="<PHI|H^SOC|PHI> in atomic orital basis",&
                        starter="<PHI|H^SOC|PHI>")
      END DO
      CALL soc_atom_release(soc_atom_env)

      !! Point at H^SOC and MOs for better readablity of the code
      sg_coeffs => soc_env%a_coeff
      tp_coeffs => soc_env%b_coeff

      soc_env%evals_a => evals_sing
      soc_env%evals_b => evals_trip

      orb_soc_x => soc_env%orb_soc(1)%matrix
      orb_soc_y => soc_env%orb_soc(2)%matrix
      orb_soc_z => soc_env%orb_soc(3)%matrix

   !   !!DEBUG::
   !   CALL print_soc(tddfpt_print_section, "SOC_PRINT", logger, ".sgwf", &
   !                  fm_matrix=sg_coeffs ,&
   !                  header="Excitation Vectors Singels",&
   !                  starter="Singels")
   !   CALL print_soc(tddfpt_print_section, "SOC_PRINT", logger, ".tpwf", &
   !                  fm_matrix=tp_coeffs ,&
   !                     header="Excitation Vectors Tripletts")



      !! Create a matrix-structure, which links all states in this calculation
      CALL cp_fm_struct_create(full_struct,&
                               context=blacs_env,&
                               para_env=para_env, &
                               nrow_global=ntot,&
                               ncol_global=ntot)
      CALL cp_fm_create(real_fm, full_struct)
      CALL cp_fm_create(img_fm, full_struct)

      !  Put the excitation energies on the diagonal of the real matrix
      DO isg = 1, nsg
         CALL cp_fm_set_element(real_fm, 1 + isg, 1 + isg, evals_sing(isg))
      END DO
      DO itp = 1, ntp
         ! first T^-1, then T^0, then T^+1
         CALL cp_fm_set_element(real_fm,&
            1 + itp + nsg, 1 + itp + nsg,evals_trip(itp))
         CALL cp_fm_set_element(real_fm,&
            1 + itp + ntp +  nsg, 1 + itp + ntp+ nsg,evals_trip(itp))
         CALL cp_fm_set_element(real_fm,&
            1 + itp + 2*ntp + nsg, 1 + itp + 2*ntp + nsg, evals_trip(itp))
      END DO


   !!Create the dbcsr structures for this calculations
      CALL get_qs_env(qs_env, dbcsr_dist=dbcsr_dist)
      CALL dbcsr_distribution_get(dbcsr_dist,&
                                  group=group,&
                                  row_dist=row_dist,&
                                  pgrid=pgrid, &
                                  npcols=npcols,&
                                  nprows=nprows)

      ALLOCATE (col_dist(nex), row_dist_new(nex))
      ALLOCATE (coeffs_dist, prod_dist)
     
      DO iex = 1,nex
         col_dist(iex) = MODULO(npcols-iex, npcols)
         row_dist_new(iex) = MODULO(nprows-iex, nprows)
      END DO

      CALL dbcsr_distribution_new(coeffs_dist,&
                                  group=group,&
                                  pgrid=pgrid,&
                                  row_dist=row_dist, &
                                  col_dist=col_dist)
      CALL dbcsr_distribution_new(prod_dist,&
                                  group=group,&
                                  pgrid=pgrid,&
                                  row_dist=row_dist_new, &
                                  col_dist=col_dist)

     !! Create the matrices
      ALLOCATE (col_blk_size(nex))
      col_blk_size = nmo
      CALL dbcsr_get_info(matrix_s(1)%matrix, row_blk_size=row_blk_size)


      ALLOCATE (dbcsr_sg, dbcsr_tp, dbcsr_work, dbcsr_ovlp, dbcsr_tmp,&
                dbcsr_prod)

      CALL dbcsr_create(matrix=dbcsr_sg,&
                        name="SINGLETS",&
                        matrix_type=dbcsr_type_no_symmetry, &
                        dist=coeffs_dist,&
                        row_blk_size=row_blk_size,&
                        col_blk_size=col_blk_size)
      CALL dbcsr_create(matrix=dbcsr_tp,&
                        name="TRIPLETS",&
                        matrix_type=dbcsr_type_no_symmetry, &
                        dist=coeffs_dist,&
                        row_blk_size=row_blk_size,&
                        col_blk_size=col_blk_size)
      CALL dbcsr_create(matrix=dbcsr_work,&
                        name="WORK",&
                        matrix_type=dbcsr_type_no_symmetry, &
                        dist=coeffs_dist,&
                        row_blk_size=row_blk_size,&
                        col_blk_size=col_blk_size)
      CALL dbcsr_create(matrix=dbcsr_prod,&
                        name="PROD",&
                        matrix_type=dbcsr_type_no_symmetry, &
                        dist=prod_dist,&
                        row_blk_size=col_blk_size,&
                        col_blk_size=col_blk_size)
      CALL dbcsr_create(matrix=dbcsr_ovlp,&
                        name="OVLP", matrix_type=dbcsr_type_no_symmetry, &
                        dist=prod_dist,&
                        row_blk_size=col_blk_size,&
                        col_blk_size=col_blk_size)

      col_blk_size = 1 
      CALL dbcsr_create(matrix=dbcsr_tmp,&
                        name="TMP",&
                        matrix_type=dbcsr_type_no_symmetry, &
                        dist=prod_dist,&
                        row_blk_size=col_blk_size,&
                        col_blk_size=col_blk_size)
      CALL dbcsr_reserve_all_blocks(dbcsr_tmp)

      dbcsr_soc_package%dbcsr_sg => dbcsr_sg
      dbcsr_soc_package%dbcsr_tp => dbcsr_tp
      dbcsr_soc_package%dbcsr_work => dbcsr_work
      dbcsr_soc_package%dbcsr_ovlp => dbcsr_ovlp
      dbcsr_soc_package%dbcsr_prod => dbcsr_prod
      dbcsr_soc_package%dbcsr_tmp => dbcsr_tmp


      !Filling the coeffs matrices by copying from the stored fms
      CALL copy_fm_to_dbcsr(sg_coeffs, dbcsr_sg)
      CALL copy_fm_to_dbcsr(tp_coeffs, dbcsr_tp)

!Create the work and helper fms
      CALL cp_fm_get_info(gs_coeffs, matrix_struct=vec_struct)
      CALL cp_fm_struct_create(prod_struct,&
                               context=blacs_env,&
                               para_env=para_env, &
                               nrow_global=nmo,&
                               ncol_global=nmo)
      CALL cp_fm_create(prod_fm, prod_struct)
      CALL cp_fm_create(vec_soc_x, vec_struct)
      CALL cp_fm_create(vec_soc_y, vec_struct)
      CALL cp_fm_create(vec_soc_z, vec_struct)
      CALL cp_fm_struct_create(work_struct,&
                               context=blacs_env,&
                               para_env=para_env, &
                               nrow_global=nex,&
                               ncol_global=nex)
      CALL cp_fm_create(work_fm, work_struct)
      CALL cp_fm_create(tmp_fm, work_struct)
      ALLOCATE (diag(nmo))

!! Begin with the precomputation 
      sqrt2 = SQRT(2.0_dp)

      !! Precompute the <phi_i^0|H^SOC|phi_j^0> matrix elements
      ALLOCATE (mo_soc_x(nmo,nmo), mo_soc_y(nmo,nmo), mo_soc_z(nmo,nmo))

      CALL cp_dbcsr_sm_fm_multiply(orb_soc_x,gs_coeffs,vec_soc_x,ncol=nmo)
      CALL parallel_gemm('T', 'N', nmo,&
                         nmo, nao, 1.0_dp,&
                         gs_coeffs,vec_soc_x,&
                         0.0_dp, prod_fm)
      CALL cp_fm_get_submatrix(prod_fm, mo_soc_x)

      CALL cp_dbcsr_sm_fm_multiply(orb_soc_y, gs_coeffs, vec_soc_y, ncol=nmo)
      CALL parallel_gemm('T', 'N', nmo,&
                         nmo, nao, 1.0_dp,&
                         gs_coeffs, vec_soc_y,&
                         0.0_dp, prod_fm)
      CALL cp_fm_get_submatrix(prod_fm, mo_soc_y)

      CALL cp_dbcsr_sm_fm_multiply(orb_soc_z, gs_coeffs, vec_soc_z, ncol=nmo)
      CALL parallel_gemm('T', 'N', nmo,&
                         nmo, nao, 1.0_dp,&
                         gs_coeffs, vec_soc_z,&
                         0.0_dp, prod_fm)
      CALL cp_fm_get_submatrix(prod_fm, mo_soc_z)

!  Only have SOC between singlet-triplet triplet-triplet and ground_state-triplet, the resulting
!  matrix is Hermitian i.e. the real part is symmetric and the imaginary part is anti-symmetric.
!  Can only fill upper half

      !Start with the ground state/triplet SOC, SOC*gs_coeffs already computed above
      !note: we are computing <0|H|T>, but have SOC*gs_coeffs instead of gs_coeffs*SOC in store. 
      !Since the SOC Hamiltonian is anti-symmetric, a - signs pops up in the gemms below

      CALL cp_fm_struct_create(gstp_struct,&
                               context=blacs_env,&
                               para_env=para_env, &
                               nrow_global=ntp*nmo,&
                               ncol_global=nmo)
      CALL cp_fm_create(gstp_fm, gstp_struct)
      ALLOCATE (gstp_block(nmo,nmo))

      !gs-triplet with Ms=+-1, imaginary part
      CALL parallel_gemm('T', 'N', nmo*ntp,&
                         nmo, nao, -1.0_dp,&
                         tp_coeffs, vec_soc_x,&
                         0.0_dp, gstp_fm)

      DO itp = 1, ntp
         CALL cp_fm_get_submatrix(fm=gstp_fm, target_m=gstp_block, &
                                  start_row=(itp-1)*nmo+1, start_col=1, &
                                  n_rows=nmo, n_cols=nmo)
         diag(:) = get_diag(gstp_block)
         soc_gst = SUM(diag)
         !!  <0|H_x|T^-1>
         CALL cp_fm_set_element(img_fm, 1, 1 + nsg + itp, -1_dp*soc_gst)
         !! <0|H_x|T^+1>
         CALL cp_fm_set_element(img_fm, 1, 1 + nsg + 2*ntp + itp, soc_gst)
      END DO

      !gs-triplet with Ms=+-1, real part
      CALL parallel_gemm('T', 'N', nmo*ntp,&
                         nmo, nao, -1.0_dp,&
                         tp_coeffs, vec_soc_y,&
                         0.0_dp, gstp_fm)

      DO itp = 1, ntp
         CALL cp_fm_get_submatrix(fm=gstp_fm,&
                                  target_m=gstp_block,&
                                  start_row=(itp - 1)*nmo + 1, &
                                  start_col=1,&
                                  n_rows=nmo,&
                                  n_cols=nmo)
         diag(:) = get_diag(gstp_block)
         soc_gst = SUM(diag)
         ! <0|H_y|T^-1>
         CALL cp_fm_set_element(real_fm, 1,&
                                1 + nsg + itp, -1.0_dp*soc_gst) 
         ! <0|H_y|T^+1>
         CALL cp_fm_set_element(real_fm, 1,&
                                1 + nsg + 2*ntp + itp, -1.0_dp*soc_gst) 
      END DO

      !gs-triplet with Ms=0, purely imaginary
      CALL parallel_gemm('T', 'N', nmo*ntp,&
                         nmo, nao, -1.0_dp,&
                         tp_coeffs, vec_soc_z,&
                         0.0_dp, gstp_fm)

      DO itp = 1, ntp
         CALL cp_fm_get_submatrix(fm=gstp_fm,&
                                  target_m=gstp_block,&
                                  start_row=(itp - 1)*nmo + 1, &
                                  start_col=1,&
                                  n_rows=nmo,&
                                  n_cols=nmo)
         diag(:) = get_diag(gstp_block)
         soc_gst = sqrt2*SUM(diag)
         CALL cp_fm_set_element(img_fm, 1, 1 + nsg + ntp + itp, soc_gst)
      END DO

      !gs clean-up
      CALL cp_fm_release(prod_fm)
      CALL cp_fm_release(vec_soc_x)
      CALL cp_fm_release(vec_soc_y)
      CALL cp_fm_release(vec_soc_z)
      CALL cp_fm_release(gstp_fm)
      CALL cp_fm_struct_release(gstp_struct)
      CALL cp_fm_struct_release(prod_struct)
      DEALLOCATE (gstp_block)

      !Now do the singlet-triplet SOC
      !start by computing the singlet-triplet overlap
      CALL dbcsr_multiply('N', 'N', 1.0_dp,&
                           matrix_s(1)%matrix,&
                           dbcsr_tp, 0.0_dp, &
                          dbcsr_work,&
                          filter_eps=eps_filter)
      CALL dbcsr_multiply('T', 'N', 1.0_dp,&
                          dbcsr_sg,&
                          dbcsr_work,&
                          0.0_dp, dbcsr_ovlp,&
                          filter_eps=eps_filter)

      !singlet-triplet with Ms=+-1, imaginary part
      CALL dbcsr_multiply('N', 'N', 1.0_dp,&
                           orb_soc_x,&
                           dbcsr_tp, 0.0_dp,&
                           dbcsr_work,&
                           filter_eps=eps_filter)
      CALL dbcsr_multiply('T', 'N', 1.0_dp,&
                          dbcsr_sg,&
                          dbcsr_work,0.0_dp,&
                          dbcsr_prod,&
                          filter_eps=eps_filter)

      CALL rcs_amew_soc_elements(dbcsr_tmp,&
                                 dbcsr_prod,&
                                 dbcsr_ovlp,&
                                 mo_soc_x, &
                                 pref_trace=-1.0_dp,&
                                 pref_overall=-0.5_dp*sqrt2)

      !<S|H_x|T^-1>
      CALL copy_dbcsr_to_fm(dbcsr_tmp, tmp_fm)
      CALL cp_fm_to_fm_submat(msource=tmp_fm,&
                              mtarget=img_fm,&
                              nrow=nex,&
                              ncol=nex, &
                              s_firstrow=1,&
                              s_firstcol=1,&
                              t_firstrow=2, &
                              t_firstcol=1 + nsg + 1)

      !<S|H_x|T^+1> takes a minus sign
      CALL cp_fm_scale(-1.0_dp, tmp_fm)
      CALL cp_fm_to_fm_submat(msource=tmp_fm,&
                              mtarget=img_fm,&
                              nrow=nex,&
                              ncol=nex, &
                              s_firstrow=1,&
                              s_firstcol=1,&
                              t_firstrow=2, &
                              t_firstcol=1 + nsg + 2*ntp + 1)

      !singlet-triplet with Ms=+-1, real part
      CALL dbcsr_multiply('N', 'N', 1.0_dp,&
                          orb_soc_y, dbcsr_tp,&
                          0.0_dp, dbcsr_work,&
                          filter_eps=eps_filter)
      CALL dbcsr_multiply('T', 'N', 1.0_dp,&
                          dbcsr_sg, dbcsr_work,&
                          0.0_dp, dbcsr_prod,&
                          filter_eps=eps_filter)

      CALL rcs_amew_soc_elements(dbcsr_tmp,&
                                 dbcsr_prod,&
                                 dbcsr_ovlp,&
                                 mo_soc_y, &
                                 pref_trace=-1.0_dp,&
                                 pref_overall=-0.5_dp*sqrt2)

      !<S|H_y|T^-1>
      CALL copy_dbcsr_to_fm(dbcsr_tmp, tmp_fm)
      CALL cp_fm_to_fm_submat(msource=tmp_fm,&
                              mtarget=real_fm,&
                              nrow=nex,&
                              ncol=nex, &
                              s_firstrow=1,&
                              s_firstcol=1,&
                              t_firstrow=2, &
                              t_firstcol=1 + nsg + 1)

      !<S|H_y|T^+1>
      CALL cp_fm_to_fm_submat(msource=tmp_fm,&
                              mtarget=real_fm,&
                              nrow=nex,&
                              ncol=nex, &
                              s_firstrow=1,&
                              s_firstcol=1,&
                              t_firstrow=2, &
                              t_firstcol=1 + nsg + 2*ntp + 1)

      !singlet-triplet with Ms=0, purely imaginary
      CALL dbcsr_multiply('N', 'N', 1.0_dp,&
                          orb_soc_z, dbcsr_tp,&
                          0.0_dp, dbcsr_work,&
                          filter_eps=eps_filter)
      CALL dbcsr_multiply('T', 'N', 1.0_dp,&
                          dbcsr_sg, dbcsr_work,&
                          0.0_dp, dbcsr_prod,&
                          filter_eps=eps_filter)

      CALL rcs_amew_soc_elements(dbcsr_tmp,&
                                 dbcsr_prod,&
                                 dbcsr_ovlp,&
                                 mo_soc_z, &
                                 pref_trace=-1.0_dp,&
                                 pref_overall=1.0_dp)

      !<S|H_z|T^0>
      CALL copy_dbcsr_to_fm(dbcsr_tmp, tmp_fm)
      CALL cp_fm_to_fm_submat(msource=tmp_fm,&
                              mtarget=img_fm,&
                              nrow=nex,&
                              ncol=nex, &
                              s_firstrow=1,&
                              s_firstcol=1,&
                              t_firstrow=2, &
                              t_firstcol=1 + nsg + ntp + 1)

      !Now the triplet-triplet SOC
      !start by computing the overlap
      CALL dbcsr_multiply('N', 'N', 1.0_dp,&
                          matrix_s(1)%matrix,&
                          dbcsr_tp, 0.0_dp, &
                          dbcsr_work,&
                          filter_eps=eps_filter)
      CALL dbcsr_multiply('T', 'N', 1.0_dp,&
                          dbcsr_tp, dbcsr_work,&
                          0.0_dp, dbcsr_ovlp,&
                          filter_eps=eps_filter)

      !Ms=0 to Ms=+-1 SOC, imaginary part
      CALL dbcsr_multiply('N', 'N', 1.0_dp,&
                          orb_soc_x, dbcsr_tp,&
                          0.0_dp, dbcsr_work,&
                          filter_eps=eps_filter)
      CALL dbcsr_multiply('T', 'N', 1.0_dp,&
                          dbcsr_tp, dbcsr_work,&
                          0.0_dp, dbcsr_prod,&
                          filter_eps=eps_filter)

      CALL rcs_amew_soc_elements(dbcsr_tmp,&
                                 dbcsr_prod,&
                                 dbcsr_ovlp,&
                                 mo_soc_x, &
                                 pref_trace=1.0_dp,&
                                 pref_overall=-0.5_dp*sqrt2)

      !<T^0|H_x|T^+1>
      CALL copy_dbcsr_to_fm(dbcsr_tmp, tmp_fm)
      CALL cp_fm_to_fm_submat(msource=tmp_fm,&
                              mtarget=img_fm,&
                              nrow=nex,&
                              ncol=nex, &
                              s_firstrow=1,&
                              s_firstcol=1,&
                              t_firstrow=1 + nsg + ntp + 1, &
                              t_firstcol=1 + nsg + 2*ntp + 1)

      !<T^-1|H_x|T^0>, takes a minus sign and a transpose (because computed <T^0|H_x|T^-1>)
      CALL cp_fm_transpose(tmp_fm, work_fm)
      CALL cp_fm_scale(-1.0_dp, work_fm)
      CALL cp_fm_to_fm_submat(msource=work_fm,&
                              mtarget=img_fm,&
                              nrow=nex,&
                              ncol=nex, &
                              s_firstrow=1,&
                              s_firstcol=1,&
                              t_firstrow=1 + nsg + 1, &
                              t_firstcol=1 + nsg + ntp + 1)

      !Ms=0 to Ms=+-1 SOC, real part
      CALL dbcsr_multiply('N', 'N', 1.0_dp,&
                          orb_soc_y, dbcsr_tp,&
                          0.0_dp, dbcsr_work,&
                          filter_eps=eps_filter)
      CALL dbcsr_multiply('T', 'N', 1.0_dp,&
                          dbcsr_tp, dbcsr_work,&
                          0.0_dp, dbcsr_prod,&
                          filter_eps=eps_filter)

      CALL rcs_amew_soc_elements(dbcsr_tmp,&
                                 dbcsr_prod,&
                                 dbcsr_ovlp,&
                                 mo_soc_y, &
                                 pref_trace=1.0_dp,&
                                 pref_overall=0.5_dp*sqrt2)

      !<T^0|H_y|T^+1>
      CALL copy_dbcsr_to_fm(dbcsr_tmp, tmp_fm)
      CALL cp_fm_to_fm_submat(msource=tmp_fm,&
                              mtarget=real_fm,&
                              nrow=nex,&
                              ncol=nex, &
                              s_firstrow=1,&
                              s_firstcol=1, t_firstrow=1 + nsg + ntp + 1, &
                              t_firstcol=1 + nsg + 2*ntp + 1)

      !<T^-1|H_y|T^0>, takes a minus sign and a transpose
      CALL cp_fm_transpose(tmp_fm, work_fm)
      CALL cp_fm_scale(-1.0_dp, work_fm)
      CALL cp_fm_to_fm_submat(msource=work_fm,&
                              mtarget=real_fm,&
                              nrow=nex,&
                              ncol=nex, &
                              s_firstrow=1,&
                              s_firstcol=1, t_firstrow=1 + nsg + 1, &
                              t_firstcol=1 + nsg + ntp + 1)

      !Ms=1 to Ms=1 and Ms=-1 to Ms=-1 SOC, purely imaginary
      CALL dbcsr_multiply('N', 'N', 1.0_dp,&
                          orb_soc_z, dbcsr_tp,&
                          0.0_dp, dbcsr_work,&
                          filter_eps=eps_filter)
      CALL dbcsr_multiply('T', 'N', 1.0_dp,&
                          dbcsr_tp, dbcsr_work,&
                          0.0_dp, dbcsr_prod,&
                          filter_eps=eps_filter)

      CALL rcs_amew_soc_elements(dbcsr_tmp,&
                                 dbcsr_prod,&
                                 dbcsr_ovlp,&
                                 mo_soc_z, &
                                 pref_trace=1.0_dp,&
                                 pref_overall=1.0_dp)

      !<T^+1|H_z|T^+1>
      CALL copy_dbcsr_to_fm(dbcsr_tmp, tmp_fm)
      CALL cp_fm_to_fm_submat(msource=tmp_fm,&
                              mtarget=img_fm,&
                              nrow=nex,&
                              ncol=nex, &
                              s_firstrow=1,&
                              s_firstcol=1,&
                              t_firstrow=1 + nsg + 2*ntp + 1, &
                              t_firstcol=1 + nsg + 2*ntp + 1)

      !<T^-1|H_z|T^-1>, takes a minus sign
      CALL cp_fm_scale(-1.0_dp, tmp_fm)
      CALL cp_fm_to_fm_submat(msource=tmp_fm,&
                              mtarget=img_fm,&
                              nrow=nex,&
                              ncol=nex, &
                              s_firstrow=1,&
                              s_firstcol=1,&
                              t_firstrow=1 + nsg + 1, &
                              t_firstcol=1 + nsg + 1)

!  Intermediate clean-up
      CALL cp_fm_struct_release(work_struct)
      CALL cp_fm_release(work_fm)
      CALL cp_fm_release(tmp_fm)
      DEALLOCATE (diag, mo_soc_x, mo_soc_y, mo_soc_z)

!  Set-up the complex hermitian perturbation matrix
      CALL cp_cfm_create(hami_cfm, full_struct)
      CALL cp_fm_to_cfm(real_fm, img_fm, hami_cfm)

      !!DEBUG::
      CALL print_soc(tddfpt_print_section, "SOC_PRINT", logger,".compsoc",&
                        fm_matrix=real_fm, header="SOC terms real",&
                        starter="fm real matrix of SOC terms")

      !!DEBUG::
      CALL print_soc(tddfpt_print_section, "SOC_PRINT", logger,".compsoc",&
                        fm_matrix=img_fm, header="SOC terms img",&
                        starter="fm img matrix of SOC terms")

      CALL cp_fm_release(real_fm)
      CALL cp_fm_release(img_fm)

!  Diagonalize the Hamiltonian
      ALLOCATE (tmp_evals(ntot))
      CALL cp_cfm_create(evecs_cfm, full_struct)
      CALL cp_cfm_heevd(hami_cfm, evecs_cfm, tmp_evals)

      !  Adjust the energies so the GS has zero, and store in the donor_state (without the GS)
      ALLOCATE (soc_env%soc_evals(ntot - 1))
      soc_env%soc_evals(:) = tmp_evals(2:ntot) - tmp_evals(1)

!  Compute the dipole oscillator strengths
      CALL soc_dipol(qs_env, dbcsr_soc_package, soc_env,&
                     evecs_cfm, eps_filter, dipol_form,.FALSE.)

      !!DEBUG::
      CALL print_soc(tddfpt_print_section, "SOC_PRINT",logger,".spectrum",&
                        real_vec=soc_env%soc_evals,&
                        real_vec2=soc_env%soc_osc,&
                        real1=27.21139613179_dp,&
                        header="SOC OSCILLATOR STRENGTH")

      !!DEBUG::
  !    DO i = 1, SIZE(soc_env%dipmat)
  !    CALL print_soc(tddfpt_print_section, "SOC_PRINT", logger, ".debug", &
  !                   dbcsr_matrix=soc_env%dipmat(i),&
  !                   !dbcsr2fm=.TRUE.,&
  !                   !blacs_env=blacs_env,&
  !                   header="DEBUG DIPOL OPERATOR")
  !    END DO

  !    !!DEBUG
  !    log_unit = cp_print_key_unit_nr(logger,tddfpt_print_section, &
  !                                           "PROGRAM_BANNER", &
  !                                           extension=".tddfptLog")

  !    CALL tddfpt_dipole_operator(dipole_op_mos_occ, tddfpt_control,gs_mos, qs_env)

  !    !!DEBUG::
  !    DO i = 1, SIZE(dipole_op_mos_occ,1)
  !       DO itp = 1, SIZE(dipole_op_mos_occ,2)
  !       CALL print_soc(tddfpt_print_section, "SOC_PRINT", logger, ".debug2", &
  !                   fm_matrix=dipole_op_mos_occ(i,itp),&
  !                   header="DEBUG DIPOL OPERATOR")
  !       END DO
  !    END DO

  !    CALL split_cfm_into_fm(evecs_cfm,&
  !                           soc_evects,&
  !                           ntot,&
  !                           evects_sing(1,1)%matrix_struct) 

  !    soc_env%soc_osc = 0.0_dp
  !    CALL tddfpt_print_summary(log_unit,&
  !                              soc_evects,&
  !                              soc_env%soc_evals,&
  !                              soc_env%soc_osc,&
  !                              1, &
  !                              dipole_op_mos_occ,&
  !                              dipol_form)
  ! 
  !   !!END DEBUG

!! clean up
     CALL soc_env_release(soc_env)
     CALL cp_fm_struct_release(full_struct)
     CALL cp_cfm_release(hami_cfm)
     CALL cp_cfm_release(evecs_cfm)
     CALL dbcsr_distribution_release(coeffs_dist)
     CALL dbcsr_distribution_release(prod_dist)
     CALL dbcsr_release(dbcsr_sg)
     CALL dbcsr_release(dbcsr_tp)
     CALL dbcsr_release(dbcsr_prod)
     CALL dbcsr_release(dbcsr_ovlp)
     CALL dbcsr_release(dbcsr_tmp)
     CALL dbcsr_release(dbcsr_work)
     DEALLOCATE (dbcsr_work, dbcsr_prod, dbcsr_ovlp, dbcsr_tmp)
     DEALLOCATE (coeffs_dist, prod_dist, col_dist, col_blk_size, row_dist_new)
     DEALLOCATE (dbcsr_sg, dbcsr_tp, tmp_evals)

     CALL timestop(handle)

   END SUBROUTINE tddfpt_soc_rcs

   SUBROUTINE tddfpt_soc_os(qs_env, evals_sc,&
                         evals_sf, evects_sc, evects_sf)

      TYPE(qs_environment_type), POINTER             :: qs_env
      CHARACTER(LEN=*), PARAMETER                    :: routineN ="tddfpt_soc_os"
      INTEGER                                        :: handle, &
                                                    !    i,&
                                                        ntot,&
                                                        isc,&
                                                        isf,&
                                                        nsc,&
                                                        nsf,&
                                                        group
      TYPE(cp_blacs_env_type), POINTER               :: blacs_env
      TYPE(cp_para_env_type), POINTER                :: para_env
      TYPE(cp_fm_type)                               :: real_fm,&
                                                        img_fm
      TYPE(cp_fm_struct_type), POINTER               :: full_struct
      TYPE(dft_control_type), POINTER                :: dft_control
      TYPE(cp_logger_type), POINTER                  :: logger
      TYPE(section_vals_type), POINTER               :: tddfpt_print_section
      !!Under Contsruction!!
      INTEGER                                        :: nex, iex,&
                                                        npcols,&
                                                        nso, &
                                                        nprows,&
                                                        dipol_form, &
                                                        homo,&
                                                        tas(2),&
                                                        tbs(2)
      INTEGER, POINTER                               :: nmo, &
                                                        nao
      INTEGER, DIMENSION(:), POINTER                 :: row_dist,&
                                                        row_dist_new,&
                                                        col_dist,&
                                                        col_blk_size,&
                                                        row_blk_size
      INTEGER, DIMENSION(:, :), POINTER              :: pgrid
      TYPE(cp_cfm_type)                              :: evecs_cfm,pert_cfm
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER      :: matrix_s
      TYPE(cp_fm_type), POINTER                      :: gs_coeffs,&
                                                        mo_coeff, &
                                                        sc_coeffs, &
                                                        sf_coeffs
      TYPE(cp_fm_type)                               :: work_fm, &
                                                        prod_work,&
                                                        vec_work, &
                                                        gsex_fm, &
                                                        vec_soc_x, &
                                                        vec_soc_y, &
                                                        vec_soc_z
      TYPE(cp_fm_type),DIMENSION(:,:)                :: evects_sc, &
                                                        evects_sf
      TYPE(cp_fm_struct_type), POINTER               :: vec_struct,&
                                                        work_struct,&
                                                        prod_struct,&
                                                        gsex_struct
      TYPE(mo_set_type), DIMENSION(:), POINTER       :: mos
      TYPE(dbcsr_distribution_type), POINTER         :: coeffs_dist, dbcsr_dist, prod_dist
      TYPE(soc_atom_env_type), POINTER               :: soc_atom_env
      TYPE(soc_env_type), TARGET                     :: soc_env
      TYPE(dbcsr_type), POINTER                      :: dbcsr_ovlp,&
                                                        dbcsr_prod,&
                                                        dbcsr_sc, &
                                                        dbcsr_tmp, &
                                                        dbcsr_sf, &
                                                        dbcsr_work
      TYPE(dbcsr_type), POINTER                      :: orb_soc_x, &
                                                        orb_soc_y, &
                                                        orb_soc_z
      REAL(dp)                                       :: soc, &
                                                        eps_filter, &
                                                        gs_sum
      TYPE(section_vals_type), POINTER               :: soc_section
      REAL(dp), ALLOCATABLE, DIMENSION(:)            :: diag, &
                                                        tmp_evals
      REAL(kind=dp), INTENT(in), DIMENSION(:),&
         TARGET                                      :: evals_sc, &
                                                        evals_sf
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:,:)     :: mo_soc_x, &
                                                        mo_soc_y, &
                                                        mo_soc_z, &
                                                        gsex_block
      TYPE(dbcsr_soc_package_type)                   :: dbcsr_soc_package
      LOGICAL, POINTER                               :: do_uks,&
                                                        do_roks

      CALL timeset(routineN, handle)

       NULLIFY (logger, full_struct, pgrid, row_dist, row_dist_new)
       NULLIFY (col_dist, col_blk_size, row_blk_size, matrix_s)
       NULLIFY (gsex_struct, dft_control, do_uks, do_roks)
      logger => cp_get_default_logger()

      tddfpt_print_section => section_vals_get_subs_vals(qs_env%input, &
                                       "PROPERTIES%TDDFPT%PRINT")
      soc_section => section_vals_get_subs_vals(qs_env%input, &
                                       "PROPERTIES%TDDFPT%SOC")
      CALL section_vals_val_get(soc_section, "EPS_FILTER", r_val=eps_filter)
      CALL section_vals_val_get(soc_section, "DIPOLE_FORM",i_val=dipol_form)
      nsc = SIZE(evals_sc)  !! Number of excited spin-conserving states
      nsf =  SIZE(evals_sf) !! Number of excited spin-flip states
      nex = nsc               !! Number of excited states of each multiplicity
      ntot = 1 + nsc + nsf !! number of (GS + S + T^-1 + T^0 + T^1)

      CALL inititialize_soc(qs_env, soc_atom_env, soc_env,&
                               evects_sc, & 
                               evects_sf,&
                               dipol_form)


      blacs_env => soc_env%blacs_env
      para_env => soc_env%para_env
      mos => soc_env%mos
      nmo => soc_env%nmo
      nso = 2*nmo                        !!number spin orbitals
      nao => soc_env%nao
      CALL get_qs_env(qs_env, matrix_s= matrix_s)

      !! This will create the H^SOC, but in an atomic basis !!
      !CALL integrate_soc_atoms(soc_env%orb_soc, qs_env, soc_atom_env=soc_atom_env)
      CALL soc_atom_release(soc_atom_env)

      CALL get_qs_env(qs_env, dft_control=dft_control)
      do_uks => dft_control%uks 
      do_roks => dft_control%roks

      IF (do_uks) gs_coeffs => soc_env%gs_coeff
      IF (do_roks) THEN
         CALL cp_fm_struct_create(vec_struct, context=blacs_env, para_env=para_env, &
                                  nrow_global=nao, ncol_global=nso)
         ALLOCATE (gs_coeffs)
         CALL cp_fm_create(gs_coeffs, vec_struct)

         ! only alpha donor MOs are stored, need to copy them intoboth the alpha and the beta slot
         CALL cp_fm_to_fm_submat(msource=soc_env%gs_coeff,mtarget=gs_coeffs, nrow=nao, &
                                 ncol=nmo, s_firstrow=1, s_firstcol=1, t_firstrow=1, &
                                 t_firstcol=1)
         CALL cp_fm_to_fm_submat(msource=soc_env%gs_coeff,mtarget=gs_coeffs, nrow=nao, &
                                 ncol=nmo, s_firstrow=1, s_firstcol=1, t_firstrow=1, &
                                 t_firstcol=nmo + 1)

         CALL cp_fm_struct_release(vec_struct)
      END IF

      sc_coeffs => soc_env%a_coeff
      sf_coeffs => soc_env%b_coeff

      soc_env%evals_a => evals_sc
      soc_env%evals_b => evals_sf

      orb_soc_x => soc_env%orb_soc(1)%matrix
      orb_soc_y => soc_env%orb_soc(2)%matrix
      orb_soc_z => soc_env%orb_soc(3)%matrix


      !! Create a matrix-structure, which links all states in this calculation
      CALL cp_fm_struct_create(full_struct, context=blacs_env, para_env=para_env, &
                               nrow_global=ntot, ncol_global=ntot)
      CALL cp_fm_create(real_fm, full_struct)
      CALL cp_fm_create(img_fm, full_struct)

      !  Put the excitation energies on the diagonal of the real matrix
      DO isc = 1, nsc
         CALL cp_fm_set_element(real_fm, 1 + isc, 1 + isc,evals_sc(isc))
      END DO
      DO isf = 1, nsf
         ! first T^-1, then T^0, then T^+1
         CALL cp_fm_set_element(real_fm, 1 + isf + nsc, 1 + isf +nsc,evals_sf(isf))
      END DO

      ! Create the dbcsr machinery
      CALL get_qs_env(qs_env, dbcsr_dist=dbcsr_dist)
      CALL dbcsr_distribution_get(dbcsr_dist, group=group, row_dist=row_dist, pgrid=pgrid, &
                                  npcols=npcols, nprows=nprows)
      ALLOCATE (col_dist(nex), row_dist_new(nex))
      DO iex = 1, nex
         col_dist(iex) = MODULO(npcols - iex, npcols)
         row_dist_new(iex) = MODULO(nprows - iex, nprows)
      END DO
      ALLOCATE (coeffs_dist, prod_dist)
      CALL dbcsr_distribution_new(coeffs_dist, group=group, pgrid=pgrid, row_dist=row_dist, &
                                  col_dist=col_dist)
      CALL dbcsr_distribution_new(prod_dist, group=group, pgrid=pgrid, row_dist=row_dist_new, &
                                  col_dist=col_dist)

      !Create the matrices
      ALLOCATE (col_blk_size(nex))
      col_blk_size = nso
      CALL dbcsr_get_info(matrix_s(1)%matrix, row_blk_size=row_blk_size)

      ALLOCATE (dbcsr_sc, dbcsr_sf, dbcsr_work, dbcsr_ovlp, dbcsr_tmp, dbcsr_prod)
      CALL dbcsr_create(matrix=dbcsr_sc, name="SPIN CONS", matrix_type=dbcsr_type_no_symmetry, &
                        dist=coeffs_dist, row_blk_size=row_blk_size, col_blk_size=col_blk_size)
      CALL dbcsr_create(matrix=dbcsr_sf, name="SPIN FLIP", matrix_type=dbcsr_type_no_symmetry, &
                        dist=coeffs_dist, row_blk_size=row_blk_size, col_blk_size=col_blk_size)
      CALL dbcsr_create(matrix=dbcsr_work, name="WORK", matrix_type=dbcsr_type_no_symmetry, &
                        dist=coeffs_dist, row_blk_size=row_blk_size, col_blk_size=col_blk_size)
      CALL dbcsr_create(matrix=dbcsr_prod, name="PROD", matrix_type=dbcsr_type_no_symmetry, &
                        dist=prod_dist, row_blk_size=col_blk_size, col_blk_size=col_blk_size)
      CALL dbcsr_create(matrix=dbcsr_ovlp, name="OVLP", matrix_type=dbcsr_type_no_symmetry, &
                        dist=prod_dist, row_blk_size=col_blk_size, col_blk_size=col_blk_size)

      col_blk_size = 1
      CALL dbcsr_create(matrix=dbcsr_tmp, name="TMP", matrix_type=dbcsr_type_no_symmetry, &
                        dist=prod_dist, row_blk_size=col_blk_size, col_blk_size=col_blk_size)
      CALL dbcsr_reserve_all_blocks(dbcsr_tmp)

      dbcsr_soc_package%dbcsr_sc => dbcsr_sc
      dbcsr_soc_package%dbcsr_sf => dbcsr_sf
      dbcsr_soc_package%dbcsr_work => dbcsr_work
      dbcsr_soc_package%dbcsr_ovlp => dbcsr_ovlp
      dbcsr_soc_package%dbcsr_prod => dbcsr_prod
      dbcsr_soc_package%dbcsr_tmp => dbcsr_tmp

      !Filling the coeffs matrices by copying from the stored fms
      CALL copy_fm_to_dbcsr(sc_coeffs, dbcsr_sc)
      CALL copy_fm_to_dbcsr(sf_coeffs, dbcsr_sf)

      !!09.03

      ! Precompute what we can before looping over excited states.
      ! Need to compute the scalar: sum_i sum_sigma <phi^0_i,sigma|SOC|phi^0_i,sigma>, where all
      ! occupied MOs are taken into account

      !start with the alpha MOs
      CALL get_mo_set(mos(1), mo_coeff=mo_coeff, homo=homo)
      ALLOCATE (diag(homo))
      CALL cp_fm_get_info(mo_coeff, matrix_struct=vec_struct)
      CALL cp_fm_struct_create(prod_struct, context=blacs_env, para_env=para_env, &
                               nrow_global=homo, ncol_global=homo)
      CALL cp_fm_create(vec_work, vec_struct)
      CALL cp_fm_create(prod_work, prod_struct)

      ! <alpha|SOC_z|alpha> => spin integration yields +1
      CALL cp_dbcsr_sm_fm_multiply(orb_soc_z, mo_coeff, vec_work, ncol=homo)
      CALL parallel_gemm('T', 'N', homo, homo, nao, 1.0_dp, mo_coeff, vec_work, 0.0_dp, prod_work)
      CALL cp_fm_get_diag(prod_work, diag)
      gs_sum = SUM(diag)

      CALL cp_fm_release(vec_work)
      CALL cp_fm_release(prod_work)
      CALL cp_fm_struct_release(prod_struct)
      DEALLOCATE (diag)
      NULLIFY (vec_struct)

      ! Now do the same with the beta gs coeffs
      CALL get_mo_set(mos(2), mo_coeff=mo_coeff, homo=homo)
      ALLOCATE (diag(homo))
      CALL cp_fm_get_info(mo_coeff, matrix_struct=vec_struct)
      CALL cp_fm_struct_create(prod_struct, context=blacs_env, para_env=para_env, &
                               nrow_global=homo, ncol_global=homo)
      CALL cp_fm_create(vec_work, vec_struct)
      CALL cp_fm_create(prod_work, prod_struct)

      ! <beta|SOC_z|beta> => spin integration yields -1
      CALL cp_dbcsr_sm_fm_multiply(orb_soc_z, mo_coeff, vec_work, ncol=homo)
      CALL parallel_gemm('T', 'N', homo, homo, nao, 1.0_dp, mo_coeff, vec_work, 0.0_dp, prod_work)
      CALL cp_fm_get_diag(prod_work, diag)
      gs_sum = gs_sum - SUM(diag) ! -1 because of spin integration

      CALL cp_fm_release(vec_work)
      CALL cp_fm_release(prod_work)
      CALL cp_fm_struct_release(prod_struct)
      DEALLOCATE (diag)

      ! Need to compute: <phi^0_Isigma|SOC|phi^0_Jtau> for the donor MOs

      CALL cp_fm_struct_create(vec_struct, context=blacs_env, para_env=para_env, &
                               nrow_global=nao, ncol_global=nso)
      CALL cp_fm_struct_create(prod_struct, context=blacs_env, para_env=para_env, &
                               nrow_global=nso, ncol_global=nso)
      CALL cp_fm_create(vec_soc_x, vec_struct) ! for SOC_x*gs_coeffs
      CALL cp_fm_create(vec_soc_y, vec_struct) ! for SOC_y*gs_coeffs
      CALL cp_fm_create(vec_soc_z, vec_struct) ! for SOC_z*gs_coeffs
      CALL cp_fm_create(prod_work, prod_struct)
      ALLOCATE (diag(nso))

      ALLOCATE (mo_soc_x(nso, nso), mo_soc_y(nso, nso), mo_soc_z(nso, nso))

      CALL cp_dbcsr_sm_fm_multiply(orb_soc_x, gs_coeffs, vec_soc_x, ncol=nso)
      CALL parallel_gemm('T', 'N', nso, nso, nao, 1.0_dp, gs_coeffs, vec_soc_x, 0.0_dp, prod_work)
      CALL cp_fm_get_submatrix(prod_work, mo_soc_x)

      CALL cp_dbcsr_sm_fm_multiply(orb_soc_y, gs_coeffs, vec_soc_y, ncol=nso)
      CALL parallel_gemm('T', 'N', nso, nso, nao, 1.0_dp, gs_coeffs, vec_soc_y, 0.0_dp, prod_work)
      CALL cp_fm_get_submatrix(prod_work, mo_soc_y)

      CALL cp_dbcsr_sm_fm_multiply(orb_soc_z, gs_coeffs, vec_soc_z, ncol=nso)
      CALL parallel_gemm('T', 'N', nso, nso, nao, 1.0_dp, gs_coeffs, vec_soc_z, 0.0_dp, prod_work)
      CALL cp_fm_get_submatrix(prod_work, mo_soc_z)

      ! some more useful work matrices
      CALL cp_fm_struct_create(work_struct, context=blacs_env, para_env=para_env, &
                               nrow_global=nex, ncol_global=nex)
      CALL cp_fm_create(work_fm, work_struct)

!  Looping over the excited states, computing the SOC and filling the perturbation matrix
!  There are 3 loops to do: sc-sc, sc-sf and sf-sf
!  The final perturbation matrix is Hermitian, only the upper diagonal is needed

      !need some work matrices for the GS stuff
      CALL cp_fm_struct_create(gsex_struct, context=blacs_env, para_env=para_env, &
                               nrow_global=nex*nso, ncol_global=nso)
      CALL cp_fm_create(gsex_fm, gsex_struct)
      ALLOCATE (gsex_block(nso, nso))

!  Start with ground-state/spin-conserving SOC:
      !  <Psi_0|SOC|Psi_Jsc> = sum_k,sigma <phi^0_k,sigma|SOC|phi^Jsc_k,sigma>

      !compute -sc_coeffs*SOC_Z*gs_coeffs, minus sign because SOC_z antisymmetric
      CALL parallel_gemm('T', 'N', nex*nso, nso, nao, -1.0_dp, sc_coeffs, vec_soc_z, 0.0_dp, gsex_fm)

      DO isc = 1, nsc
         CALL cp_fm_get_submatrix(fm=gsex_fm, target_m=gsex_block, start_row=(isc - 1)*nso + 1, &
                                  start_col=1, n_rows=nso, n_cols=nso)
         diag(:) = get_diag(gsex_block)
         soc = SUM(diag(1:nmo)) - SUM(diag(nmo + 1:nso)) !minus sign because of spin integration

         !purely imaginary contribution
         CALL cp_fm_set_element(img_fm, 1, 1 + isc, soc)
      END DO !isc

!  Then ground-state/spin-flip SOC:
      !<Psi_0|SOC|Psi_Jsf> = sum_k,sigma <phi^0_k,sigma|SOC|phi^Jsc_k,tau>   sigma != tau

      !compute  -sc_coeffs*SOC_x*gs_coeffs, imaginary contribution
      CALL parallel_gemm('T', 'N', nex*nso, nso, nao, -1.0_dp, sc_coeffs, vec_soc_x, 0.0_dp, gsex_fm)

      DO isf = 1, nsf
         CALL cp_fm_get_submatrix(fm=gsex_fm, target_m=gsex_block, start_row=(isf - 1)*nso + 1, &
                                  start_col=1, n_rows=nso, n_cols=nso)
         diag(:) = get_diag(gsex_block)
         soc = SUM(diag) !alpha and beta parts are simply added due to spin integration
         CALL cp_fm_set_element(img_fm, 1, 1 + nsc + isf, soc)
      END DO !isf

      !compute -sc_coeffs*SOC_y*gs_coeffs, real contribution
      CALL parallel_gemm('T', 'N', nex*nso, nso, nao, -1.0_dp, sc_coeffs, vec_soc_y, 0.0_dp, gsex_fm)

      DO isf = 1, nsf
         CALL cp_fm_get_submatrix(fm=gsex_fm, target_m=gsex_block, start_row=(isf - 1)*nso + 1, &
                                  start_col=1, n_rows=nso, n_cols=nso)
         diag(:) = get_diag(gsex_block)
         soc = SUM(diag(1:nmo)) ! alpha-beta
         soc = soc - SUM(diag(nmo + 1:nso)) !beta-alpha
         CALL cp_fm_set_element(real_fm, 1, 1 + nsc + isf, soc)
      END DO !isf

      !ground-state cleanup
      CALL cp_fm_release(gsex_fm)
      CALL cp_fm_release(vec_soc_x)
      CALL cp_fm_release(vec_soc_y)
      CALL cp_fm_release(vec_soc_z)
      CALL cp_fm_release(prod_work)
      CALL cp_fm_struct_release(gsex_struct)
      CALL cp_fm_struct_release(prod_struct)
      CALL cp_fm_struct_release(vec_struct)
      DEALLOCATE (gsex_block)

!  Then spin-conserving/spin-conserving SOC
!  <Psi_Isc|SOC|Psi_Jsc> =
!  sum_k,sigma [<psi^Isc_k,sigma|SOC|psi^Jsc_k,sigma> + <psi^Isc_k,sigma|psi^Jsc_k,sigma> * gs_sum]
!  - sum_k,l,sigma <psi^0_k,sigma|SOC|psi^0_l,sigma> * <psi^Isc_l,sigma|psi^Jsc_k,sigma>

      !Same spin integration => only SOC_z matters, and the contribution is purely imaginary
      CALL dbcsr_multiply('N', 'N', 1.0_dp, orb_soc_z, dbcsr_sc, 0.0_dp, dbcsr_work, filter_eps=eps_filter)
      CALL dbcsr_multiply('T', 'N', 1.0_dp, dbcsr_sc, dbcsr_work, 0.0_dp, dbcsr_prod, filter_eps=eps_filter)

      !the overlap as well
      CALL dbcsr_multiply('N', 'N', 1.0_dp, matrix_s(1)%matrix, dbcsr_sc, 0.0_dp, dbcsr_work, &
                          filter_eps=eps_filter)
      CALL dbcsr_multiply('T', 'N', 1.0_dp, dbcsr_sc, dbcsr_work, 0.0_dp, dbcsr_ovlp, filter_eps=eps_filter)

      CALL os_amew_soc_elements(dbcsr_tmp, dbcsr_prod, dbcsr_ovlp, mo_soc_z, pref_diaga=1.0_dp, &
                                pref_diagb=-1.0_dp, pref_tracea=-1.0_dp, pref_traceb=1.0_dp, &
                                pref_diags=gs_sum, symmetric=.TRUE.)

      CALL copy_dbcsr_to_fm(dbcsr_tmp, work_fm)
      CALL cp_fm_to_fm_submat(msource=work_fm, mtarget=img_fm, nrow=nex, ncol=nex, s_firstrow=1, &
                              s_firstcol=1, t_firstrow=2, t_firstcol=2)

!  Then spin-flip/spin-flip SOC
!  <Psi_Isf|SOC|Psi_Jsf> =
!  sum_k,sigma [<psi^Isf_k,tau|SOC|psi^Jsf_k,tau> + <psi^Isf_k,tau|psi^Jsf_k,tau> * gs_sum]
!  - sum_k,l,sigma <psi^0_k,sigma|SOC|psi^0_l,sigma> * <psi^Isf_l,tau|psi^Jsf_k,tau> , tau != sigma

      !Same spin integration => only SOC_z matters, and the contribution is purely imaginary
      CALL dbcsr_multiply('N', 'N', 1.0_dp, orb_soc_z, dbcsr_sf, 0.0_dp, dbcsr_work, filter_eps=eps_filter)
      CALL dbcsr_multiply('T', 'N', 1.0_dp, dbcsr_sf, dbcsr_work, 0.0_dp, dbcsr_prod, filter_eps=eps_filter)

      !the overlap as well
      CALL dbcsr_multiply('N', 'N', 1.0_dp, matrix_s(1)%matrix, dbcsr_sf, 0.0_dp, &
                          dbcsr_work, filter_eps=eps_filter)
      CALL dbcsr_multiply('T', 'N', 1.0_dp, dbcsr_sf, dbcsr_work, 0.0_dp, dbcsr_ovlp, filter_eps=eps_filter)

      !note: the different prefactors are derived from the fact that because of spin-flip, we have
      !alpha-gs and beta-lr interaction
      CALL os_amew_soc_elements(dbcsr_tmp, dbcsr_prod, dbcsr_ovlp, mo_soc_z, pref_diaga=-1.0_dp, &
                                pref_diagb=1.0_dp, pref_tracea=-1.0_dp, pref_traceb=1.0_dp, &
                                pref_diags=gs_sum, symmetric=.TRUE.)

      CALL copy_dbcsr_to_fm(dbcsr_tmp, work_fm)
      CALL cp_fm_to_fm_submat(msource=work_fm, mtarget=img_fm, nrow=nex, ncol=nex, s_firstrow=1, &
                              s_firstcol=1, t_firstrow=1 + nsc + 1, t_firstcol=1 + nsc + 1)

!  Finally the spin-conserving/spin-flip interaction
! <Psi_Isc|SOC|Psi_Jsf> =   sum_k,sigma <psi^Isc_k,sigma|SOC|psi^Isf_k,tau>
!                           - sum_k,l,sigma <psi^0_k,tau|SOC|psi^0_l,sigma

      tas(1) = nmo + 1; tbs(1) = 1
      tas(2) = 1; tbs(2) = nmo + 1

      !the overlap
      CALL dbcsr_multiply('N', 'N', 1.0_dp, matrix_s(1)%matrix, dbcsr_sf, 0.0_dp, &
                          dbcsr_work, filter_eps=eps_filter)
      CALL dbcsr_multiply('T', 'N', 1.0_dp, dbcsr_sc, dbcsr_work, 0.0_dp, dbcsr_ovlp, filter_eps=eps_filter)

      !start with the imaginary contribution
      CALL dbcsr_multiply('N', 'N', 1.0_dp, orb_soc_x, dbcsr_sc, 0.0_dp, dbcsr_work, filter_eps=eps_filter)
      CALL dbcsr_multiply('T', 'N', 1.0_dp, dbcsr_sf, dbcsr_work, 0.0_dp, dbcsr_prod, filter_eps=eps_filter)

      CALL os_amew_soc_elements(dbcsr_tmp, dbcsr_prod, dbcsr_ovlp, mo_soc_x, pref_diaga=1.0_dp, &
                                pref_diagb=1.0_dp, pref_tracea=-1.0_dp, pref_traceb=-1.0_dp, &
                                tracea_start=tas, traceb_start=tbs)

      CALL copy_dbcsr_to_fm(dbcsr_tmp, work_fm)
      CALL cp_fm_to_fm_submat(msource=work_fm, mtarget=img_fm, nrow=nex, ncol=nex, s_firstrow=1, &
                              s_firstcol=1, t_firstrow=2, t_firstcol=1 + nsc + 1)

      !follow up with the real contribution
      CALL dbcsr_multiply('N', 'N', 1.0_dp, orb_soc_y, dbcsr_sf, 0.0_dp, dbcsr_work, filter_eps=eps_filter)
      CALL dbcsr_multiply('T', 'N', 1.0_dp, dbcsr_sc, dbcsr_work, 0.0_dp, dbcsr_prod, filter_eps=eps_filter)

      CALL os_amew_soc_elements(dbcsr_tmp, dbcsr_prod, dbcsr_ovlp, mo_soc_y, pref_diaga=1.0_dp, &
                                pref_diagb=-1.0_dp, pref_tracea=1.0_dp, pref_traceb=-1.0_dp, &
                                tracea_start=tas, traceb_start=tbs)

      CALL copy_dbcsr_to_fm(dbcsr_tmp, work_fm)
      CALL cp_fm_to_fm_submat(msource=work_fm, mtarget=real_fm, nrow=nex, ncol=nex, s_firstrow=1, &
                              s_firstcol=1, t_firstrow=2, t_firstcol=1 + nsc + 1)

!  Setting up the complex Hermitian perturbed matrix
      CALL cp_cfm_create(pert_cfm, full_struct)
      CALL cp_fm_to_cfm(real_fm, img_fm, pert_cfm)

      CALL cp_fm_release(real_fm)
      CALL cp_fm_release(img_fm)

!  Diagonalize the perturbed matrix
      ALLOCATE (tmp_evals(ntot))
      CALL cp_cfm_create(evecs_cfm, full_struct)
      CALL cp_cfm_heevd(pert_cfm, evecs_cfm, tmp_evals)

      !shift the energies such that the GS has zero and store all that in soc_evals (\wo the GS)
      ALLOCATE (soc_env%soc_evals(ntot - 1))
      soc_env%soc_evals(:) = tmp_evals(2:ntot) - tmp_evals(1)

!  The SOC dipole oscillator strengths
      CALL soc_dipol(qs_env, dbcsr_soc_package, soc_env,evecs_cfm,&
                     eps_filter, dipol_form, .TRUE.)

!  And quadrupole
!      IF (xas_tdp_control%do_quad) THEN
!         CALL compute_soc_quadrupole_fosc(evecs_cfm, dbcsr_soc_package, donor_state, xas_tdp_env, &
!                                          xas_tdp_control, qs_env, gs_coeffs=gs_coeffs)
!      END IF

! Clean-up
      CALL cp_cfm_release(pert_cfm)
      CALL cp_cfm_release(evecs_cfm)
      CALL cp_fm_struct_release(full_struct)
      IF (do_roks) THEN
         CALL cp_fm_release(gs_coeffs)
         DEALLOCATE (gs_coeffs)
      END IF
      CALL dbcsr_distribution_release(coeffs_dist)
      CALL dbcsr_distribution_release(prod_dist)
      CALL dbcsr_release(dbcsr_sc)
      CALL dbcsr_release(dbcsr_sf)
      CALL dbcsr_release(dbcsr_prod)
      CALL dbcsr_release(dbcsr_ovlp)
      CALL dbcsr_release(dbcsr_tmp)
      CALL dbcsr_release(dbcsr_work)
      CALL cp_fm_release(work_fm)
      CALL cp_fm_struct_release(work_struct)
      DEALLOCATE (coeffs_dist, prod_dist, col_dist, col_blk_size, row_dist_new)
      DEALLOCATE (dbcsr_sc, dbcsr_sf, dbcsr_work, dbcsr_prod, dbcsr_ovlp, dbcsr_tmp)

      CALL timestop(handle)


   END SUBROUTINE tddfpt_soc_os

!******************************************************************************
! \brief This routine will be handling nessesary initalizations for the
! soc calcualtions
!\param
!
!******************************************************************************

   SUBROUTINE inititialize_soc(qs_env, soc_atom_env, soc_env,&
                               evects_a, evects_b,&
                               dipole_form)
         !! Different Types 
         TYPE(cp_fm_type),DIMENSION(:,:)                  ::evects_a,&
                                                            evects_b
         TYPE(qs_environment_type), POINTER, INTENT(IN)   :: qs_env
         TYPE(qs_kind_type), DIMENSION(:), POINTER        :: qs_kind_set
         TYPE(particle_type), DIMENSION(:), POINTER       :: particle_set
         TYPE(cp_blacs_env_type), POINTER                 :: blacs_env
         TYPE(cp_para_env_type), POINTER                  :: para_env
         TYPE(soc_atom_env_type), POINTER, INTENT(OUT)    :: soc_atom_env
         TYPE(soc_env_type), INTENT(OUT), TARGET          :: soc_env
         TYPE(dbcsr_p_type), DIMENSION(:), POINTER        :: matrix_s,work
         INTEGER                                          :: ikind, nkind, & 
                                                             natom, i, j,&
                                                             dipole_form
         REAL(dp), DIMENSION(3)                           :: rc
         INTEGER                                          :: nrep,&
                                                             irep
         TYPE(cp_fm_type), POINTER                        :: a_coeff, &
                                                             b_coeff
         TYPE(section_vals_type), POINTER                 :: soc_section
         CHARACTER(len=default_string_length), &
            DIMENSION(:,:), POINTER                       :: grid_info
         CHARACTER(len=default_string_length), &
            DIMENSION(:), POINTER                         :: k_list
         CHARACTER(len=default_string_length)             :: message
         CHARACTER(len=10)                                :: str_nkinds,&
                                                             str_nrep
         LOGICAL                                          :: do_os

         NULLIFY(qs_kind_set, particle_set, blacs_env, para_env,&
                 a_coeff, b_coeff, grid_info)

         do_os = .FALSE.
         !IF (SIZE(evects_a,1) > 1) do_os=.TRUE.

         CALL soc_env_create(soc_env)
         CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set,& 
                        mos=soc_env%mos, natom=natom, & 
                        particle_set=particle_set,&
                        blacs_env=soc_env%blacs_env,& 
                        para_env=soc_env%para_env, &
                        matrix_s=matrix_s)

         para_env => soc_env%para_env
         blacs_env => soc_env%blacs_env
         a_coeff => soc_env%a_coeff
         b_coeff => soc_env%b_coeff

         soc_section => section_vals_get_subs_vals(qs_env%input, &
                                       "PROPERTIES%TDDFPT%SOC")

         CALL section_vals_val_get(soc_section,"GRID", n_rep_val=nrep)
         ALLOCATE(grid_info(nrep,3))
         DO irep = 1,nrep
            CALL section_vals_val_get(soc_section,&
                                      "GRID",&
                                      i_rep_val=irep,&
                                      c_vals=k_list)
            grid_info(irep,:) = k_list
         END DO


         ALLOCATE (soc_env%dipmat(3))
         DO i = 1, 3
            ALLOCATE (soc_env%dipmat(i)%matrix)
            CALL dbcsr_copy(soc_env%dipmat(i)%matrix,&
                             matrix_s(1)%matrix,&
                             name="dipole matrix")
            CALL dbcsr_set(soc_env%dipmat(i)%matrix, 0.0_dp)
         END DO

         !     Precompute it in the velocity representation, if so chosen
         IF (dipole_form == tddfpt_soc_dipol_vel) THEN
            !enforce minimum image to avoid any PBCs related issues. Ok because very localized densities
            CALL p_xyz_ao(soc_env%dipmat, qs_env, minimum_image=.TRUE.)
         ELSE
            DO i = 1, SIZE(particle_set)
               NULLIFY (work)
               rc = particle_set(i)%r
               ALLOCATE (work(9))
               DO j = 1,3
                  CALL dbcsr_set(soc_env%dipmat(j)%matrix, 0.0_dp)
                  work(j)%matrix => soc_env%dipmat(j)%matrix
               END DO
               CALL rRc_xyz_ao(work, qs_env, rc, order=1, minimum_image=.TRUE.)
            END DO
         END IF

         CALL get_mo_set(soc_env%mos(1), mo_coeff=soc_env%gs_coeff,&
         nmo=soc_env%nmo, nao=soc_env%nao)

         !!DEBUG
         !a_coeff = evects_a
         !b_coeff = evects_b

         CALL soc_contract_evect(evects_a,&
                                 a_coeff,&
                                 blacs_env,&
                                 para_env)
         CALL soc_contract_evect(evects_b,&
                                 b_coeff,&
                                 blacs_env,&
                                 para_env)

         ALLOCATE (soc_env%orb_soc(3))
         DO i = 1,3
            ALLOCATE(soc_env%orb_soc(i)%matrix)
         end DO


!! Make soc_atom_env for H^SOC calculations
         nkind = SIZE(qs_kind_set)
         !  Compute the contraction coefficients for spherical orbitals
         CALL soc_atom_create(soc_atom_env)
         IF (do_os) THEN
            soc_atom_env%nspins=2
         ELSE
            soc_atom_env%nspins=1
         END IF

         soc_atom_env%ri_radius = 8.0_dp

         ALLOCATE(soc_atom_env%grid_atom_set(nkind))
         ALLOCATE(soc_atom_env%harmonics_atom_set(nkind))  
         ALLOCATE (soc_atom_env%orb_sphi_so(nkind))
         ALLOCATE (soc_atom_env%ri_sphi_so(nkind))
  !       IF (do_xc) THEN
         ALLOCATE (soc_atom_env%gr(nkind))
         ALLOCATE (soc_atom_env%ga(nkind))
         ALLOCATE (soc_atom_env%dgr1(nkind))
         ALLOCATE (soc_atom_env%dgr2(nkind))
         ALLOCATE (soc_atom_env%dga1(nkind))
         ALLOCATE (soc_atom_env%dga2(nkind))
!         END IF

         ALLOCATE(soc_env%all_atoms(natom))
         ALLOCATE(soc_env%all_kinds(nkind))

         DO i = 1,natom
            soc_env%all_atoms(i) = i
         END DO

         DO i = 1,nkind
            soc_env%all_kinds(i) = i
         END DO

         soc_atom_env%excited_atoms => soc_env%all_atoms
         soc_atom_env%excited_kinds => soc_env%all_kinds
         soc_atom_env%qs_kind_set => qs_kind_set

         !! Possible Error detection
         IF (.NOT. (nrep == nkind)) THEN
            CALL integer_to_string(nkind,str_nkinds)
            CALL integer_to_string(nrep,str_nrep)
            message = "Every kind needs a Grid. There are "//trim(str_nrep)
            message = trim(message)//"and "//trim(str_nkinds)//" are nedded"
            !CPABORT(message)
         END IF
         !! Prepare the atomic grid we need to calculate the SOC Operator in an Atomic basis
         
         CALL init_atom_grid(soc_atom_env, grid_info, .TRUE., qs_env)

         DO ikind = 1, nkind
            NULLIFY (soc_atom_env%orb_sphi_so(ikind)%array)!, soc_atom_env%ri_sphi_so(ikind)%array)
            CALL compute_sphi_so(ikind,&
                                 "ORB",&
                                 soc_atom_env%orb_sphi_so(ikind)%array,&
                                 qs_env)
            CALL compute_sphi_so(ikind, "RI_XAS", soc_atom_env%ri_sphi_so(ikind)%array, qs_env)
         END DO !ikind

         !IF (do_xc) THEN
         !CALL calculate_density_coeffs(soc_atom_env=soc_atom_env,&
         !                              qs_env=qs_env)
         !END IF

         NULLIFY (para_env, blacs_env, a_coeff, b_coeff)
   END SUBROUTINE inititialize_soc

! **************************************************************************************************
!> \brief Initializes the atomic grids and harmonics for the RI atomic calculations
!> \param xas_atom_env ...
!> \param grid_info ...
!> \param do_xc Whether the xc kernel will ne computed on the atomic grids. If not, the harmonics
!>        are built for the orbital basis for all kinds.
!> \param qs_env ...
!> \note Copied and modified from init_xas_atom_grid_harmo
! **************************************************************************************************
   SUBROUTINE init_atom_grid(soc_atom_env, grid_info, do_xc, qs_env)

         TYPE(soc_atom_env_type)                          :: soc_atom_env
         TYPE(qs_environment_type)                        :: qs_env
         CHARACTER(len=default_string_length), &
            DIMENSION(:,:), POINTER                       :: grid_info
         LOGICAL,INTENT(in)                               :: do_xc

         TYPE(dft_control_type), POINTER                  :: dft_control
         TYPE(grid_atom_type), POINTER                    :: grid_atom
         TYPE(gto_basis_set_type), POINTER                :: tmp_basis
         TYPE(harmonics_atom_type), POINTER               :: harmonics
         TYPE(qs_control_type), POINTER                   :: qs_control
         TYPE(qs_kind_type), DIMENSION(:), POINTER        :: qs_kind_set
         CHARACTER(LEN=default_string_length)             :: kind_name
         INTEGER                                          :: ikind, & 
                                                             igrid, &
                                                             il, &
                                                             iso, &
                                                             iso1, &
                                                             iso2, &
                                                             l1, &
                                                             l1l2, &
                                                             l2, &
                                                             la, &
                                                             lc1, &
                                                             lc2, &
                                                             lcleb, &
                                                             ll, &
                                                             llmax, &
                                                             lp, &
                                                             m1, &
                                                             m2, &
                                                             max_s_harm,&
                                                             max_s_set, &
                                                             maxl, &
                                                             maxlgto, &
                                                             maxs, &
                                                             mm, &
                                                             mp, &
                                                             na, &
                                                             nr, &
                                                             quadrature,&
                                                             stat
         REAL(dp)                                         :: kind_radius
         REAL(dp), ALLOCATABLE, DIMENSION(:, :)           :: rga
         REAL(dp), DIMENSION(:, :, :), POINTER            :: my_CG

      NULLIFY (my_CG, qs_kind_set, dft_control)
      NULLIFY (tmp_basis, grid_atom, harmonics, qs_control) 

!!  Initialization of some integer for the CG coeff generation
      CALL get_qs_env(qs_env,&
         qs_kind_set=qs_kind_set,&
         dft_control=dft_control)

      qs_control => dft_control%qs_control

      DO ikind = 1, SIZE(qs_kind_set)
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=tmp_basis, basis_type="ORB")
         CALL init_interaction_radii_orb_basis(tmp_basis, 10.0_dp) !eps_pgf_orb=qs_control%eps_pgf_orb)
         IF (.NOT. ASSOCIATED(tmp_basis)) CPABORT("tmp_basis ORB is not Present!")
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=tmp_basis, basis_type="RI_XAS")
         CALL init_interaction_radii_orb_basis(tmp_basis, 10.0_dp) !eps_pgf_orb=qs_control%eps_pgf)
         IF (.NOT. ASSOCIATED(tmp_basis)) CPABORT("tmp_basis RI_XAS is not Present!")
       END DO


      IF (do_xc) THEN
         CALL get_qs_kind_set(qs_kind_set, maxlgto=maxlgto,basis_type="RI_XAS")
      ELSE
         CALL get_qs_kind_set(qs_kind_set, maxlgto=maxlgto, basis_type="ORB")
      END IF

      IF (maxlgto < 0) THEN
         CPABORT("tmp_basis is not Present!")
      END IF

      !maximum expansion
      llmax = 2*maxlgto
      max_s_harm = nsoset(llmax)
      max_s_set = nsoset(maxlgto)
      lcleb = llmax

!  Allocate and compute the CG coeffs (copied from init_rho_atom)
      CALL clebsch_gordon_init(lcleb)
      CALL reallocate(my_CG, 1, max_s_set, 1, max_s_set, 1, max_s_harm)

      ALLOCATE (rga(lcleb, 2))
      DO lc1 = 0, maxlgto
         DO iso1 = nsoset(lc1 - 1) + 1, nsoset(lc1)
            l1 = indso(1, iso1)
            m1 = indso(2, iso1)
            DO lc2 = 0, maxlgto
               DO iso2 = nsoset(lc2 - 1) + 1, nsoset(lc2)
                  l2 = indso(1, iso2)
                  m2 = indso(2, iso2)
                  CALL clebsch_gordon(l1, m1, l2, m2, rga)
                  IF (l1 + l2 > llmax) THEN
                     l1l2 = llmax
                  ELSE
                     l1l2 = l1 + l2
                  END IF
                  mp = m1 + m2
                  mm = m1 - m2
                  IF (m1*m2<0 .OR. (m1*m2==0 .AND. (m1<0 .OR. m2<0)))THEN
                     mp = -ABS(mp)
                     mm = -ABS(mm)
                  ELSE
                     mp = ABS(mp)
                     mm = ABS(mm)
                  END IF
                  DO lp = MOD(l1 + l2, 2), l1l2, 2
                     il = lp/2 + 1
                     IF (ABS(mp) <= lp) THEN
                     IF (mp >= 0) THEN
                        iso = nsoset(lp - 1) + lp + 1 + mp
                     ELSE
                        iso = nsoset(lp - 1) + lp + 1 - ABS(mp)
                     END IF
                     my_CG(iso1, iso2, iso) = rga(il, 1)
                     END IF
                     IF (mp /= mm .AND. ABS(mm) <= lp) THEN
                     IF (mm >= 0) THEN
                        iso = nsoset(lp - 1) + lp + 1 + mm
                     ELSE
                        iso = nsoset(lp - 1) + lp + 1 - ABS(mm)
                     END IF
                     my_CG(iso1, iso2, iso) = rga(il, 2)
                     END IF
                  END DO
               END DO ! iso2
            END DO ! lc2
         END DO ! iso1
      END DO ! lc1
      DEALLOCATE (rga)
      CALL clebsch_gordon_deallocate()

!  Create the Lebedev grids and compute the spherical harmonics
      CALL init_lebedev_grids()
      quadrature = qs_control%gapw_control%quadrature

      DO ikind = 1, SIZE(soc_atom_env%grid_atom_set)

!        Allocate the grid and the harmonics for this kind
         NULLIFY (soc_atom_env%grid_atom_set(ikind)%grid_atom)
         NULLIFY (soc_atom_env%harmonics_atom_set(ikind)%harmonics_atom)
         CALL allocate_grid_atom(soc_atom_env%grid_atom_set(ikind)%grid_atom)
         CALL allocate_harmonics_atom(&
            soc_atom_env%harmonics_atom_set(ikind)%harmonics_atom)

         NULLIFY (grid_atom, harmonics)
         grid_atom => soc_atom_env%grid_atom_set(ikind)%grid_atom
         harmonics => soc_atom_env%harmonics_atom_set(ikind)%harmonics_atom

!        Initialize some integers
         CALL get_qs_kind(qs_kind_set(ikind),&
                          ngrid_rad=nr,&
                          ngrid_ang=na,&
                          name=kind_name)

         !take the grid dimension given as input, if none, take the GAPW ones above
               !hack to convert string into integer
         DO igrid = 1, SIZE(grid_info, 1)
             IF (grid_info(igrid, 1) == kind_name) THEN
                 READ (grid_info(igrid,2), *, iostat=stat) na
                 IF (stat .NE. 0) CPABORT("The 'na' value for the GRID keyword must be an integer")
                 READ (grid_info(igrid,3), *, iostat=stat) nr
                 IF (stat .NE. 0) CPABORT("The 'nr' value for the GRID keyword must be an integer")
                 EXIT
             END IF
         END DO

         ll = get_number_of_lebedev_grid(n=na)
         na = lebedev_grid(ll)%n
         la = lebedev_grid(ll)%l
         grid_atom%ng_sphere = na
         grid_atom%nr = nr

        !Create the harmonics with the RI_SOC basis, otherwise the ORB
         IF (do_xc) THEN
            CALL get_qs_kind(qs_kind_set(ikind),&
                             basis_set=tmp_basis,&
                             basis_type="RI_XAS")
         ELSE
            CALL get_qs_kind(qs_kind_set(ikind),&
                             basis_set=tmp_basis,&
                             basis_type="ORB")
         END IF
         CALL get_gto_basis_set(gto_basis_set=tmp_basis,&
                                maxl=maxl,&
                                kind_radius=kind_radius)

         CALL create_grid_atom(grid_atom, nr, na, llmax, ll, quadrature)
         CALL truncate_radial_grid(grid_atom, kind_radius)

         maxs = nsoset(maxl)
         CALL create_harmonics_atom(harmonics, &
                                    my_CG, na,&
                                    llmax, maxs,&
                                    max_s_harm, ll,&
                                    grid_atom%wa, &
                                    grid_atom%azi,&
                                    grid_atom%pol)
        CALL get_maxl_CG(harmonics, tmp_basis, llmax, max_s_harm)
      END DO

      CALL deallocate_lebedev_grids()
      DEALLOCATE (my_CG)

   END SUBROUTINE init_atom_grid

!******************************************************************************
! \brief This routine will write a matrix into an output file. it can
! handle fm as well as dbcsr. It is able to transform an dbcsr to an fm
! befor writing the output
! \param tddfpt_print_section ...
! \param key : Key of the used section within the tddfpt_print_section
! \param logger ...
! \param dbcsr_matrix: dbcsr_p_type-Matrix for printing, optional
! \param fm_matrix: cp_fm_type-Matrix for printing, optional
! If dbcsr_matrix and fm_matrix are present, only dbcsr_matrix will be
! printed!
! \param dbcsr2fm: Should a dbcsr be transformed into an fm before
! printing?
! \param blacs_env ...
! \param header : Will be printed in Front of the matrix
! \param starter: Additional Print in Front of the matrix, specaial to
! fm-Print!
!******************************************************************************

   SUBROUTINE print_soc(tddfpt_print_section, key, logger, ext, dbcsr_matrix,&
                        fm_matrix, integ,real_vec, real_vec2, real1, real2, &
                        dbcsr2fm, blacs_env, header,&
                        starter)

         TYPE(section_vals_type), INTENT(in),  POINTER :: tddfpt_print_section
         CHARACTER(len=*), INTENT(in)                  :: key, ext
         CHARACTER(len=*), INTENT(in), OPTIONAL        :: header, &
                                                          starter
         INTEGER, INTENT(in), OPTIONAL                 :: integ
         TYPE(dbcsr_p_type), INTENT(in), OPTIONAL      :: dbcsr_matrix
         TYPE(cp_fm_type), INTENT(in), OPTIONAL        :: fm_matrix
         REAL(dp), DIMENSION(:), OPTIONAL              :: real_vec,&
                                                          real_vec2
         REAL(dp), OPTIONAL                            :: real1, real2
         REAL(dp)                                      :: real_1, real_2
         LOGICAL,  OPTIONAL                            :: dbcsr2fm
         TYPE(cp_blacs_env_type), POINTER, OPTIONAL    :: blacs_env
         TYPE(cp_fm_type), ALLOCATABLE                 :: soc_fm
         TYPE(cp_fm_struct_type), POINTER              :: soc_struc
         TYPE(cp_logger_type), INTENT(in), POINTER     :: logger
         LOGICAL                                       :: print_dbcsr,&
                                                          print_fm,&
                                                          db2fm, &
                                                          print_int,&
                                                          print_realvec,&
                                                          print_2realvec
         INTEGER                                       :: ounit, &
                                                          i, &
                                                          nfullrows_total, &
                                                          nfullcols_total
         !! Inizalize logical variables
         print_dbcsr = .FALSE.
         print_fm = .FALSE.
         db2fm =  .FALSE.
         print_realvec = .FALSE.
         print_2realvec = .FALSE.
         print_int = .FALSE.

         !! Interpret the output to find out, what should be printed
         IF (PRESENT(dbcsr2fm)) db2fm = dbcsr2fm
         IF (PRESENT(integ)) print_int = .TRUE. 
         IF (PRESENT(dbcsr_matrix) .AND. .NOT. db2fm) print_dbcsr = .TRUE.
         IF (PRESENT(fm_matrix)) print_fm = .TRUE.
         IF (PRESENT(dbcsr_matrix) .AND. db2fm) print_fm=.TRUE.
         IF (PRESENT(real_vec)) THEN 
            IF (PRESENT(real_vec2)) THEN
               print_2realvec = .TRUE.
               IF (.NOT. PRESENT(real1)) THEN
                  real_1 = 1.0_dp
               ELSE
                  real_1 = real1
               END IF
               IF (.NOT. PRESENT(real2)) THEN 
                  real_2 = 1.0_dp
               ELSE
                  real_2 = real2
               END IF
            ELSE
               print_realvec=.TRUE.
            IF (.NOT. PRESENT(real1)) THEN
               real_1 = 1.0_dp
            ELSE 
               real_1 = real1
            END IF
         END IF
      END IF

      !! Make the output unit for the project-name-SOC.soc-File
      ounit = cp_print_key_unit_nr(logger, tddfpt_print_section, key, &
                                   extension=ext, &
                                   file_form="FORMATTED", & 
                                   file_action="WRITE", &
                                   file_status="UNKNOWN")

      !! If the Unit is valid, print the matrix
      IF (ounit > 0) THEN
         WRITE(ounit, '(A)') header
         IF (db2fm) THEN
            IF (.NOT. PRESENT(blacs_env)) THEN
               WRITE(ounit, '(/,A)') "blacs_env is needed for dbcsr conversion!"
               RETURN
            END IF
             !! First ew need to convert the dbcsr_type to a fm-type
             ALLOCATE(soc_fm)
             CALL dbcsr_get_info(dbcsr_matrix%matrix, &
                                 nfullrows_total=nfullrows_total, &
                                nfullcols_total=nfullcols_total)
             CALL cp_fm_struct_create(soc_struc,&
                                      context=blacs_env,&
                                      ncol_global=nfullcols_total, &
                                       nrow_global=nfullrows_total)
            CALL cp_fm_create(soc_fm, soc_struc )
            CALL copy_dbcsr_to_fm(dbcsr_matrix%matrix, soc_fm)
            CALL cp_fm_write_formatted(soc_fm, ounit, starter)
           DEALLOCATE(soc_struc)
           DEALLOCATE(soc_fm)
         ELSE IF (print_dbcsr) THEN
           CALL dbcsr_print(dbcsr_matrix%matrix, unit_nr=ounit)
         ELSE IF (print_fm) THEN
           CALL cp_fm_write_formatted(fm_matrix, ounit, starter)
         ELSE IF (print_2realvec) THEN
           DO i = 1,SIZE(real_vec)
              write(ounit,&
            '(I3,5x,F10.6,5x,F10.6)') i,real_vec(i)*real_1, real_vec2(i)*real_2
           END DO
         ELSE IF (print_realvec) THEN
           DO i = 1,SIZE(real_vec)
              write(ounit, '(I3,5x,F10.6)')  i, real_vec(i) * real_1
           END DO
         ELSE IF (print_int) THEN
           write(ounit, '(I10,/)') integ
         ELSE
            CPABORT('ERROR: There is nothing to print!')
         END IF
   ELSE
         CPWARN('ERROR: No log-unit could be generated!')
   END IF

   CALL cp_print_key_finished_output(ounit, logger,tddfpt_print_section, "SOC_PRINT")

END SUBROUTINE print_soc

!******************************************************************************
! \brief : Will write the tddft Excitation vectors into a single fm
!          It was optaint as a nspin x nstat matrix of fms
! \param fm_start: a nspin x nstat Matrix of matries
! 

SUBROUTINE soc_contract_evect(fm_start, fm_res, blacs_env,para_env)
   TYPE(cp_blacs_env_type), POINTER            :: blacs_env
   TYPE(cp_para_env_type), POINTER             :: para_env
   TYPE(cp_fm_type), DIMENSION(:,:), INTENT(in):: fm_start
   TYPE(cp_fm_type), INTENT(out)               :: fm_res
   TYPE(cp_fm_struct_type), POINTER            :: res_struct
   CHARACTER(len=*), PARAMETER                 :: routineN="soc_contract_evect"
   INTEGER                                     :: nao,&
                                                  handle,&
                                                  nstat,&
                                                  nocc,&
                                                  i,&
                                                  j

   CALL timeset(routineN, handle)
   NULLIFY (res_struct)

   nstat = SIZE(fm_start,2)
   CALL cp_fm_get_info(fm_start(1,1), nrow_global=nao, ncol_global=nocc)
   CALL cp_fm_struct_create(res_struct, context=blacs_env,&
                            para_env=para_env, nrow_global=nao,&
                            ncol_global=SIZE(fm_start,1)*nocc*nstat)

   CALL cp_fm_create(fm_res,res_struct)

   DO j = 1,SIZE(fm_start,1)
       DO i = 1,nstat
         CALL cp_fm_to_fm_submat(fm_start(j,i),&
                                 fm_res,&
                                 nao,&
                                 nocc,&
                                 1,&
                                 1,&
                                 1,&
                                 1+nocc*(i-1)+(j-1)*nocc*nstat)
      END DO
   END DO

   CALL timestop(handle)

END SUBROUTINE  soc_contract_evect

SUBROUTINE split_cfm_into_fm(cfm, fm, nstat, struct)
       TYPE(cp_fm_type), INTENT(out), DIMENSION(:,:),&
       ALLOCATABLE                                    :: fm
       TYPE(cp_fm_struct_type)                        :: struct
       TYPE(cp_cfm_type), INTENT(in)                  :: cfm
       TYPE(cp_cfm_type)                              :: help
       INTEGER                                        :: nstat,&
                                                         ncol,&
                                                         nrow,&
                                                         i
      !CALL cp_cfm_create(help,struct)
      CALL cp_fm_struct_get(struct,&
                            nrow_global=nrow, &
                            ncol_global=ncol)
      ALLOCATE(fm(1,nstat))

      DO i=1,nstat
         CPWARN("create help")
         CALL cp_cfm_create(help,struct)
         CPWARN("copy to help")
         CALL cp_cfm_to_cfm(cfm, help,ncol=ncol,source_start=(1+ncol*(i-1))) !_columns(cfm, help, ncol, 1+ncol*(i-1))
         CPWARN("create fm")
         CALL cp_fm_create(fm(1,i), struct)
         CPWARN("copy from help")
         CALL cp_cfm_to_fm(help, mtargetr=fm(1,i))
         CALL cp_cfm_release(help)
      END DO

      !CALL cp_cfm_release(help)

END SUBROUTINE split_cfm_into_fm


SUBROUTINE soc_dipol(qs_env, dbcsr_soc_package, soc_env,&
                     soc_evecs_cfm, eps_filter, dipol_form,&
                     do_os)
   TYPE(qs_environment_type)                     :: qs_env
   TYPE(cp_blacs_env_type), POINTER              :: blacs_env
   TYPE(cp_para_env_type), POINTER               :: para_env
   TYPE(soc_env_type), TARGET                    :: soc_env
   TYPE(cp_cfm_type), INTENT(in)                 :: soc_evecs_cfm
   REAL(dp), INTENT(IN)                          :: eps_filter

   CHARACTER(len=*), PARAMETER                   :: routineN = 'soc_dipole'
   COMPLEX(dp), ALLOCATABLE, DIMENSION(:, :)     :: transdip
   INTEGER                                       :: handle, i, nosc, ntot,&
                                                    dipol_form
   REAL(dp), DIMENSION(:), POINTER               :: osc_str, & 
                                                    soc_evals
   TYPE(cp_cfm_type)                             :: dip_cfm, work1_cfm, work2_cfm
   TYPE(cp_fm_struct_type), POINTER              :: dip_struct, full_struct
   TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)   :: amew_dip
   TYPE(dbcsr_soc_package_type)                  :: dbcsr_soc_package
   LOGICAL                                       :: do_os
   !!DEBUG
   !INTEGER                                        :: output_unit

   NULLIFY (para_env, blacs_env, dip_struct, full_struct, osc_str)
   NULLIFY (soc_evals)

   CALL timeset(routineN, handle)

   !!DEBUG
   !output_unit = cp_logger_get_default_io_unit()

   blacs_env => soc_env%blacs_env
   para_env => soc_env%para_env
   soc_evals => soc_env%soc_evals
   nosc = SIZE(soc_evals)
   ntot = nosc + 1
   ALLOCATE(soc_env%soc_osc(nosc))
   osc_str => soc_env%soc_osc
   osc_str(:)  = 0.0_dp

   !get some work arrays/matrix
   CALL cp_fm_struct_create(dip_struct,&
                            context=blacs_env,&
                            para_env=para_env, &
                            nrow_global=ntot,&
                            ncol_global=1)
   CALL cp_cfm_get_info(soc_evecs_cfm, matrix_struct=full_struct)
   CALL cp_cfm_create(dip_cfm, dip_struct)
   CALL cp_cfm_create(work1_cfm, full_struct)
   CALL cp_cfm_create(work2_cfm, full_struct)
   ALLOCATE (transdip(ntot, 1))

   IF (do_os) THEN
      CALL get_os_amew_op(amew_dip,&
                          soc_env,&
                          dbcsr_soc_package,&
                          eps_filter,&
                          qs_env)
   ELSE
      CALL get_rcs_amew_op(amew_dip,&
                           soc_env,&
                           dbcsr_soc_package,&
                           eps_filter,&
                           qs_env)
   END IF

   !!DEBUG
   !IF (output_unit > 0) THEN
   !   DO i = 1, SIZE(amew_dip)
   !      CALL cp_fm_write_formatted(amew_dip(i), output_unit, "DEBUG::amew_dip")
   !   END DO
   !END IF

   DO i = 1, 3 !cartesian coord x, y, z

      !Convert the real dipole into the cfm format for calculations
      CALL cp_fm_to_cfm(msourcer=amew_dip(i), mtarget=work1_cfm)

      !compute amew_coeffs^dagger * amew_dip * amew_gs to get the transition moments
      CALL parallel_gemm('C', 'N', ntot,&
                         ntot, ntot,&
                         (1.0_dp, 0.0_dp),&
                         soc_evecs_cfm,&
                         work1_cfm, &
                         (0.0_dp, 0.0_dp),&
                         work2_cfm)
      CALL parallel_gemm('N', 'N', ntot,&
                          1, ntot,&
                          (1.0_dp, 0.0_dp),&
                          work2_cfm,&
                          soc_evecs_cfm, &
                         (0.0_dp, 0.0_dp),&
                         dip_cfm)

      CALL cp_cfm_get_submatrix(dip_cfm, transdip)

      !transition dipoles are real numbers
      osc_str(:) = osc_str(:) + REAL(transdip(2:ntot, 1))**2 
      osc_str(:) = osc_str(:) + AIMAG(transdip(2:ntot, 1))**2

   END DO !i

   !multiply with appropriate prefac depending in the rep
   IF (dipol_form == tddfpt_soc_dipol_len) THEN
      osc_str(:) = 2.0_dp/3.0_dp*soc_evals(:)*osc_str(:)
   ELSE
      osc_str(:) = 2.0_dp/3.0_dp/soc_evals(:)*osc_str(:)
   END IF

   !clean-up
   CALL cp_fm_struct_release(dip_struct)
   CALL cp_cfm_release(work1_cfm)
   CALL cp_cfm_release(work2_cfm)
   CALL cp_cfm_release(dip_cfm)
   DO i = 1, 3
      CALL cp_fm_release(amew_dip(i))
   END DO
   DEALLOCATE (amew_dip, transdip)
   

   CALL timestop(handle)

END SUBROUTINE soc_dipol

SUBROUTINE get_rcs_amew_op(amew_op, soc_env, dbcsr_soc_package, eps_filter, qs_env)

   TYPE(cp_fm_type), ALLOCATABLE,&
   DIMENSION(:), INTENT(OUT)                  :: amew_op
   TYPE(dbcsr_p_type), DIMENSION(:), POINTER  :: ao_op
   TYPE(qs_environment_type)                  :: qs_env
   TYPE(soc_env_type), TARGET                 :: soc_env
   TYPE(cp_fm_struct_type), POINTER           :: tmp_struct, &
                                                 full_struct,&
                                                 gsgs_struct,&
                                                 prod_struct, &
                                                 sggs_struct,&
                                                 std_struct,&
                                                 vec_struct
   TYPE(cp_fm_type)                           :: gs_fm,&
                                                 prod_fm,&
                                                 sggs_fm,&
                                                 tmp_fm,&
                                                 vec_op, &
                                                 work_fm
   TYPE(cp_fm_type), POINTER                  :: gs_coeffs, &
                                                 mo_coeff, &
                                                 sg_coeffs
   TYPE(dbcsr_p_type), DIMENSION(:), POINTER  :: matrix_s
   TYPE(dbcsr_type), POINTER                  :: ao_op_i,&
                                                 dbcsr_ovlp,&
                                                 dbcsr_prod, &
                                                 dbcsr_sg,&
                                                 dbcsr_tmp,&
                                                 dbcsr_tp, &
                                                 dbcsr_work
   TYPE(dbcsr_soc_package_type)               :: dbcsr_soc_package
   TYPE(mo_set_type), DIMENSION(:), POINTER   :: mos
   INTEGER                                    :: dim_op,&
                                                 homo,&
                                                 i,&
                                                 isg,&
                                                 nao,&
                                                 nmo,&
                                                 nex, &
                                                 nsg,&
                                                 ntot,&
                                                 ntp
   REAL(dp)                                   :: op, sqrt2
   REAL(dp), ALLOCATABLE, DIMENSION(:)        :: diag,&
                                                 gs_diag,&
                                                 gsgs_op
   REAL(dp), ALLOCATABLE, DIMENSION(:, :)     :: mo_op,&
                                                 sggs_block
   REAL(dp), INTENT(IN)                       :: eps_filter
   TYPE(cp_blacs_env_type), POINTER           :: blacs_env
   TYPE(cp_para_env_type), POINTER            :: para_env
   INTEGER                                    :: handle
   CHARACTER(len=*), PARAMETER                :: routineN = "get_rcs_amew_op"
   !!DEBUG
   !INTEGER                                    :: output_unit
   

   CALL timeset(routineN, handle)
   !!DEBUG
   !output_unit = cp_logger_get_default_io_unit()


   NULLIFY (gs_coeffs, sg_coeffs, matrix_s, full_struct, prod_struct)
   NULLIFY (para_env, mo_coeff, mos, gsgs_struct, std_struct)
   NULLIFY (ao_op_i, dbcsr_tp, dbcsr_sg, dbcsr_ovlp, dbcsr_work)
   NULLIFY (ao_op, dbcsr_tmp, dbcsr_prod, tmp_struct, sggs_struct)
   NULLIFY (vec_struct, blacs_env)

   !  Initialization
      ao_op => soc_env%dipmat
      gs_coeffs => soc_env%gs_coeff
      sg_coeffs => soc_env%a_coeff
      nsg = SIZE(soc_env%evals_a)
      ntp = nsg; nex = nsg !all the same by construction, keep them separate for clarity
      ntot = 1 + nsg + 3*ntp
      nmo = soc_env%nmo
      blacs_env => soc_env%blacs_env
      para_env => soc_env%para_env
      mos => soc_env%mos
      CALL get_qs_env(qs_env, matrix_s=matrix_s)
      sqrt2 = SQRT(2.0_dp)
      dim_op = SIZE(ao_op)

      dbcsr_sg => dbcsr_soc_package%dbcsr_sg
      dbcsr_tp => dbcsr_soc_package%dbcsr_tp
      dbcsr_work => dbcsr_soc_package%dbcsr_work
      dbcsr_prod => dbcsr_soc_package%dbcsr_prod
      dbcsr_ovlp => dbcsr_soc_package%dbcsr_ovlp
      dbcsr_tmp => dbcsr_soc_package%dbcsr_tmp

      !  Create the amew_op matrix
      CALL cp_fm_struct_create(full_struct,&
                               context=blacs_env,&
                               para_env=para_env, &
                               nrow_global=ntot,&
                               ncol_global=ntot)
      ALLOCATE (amew_op(dim_op))
      DO i = 1, dim_op
         CALL cp_fm_create(amew_op(i), full_struct)
      END DO !i

!  Deal with the GS-GS contribution <0|0> = 2*sum_j <phi_j|op|phi_j>
      CALL get_mo_set(mos(1), mo_coeff=mo_coeff, nao=nao, homo=homo)

      CALL cp_fm_struct_create(gsgs_struct,&
                               context=blacs_env,&
                               para_env=para_env, &
                               nrow_global=homo,&
                               ncol_global=homo)
      CALL cp_fm_get_info(mo_coeff,&
                          matrix_struct=std_struct)
      CALL cp_fm_create(gs_fm, gsgs_struct)
      CALL cp_fm_create(work_fm, std_struct)
      ALLOCATE (gsgs_op(dim_op))
      ALLOCATE (gs_diag(homo))

      DO i = 1, dim_op

         ao_op_i => ao_op(i)%matrix

         CALL cp_dbcsr_sm_fm_multiply(ao_op_i, mo_coeff, work_fm,ncol=homo)
         CALL parallel_gemm('T', 'N', homo,&
                            homo, nao,&
                            1.0_dp, mo_coeff,&
                            work_fm, 0.0_dp, gs_fm)
         CALL cp_fm_get_diag(gs_fm, gs_diag)
         gsgs_op(i) = 2.0_dp*SUM(gs_diag)

      END DO !i

      CALL cp_fm_release(gs_fm)
      CALL cp_fm_release(work_fm)
      CALL cp_fm_struct_release(gsgs_struct)
      DEALLOCATE (gs_diag)

      !  Create the work and helper fms
      CALL cp_fm_get_info(gs_coeffs, matrix_struct=vec_struct)
      CALL cp_fm_struct_create(prod_struct,&
                               context=blacs_env,&
                               para_env=para_env, &
                               nrow_global=nmo,&
                               ncol_global=nmo)
      CALL cp_fm_create(prod_fm, prod_struct)
      CALL cp_fm_create(vec_op, vec_struct)
      CALL cp_fm_struct_create(tmp_struct,&
                               context=blacs_env,&
                               para_env=para_env, &
                               nrow_global=nex,&
                               ncol_global=nex)
      CALL cp_fm_struct_create(sggs_struct,&
                               context=blacs_env,&
                               para_env=para_env, &
                               nrow_global=nmo*nsg,&
                               ncol_global=nmo)
      CALL cp_fm_create(tmp_fm, tmp_struct)
      CALL cp_fm_create(work_fm, full_struct)
      CALL cp_fm_create(sggs_fm, sggs_struct)
      ALLOCATE (diag(nmo))
      ALLOCATE (mo_op(nmo, nmo))
      ALLOCATE (sggs_block(nmo, nmo))

! Iterate over the dimensions of the operator
! Note: operator matrices are asusmed symmetric, can only do upper half
      DO i = 1, dim_op

         ao_op_i => ao_op(i)%matrix

         ! The GS-GS contribution
         CALL cp_fm_set_element(amew_op(i), 1, 1, gsgs_op(i))

         ! Compute the operator in the donor MOs basis
         CALL cp_dbcsr_sm_fm_multiply(ao_op_i,&
                                      gs_coeffs,&
                                      vec_op,&
                                      ncol=nmo)
         CALL parallel_gemm('T', 'N', nmo,&
                            nmo, nao,&
                            1.0_dp, gs_coeffs,&
                            vec_op, 0.0_dp,&
                            prod_fm)
         CALL cp_fm_get_submatrix(prod_fm, mo_op)

         ! Compute the ground-state/singlet components. ao_op*gs_coeffs already stored in vec_op
         CALL parallel_gemm('T', 'N', nmo*nsg,&
                            nmo, nao,&
                            1.0_dp, sg_coeffs,&
                            vec_op, 0.0_dp,&
                            sggs_fm)
         DO isg = 1, nsg
            CALL cp_fm_get_submatrix(fm=sggs_fm,&
                                     target_m=sggs_block,&
                                     start_row=(isg - 1)*nmo + 1, &
                                     start_col=1,&
                                     n_rows=nmo,&
                                     n_cols=nmo)
            diag(:) = get_diag(sggs_block)
            op = sqrt2*SUM(diag)
            CALL cp_fm_set_element(amew_op(i), 1, 1 + isg, op)
         END DO

         ! do the singlet-singlet components
         !start with the overlap
         CALL dbcsr_multiply('N', 'N', 1.0_dp,&
                             matrix_s(1)%matrix, dbcsr_sg,&
                             0.0_dp, dbcsr_work,&
                             filter_eps=eps_filter)
         CALL dbcsr_multiply('T', 'N', 1.0_dp,&
                             dbcsr_sg, dbcsr_work,&
                             0.0_dp, dbcsr_ovlp,&
                             filter_eps=eps_filter)

         !then the operator in the LR orbital basis
         CALL dbcsr_multiply('N', 'N', 1.0_dp,&
                             ao_op_i, dbcsr_sg,&
                             0.0_dp, dbcsr_work,&
                             filter_eps=eps_filter)
         CALL dbcsr_multiply('T', 'N', 1.0_dp,&
                             dbcsr_sg, dbcsr_work,&
                             0.0_dp, dbcsr_prod,&
                             filter_eps=eps_filter)

         !use the soc routine, it is compatible
         CALL rcs_amew_soc_elements(dbcsr_tmp,&
                                    dbcsr_prod,&
                                    dbcsr_ovlp,&
                                    mo_op,&
                                    pref_trace=-1.0_dp, &
                                    pref_overall=1.0_dp,&
                                    pref_diags=gsgs_op(i),&
                                    symmetric=.TRUE.)

         CALL copy_dbcsr_to_fm(dbcsr_tmp, tmp_fm)
         CALL cp_fm_to_fm_submat(msource=tmp_fm,&
                                 mtarget=amew_op(i),&
                                 nrow=nex,&
                                 ncol=nex, &
                                 s_firstrow=1,&
                                 s_firstcol=1,&
                                 t_firstrow=2,&
                                 t_firstcol=2)

         ! compute the triplet-triplet components
         !the overlap
         CALL dbcsr_multiply('N', 'N', 1.0_dp,&
                             matrix_s(1)%matrix,&
                             dbcsr_tp, 0.0_dp, &
                             dbcsr_work,&
                             filter_eps=eps_filter)
         CALL dbcsr_multiply('T', 'N', 1.0_dp,&
                             dbcsr_tp, dbcsr_work,&
                             0.0_dp, dbcsr_ovlp,&
                             filter_eps=eps_filter)

         !the operator in the LR orbital basis
         CALL dbcsr_multiply('N', 'N', 1.0_dp,&
                             ao_op_i, dbcsr_sg,&
                             0.0_dp, dbcsr_work,&
                             filter_eps=eps_filter)
         CALL dbcsr_multiply('T', 'N', 1.0_dp,&
                             dbcsr_sg, dbcsr_work,&
                             0.0_dp, dbcsr_prod,&
                             filter_eps=eps_filter)

         CALL rcs_amew_soc_elements(dbcsr_tmp,&
                                    dbcsr_prod,&
                                    dbcsr_ovlp,&
                                    mo_op,&
                                    pref_trace=-1.0_dp, &
                                    pref_overall=1.0_dp,&
                                    pref_diags=gsgs_op(i),&
                                    symmetric=.TRUE.)

         CALL copy_dbcsr_to_fm(dbcsr_tmp, tmp_fm)
         !<T^-1|op|T^-1>
         CALL cp_fm_to_fm_submat(msource=tmp_fm,&
                                 mtarget=amew_op(i),&
                                 nrow=nex,&
                                 ncol=nex, &
                                 s_firstrow=1,&
                                 s_firstcol=1,&
                                 t_firstrow=1 + nsg + 1,&
                                 t_firstcol=1 + nsg + 1)

         !<T^0|op|T^0>
         CALL cp_fm_to_fm_submat(msource=tmp_fm,&
                                 mtarget=amew_op(i),&
                                 nrow=nex,&
                                 ncol=nex, &
                                 s_firstrow=1,&
                                 s_firstcol=1,&
                                 t_firstrow=1 + nsg + ntp + 1, &
                                 t_firstcol=1 + nsg + ntp + 1)
         !<T^-1|op|T^-1>
         CALL cp_fm_to_fm_submat(msource=tmp_fm,&
                                 mtarget=amew_op(i),&
                                 nrow=nex,&
                                 ncol=nex, &
                                 s_firstrow=1,&
                                 s_firstcol=1,&
                                 t_firstrow=1 + nsg + 2*ntp + 1, &
                                 t_firstcol=1 + nsg + 2*ntp + 1)

         ! Symmetrize the matrix (only upper triangle built)
         CALL cp_fm_upper_to_full(amew_op(i), work_fm)

      END DO !i

!  Clean-up
      CALL cp_fm_release(prod_fm)
      CALL cp_fm_release(work_fm)
      CALL cp_fm_release(tmp_fm)
      CALL cp_fm_release(vec_op)
      CALL cp_fm_release(sggs_fm)
      CALL cp_fm_struct_release(prod_struct)
      CALL cp_fm_struct_release(full_struct)
      CALL cp_fm_struct_release(tmp_struct)
      CALL cp_fm_struct_release(sggs_struct)

      DEALLOCATE(diag, gsgs_op)
      CALL timestop(handle)

   END SUBROUTINE get_rcs_amew_op

   ! **************************************************************************************************
!> \brief Computes the os SOC matrix elements between excited states AMEWs based on the LR orbitals
!> \param amew_soc output dbcsr matrix with the SOC in the AMEW basis (needs to be fully resereved)
!> \param lr_soc dbcsr matrix with the SOC wrt the LR orbitals
!> \param lr_overlap dbcsr matrix with the excited states LR orbital overlap
!> \param mo_soc the SOC in the basis of the MOs
!> \param pref_diaga ...
!> \param pref_diagb ...
!> \param pref_tracea ...
!> \param pref_traceb ...
!> \param pref_diags see notes
!> \param symmetric if the outcome is known to be symmetric, only elements with iex <= jex are done
!> \param tracea_start the indices where to start in the trace part for alpha
!> \param traceb_start the indices where to start in the trace part for beta
!> \note For an excited states pair i,j, the AMEW SOC matrix element is:
!>       soc_ij =   pref_diaga*SUM(alpha part of diag of lr_soc_ij)
!>                + pref_diagb*SUM(beta part of diag of lr_soc_ij)
!>                + pref_tracea*SUM(alpha part of lr_ovlp_ij*TRANSPOSE(mo_soc))
!>                + pref_traceb*SUM(beta part of lr_ovlp_ij*TRANSPOSE(mo_soc))
!>       optinally, one can add pref_diags*SUM(diag lr_ovlp_ij)
! **************************************************************************************************
   SUBROUTINE os_amew_soc_elements(amew_soc, lr_soc, lr_overlap, mo_soc, pref_diaga, &
                                   pref_diagb, pref_tracea, pref_traceb, pref_diags, &
                                   symmetric, tracea_start, traceb_start)

      TYPE(dbcsr_type)                                        :: amew_soc, lr_soc, lr_overlap
      REAL(dp), DIMENSION(:, :)                               :: mo_soc
      REAL(dp)                                                :: pref_diaga, pref_diagb, pref_tracea, &
                                                                 pref_traceb
      REAL(dp), OPTIONAL                                      :: pref_diags
      LOGICAL, OPTIONAL                                       :: symmetric
      INTEGER, DIMENSION(2), OPTIONAL                         :: tracea_start, traceb_start

      INTEGER                                                 :: blk, iex, jex, nmo, nso
      INTEGER, DIMENSION(2)                                   :: tas, tbs
      LOGICAL                                                 :: do_diags, found, my_symm
      REAL(dp)                                                :: soc_elem
      REAL(dp), ALLOCATABLE, DIMENSION(:)                     :: diag
      REAL(dp), DIMENSION(:, :), POINTER                      :: pblock
      TYPE(dbcsr_iterator_type)                               :: iter

      nso = SIZE(mo_soc, 1)
      nmo = nso/2
      ALLOCATE (diag(nso))
      my_symm = .FALSE.
      IF (PRESENT(symmetric)) my_symm = symmetric
      do_diags = .FALSE.
      IF (PRESENT(pref_diags)) do_diags = .TRUE.

      !by default, alpha part is (1:ndo_mo,1:ndo_mo) and beta is (ndo_mo+1:ndo_so,ndo_mo+1:ndo_so)
      !note: in some SF cases, that might change, mainly because the spin-flip LR-coeffs have
      !inverse order, that is: the beta-coeffs in the alpha spot and the alpha coeffs in the
      !beta spot
      tas = 1
      tbs = nmo + 1
      IF (PRESENT(tracea_start)) tas = tracea_start
      IF (PRESENT(traceb_start)) tbs = traceb_start

      CALL dbcsr_set(amew_soc, 0.0_dp)
      !loop over the excited states pairs as the block of amew_soc (which are all reserved)
      CALL dbcsr_iterator_start(iter, amew_soc)
      DO WHILE (dbcsr_iterator_blocks_left(iter))

         CALL dbcsr_iterator_next_block(iter, row=iex, column=jex, blk=blk)

         IF (my_symm .AND. iex > jex) CYCLE

         !compute the soc matrix element
         soc_elem = 0.0_dp
         CALL dbcsr_get_block_p(lr_soc, iex, jex, pblock, found)
         IF (found) THEN
            diag(:) = get_diag(pblock)
            soc_elem = soc_elem + pref_diaga*SUM(diag(1:nmo)) + pref_diagb*(SUM(diag(nmo + 1:nso)))
         END IF

         CALL dbcsr_get_block_p(lr_overlap, iex, jex, pblock, found)
         IF (found) THEN
            soc_elem = soc_elem &
                       + pref_tracea*SUM(pblock(tas(1):tas(1) + nmo - 1, tas(2):tas(2) + nmo - 1)* &
                                         mo_soc(tas(1):tas(1) + nmo - 1, tas(2):tas(2) + nmo - 1)) &
                       + pref_traceb*SUM(pblock(tbs(1):tbs(1) + nmo - 1, tbs(2):tbs(2) + nmo - 1)* &
                                         mo_soc(tbs(1):tbs(1) + nmo - 1, tbs(2):tbs(2) + nmo - 1))

            IF (do_diags) THEN
               diag(:) = get_diag(pblock)
               soc_elem = soc_elem + pref_diags*SUM(diag)
            END IF
         END IF

         CALL dbcsr_get_block_p(amew_soc, iex, jex, pblock, found)
         pblock = soc_elem

      END DO
      CALL dbcsr_iterator_stop(iter)

   END SUBROUTINE os_amew_soc_elements

! **************************************************************************************************
!> \brief Computes the matrix elements of a one-body operator (given wrt AOs) in the basis of the
!>        excited state AMEWs with ground state, for the open-shell case
!> \param amew_op the operator in the basis of the AMEWs (array because could have x,y,z components)
!> \param ao_op the operator in the basis of the atomic orbitals
!> \param gs_coeffs the coefficient of the GS donor MOs. Ecplicitely passed because of special
!>                  format in the ROKS case (see include_os_soc routine)
!> \param dbcsr_soc_package inhertited from the main SOC routine
!> \param donor_state ...
!> \param eps_filter ...
!> \param qs_env ...
!> \note The ordering of the AMEWs is consistent with SOC and is gs, sc, sf
!>       We assume that the operator is spin-independent => only <0|0>, <0|sc>, <sc|sc> and <sf|sf>
!>       yield non-zero matrix elements
!>       Only for open-shell calculations
! **************************************************************************************************
   SUBROUTINE get_os_amew_op(amew_op, soc_env, dbcsr_soc_package, eps_filter, qs_env)

      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:), &
         INTENT(OUT)                                          :: amew_op
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER               :: ao_op
      TYPE(qs_environment_type)                               :: qs_env
      TYPE(soc_env_type), TARGET                              :: soc_env
      TYPE(cp_fm_struct_type), POINTER                        :: full_struct, gsex_struct, prod_struct, &
                                                                 tmp_struct, vec_struct
      TYPE(cp_fm_type)                                        :: gsex_fm, tmp_fm,&
                                                                 prod_work, vec_work, work_fm
      TYPE(cp_fm_type), POINTER                               :: sc_coeffs, &
                                                                 gs_coeff,&
                                                                 mo_coeff, &
                                                                 sf_coeffs
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER               :: matrix_s
      TYPE(dbcsr_type), POINTER                               :: ao_op_i, dbcsr_ovlp, dbcsr_prod, &
                                                                 dbcsr_sf, dbcsr_tmp, dbcsr_sc, &
                                                                 dbcsr_work
      TYPE(dbcsr_soc_package_type)                            :: dbcsr_soc_package
      TYPE(mo_set_type), DIMENSION(:), POINTER                :: mos
      INTEGER                                                 :: dim_op, homo, i, nao, nex, &
                                                                 ntot
      REAL(dp)                                                :: op
      REAL(dp), ALLOCATABLE, DIMENSION(:)                     :: diag, gsgs_op
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)                  :: mo_op, gsex_block, tmp
      REAL(dp), INTENT(IN)                                    :: eps_filter
      TYPE(cp_blacs_env_type), POINTER                        :: blacs_env
      TYPE(cp_para_env_type), POINTER                         :: para_env
      INTEGER                                                 :: handle,&
                                                                 nmo, &
                                                                 nso, &
                                                                 nsc, &
                                                                 nsf, &
                                                                 isc
      CHARACTER(len=*), PARAMETER                             :: routineN = "get_rcs_amew_op"

      CALL timeset(routineN, handle)


      NULLIFY (matrix_s, para_env, blacs_env, full_struct, vec_struct, prod_struct, mos)
      NULLIFY (mo_coeff, ao_op_i, tmp_struct)
      NULLIFY (dbcsr_sc, dbcsr_sf, dbcsr_ovlp, dbcsr_work, dbcsr_tmp, dbcsr_prod)

   !  Initialization
       ao_op => soc_env%dipmat
      gs_coeff => soc_env%gs_coeff
      sc_coeffs => soc_env%a_coeff
      sf_coeffs => soc_env%b_coeff
      nsc = SIZE(soc_env%evals_a)
      nsf = nsc; nex = nsc        !all the same by construction, keep them separate for clarity
      ntot = 1 + nsc + nsf
      nmo = soc_env%nmo
      nso = 2 * nmo
      nex = nsc
      blacs_env => soc_env%blacs_env
      para_env => soc_env%para_env
      mos => soc_env%mos
      CALL get_qs_env(qs_env, matrix_s=matrix_s)
      dim_op = SIZE(ao_op)

      dbcsr_sc => dbcsr_soc_package%dbcsr_sc
      dbcsr_sf => dbcsr_soc_package%dbcsr_sf
      dbcsr_work => dbcsr_soc_package%dbcsr_work
      dbcsr_prod => dbcsr_soc_package%dbcsr_prod
      dbcsr_ovlp => dbcsr_soc_package%dbcsr_ovlp
      dbcsr_tmp => dbcsr_soc_package%dbcsr_tmp

      !  Create the amew_op matrix set
      CALL cp_fm_struct_create(full_struct, context=blacs_env, para_env=para_env, &
                               nrow_global=ntot, ncol_global=ntot)
      ALLOCATE (amew_op(dim_op))
      DO i = 1, dim_op
         CALL cp_fm_create(amew_op(i), full_struct)
      END DO

!  Before looping, need to evaluate sum_j,sigma <phi^0_j,sgima|op|phi^0_j,sigma>, for each dimension
!  of the operator
      ALLOCATE (gsgs_op(dim_op))

      !start with the alpha MOs
      CALL get_mo_set(mos(1), mo_coeff=mo_coeff, homo=homo)
      ALLOCATE (diag(homo))
      CALL cp_fm_get_info(mo_coeff, matrix_struct=vec_struct)
      CALL cp_fm_struct_create(prod_struct, context=blacs_env, para_env=para_env, &
                               nrow_global=homo, ncol_global=homo)
      CALL cp_fm_create(vec_work, vec_struct)
      CALL cp_fm_create(prod_work, prod_struct)

      DO i = 1, dim_op

         ao_op_i => ao_op(i)%matrix

         CALL cp_dbcsr_sm_fm_multiply(ao_op_i, mo_coeff, vec_work, ncol=homo)
         CALL parallel_gemm('T', 'N', homo, homo, nao, 1.0_dp, mo_coeff, vec_work, 0.0_dp, prod_work)
         CALL cp_fm_get_diag(prod_work, diag)
         gsgs_op(i) = SUM(diag)

      END DO !i

      CALL cp_fm_release(vec_work)
      CALL cp_fm_release(prod_work)
      CALL cp_fm_struct_release(prod_struct)
      DEALLOCATE (diag)
      NULLIFY (vec_struct)

      !then beta orbitals
      CALL get_mo_set(mos(2), mo_coeff=mo_coeff, homo=homo)
      ALLOCATE (diag(homo))
      CALL cp_fm_get_info(mo_coeff, matrix_struct=vec_struct)
      CALL cp_fm_struct_create(prod_struct, context=blacs_env, para_env=para_env, &
                               nrow_global=homo, ncol_global=homo)
      CALL cp_fm_create(vec_work, vec_struct)
      CALL cp_fm_create(prod_work, prod_struct)

      DO i = 1, dim_op

         ao_op_i => ao_op(i)%matrix

         CALL cp_dbcsr_sm_fm_multiply(ao_op_i, mo_coeff, vec_work, ncol=homo)
         CALL parallel_gemm('T', 'N', homo, homo, nao, 1.0_dp, mo_coeff, vec_work, 0.0_dp, prod_work)
         CALL cp_fm_get_diag(prod_work, diag)
         gsgs_op(i) = gsgs_op(i) + SUM(diag)

      END DO !i

      CALL cp_fm_release(vec_work)
      CALL cp_fm_release(prod_work)
      CALL cp_fm_struct_release(prod_struct)
      DEALLOCATE (diag)
      NULLIFY (vec_struct)

      !then beta orbitals
      CALL get_mo_set(mos(2), mo_coeff=mo_coeff, homo=homo)
      ALLOCATE (diag(homo))
      CALL cp_fm_get_info(mo_coeff, matrix_struct=vec_struct)
      CALL cp_fm_struct_create(prod_struct, context=blacs_env, para_env=para_env, &
                               nrow_global=homo, ncol_global=homo)
      CALL cp_fm_create(vec_work, vec_struct)
      CALL cp_fm_create(prod_work, prod_struct)

      DO i = 1, dim_op

         ao_op_i => ao_op(i)%matrix

         CALL cp_dbcsr_sm_fm_multiply(ao_op_i, mo_coeff, vec_work, ncol=homo)
         CALL parallel_gemm('T', 'N', homo, homo, nao, 1.0_dp, mo_coeff, vec_work, 0.0_dp, prod_work)
         CALL cp_fm_get_diag(prod_work, diag)
         gsgs_op(i) = gsgs_op(i) + SUM(diag)

      END DO !i

      CALL cp_fm_release(vec_work)
      CALL cp_fm_release(prod_work)
      CALL cp_fm_struct_release(prod_struct)
      DEALLOCATE (diag)
      NULLIFY (vec_struct)

!  Before looping over excited AMEWs, define some work matrices and structures
      CALL cp_fm_struct_create(vec_struct, context=blacs_env, para_env=para_env, &
                               nrow_global=nao, ncol_global=nso)
      CALL cp_fm_struct_create(prod_struct, context=blacs_env, para_env=para_env, &
                               nrow_global=nso, ncol_global=nso)
      CALL cp_fm_struct_create(gsex_struct, context=blacs_env, para_env=para_env, &
                               nrow_global=nso*nex, ncol_global=nso)
      CALL cp_fm_struct_create(tmp_struct, context=blacs_env, para_env=para_env, &
                               nrow_global=nex, ncol_global=nex)
      CALL cp_fm_create(vec_work, vec_struct) !for op*|phi>
      CALL cp_fm_create(prod_work, prod_struct) !for any <phi|op|phi>
      CALL cp_fm_create(work_fm, full_struct)
      CALL cp_fm_create(gsex_fm, gsex_struct)
      CALL cp_fm_create(tmp_fm, tmp_struct)
      ALLOCATE (diag(nso))
      ALLOCATE (mo_op(nso, nso))
      ALLOCATE (tmp(nso, nso))
      ALLOCATE (gsex_block(nso, nso))

!  Loop over the dimensions of the operator
      DO i = 1, dim_op

         ao_op_i => ao_op(i)%matrix
         !put the gs-gs contribution
         CALL cp_fm_set_element(amew_op(i), 1, 1, gsgs_op(i))

         !  Precompute what we can before looping over excited states
         ! Need the operator in the donor MOs basis <phi^0_I,sigma|op_i|phi^0_J,tau>
         CALL cp_dbcsr_sm_fm_multiply(ao_op_i, gs_coeff, vec_work, ncol=nso)
         CALL parallel_gemm('T', 'N', nso, nso, nao, 1.0_dp, gs_coeff, vec_work, 0.0_dp, prod_work)
         CALL cp_fm_get_submatrix(prod_work, mo_op)

         !  Do the ground-state/spin-conserving operator
         CALL parallel_gemm('T', 'N', nso*nsc, nso, nao, 1.0_dp, sc_coeffs, vec_work, 0.0_dp, gsex_fm)
         DO isc = 1, nsc
            CALL cp_fm_get_submatrix(fm=gsex_fm, target_m=gsex_block, start_row=(isc - 1)*nso + 1, &
                                     start_col=1, n_rows=nso, n_cols=nso)
            diag(:) = get_diag(gsex_block)
            op = SUM(diag)
            CALL cp_fm_set_element(amew_op(i), 1, 1 + isc, op)
         END DO !isc

         !  The spin-conserving/spin-conserving operator
         !overlap
         CALL dbcsr_multiply('N', 'N', 1.0_dp, matrix_s(1)%matrix, dbcsr_sc, 0.0_dp, &
                             dbcsr_work, filter_eps=eps_filter)
         CALL dbcsr_multiply('T', 'N', 1.0_dp, dbcsr_sc, dbcsr_work, 0.0_dp, dbcsr_ovlp, filter_eps=eps_filter)

         !operator in SC LR-orbital basis
         CALL dbcsr_multiply('N', 'N', 1.0_dp, ao_op_i, dbcsr_sc, 0.0_dp, dbcsr_work, filter_eps=eps_filter)
         CALL dbcsr_multiply('T', 'N', 1.0_dp, dbcsr_sc, dbcsr_work, 0.0_dp, dbcsr_prod, filter_eps=eps_filter)

         CALL os_amew_soc_elements(dbcsr_tmp, dbcsr_prod, dbcsr_ovlp, mo_op, pref_diaga=1.0_dp, &
                                   pref_diagb=1.0_dp, pref_tracea=-1.0_dp, pref_traceb=-1.0_dp, &
                                   pref_diags=gsgs_op(i), symmetric=.TRUE.)

         CALL copy_dbcsr_to_fm(dbcsr_tmp, tmp_fm)
         CALL cp_fm_to_fm_submat(msource=tmp_fm, mtarget=amew_op(i), nrow=nex, ncol=nex, &
                                 s_firstrow=1, s_firstcol=1, t_firstrow=2, t_firstcol=2)

         !  The spin-flip/spin-flip operator
         !overlap
         CALL dbcsr_multiply('N', 'N', 1.0_dp, matrix_s(1)%matrix, dbcsr_sf, 0.0_dp, &
                             dbcsr_work, filter_eps=eps_filter)
         CALL dbcsr_multiply('T', 'N', 1.0_dp, dbcsr_sf, dbcsr_work, 0.0_dp, dbcsr_ovlp, filter_eps=eps_filter)

         !operator in SF LR-orbital basis
         CALL dbcsr_multiply('N', 'N', 1.0_dp, ao_op_i, dbcsr_sf, 0.0_dp, dbcsr_work, filter_eps=eps_filter)
         CALL dbcsr_multiply('T', 'N', 1.0_dp, dbcsr_sf, dbcsr_work, 0.0_dp, dbcsr_prod, filter_eps=eps_filter)

         !need to reorganize the mo_op array by swapping the alpha-alpha and the beta-beta quarter
         tmp(1:nmo, 1:nmo) = mo_op(nmo + 1:nso, nmo + 1:nso)
         tmp(nmo + 1:nso, nmo + 1:nso) = mo_op(1:nmo, 1:nmo)

         CALL os_amew_soc_elements(dbcsr_tmp, dbcsr_prod, dbcsr_ovlp, tmp, pref_diaga=1.0_dp, &
                                   pref_diagb=1.0_dp, pref_tracea=-1.0_dp, pref_traceb=-1.0_dp, &
                                   pref_diags=gsgs_op(i), symmetric=.TRUE.)

         CALL copy_dbcsr_to_fm(dbcsr_tmp, tmp_fm)
         CALL cp_fm_to_fm_submat(msource=tmp_fm,&
                                 mtarget=amew_op(i),&
                                 nrow=nex,&
                                 ncol=nex, &
                                 s_firstrow=1,&
                                 s_firstcol=1,&
                                 t_firstrow=1 + nsc + 1,&
                                 t_firstcol=1 + nsc + 1)

         !Symmetry => only upper diag explicitly built
         CALL cp_fm_upper_to_full(amew_op(i), work_fm)

      END DO !i

!  Clean-up
      CALL cp_fm_struct_release(full_struct)
      CALL cp_fm_struct_release(prod_struct)
      CALL cp_fm_struct_release(vec_struct)
      CALL cp_fm_struct_release(tmp_struct)
      CALL cp_fm_struct_release(gsex_struct)
      CALL cp_fm_release(work_fm)
      CALL cp_fm_release(tmp_fm)
      CALL cp_fm_release(vec_work)
      CALL cp_fm_release(prod_work)
      CALL cp_fm_release(gsex_fm)


   END SUBROUTINE get_os_amew_op

END MODULE qs_tddfpt2_soc
