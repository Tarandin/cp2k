!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Utilities for X-ray absorption spectroscopy using TDDFPT
!> \author AB (01.2018)
!> This subroutine was directly copied from xas_tdp_utils to prevent a
!> circular dependency during the debug-Phase
! **************************************************************************************************

MODULE qs_tddfpt2_soc_utils
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_cfm_diag,                     ONLY: cp_cfm_heevd
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_get_info,&
                                              cp_cfm_get_submatrix,&
                                              cp_cfm_release,&
                                              cp_cfm_type,&
                                              cp_fm_to_cfm
   USE cp_control_types,                ONLY: tddfpt2_control_type
   USE cp_dbcsr_cholesky,               ONLY: cp_dbcsr_cholesky_decompose,&
                                              cp_dbcsr_cholesky_invert
   USE cp_dbcsr_diag,                   ONLY: cp_dbcsr_power
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr,&
                                              cp_dbcsr_sm_fm_multiply,&
                                              dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale,&
                                              cp_fm_scale,&
                                              cp_fm_transpose,&
                                              cp_fm_upper_to_full
   USE cp_fm_diag,                      ONLY: choose_eigv_solver,&
                                              cp_fm_geeig
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_diag,&
                                              cp_fm_get_info,&
                                              cp_fm_get_submatrix,&
                                              cp_fm_release,&
                                              cp_fm_set_element,&
                                              cp_fm_to_fm_submat,&
                                              cp_fm_type,&
                                              cp_fm_write_formatted
   USE cp_log_handling,                 ONLY: cp_logger_get_default_io_unit
   USE dbcsr_api,                       ONLY: &
        dbcsr_add, dbcsr_copy, dbcsr_create, dbcsr_distribution_get, dbcsr_distribution_new, &
        dbcsr_distribution_release, dbcsr_distribution_type, dbcsr_finalize, dbcsr_get_block_p, &
        dbcsr_get_info, dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, &
        dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_multiply, &
        dbcsr_p_type, dbcsr_put_block, dbcsr_release, dbcsr_reserve_all_blocks, dbcsr_set, &
        dbcsr_type, dbcsr_type_no_symmetry, dbcsr_type_symmetric, dbcsr_print, dbcsr_init_p
   USE input_constants,                 ONLY: ot_precond_full_single,&
                                              tddfpt_singlet,&
                                              tddfpt_spin_cons,&
                                              tddfpt_spin_flip,&
                                              tddfpt_triplet,&
                                              tddfpt_dipole_length,&
                                              tddfpt_dipole_velocity,&
                                              tddfpt_dipole_berry
   USE kinds,                           ONLY: dp
   USE mathlib,                         ONLY: get_diag
   USE moments_utils,                   ONLY: get_reference_point
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE physcon,                         ONLY: a_fine, wavenumbers
   USE preconditioner_types,            ONLY: destroy_preconditioner,&
                                              init_preconditioner,&
                                              preconditioner_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_mo_methods,                   ONLY: calculate_subspace_eigenvalues
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_type
   USE qs_operators_ao,                 ONLY: p_xyz_ao,&
                                              rRc_xyz_ao
   USE qs_ot_eigensolver,               ONLY: ot_eigensolver
   USE qs_tddfpt2_types,                ONLY: tddfpt_ground_state_mos
   USE xas_tdp_kernel,                  ONLY: kernel_coulomb_xc,&
                                              kernel_exchange
   USE xas_tdp_types,                   ONLY: donor_state_type,&
                                              xas_tdp_control_type,&
                                              xas_tdp_env_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                                     cp_logger_get_default_unit_nr,&
                                              cp_logger_type

!$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xas_tdp_utils'

   PUBLIC :: rcs_amew_soc_elements, soc_dipole_operator 

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

CONTAINS

! **************************************************************************************************
!> \brief Computes the rcs SOC matrix elements between excited states AMEWs based on the LR orbitals
!> \param amew_soc output dbcsr matrix with the SOC in the AMEW basis (needs to be fully resereved)
!> \param lr_soc dbcsr matrix with the SOC wrt the LR orbitals
!> \param lr_overlap dbcsr matrix with the excited states LR orbital overlap
!> \param domo_soc the SOC in the basis of the donor MOs
!> \param pref_trace see notes
!> \param pref_overall see notes
!> \param pref_diags see notes
!> \param symmetric if the outcome is known to be symmetric, only elements with iex <= jex are done
!> \note For an excited states pair i,j, the AMEW SOC matrix element is:
!>       soc_ij = pref_overall*(SUM(diag(lr_soc_ij)) + pref_trace*SUM(lr_overlap_ij*TRANSPOSE(domo_soc)))
!>       optionally, the value pref_diags*SUM(diag(lr_overlap_ij)) can be added (before pref_overall)
! **************************************************************************************************
   SUBROUTINE rcs_amew_soc_elements(amew_soc, lr_soc, lr_overlap, mo_soc, pref_trace, &
                                    pref_overall, pref_diags, symmetric)

      TYPE(dbcsr_type)                                   :: amew_soc, lr_soc, lr_overlap
      REAL(dp), DIMENSION(:, :)                          :: mo_soc
      REAL(dp)                                           :: pref_trace, pref_overall
      REAL(dp), OPTIONAL                                 :: pref_diags
      LOGICAL, OPTIONAL                                  :: symmetric
      INTEGER                                            :: blk, iex, jex
      LOGICAL                                            :: do_diags, found, my_symm
      REAL(dp)                                           :: soc_elem
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: diag
      REAL(dp), DIMENSION(:, :), POINTER                 :: pblock
      TYPE(dbcsr_iterator_type)                          :: iter
             

      ALLOCATE (diag(SIZE(mo_soc, 1)))
      my_symm = .FALSE.      
      IF (PRESENT(symmetric)) my_symm = symmetric
      do_diags = .FALSE.
      IF (PRESENT(pref_diags)) do_diags = .TRUE.


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
            soc_elem = soc_elem + SUM(diag)
         END IF

         CALL dbcsr_get_block_p(lr_overlap, iex, jex, pblock, found)
         IF (found) THEN
            soc_elem = soc_elem + pref_trace*SUM(pblock*TRANSPOSE(mo_soc))

          IF (do_diags) THEN
             diag(:) = get_diag(pblock)
             soc_elem = soc_elem + pref_diags*SUM(diag)
          END IF
         END IF

         CALL dbcsr_get_block_p(amew_soc, iex, jex, pblock, found)
         pblock = pref_overall*soc_elem
      END DO
      CALL dbcsr_iterator_stop(iter)
   
      DEALLOCATE(diag)

   END SUBROUTINE rcs_amew_soc_elements

   SUBROUTINE soc_dipole_operator(dipmat, tddfpt_control, qs_env)
           TYPE(dbcsr_p_type),DIMENSION(:),POINTER,INTENT(OUT)  :: dipmat
           TYPE(tddfpt2_control_type), POINTER, INTENT(IN)      :: tddfpt_control
           TYPE(qs_environment_type), POINTER, INTENT(IN)       :: qs_env

           CHARACTER(len=*), PARAMETER                          :: routineN= "soc_dipole_operator"

           INTEGER                                              :: handle, nao, nspin, ispin, iao, dim_op, i_dim
           TYPE(dbcsr_p_type), DIMENSION(:), POINTER            :: matrix_s
           REAL(kind=dp), DIMENSION(3)                          :: reference_point

           CALL timeset(routineN, handle)

           NULLIFY (matrix_s)

           IF (tddfpt_control%dipole_form == tddfpt_dipole_berry) CPABORT("BERRY DIPOLE FORM NOT IMPLEMENTED FOR SOC")

           !! ONLY RCS have been implemented, Therefore, nspin sould always be 1! 
           nspin = 1
           !! Number of dimensions should be 3, unless multipole is implemented in the future
           dim_op = 3


           !! Initzilize the dipmat structure
           CALL get_qs_env(qs_env, matrix_s=matrix_s)
           CALL dbcsr_get_info(matrix_s(1)%matrix,nfullrows_total=nao)

           NULLIFY(dipmat)
           ALLOCATE(dipmat(dim_op))
           DO i_dim = 1,dim_op
                ALLOCATE(dipmat(i_dim)%matrix)
                CALL dbcsr_copy(dipmat(i_dim)%matrix,&
                                matrix_s(1)%matrix,&
                                name="dipole operator matrix")
                CALL dbcsr_set(dipmat(i_dim)%matrix, 0.0_dp) !!This line needs to go for qs_tddfpt_properties like ao_op
           END DO

           SELECT CASE (tddfpt_control%dipole_form)
           CASE (tddfpt_dipole_length)
                   !!This routine is analog to qs_tddfpt_prperties but only until the rRc_xyz_ao routine
                   !! This will lead to an operator within the nao x nao basis
                   !! qs_tddpft_properies uses nvirt x nocc
                   CALL get_reference_point(reference_point, qs_env=qs_env, &
                                  reference=tddfpt_control%dipole_reference, &
                                  ref_point=tddfpt_control%dipole_ref_point)

                   CALL rRc_xyz_ao(op=dipmat, qs_env=qs_env, rc=reference_point, order=1, &
                         minimum_image=.FALSE., soft=.FALSE.)
          CASE (tddfpt_dipole_velocity)
                  !!This Routine calcluates the dipole OPerator within the velocity-form within the ao basis
                  !!This Operation is only used in xas_tdp and qs_tddfpt_soc, lines uses rmc_x_p_xyz_ao
                  CALL p_xyz_ao(dipmat, qs_env, minimum_image=.FALSE.)
          CASE DEFAULT
                  CPABORT("Unimplemented form of the dipole operator")
          END SELECT


          CALL timestop(handle)
   
   END SUBROUTINE soc_dipole_operator
END MODULE qs_tddfpt2_soc_utils

