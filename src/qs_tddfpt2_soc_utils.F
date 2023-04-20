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
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_add, dbcsr_copy, dbcsr_create, dbcsr_distribution_get, dbcsr_distribution_new, &
        dbcsr_distribution_release, dbcsr_distribution_type, dbcsr_finalize, dbcsr_get_block_p, &
        dbcsr_get_info, dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, &
        dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_multiply, &
        dbcsr_p_type, dbcsr_put_block, dbcsr_release, dbcsr_reserve_all_blocks, dbcsr_set, &
        dbcsr_type, dbcsr_type_no_symmetry, dbcsr_type_symmetric, dbcsr_print
   USE input_constants,                 ONLY: ot_precond_full_single,&
                                              tddfpt_singlet,&
                                              tddfpt_spin_cons,&
                                              tddfpt_spin_flip,&
                                              tddfpt_triplet,&
                                              xas_dip_len
   USE kinds,                           ONLY: dp
   USE mathlib,                         ONLY: get_diag
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
   USE qs_ot_eigensolver,               ONLY: ot_eigensolver
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

   PUBLIC :: rcs_amew_soc_elements 

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
      my_symm = .FALSE.


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

   END SUBROUTINE rcs_amew_soc_elements
END MODULE qs_tddfpt2_soc_utils

