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
   USE cp_control_types,                ONLY: tddfpt2_control_type
   USE cp_cfm_types,                    ONLY: cp_cfm_get_submatrix,&
                                              cp_cfm_type,&
                                              cp_cfm_set_submatrix,&
                                              cp_cfm_get_info

   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr,&
                                              cp_dbcsr_sm_fm_multiply,&
                                              dbcsr_allocate_matrix_set
   USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm_submat,&
                                              cp_fm_create,&
                                              cp_fm_release,&
                                              cp_fm_type
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE dbcsr_api,                       ONLY: dbcsr_copy,&
                                              dbcsr_get_info,&
                                              dbcsr_p_type,&
                                              dbcsr_create,&
                                              dbcsr_release,&
                                              dbcsr_desymmetrize,&
                                              dbcsr_type
   USE input_constants,                 ONLY: tddfpt_dipole_berry,&
                                              tddfpt_dipole_length,&
                                              tddfpt_dipole_velocity
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE moments_utils,                   ONLY: get_reference_point
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_operators_ao,                 ONLY: p_xyz_ao,&
                                              rRc_xyz_ao
   USE qs_tddfpt2_types,                ONLY: tddfpt_ground_state_mos
   USE qs_tddfpt2_soc_types,            ONLY: soc_env_type

!$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_tddfpt2_soc_utils'

   PUBLIC :: soc_dipole_operator, soc_contract_evect, resort_evects

   !A helper type for SOC
   TYPE dbcsr_soc_package_type
      TYPE(dbcsr_type), POINTER     :: dbcsr_sg => Null()
      TYPE(dbcsr_type), POINTER     :: dbcsr_tp => Null()
      TYPE(dbcsr_type), POINTER     :: dbcsr_sc => Null()
      TYPE(dbcsr_type), POINTER     :: dbcsr_sf => Null()
      TYPE(dbcsr_type), POINTER     :: dbcsr_prod => Null()
      TYPE(dbcsr_type), POINTER     :: dbcsr_ovlp => Null()
      TYPE(dbcsr_type), POINTER     :: dbcsr_tmp => Null()
      TYPE(dbcsr_type), POINTER     :: dbcsr_work => Null()
   END TYPE dbcsr_soc_package_type

CONTAINS

! **************************************************************************************************
!> \brief Build the atomic dipole operator
!> \param dipmat dipole operator
!> \param tddfpt_control informations on how to build the operaot
!> \param qs_env Qucikstep environment
! **************************************************************************************************
   SUBROUTINE soc_dipole_operator(soc_env, tddfpt_control, qs_env, gs_mos)
      TYPE(soc_env_type), TARGET                         :: soc_env
      !TYPE(dbcsr_p_type), DIMENSION(:),POINTER           :: dipmat
      TYPE(tddfpt2_control_type), POINTER                :: tddfpt_control
      TYPE(qs_environment_type), INTENT(IN), POINTER     :: qs_env
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         INTENT(in)                                      :: gs_mos

      CHARACTER(len=*), PARAMETER :: routineN = 'soc_dipole_operator'

      INTEGER                                            :: dim_op, handle, i_dim, nao, nspin
      REAL(kind=dp), DIMENSION(3)                        :: reference_point
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s

      CALL timeset(routineN, handle)

      NULLIFY (matrix_s)

      IF (tddfpt_control%dipole_form == tddfpt_dipole_berry) THEN
         CPABORT("BERRY DIPOLE FORM NOT IMPLEMENTED FOR SOC")
      END IF
           !! ONLY RCS have been implemented, Therefore, nspin sould always be 1!
      nspin = 1
           !! Number of dimensions should be 3, unless multipole is implemented in the future
      dim_op = 3
      !dipmat => soc_env%dipmat_ao

           !! Initzilize the dipmat structure
      CALL get_qs_env(qs_env, matrix_s=matrix_s)
      CALL dbcsr_get_info(matrix_s(1)%matrix, nfullrows_total=nao)

      !NULLIFY (dipmat)
      ALLOCATE (soc_env%dipmat_ao(dim_op))
      DO i_dim = 1, dim_op
         ALLOCATE (soc_env%dipmat_ao(i_dim)%matrix)
         CALL dbcsr_copy(soc_env%dipmat_ao(i_dim)%matrix, &
                         matrix_s(1)%matrix, &
                         name="dipole operator matrix")
         !CALL dbcsr_set(dipmat(i_dim)%matrix, 0.0_dp) !!This line needs to go for qs_tddfpt_properties like ao_op
      END DO

      SELECT CASE (tddfpt_control%dipole_form)
      CASE (tddfpt_dipole_length)
                   !!This routine is analog to qs_tddfpt_prperties but only until the rRc_xyz_ao routine
                   !! This will lead to an operator within the nao x nao basis
                   !! qs_tddpft_properies uses nvirt x nocc
         CALL get_reference_point(reference_point, qs_env=qs_env, &
                                  reference=tddfpt_control%dipole_reference, &
                                  ref_point=tddfpt_control%dipole_ref_point)

         CALL rRc_xyz_ao(op=soc_env%dipmat_ao, qs_env=qs_env, rc=reference_point, order=1, &
                         minimum_image=.FALSE., soft=.FALSE.)
         !! This will lead to S C^virt C^virt,T Q_q (vgl Strand et al., J. Chem Phys. 150, 044702, 2019)
         CALL length_rep(qs_env, gs_mos, soc_env)
      CASE (tddfpt_dipole_velocity)
                  !!This Routine calcluates the dipole OPerator within the velocity-form within the ao basis
                  !!This Operation is only used in xas_tdp and qs_tddfpt_soc, lines uses rmc_x_p_xyz_ao
         CALL p_xyz_ao(soc_env%dipmat_ao, qs_env, minimum_image=.FALSE.)
      CASE DEFAULT
         CPABORT("Unimplemented form of the dipole operator")
      END SELECT

      CALL timestop(handle)

   END SUBROUTINE soc_dipole_operator

   SUBROUTINE length_rep(qs_env, gs_mos, soc_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         INTENT(in)                                      :: gs_mos
      TYPE(soc_env_type), TARGET                         :: soc_env
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: S_mos_virt
      TYPE(cp_fm_type), POINTER                          :: dipmat_tmp, wfm_ao_ao
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: dipole_op_mos_occ
      TYPE(cp_fm_struct_type), POINTER                   :: dip_struct, fm_struct
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(dbcsr_type), POINTER                          :: symm_tmp
      TYPE(mp_para_env_type), POINTER                    :: para_env
      INTEGER                                            :: nao, nderivs, ideriv, ispin, nspins
      INTEGER,DIMENSION(:), ALLOCATABLE                  :: nmo_occ, nmo_virt

      NULLIFY(matrix_s, para_env, blacs_env)

      CALL get_qs_env(qs_env,matrix_s=matrix_s, blacs_env=blacs_env, para_env=para_env)

      nderivs = 3
      nspins = 1  !!We only account for rcs, will be changed in the future
      CALL dbcsr_get_info(matrix_s(1)%matrix, nfullrows_total=nao)
      ALLOCATE (S_mos_virt(nspins),dipole_op_mos_occ(3, nspins), nmo_occ(nspins),&
                wfm_ao_ao, nmo_virt(nspins), symm_tmp, dipmat_tmp)

      CALL cp_fm_struct_create(dip_struct, context=blacs_env,  ncol_global=nao, nrow_global=nao, para_env=para_env)
      
      CALL dbcsr_allocate_matrix_set(soc_env%dipmat, 3)
      CALL dbcsr_desymmetrize(matrix_s(1)%matrix, symm_tmp)
      DO ideriv=1,nderivs
         ALLOCATE(soc_env%dipmat(ideriv)%matrix)
         CALL dbcsr_create(soc_env%dipmat(ideriv)%matrix, template=symm_tmp,&
                 name="contracted operator", matrix_type="N")
         DO ispin = 1, nspins
            CALL cp_fm_create(dipole_op_mos_occ(ideriv, ispin), matrix_struct=dip_struct)
         END DO
      END DO

      CALL dbcsr_release(symm_tmp)
      DEALLOCATE(symm_tmp)

      DO ispin = 1, nspins
         nmo_virt(ispin) = SIZE(gs_mos(ispin)%evals_virt)
         CALL cp_fm_get_info(gs_mos(ispin)%mos_virt, matrix_struct=fm_struct)
         CALL cp_fm_create(wfm_ao_ao, dip_struct)
         CALL cp_fm_create(S_mos_virt(ispin), fm_struct)

         CALL cp_dbcsr_sm_fm_multiply(matrix_s(1)%matrix, &
                                      gs_mos(ispin)%mos_virt, &
                                      S_mos_virt(ispin), &
                                      ncol=nmo_virt(ispin), alpha=1.0_dp, beta=0.0_dp)
         CALL parallel_gemm('N', 'T', nao, nao, nmo_virt(ispin), &
                               1.0_dp, S_mos_virt(ispin), gs_mos(ispin)%mos_virt, &
                               0.0_dp, wfm_ao_ao)

         DO ideriv = 1,nderivs
                CALL cp_fm_create(dipmat_tmp,dip_struct)
                CALL copy_dbcsr_to_fm(soc_env%dipmat_ao(ideriv)%matrix,dipmat_tmp)
                CALL parallel_gemm('N', 'T', nao, nao, nao, &
                                  1.0_dp, wfm_ao_ao, dipmat_tmp, &
                                  0.0_dp, dipole_op_mos_occ(ideriv, ispin))
                CALL copy_fm_to_dbcsr(dipole_op_mos_occ(ideriv, ispin),soc_env%dipmat(ideriv)%matrix)
                CALL cp_fm_release(dipmat_tmp)
         END DO
         CALL cp_fm_release(wfm_ao_ao)
         DEALLOCATE(wfm_ao_ao)
      END DO

      CALL cp_fm_struct_release(dip_struct)
      DO ispin = 1,nspins
        CALL cp_fm_release(S_mos_virt(ispin))
        DO ideriv = 1, nderivs
                CALL cp_fm_release(dipole_op_mos_occ(ideriv, ispin))
        END DO
      END DOveloclity_rep
      DEALLOCATE(S_mos_virt, dipole_op_mos_occ, nmo_virt, dipmat_tmp)


   END SUBROUTINE length_rep

   SUBROUTINE veloclity_rep(qs_env, gs_mos, soc_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         INTENT(in)                                      :: gs_mos
      TYPE(soc_env_type), TARGET                         :: soc_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: scrm
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      INTEGER                                            :: ispin, nspins, ideriv, nderivs

      nspins = 1
      nderivs = 3

      CALL get_qs_env(qs_env, ks_env=ks_env, sab_orb=sab_orb)
      NULLIFY (scrm)
      CALL build_overlap_matrix(ks_env, matrix_s=scrm, nderivative=1, &
                                   basis_type_a="ORB", basis_type_b="ORB", &
                                   sab_nl=sab_orb)

      DO ispin= 1,nspins
      END DO

   END SUBROUTINE veloclity_rep

! **************************************************************************************************
!> \brief ...
!> \param fm_start ...
!> \param fm_res ...
! **************************************************************************************************
   SUBROUTINE soc_contract_evect(fm_start, fm_res)

      TYPE(cp_fm_type), DIMENSION(:, :), INTENT(in)      :: fm_start
      TYPE(cp_fm_type), INTENT(inout)                    :: fm_res

      CHARACTER(len=*), PARAMETER :: routineN = 'soc_contract_evect'

      INTEGER                                            :: handle, ii, jj, nactive, nao, nspins, &
                                                            nstates, ntmp1, ntmp2

      CALL timeset(routineN, handle)

      nstates = SIZE(fm_start, 2)
      nspins = SIZE(fm_start, 1)

      CALL cp_fm_set_all(fm_res, 0.0_dp)
         !! Evects are written into one matrix.
      DO ii = 1, nstates
         DO jj = 1, nspins
            CALL cp_fm_get_info(fm_start(jj, ii), nrow_global=nao, ncol_global=nactive)
            CALL cp_fm_get_info(fm_res, nrow_global=ntmp1, ncol_global=ntmp2)
            CALL cp_fm_to_fm_submat(fm_start(jj, ii), &
                                    fm_res, &
                                    nao, nactive, &
                                    1, 1, 1, &
                                    1 + nactive*(ii - 1) + (jj - 1)*nao*nstates)
         END DO !nspins
      END DO !nsstates

      CALL timestop(handle)

   END SUBROUTINE soc_contract_evect

   SUBROUTINE test_repetition(vec, new_entry, res, res_int)
     INTEGER, DIMENSION(:), INTENT(IN)              :: vec
     INTEGER, INTENT(IN)                            :: new_entry
     LOGICAL, INTENT(OUT)                           :: res
     INTEGER, INTENT(OUT), OPTIONAL                 :: res_int
     INTEGER                                        :: i

     res = .TRUE.
     IF (PRESENT(res_int)) res_int = -1
 
     i = 1
     DO while (vec(i) > 0 .and. i .le. SIZE(vec)) 
        IF (vec(i) == new_entry) THEN 
                res = .FALSE. 
                IF (PRESENT(res_int)) res_int = i
                i = SIZE(vec)
        END IF
        i = i + 1
     END DO

   END SUBROUTINE test_repetition

   SUBROUTINE resort_evects(evects_cfm, sort, cpl_tmp)
           TYPE(cp_cfm_type), INTENT(INOUT)         :: evects_cfm
           INTEGER, DIMENSION(:), ALLOCATABLE, &
                   INTENT(OUT)                      :: sort
           COMPLEX(dp), DIMENSION(:,:), ALLOCATABLE :: cpl_tmp
           REAL(dp)                                 :: wfn_sq, max_wfn, max_dev
           INTEGER                                  :: ii, jj, tmp, ntot, i_rep
           INTEGER, DIMENSION(:), ALLOCATABLE       :: rep_int
           LOGICAL                                  :: rep

           CALL cp_cfm_get_info(evects_cfm, nrow_global=ntot)
           ALLOCATE (cpl_tmp(ntot,ntot))
           ALLOCATE (sort(ntot), rep_int(ntot))
           cpl_tmp = 0_dp
           sort = 0
           max_dev = 0.5
           CALL cp_cfm_get_submatrix(evects_cfm, cpl_tmp)

           DO jj = 1,ntot
             rep_int= 0
             tmp = 0
             max_wfn = 0_dp
             DO ii = 1, ntot
                wfn_sq = ABS(REAL(cpl_tmp(ii,jj)**2-AIMAG(cpl_tmp(ii,jj)**2)))
                IF (max_wfn .le. wfn_sq) THEN
                        CALL test_repetition(sort,ii, rep, rep_int(ii))
                        IF (rep) THEN 
                                max_wfn = wfn_sq
                                tmp = ii
                        END IF
                END IF
             END DO
             IF (tmp > 0) THEN
                     sort(jj) = tmp
             ELSE
                     DO i_rep = 1,ntot
                        IF (rep_int(i_rep) > 0) THEN
                           max_wfn = ABS(REAL(cpl_tmp(sort(i_rep),jj)**2-AIMAG(cpl_tmp(sort(i_rep),jj)**2)))-max_dev
                           DO ii = 1,ntot
                              wfn_sq = ABS(REAL(cpl_tmp(ii,jj)**2-AIMAG(cpl_tmp(ii,jj)**2)))
                              IF ((max_wfn - wfn_sq)/max_wfn .le. max_dev) THEN
                                      CALL test_repetition(sort,ii, rep)
                                      IF (rep .and. ii /= i_rep) THEN
                                              sort(jj) = sort(i_rep)
                                              sort(i_rep) = ii
                                      END IF
                              END IF
                           END DO
                        ENDIF
                     END DO
             END IF
           END DO

           DO ii=1,ntot
             CALL cp_cfm_set_submatrix(evects_cfm, cpl_tmp(1:ntot,sort(ii):sort(ii)),1,sort(ii), ntot, 1)
           END DO
           CALL cp_cfm_get_submatrix(evects_cfm, cpl_tmp)
           !DEALLOCATE(cpl_tmp)

   END SUBROUTINE resort_evects
END MODULE qs_tddfpt2_soc_utils
