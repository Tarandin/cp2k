!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
MODULE qs_tddfpt2_soc_types
   USE cp_array_utils,                  ONLY: cp_1d_i_p_type,&
                                              cp_1d_r_p_type,&
                                              cp_2d_i_p_type,&
                                              cp_2d_r_p_type,&
                                              cp_3d_r_p_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_fm_pool_types,                ONLY: fm_pool_create_fm
   USE cp_fm_struct,                    ONLY: cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: dbcsr_copy,&
                                              dbcsr_create,&
                                              dbcsr_p_type,&
                                              dbcsr_set, &
                                              dbcsr_release
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE libint_wrapper,                  ONLY: cp_libint_static_cleanup
   USE qs_mo_types,                     ONLY: mo_set_type
   USE qs_grid_atom,                    ONLY: deallocate_grid_atom,&
                                              grid_atom_type
   USE qs_harmonics_atom,               ONLY: deallocate_harmonics_atom,&
                                              harmonics_atom_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              get_qs_kind_set,&
                                              qs_kind_type

   IMPLICIT NONE

   PUBLIC  :: soc_control_type, soc_env_type, soc_env_create, &
               soc_env_release,  soc_atom_env_type, &
               soc_atom_create, soc_atom_release

   PRIVATE

      TYPE soc_control_type
         LOGICAL                                         :: do_soc
      END TYPE
!! Helper copied from xas_tdp_types
   TYPE grid_atom_p_type
      TYPE(grid_atom_type), POINTER                   :: grid_atom
   END TYPE grid_atom_p_type

   TYPE harmonics_atom_p_type
      TYPE(harmonics_atom_type), POINTER              :: harmonics_atom
   END TYPE harmonics_atom_p_type



!*************************************************************************************************
! \brief
! \param
!*************************************************************************************************
      TYPE soc_env_type
         !! a :: singlett or spin-conservm b :: triplett or spin flip
         TYPE(dbcsr_p_type), DIMENSION(:), POINTER       :: orb_soc
         TYPE(dbcsr_p_type), DIMENSION(:), POINTER       :: dipmat
         REAL(dp), POINTER, DIMENSION(:)                 :: evals_a, &
                                                            evals_b
         TYPE(cp_fm_type)                                :: a_coeff, &
                                                            b_coeff
         REAL(dp), ALLOCATABLE, DIMENSION(:)             :: soc_evals, &
                                                            soc_osc
      END TYPE

! ************************************************************************************************
! \bief an environment type analog to the xas_atom_env in
! xas_tdp_types:343  All unused parameters have been droped
! \param orb_sphi_so contains the coefficient for direct contraction from so to sgf, for the orb basis
! \param ri_sphi_so contains the coefficient for direct contraction from
! so to sgf, for the ri basis. NOT IN USE YET!!! 
!*************************************************************************************************
      TYPE soc_atom_env_type
         INTEGER                                           :: nspins
         TYPE(grid_atom_p_type), DIMENSION(:), POINTER     :: grid_atom_set
         TYPE(harmonics_atom_p_type), DIMENSION(:), &
         POINTER                                           :: harmonics_atom_set
         TYPE(cp_2d_r_p_type), DIMENSION(:), POINTER       :: orb_sphi_so
      END TYPE

CONTAINS

      SUBROUTINE soc_env_create(soc_env)
         TYPE(soc_env_type), TARGET                      :: soc_env

         NULLIFY (soc_env%orb_soc)
         NULLIFY (soc_env%dipmat)
         NULLIFY (soc_env%evals_a)
         NULLIFY (soc_env%evals_b)

      END SUBROUTINE soc_env_create

      SUBROUTINE soc_env_release(soc_env)
         TYPE(soc_env_type), TARGET                     :: soc_env
         INTEGER                                        :: i

         IF (ASSOCIATED(soc_env%orb_soc)) THEN
            DO i = 1, SIZE(soc_env%orb_soc)
               CALL dbcsr_release(soc_env%orb_soc(i)%matrix)
               DEALLOCATE (soc_env%orb_soc(i)%matrix)
            END DO
            DEALLOCATE (soc_env%orb_soc)
         END IF

         CALL cp_fm_release(soc_env%a_coeff)
         CALL cp_fm_release(soc_env%b_coeff)
         IF (ASSOCIATED(soc_env%evals_a)) NULLIFY(soc_env%evals_a)
         IF (ASSOCIATED(soc_env%evals_b)) NULLIFY(soc_env%evals_b)
         IF (ASSOCIATED(soc_env%dipmat)) NULLIFY(soc_env%dipmat)
         IF (ALLOCATED(soc_env%soc_evals)) DEALLOCATE(soc_env%soc_evals)
         IF (ALLOCATED(soc_env%soc_osc))  DEALLOCATE(soc_env%soc_osc) 

      END SUBROUTINE soc_env_release

! ************************************************************************************************
!*************************************************************************************************


      SUBROUTINE  soc_atom_create(soc_atom_env)
         TYPE (soc_atom_env_type), POINTER               :: soc_atom_env

         ALLOCATE(soc_atom_env)

         soc_atom_env%nspins = 1
         NULLIFY(soc_atom_env%grid_atom_set)
         NULLIFY(soc_atom_env%harmonics_atom_set)
         NULLIFY(soc_atom_env%orb_sphi_so)

      END SUBROUTINE soc_atom_create

      SUBROUTINE soc_atom_release(soc_atom_env)
         TYPE(soc_atom_env_type), POINTER                 :: soc_atom_env
         INTEGER                                          :: i, j, k

         IF (ASSOCIATED(soc_atom_env%grid_atom_set)) THEN
            DO i = 1, SIZE(soc_atom_env%grid_atom_set)
               IF (ASSOCIATED(soc_atom_env%grid_atom_set(i)%grid_atom)) THEN
                  CALL deallocate_grid_atom(soc_atom_env%grid_atom_set(i)%grid_atom)
               END IF
            END DO
            DEALLOCATE (soc_atom_env%grid_atom_set)
         END IF

         IF (ASSOCIATED(soc_atom_env%harmonics_atom_set)) THEN
            DO i = 1, SIZE(soc_atom_env%harmonics_atom_set)
               IF (ASSOCIATED(soc_atom_env%harmonics_atom_set(i)%harmonics_atom)) THEN
                  CALL deallocate_harmonics_atom(soc_atom_env%harmonics_atom_set(i)%harmonics_atom)
               END IF
            END DO
            DEALLOCATE (soc_atom_env%harmonics_atom_set)
         END IF

         IF (ASSOCIATED(soc_atom_env%orb_sphi_so)) THEN
            DO i = 1, SIZE(soc_atom_env%orb_sphi_so)
               IF (ASSOCIATED(soc_atom_env%orb_sphi_so(i)%array)) THEN
                  DEALLOCATE (soc_atom_env%orb_sphi_so(i)%array)
               END IF
            END DO
            DEALLOCATE (soc_atom_env%orb_sphi_so)
         END IF

         !Clean-up libint
         CALL cp_libint_static_cleanup()

         DEALLOCATE (soc_atom_env)

      END SUBROUTINE soc_atom_release

END MODULE qs_tddfpt2_soc_types
