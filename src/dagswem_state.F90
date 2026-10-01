MODULE DAGSWEM_STATE
  USE SIZES
  USE GLOBAL
  USE DG
  USE NodalAttributes, ONLY: SwanWaveRefrac, STARTDRY, FRIC, TAU0VAR, TAU0BASE, &
                             z0land, vcanopy, BridgePilings, Chezy, ManningsN, &
                             GeoidOffset, EVM, EVC
  USE FSTARPU_MOD
  USE ISO_C_BINDING
  IMPLICIT NONE

  INTEGER, PARAMETER :: NUM_STATE_HANDLES = 516
  INTEGER, PARAMETER :: NUM_INT_SCALARS = 847
  INTEGER, PARAMETER :: NUM_REAL_SCALARS = 604
  REAL(SZ), TARGET, SAVE :: DUMMY_BUF(1) = 0.0_SZ

  CONTAINS

  SUBROUTINE DGSWEM_STATE_REGISTER(handles)
    TYPE(C_PTR), INTENT(OUT) :: handles(NUM_STATE_HANDLES)
    INTEGER, POINTER :: SCALAR_INT_BUF(:)
    REAL(SZ), POINTER :: SCALAR_REAL_BUF(:)
    ALLOCATE(SCALAR_INT_BUF(NUM_INT_SCALARS))
    ALLOCATE(SCALAR_REAL_BUF(NUM_REAL_SCALARS))
    IF (ASSOCIATED(WDFLG)) THEN
      SCALAR_INT_BUF(1) = 1
      CALL fstarpu_vector_data_register(handles(1), 0, C_LOC(WDFLG(LBOUND(WDFLG,1))), SIZE(WDFLG,1), C_SIZEOF(WDFLG(LBOUND(WDFLG,1))))
    ELSE
      SCALAR_INT_BUF(1) = 0
      CALL fstarpu_vector_data_register(handles(1), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(WDFLG)
    IF (ASSOCIATED(WDFLG_TMP)) THEN
      SCALAR_INT_BUF(2) = 1
      CALL fstarpu_vector_data_register(handles(2), 0, C_LOC(WDFLG_TMP(LBOUND(WDFLG_TMP,1))), SIZE(WDFLG_TMP,1), C_SIZEOF(WDFLG_TMP(LBOUND(WDFLG_TMP,1))))
    ELSE
      SCALAR_INT_BUF(2) = 0
      CALL fstarpu_vector_data_register(handles(2), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(WDFLG_TMP)
    IF (ASSOCIATED(DOFW)) THEN
      SCALAR_INT_BUF(3) = 1
      CALL fstarpu_vector_data_register(handles(3), 0, C_LOC(DOFW(LBOUND(DOFW,1))), SIZE(DOFW,1), C_SIZEOF(DOFW(LBOUND(DOFW,1))))
    ELSE
      SCALAR_INT_BUF(3) = 0
      CALL fstarpu_vector_data_register(handles(3), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(DOFW)
    IF (ASSOCIATED(EL_UPDATED)) THEN
      SCALAR_INT_BUF(4) = 1
      CALL fstarpu_vector_data_register(handles(4), 0, C_LOC(EL_UPDATED(LBOUND(EL_UPDATED,1))), SIZE(EL_UPDATED,1), C_SIZEOF(EL_UPDATED(LBOUND(EL_UPDATED,1))))
    ELSE
      SCALAR_INT_BUF(4) = 0
      CALL fstarpu_vector_data_register(handles(4), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(EL_UPDATED)
    IF (ASSOCIATED(LEDGE_NVEC)) THEN
      SCALAR_INT_BUF(5) = 1
      CALL fstarpu_block_data_register(handles(5), 0, C_LOC(LEDGE_NVEC(LBOUND(LEDGE_NVEC,1),LBOUND(LEDGE_NVEC,2),LBOUND(LEDGE_NVEC,3))), SIZE(LEDGE_NVEC,1), SIZE(LEDGE_NVEC,1)*SIZE(LEDGE_NVEC,2), SIZE(LEDGE_NVEC,1), SIZE(LEDGE_NVEC,2), SIZE(LEDGE_NVEC,3), C_SIZEOF(LEDGE_NVEC(LBOUND(LEDGE_NVEC,1),LBOUND(LEDGE_NVEC,2),LBOUND(LEDGE_NVEC,3))))
    ELSE
      SCALAR_INT_BUF(5) = 0
      CALL fstarpu_block_data_register(handles(5), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(LEDGE_NVEC)
    IF (ASSOCIATED(DOFS)) THEN
      SCALAR_INT_BUF(6) = 1
      CALL fstarpu_vector_data_register(handles(6), 0, C_LOC(DOFS(LBOUND(DOFS,1))), SIZE(DOFS,1), C_SIZEOF(DOFS(LBOUND(DOFS,1))))
    ELSE
      SCALAR_INT_BUF(6) = 0
      CALL fstarpu_vector_data_register(handles(6), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(DOFS)
    IF (ASSOCIATED(PCOUNT)) THEN
      SCALAR_INT_BUF(7) = 1
      CALL fstarpu_vector_data_register(handles(7), 0, C_LOC(PCOUNT(LBOUND(PCOUNT,1))), SIZE(PCOUNT,1), C_SIZEOF(PCOUNT(LBOUND(PCOUNT,1))))
    ELSE
      SCALAR_INT_BUF(7) = 0
      CALL fstarpu_vector_data_register(handles(7), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(PCOUNT)
    IF (ASSOCIATED(PDG)) THEN
      SCALAR_INT_BUF(8) = 1
      CALL fstarpu_vector_data_register(handles(8), 0, C_LOC(PDG(LBOUND(PDG,1))), SIZE(PDG,1), C_SIZEOF(PDG(LBOUND(PDG,1))))
    ELSE
      SCALAR_INT_BUF(8) = 0
      CALL fstarpu_vector_data_register(handles(8), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(PDG)
    IF (ASSOCIATED(NCOUNT)) THEN
      SCALAR_INT_BUF(9) = 1
      CALL fstarpu_vector_data_register(handles(9), 0, C_LOC(NCOUNT(LBOUND(NCOUNT,1))), SIZE(NCOUNT,1), C_SIZEOF(NCOUNT(LBOUND(NCOUNT,1))))
    ELSE
      SCALAR_INT_BUF(9) = 0
      CALL fstarpu_vector_data_register(handles(9), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(NCOUNT)
    IF (ASSOCIATED(NEDEL)) THEN
      SCALAR_INT_BUF(10) = 1
      CALL fstarpu_matrix_data_register(handles(10), 0, C_LOC(NEDEL(LBOUND(NEDEL,1),LBOUND(NEDEL,2))), SIZE(NEDEL,1), SIZE(NEDEL,1), SIZE(NEDEL,2), C_SIZEOF(NEDEL(LBOUND(NEDEL,1),LBOUND(NEDEL,2))))
    ELSE
      SCALAR_INT_BUF(10) = 0
      CALL fstarpu_matrix_data_register(handles(10), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(NEDEL)
    IF (ASSOCIATED(NEDSD)) THEN
      SCALAR_INT_BUF(11) = 1
      CALL fstarpu_matrix_data_register(handles(11), 0, C_LOC(NEDSD(LBOUND(NEDSD,1),LBOUND(NEDSD,2))), SIZE(NEDSD,1), SIZE(NEDSD,1), SIZE(NEDSD,2), C_SIZEOF(NEDSD(LBOUND(NEDSD,1),LBOUND(NEDSD,2))))
    ELSE
      SCALAR_INT_BUF(11) = 0
      CALL fstarpu_matrix_data_register(handles(11), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(NEDSD)
    IF (ASSOCIATED(NEDNO)) THEN
      SCALAR_INT_BUF(12) = 1
      CALL fstarpu_matrix_data_register(handles(12), 0, C_LOC(NEDNO(LBOUND(NEDNO,1),LBOUND(NEDNO,2))), SIZE(NEDNO,1), SIZE(NEDNO,1), SIZE(NEDNO,2), C_SIZEOF(NEDNO(LBOUND(NEDNO,1),LBOUND(NEDNO,2))))
    ELSE
      SCALAR_INT_BUF(12) = 0
      CALL fstarpu_matrix_data_register(handles(12), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(NEDNO)
    IF (ASSOCIATED(NEDNO1)) THEN
      SCALAR_INT_BUF(13) = 1
      CALL fstarpu_vector_data_register(handles(13), 0, C_LOC(NEDNO1(LBOUND(NEDNO1,1))), SIZE(NEDNO1,1), C_SIZEOF(NEDNO1(LBOUND(NEDNO1,1))))
    ELSE
      SCALAR_INT_BUF(13) = 0
      CALL fstarpu_vector_data_register(handles(13), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(NEDNO1)
    IF (ASSOCIATED(NEDNO2)) THEN
      SCALAR_INT_BUF(14) = 1
      CALL fstarpu_vector_data_register(handles(14), 0, C_LOC(NEDNO2(LBOUND(NEDNO2,1))), SIZE(NEDNO2,1), C_SIZEOF(NEDNO2(LBOUND(NEDNO2,1))))
    ELSE
      SCALAR_INT_BUF(14) = 0
      CALL fstarpu_vector_data_register(handles(14), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(NEDNO2)
    IF (ASSOCIATED(NIEDN)) THEN
      SCALAR_INT_BUF(15) = 1
      CALL fstarpu_vector_data_register(handles(15), 0, C_LOC(NIEDN(LBOUND(NIEDN,1))), SIZE(NIEDN,1), C_SIZEOF(NIEDN(LBOUND(NIEDN,1))))
    ELSE
      SCALAR_INT_BUF(15) = 0
      CALL fstarpu_vector_data_register(handles(15), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(NIEDN)
    IF (ASSOCIATED(NLEDN)) THEN
      SCALAR_INT_BUF(16) = 1
      CALL fstarpu_vector_data_register(handles(16), 0, C_LOC(NLEDN(LBOUND(NLEDN,1))), SIZE(NLEDN,1), C_SIZEOF(NLEDN(LBOUND(NLEDN,1))))
    ELSE
      SCALAR_INT_BUF(16) = 0
      CALL fstarpu_vector_data_register(handles(16), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(NLEDN)
    IF (ASSOCIATED(NEEDN)) THEN
      SCALAR_INT_BUF(17) = 1
      CALL fstarpu_vector_data_register(handles(17), 0, C_LOC(NEEDN(LBOUND(NEEDN,1))), SIZE(NEEDN,1), C_SIZEOF(NEEDN(LBOUND(NEEDN,1))))
    ELSE
      SCALAR_INT_BUF(17) = 0
      CALL fstarpu_vector_data_register(handles(17), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(NEEDN)
    IF (ASSOCIATED(NFEDN)) THEN
      SCALAR_INT_BUF(18) = 1
      CALL fstarpu_vector_data_register(handles(18), 0, C_LOC(NFEDN(LBOUND(NFEDN,1))), SIZE(NFEDN,1), C_SIZEOF(NFEDN(LBOUND(NFEDN,1))))
    ELSE
      SCALAR_INT_BUF(18) = 0
      CALL fstarpu_vector_data_register(handles(18), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(NFEDN)
    IF (ASSOCIATED(NREDN)) THEN
      SCALAR_INT_BUF(19) = 1
      CALL fstarpu_vector_data_register(handles(19), 0, C_LOC(NREDN(LBOUND(NREDN,1))), SIZE(NREDN,1), C_SIZEOF(NREDN(LBOUND(NREDN,1))))
    ELSE
      SCALAR_INT_BUF(19) = 0
      CALL fstarpu_vector_data_register(handles(19), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(NREDN)
    IF (ASSOCIATED(NEBEDN)) THEN
      SCALAR_INT_BUF(20) = 1
      CALL fstarpu_vector_data_register(handles(20), 0, C_LOC(NEBEDN(LBOUND(NEBEDN,1))), SIZE(NEBEDN,1), C_SIZEOF(NEBEDN(LBOUND(NEBEDN,1))))
    ELSE
      SCALAR_INT_BUF(20) = 0
      CALL fstarpu_vector_data_register(handles(20), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(NEBEDN)
    IF (ASSOCIATED(NIBEDN)) THEN
      SCALAR_INT_BUF(21) = 1
      CALL fstarpu_vector_data_register(handles(21), 0, C_LOC(NIBEDN(LBOUND(NIBEDN,1))), SIZE(NIBEDN,1), C_SIZEOF(NIBEDN(LBOUND(NIBEDN,1))))
    ELSE
      SCALAR_INT_BUF(21) = 0
      CALL fstarpu_vector_data_register(handles(21), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(NIBEDN)
    IF (ASSOCIATED(NIBSEGN)) THEN
      SCALAR_INT_BUF(22) = 1
      CALL fstarpu_matrix_data_register(handles(22), 0, C_LOC(NIBSEGN(LBOUND(NIBSEGN,1),LBOUND(NIBSEGN,2))), SIZE(NIBSEGN,1), SIZE(NIBSEGN,1), SIZE(NIBSEGN,2), C_SIZEOF(NIBSEGN(LBOUND(NIBSEGN,1),LBOUND(NIBSEGN,2))))
    ELSE
      SCALAR_INT_BUF(22) = 0
      CALL fstarpu_matrix_data_register(handles(22), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(NIBSEGN)
    IF (ASSOCIATED(NEBSEGN)) THEN
      SCALAR_INT_BUF(23) = 1
      CALL fstarpu_vector_data_register(handles(23), 0, C_LOC(NEBSEGN(LBOUND(NEBSEGN,1))), SIZE(NEBSEGN,1), C_SIZEOF(NEBSEGN(LBOUND(NEBSEGN,1))))
    ELSE
      SCALAR_INT_BUF(23) = 0
      CALL fstarpu_vector_data_register(handles(23), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(NEBSEGN)
    IF (ASSOCIATED(EL_NBORS)) THEN
      SCALAR_INT_BUF(24) = 1
      CALL fstarpu_matrix_data_register(handles(24), 0, C_LOC(EL_NBORS(LBOUND(EL_NBORS,1),LBOUND(EL_NBORS,2))), SIZE(EL_NBORS,1), SIZE(EL_NBORS,1), SIZE(EL_NBORS,2), C_SIZEOF(EL_NBORS(LBOUND(EL_NBORS,1),LBOUND(EL_NBORS,2))))
    ELSE
      SCALAR_INT_BUF(24) = 0
      CALL fstarpu_matrix_data_register(handles(24), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(EL_NBORS)
    IF (ASSOCIATED(BACKNODES)) THEN
      SCALAR_INT_BUF(25) = 1
      CALL fstarpu_matrix_data_register(handles(25), 0, C_LOC(BACKNODES(LBOUND(BACKNODES,1),LBOUND(BACKNODES,2))), SIZE(BACKNODES,1), SIZE(BACKNODES,1), SIZE(BACKNODES,2), C_SIZEOF(BACKNODES(LBOUND(BACKNODES,1),LBOUND(BACKNODES,2))))
    ELSE
      SCALAR_INT_BUF(25) = 0
      CALL fstarpu_matrix_data_register(handles(25), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(BACKNODES)
    IF (ASSOCIATED(MARK)) THEN
      SCALAR_INT_BUF(26) = 1
      CALL fstarpu_vector_data_register(handles(26), 0, C_LOC(MARK(LBOUND(MARK,1))), SIZE(MARK,1), C_SIZEOF(MARK(LBOUND(MARK,1))))
    ELSE
      SCALAR_INT_BUF(26) = 0
      CALL fstarpu_vector_data_register(handles(26), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(MARK)
    IF (ASSOCIATED(ATVD)) THEN
      SCALAR_INT_BUF(27) = 1
      CALL fstarpu_matrix_data_register(handles(27), 0, C_LOC(ATVD(LBOUND(ATVD,1),LBOUND(ATVD,2))), SIZE(ATVD,1), SIZE(ATVD,1), SIZE(ATVD,2), C_SIZEOF(ATVD(LBOUND(ATVD,1),LBOUND(ATVD,2))))
    ELSE
      SCALAR_INT_BUF(27) = 0
      CALL fstarpu_matrix_data_register(handles(27), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(ATVD)
    IF (ASSOCIATED(BTVD)) THEN
      SCALAR_INT_BUF(28) = 1
      CALL fstarpu_matrix_data_register(handles(28), 0, C_LOC(BTVD(LBOUND(BTVD,1),LBOUND(BTVD,2))), SIZE(BTVD,1), SIZE(BTVD,1), SIZE(BTVD,2), C_SIZEOF(BTVD(LBOUND(BTVD,1),LBOUND(BTVD,2))))
    ELSE
      SCALAR_INT_BUF(28) = 0
      CALL fstarpu_matrix_data_register(handles(28), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(BTVD)
    IF (ASSOCIATED(CTVD)) THEN
      SCALAR_INT_BUF(29) = 1
      CALL fstarpu_matrix_data_register(handles(29), 0, C_LOC(CTVD(LBOUND(CTVD,1),LBOUND(CTVD,2))), SIZE(CTVD,1), SIZE(CTVD,1), SIZE(CTVD,2), C_SIZEOF(CTVD(LBOUND(CTVD,1),LBOUND(CTVD,2))))
    ELSE
      SCALAR_INT_BUF(29) = 0
      CALL fstarpu_matrix_data_register(handles(29), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(CTVD)
    IF (ASSOCIATED(DTVD)) THEN
      SCALAR_INT_BUF(30) = 1
      CALL fstarpu_vector_data_register(handles(30), 0, C_LOC(DTVD(LBOUND(DTVD,1))), SIZE(DTVD,1), C_SIZEOF(DTVD(LBOUND(DTVD,1))))
    ELSE
      SCALAR_INT_BUF(30) = 0
      CALL fstarpu_vector_data_register(handles(30), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(DTVD)
    IF (ASSOCIATED(MAX_BOA_DT)) THEN
      SCALAR_INT_BUF(31) = 1
      CALL fstarpu_vector_data_register(handles(31), 0, C_LOC(MAX_BOA_DT(LBOUND(MAX_BOA_DT,1))), SIZE(MAX_BOA_DT,1), C_SIZEOF(MAX_BOA_DT(LBOUND(MAX_BOA_DT,1))))
    ELSE
      SCALAR_INT_BUF(31) = 0
      CALL fstarpu_vector_data_register(handles(31), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(MAX_BOA_DT)
    IF (ASSOCIATED(e1)) THEN
      SCALAR_INT_BUF(32) = 1
      CALL fstarpu_vector_data_register(handles(32), 0, C_LOC(e1(LBOUND(e1,1))), SIZE(e1,1), C_SIZEOF(e1(LBOUND(e1,1))))
    ELSE
      SCALAR_INT_BUF(32) = 0
      CALL fstarpu_vector_data_register(handles(32), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(e1)
    IF (ASSOCIATED(balance)) THEN
      SCALAR_INT_BUF(33) = 1
      CALL fstarpu_vector_data_register(handles(33), 0, C_LOC(balance(LBOUND(balance,1))), SIZE(balance,1), C_SIZEOF(balance(LBOUND(balance,1))))
    ELSE
      SCALAR_INT_BUF(33) = 0
      CALL fstarpu_vector_data_register(handles(33), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(balance)
    IF (ASSOCIATED(RKC_T)) THEN
      SCALAR_INT_BUF(34) = 1
      CALL fstarpu_vector_data_register(handles(34), 0, C_LOC(RKC_T(LBOUND(RKC_T,1))), SIZE(RKC_T,1), C_SIZEOF(RKC_T(LBOUND(RKC_T,1))))
    ELSE
      SCALAR_INT_BUF(34) = 0
      CALL fstarpu_vector_data_register(handles(34), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(RKC_T)
    IF (ASSOCIATED(RKC_U)) THEN
      SCALAR_INT_BUF(35) = 1
      CALL fstarpu_vector_data_register(handles(35), 0, C_LOC(RKC_U(LBOUND(RKC_U,1))), SIZE(RKC_U,1), C_SIZEOF(RKC_U(LBOUND(RKC_U,1))))
    ELSE
      SCALAR_INT_BUF(35) = 0
      CALL fstarpu_vector_data_register(handles(35), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(RKC_U)
    IF (ASSOCIATED(RKC_Tprime)) THEN
      SCALAR_INT_BUF(36) = 1
      CALL fstarpu_vector_data_register(handles(36), 0, C_LOC(RKC_Tprime(LBOUND(RKC_Tprime,1))), SIZE(RKC_Tprime,1), C_SIZEOF(RKC_Tprime(LBOUND(RKC_Tprime,1))))
    ELSE
      SCALAR_INT_BUF(36) = 0
      CALL fstarpu_vector_data_register(handles(36), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(RKC_Tprime)
    IF (ASSOCIATED(RKC_Tdprime)) THEN
      SCALAR_INT_BUF(37) = 1
      CALL fstarpu_vector_data_register(handles(37), 0, C_LOC(RKC_Tdprime(LBOUND(RKC_Tdprime,1))), SIZE(RKC_Tdprime,1), C_SIZEOF(RKC_Tdprime(LBOUND(RKC_Tdprime,1))))
    ELSE
      SCALAR_INT_BUF(37) = 0
      CALL fstarpu_vector_data_register(handles(37), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(RKC_Tdprime)
    IF (ASSOCIATED(RKC_a)) THEN
      SCALAR_INT_BUF(38) = 1
      CALL fstarpu_vector_data_register(handles(38), 0, C_LOC(RKC_a(LBOUND(RKC_a,1))), SIZE(RKC_a,1), C_SIZEOF(RKC_a(LBOUND(RKC_a,1))))
    ELSE
      SCALAR_INT_BUF(38) = 0
      CALL fstarpu_vector_data_register(handles(38), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(RKC_a)
    IF (ASSOCIATED(RKC_b)) THEN
      SCALAR_INT_BUF(39) = 1
      CALL fstarpu_vector_data_register(handles(39), 0, C_LOC(RKC_b(LBOUND(RKC_b,1))), SIZE(RKC_b,1), C_SIZEOF(RKC_b(LBOUND(RKC_b,1))))
    ELSE
      SCALAR_INT_BUF(39) = 0
      CALL fstarpu_vector_data_register(handles(39), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(RKC_b)
    IF (ASSOCIATED(RKC_c)) THEN
      SCALAR_INT_BUF(40) = 1
      CALL fstarpu_vector_data_register(handles(40), 0, C_LOC(RKC_c(LBOUND(RKC_c,1))), SIZE(RKC_c,1), C_SIZEOF(RKC_c(LBOUND(RKC_c,1))))
    ELSE
      SCALAR_INT_BUF(40) = 0
      CALL fstarpu_vector_data_register(handles(40), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(RKC_c)
    IF (ASSOCIATED(RKC_mu)) THEN
      SCALAR_INT_BUF(41) = 1
      CALL fstarpu_vector_data_register(handles(41), 0, C_LOC(RKC_mu(LBOUND(RKC_mu,1))), SIZE(RKC_mu,1), C_SIZEOF(RKC_mu(LBOUND(RKC_mu,1))))
    ELSE
      SCALAR_INT_BUF(41) = 0
      CALL fstarpu_vector_data_register(handles(41), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(RKC_mu)
    IF (ASSOCIATED(RKC_tildemu)) THEN
      SCALAR_INT_BUF(42) = 1
      CALL fstarpu_vector_data_register(handles(42), 0, C_LOC(RKC_tildemu(LBOUND(RKC_tildemu,1))), SIZE(RKC_tildemu,1), C_SIZEOF(RKC_tildemu(LBOUND(RKC_tildemu,1))))
    ELSE
      SCALAR_INT_BUF(42) = 0
      CALL fstarpu_vector_data_register(handles(42), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(RKC_tildemu)
    IF (ASSOCIATED(RKC_nu)) THEN
      SCALAR_INT_BUF(43) = 1
      CALL fstarpu_vector_data_register(handles(43), 0, C_LOC(RKC_nu(LBOUND(RKC_nu,1))), SIZE(RKC_nu,1), C_SIZEOF(RKC_nu(LBOUND(RKC_nu,1))))
    ELSE
      SCALAR_INT_BUF(43) = 0
      CALL fstarpu_vector_data_register(handles(43), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(RKC_nu)
    IF (ASSOCIATED(RKC_gamma)) THEN
      SCALAR_INT_BUF(44) = 1
      CALL fstarpu_vector_data_register(handles(44), 0, C_LOC(RKC_gamma(LBOUND(RKC_gamma,1))), SIZE(RKC_gamma,1), C_SIZEOF(RKC_gamma(LBOUND(RKC_gamma,1))))
    ELSE
      SCALAR_INT_BUF(44) = 0
      CALL fstarpu_vector_data_register(handles(44), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(RKC_gamma)
    IF (ASSOCIATED(BATH)) THEN
      SCALAR_INT_BUF(45) = 1
      CALL fstarpu_block_data_register(handles(45), 0, C_LOC(BATH(LBOUND(BATH,1),LBOUND(BATH,2),LBOUND(BATH,3))), SIZE(BATH,1), SIZE(BATH,1)*SIZE(BATH,2), SIZE(BATH,1), SIZE(BATH,2), SIZE(BATH,3), C_SIZEOF(BATH(LBOUND(BATH,1),LBOUND(BATH,2),LBOUND(BATH,3))))
    ELSE
      SCALAR_INT_BUF(45) = 0
      CALL fstarpu_block_data_register(handles(45), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(BATH)
    IF (ASSOCIATED(DBATHDX)) THEN
      SCALAR_INT_BUF(46) = 1
      CALL fstarpu_block_data_register(handles(46), 0, C_LOC(DBATHDX(LBOUND(DBATHDX,1),LBOUND(DBATHDX,2),LBOUND(DBATHDX,3))), SIZE(DBATHDX,1), SIZE(DBATHDX,1)*SIZE(DBATHDX,2), SIZE(DBATHDX,1), SIZE(DBATHDX,2), SIZE(DBATHDX,3), C_SIZEOF(DBATHDX(LBOUND(DBATHDX,1),LBOUND(DBATHDX,2),LBOUND(DBATHDX,3))))
    ELSE
      SCALAR_INT_BUF(46) = 0
      CALL fstarpu_block_data_register(handles(46), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(DBATHDX)
    IF (ASSOCIATED(DBATHDY)) THEN
      SCALAR_INT_BUF(47) = 1
      CALL fstarpu_block_data_register(handles(47), 0, C_LOC(DBATHDY(LBOUND(DBATHDY,1),LBOUND(DBATHDY,2),LBOUND(DBATHDY,3))), SIZE(DBATHDY,1), SIZE(DBATHDY,1)*SIZE(DBATHDY,2), SIZE(DBATHDY,1), SIZE(DBATHDY,2), SIZE(DBATHDY,3), C_SIZEOF(DBATHDY(LBOUND(DBATHDY,1),LBOUND(DBATHDY,2),LBOUND(DBATHDY,3))))
    ELSE
      SCALAR_INT_BUF(47) = 0
      CALL fstarpu_block_data_register(handles(47), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(DBATHDY)
    IF (ASSOCIATED(SFAC_ELEM)) THEN
      SCALAR_INT_BUF(48) = 1
      CALL fstarpu_block_data_register(handles(48), 0, C_LOC(SFAC_ELEM(LBOUND(SFAC_ELEM,1),LBOUND(SFAC_ELEM,2),LBOUND(SFAC_ELEM,3))), SIZE(SFAC_ELEM,1), SIZE(SFAC_ELEM,1)*SIZE(SFAC_ELEM,2), SIZE(SFAC_ELEM,1), SIZE(SFAC_ELEM,2), SIZE(SFAC_ELEM,3), C_SIZEOF(SFAC_ELEM(LBOUND(SFAC_ELEM,1),LBOUND(SFAC_ELEM,2),LBOUND(SFAC_ELEM,3))))
    ELSE
      SCALAR_INT_BUF(48) = 0
      CALL fstarpu_block_data_register(handles(48), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(SFAC_ELEM)
    IF (ASSOCIATED(BATHED)) THEN
      SCALAR_INT_BUF(49) = 1
      CALL fstarpu_tensor_data_register(handles(49), 0, C_LOC(BATHED(LBOUND(BATHED,1),LBOUND(BATHED,2),LBOUND(BATHED,3),LBOUND(BATHED,4))), SIZE(BATHED,1), SIZE(BATHED,1)*SIZE(BATHED,2), SIZE(BATHED,1)*SIZE(BATHED,2)*SIZE(BATHED,3), SIZE(BATHED,1), SIZE(BATHED,2), SIZE(BATHED,3), SIZE(BATHED,4), C_SIZEOF(BATHED(LBOUND(BATHED,1),LBOUND(BATHED,2),LBOUND(BATHED,3),LBOUND(BATHED,4))))
    ELSE
      SCALAR_INT_BUF(49) = 0
      CALL fstarpu_tensor_data_register(handles(49), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(BATHED)
    IF (ASSOCIATED(SFACED)) THEN
      SCALAR_INT_BUF(50) = 1
      CALL fstarpu_tensor_data_register(handles(50), 0, C_LOC(SFACED(LBOUND(SFACED,1),LBOUND(SFACED,2),LBOUND(SFACED,3),LBOUND(SFACED,4))), SIZE(SFACED,1), SIZE(SFACED,1)*SIZE(SFACED,2), SIZE(SFACED,1)*SIZE(SFACED,2)*SIZE(SFACED,3), SIZE(SFACED,1), SIZE(SFACED,2), SIZE(SFACED,3), SIZE(SFACED,4), C_SIZEOF(SFACED(LBOUND(SFACED,1),LBOUND(SFACED,2),LBOUND(SFACED,3),LBOUND(SFACED,4))))
    ELSE
      SCALAR_INT_BUF(50) = 0
      CALL fstarpu_tensor_data_register(handles(50), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(SFACED)
    IF (ASSOCIATED(COSNX)) THEN
      SCALAR_INT_BUF(51) = 1
      CALL fstarpu_vector_data_register(handles(51), 0, C_LOC(COSNX(LBOUND(COSNX,1))), SIZE(COSNX,1), C_SIZEOF(COSNX(LBOUND(COSNX,1))))
    ELSE
      SCALAR_INT_BUF(51) = 0
      CALL fstarpu_vector_data_register(handles(51), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(COSNX)
    IF (ASSOCIATED(SINNX)) THEN
      SCALAR_INT_BUF(52) = 1
      CALL fstarpu_vector_data_register(handles(52), 0, C_LOC(SINNX(LBOUND(SINNX,1))), SIZE(SINNX,1), C_SIZEOF(SINNX(LBOUND(SINNX,1))))
    ELSE
      SCALAR_INT_BUF(52) = 0
      CALL fstarpu_vector_data_register(handles(52), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(SINNX)
    IF (ASSOCIATED(DP_NODE)) THEN
      SCALAR_INT_BUF(53) = 1
      CALL fstarpu_block_data_register(handles(53), 0, C_LOC(DP_NODE(LBOUND(DP_NODE,1),LBOUND(DP_NODE,2),LBOUND(DP_NODE,3))), SIZE(DP_NODE,1), SIZE(DP_NODE,1)*SIZE(DP_NODE,2), SIZE(DP_NODE,1), SIZE(DP_NODE,2), SIZE(DP_NODE,3), C_SIZEOF(DP_NODE(LBOUND(DP_NODE,1),LBOUND(DP_NODE,2),LBOUND(DP_NODE,3))))
    ELSE
      SCALAR_INT_BUF(53) = 0
      CALL fstarpu_block_data_register(handles(53), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(DP_NODE)
    IF (ASSOCIATED(DP_VOL)) THEN
      SCALAR_INT_BUF(54) = 1
      CALL fstarpu_matrix_data_register(handles(54), 0, C_LOC(DP_VOL(LBOUND(DP_VOL,1),LBOUND(DP_VOL,2))), SIZE(DP_VOL,1), SIZE(DP_VOL,1), SIZE(DP_VOL,2), C_SIZEOF(DP_VOL(LBOUND(DP_VOL,1),LBOUND(DP_VOL,2))))
    ELSE
      SCALAR_INT_BUF(54) = 0
      CALL fstarpu_matrix_data_register(handles(54), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(DP_VOL)
    IF (ASSOCIATED(DRPHI)) THEN
      SCALAR_INT_BUF(55) = 1
      CALL fstarpu_block_data_register(handles(55), 0, C_LOC(DRPHI(LBOUND(DRPHI,1),LBOUND(DRPHI,2),LBOUND(DRPHI,3))), SIZE(DRPHI,1), SIZE(DRPHI,1)*SIZE(DRPHI,2), SIZE(DRPHI,1), SIZE(DRPHI,2), SIZE(DRPHI,3), C_SIZEOF(DRPHI(LBOUND(DRPHI,1),LBOUND(DRPHI,2),LBOUND(DRPHI,3))))
    ELSE
      SCALAR_INT_BUF(55) = 0
      CALL fstarpu_block_data_register(handles(55), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(DRPHI)
    IF (ASSOCIATED(DSPHI)) THEN
      SCALAR_INT_BUF(56) = 1
      CALL fstarpu_block_data_register(handles(56), 0, C_LOC(DSPHI(LBOUND(DSPHI,1),LBOUND(DSPHI,2),LBOUND(DSPHI,3))), SIZE(DSPHI,1), SIZE(DSPHI,1)*SIZE(DSPHI,2), SIZE(DSPHI,1), SIZE(DSPHI,2), SIZE(DSPHI,3), C_SIZEOF(DSPHI(LBOUND(DSPHI,1),LBOUND(DSPHI,2),LBOUND(DSPHI,3))))
    ELSE
      SCALAR_INT_BUF(56) = 0
      CALL fstarpu_block_data_register(handles(56), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(DSPHI)
    IF (ASSOCIATED(DRDX)) THEN
      SCALAR_INT_BUF(57) = 1
      CALL fstarpu_vector_data_register(handles(57), 0, C_LOC(DRDX(LBOUND(DRDX,1))), SIZE(DRDX,1), C_SIZEOF(DRDX(LBOUND(DRDX,1))))
    ELSE
      SCALAR_INT_BUF(57) = 0
      CALL fstarpu_vector_data_register(handles(57), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(DRDX)
    IF (ASSOCIATED(DSDX)) THEN
      SCALAR_INT_BUF(58) = 1
      CALL fstarpu_vector_data_register(handles(58), 0, C_LOC(DSDX(LBOUND(DSDX,1))), SIZE(DSDX,1), C_SIZEOF(DSDX(LBOUND(DSDX,1))))
    ELSE
      SCALAR_INT_BUF(58) = 0
      CALL fstarpu_vector_data_register(handles(58), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(DSDX)
    IF (ASSOCIATED(DRDY)) THEN
      SCALAR_INT_BUF(59) = 1
      CALL fstarpu_vector_data_register(handles(59), 0, C_LOC(DRDY(LBOUND(DRDY,1))), SIZE(DRDY,1), C_SIZEOF(DRDY(LBOUND(DRDY,1))))
    ELSE
      SCALAR_INT_BUF(59) = 0
      CALL fstarpu_vector_data_register(handles(59), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(DRDY)
    IF (ASSOCIATED(DSDY)) THEN
      SCALAR_INT_BUF(60) = 1
      CALL fstarpu_vector_data_register(handles(60), 0, C_LOC(DSDY(LBOUND(DSDY,1))), SIZE(DSDY,1), C_SIZEOF(DSDY(LBOUND(DSDY,1))))
    ELSE
      SCALAR_INT_BUF(60) = 0
      CALL fstarpu_vector_data_register(handles(60), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(DSDY)
    IF (ASSOCIATED(DXPHI2)) THEN
      SCALAR_INT_BUF(61) = 1
      CALL fstarpu_block_data_register(handles(61), 0, C_LOC(DXPHI2(LBOUND(DXPHI2,1),LBOUND(DXPHI2,2),LBOUND(DXPHI2,3))), SIZE(DXPHI2,1), SIZE(DXPHI2,1)*SIZE(DXPHI2,2), SIZE(DXPHI2,1), SIZE(DXPHI2,2), SIZE(DXPHI2,3), C_SIZEOF(DXPHI2(LBOUND(DXPHI2,1),LBOUND(DXPHI2,2),LBOUND(DXPHI2,3))))
    ELSE
      SCALAR_INT_BUF(61) = 0
      CALL fstarpu_block_data_register(handles(61), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(DXPHI2)
    IF (ASSOCIATED(DYPHI2)) THEN
      SCALAR_INT_BUF(62) = 1
      CALL fstarpu_block_data_register(handles(62), 0, C_LOC(DYPHI2(LBOUND(DYPHI2,1),LBOUND(DYPHI2,2),LBOUND(DYPHI2,3))), SIZE(DYPHI2,1), SIZE(DYPHI2,1)*SIZE(DYPHI2,2), SIZE(DYPHI2,1), SIZE(DYPHI2,2), SIZE(DYPHI2,3), C_SIZEOF(DYPHI2(LBOUND(DYPHI2,1),LBOUND(DYPHI2,2),LBOUND(DYPHI2,3))))
    ELSE
      SCALAR_INT_BUF(62) = 0
      CALL fstarpu_block_data_register(handles(62), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(DYPHI2)
    IF (ASSOCIATED(PHI2)) THEN
      SCALAR_INT_BUF(63) = 1
      CALL fstarpu_block_data_register(handles(63), 0, C_LOC(PHI2(LBOUND(PHI2,1),LBOUND(PHI2,2),LBOUND(PHI2,3))), SIZE(PHI2,1), SIZE(PHI2,1)*SIZE(PHI2,2), SIZE(PHI2,1), SIZE(PHI2,2), SIZE(PHI2,3), C_SIZEOF(PHI2(LBOUND(PHI2,1),LBOUND(PHI2,2),LBOUND(PHI2,3))))
    ELSE
      SCALAR_INT_BUF(63) = 0
      CALL fstarpu_block_data_register(handles(63), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(PHI2)
    IF (ASSOCIATED(EFA_DG)) THEN
      SCALAR_INT_BUF(64) = 1
      CALL fstarpu_block_data_register(handles(64), 0, C_LOC(EFA_DG(LBOUND(EFA_DG,1),LBOUND(EFA_DG,2),LBOUND(EFA_DG,3))), SIZE(EFA_DG,1), SIZE(EFA_DG,1)*SIZE(EFA_DG,2), SIZE(EFA_DG,1), SIZE(EFA_DG,2), SIZE(EFA_DG,3), C_SIZEOF(EFA_DG(LBOUND(EFA_DG,1),LBOUND(EFA_DG,2),LBOUND(EFA_DG,3))))
    ELSE
      SCALAR_INT_BUF(64) = 0
      CALL fstarpu_block_data_register(handles(64), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(EFA_DG)
    IF (ASSOCIATED(EMO_DG)) THEN
      SCALAR_INT_BUF(65) = 1
      CALL fstarpu_block_data_register(handles(65), 0, C_LOC(EMO_DG(LBOUND(EMO_DG,1),LBOUND(EMO_DG,2),LBOUND(EMO_DG,3))), SIZE(EMO_DG,1), SIZE(EMO_DG,1)*SIZE(EMO_DG,2), SIZE(EMO_DG,1), SIZE(EMO_DG,2), SIZE(EMO_DG,3), C_SIZEOF(EMO_DG(LBOUND(EMO_DG,1),LBOUND(EMO_DG,2),LBOUND(EMO_DG,3))))
    ELSE
      SCALAR_INT_BUF(65) = 0
      CALL fstarpu_block_data_register(handles(65), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(EMO_DG)
    IF (ASSOCIATED(UFA_DG)) THEN
      SCALAR_INT_BUF(66) = 1
      CALL fstarpu_block_data_register(handles(66), 0, C_LOC(UFA_DG(LBOUND(UFA_DG,1),LBOUND(UFA_DG,2),LBOUND(UFA_DG,3))), SIZE(UFA_DG,1), SIZE(UFA_DG,1)*SIZE(UFA_DG,2), SIZE(UFA_DG,1), SIZE(UFA_DG,2), SIZE(UFA_DG,3), C_SIZEOF(UFA_DG(LBOUND(UFA_DG,1),LBOUND(UFA_DG,2),LBOUND(UFA_DG,3))))
    ELSE
      SCALAR_INT_BUF(66) = 0
      CALL fstarpu_block_data_register(handles(66), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(UFA_DG)
    IF (ASSOCIATED(UMO_DG)) THEN
      SCALAR_INT_BUF(67) = 1
      CALL fstarpu_block_data_register(handles(67), 0, C_LOC(UMO_DG(LBOUND(UMO_DG,1),LBOUND(UMO_DG,2),LBOUND(UMO_DG,3))), SIZE(UMO_DG,1), SIZE(UMO_DG,1)*SIZE(UMO_DG,2), SIZE(UMO_DG,1), SIZE(UMO_DG,2), SIZE(UMO_DG,3), C_SIZEOF(UMO_DG(LBOUND(UMO_DG,1),LBOUND(UMO_DG,2),LBOUND(UMO_DG,3))))
    ELSE
      SCALAR_INT_BUF(67) = 0
      CALL fstarpu_block_data_register(handles(67), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(UMO_DG)
    IF (ASSOCIATED(VFA_DG)) THEN
      SCALAR_INT_BUF(68) = 1
      CALL fstarpu_block_data_register(handles(68), 0, C_LOC(VFA_DG(LBOUND(VFA_DG,1),LBOUND(VFA_DG,2),LBOUND(VFA_DG,3))), SIZE(VFA_DG,1), SIZE(VFA_DG,1)*SIZE(VFA_DG,2), SIZE(VFA_DG,1), SIZE(VFA_DG,2), SIZE(VFA_DG,3), C_SIZEOF(VFA_DG(LBOUND(VFA_DG,1),LBOUND(VFA_DG,2),LBOUND(VFA_DG,3))))
    ELSE
      SCALAR_INT_BUF(68) = 0
      CALL fstarpu_block_data_register(handles(68), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(VFA_DG)
    IF (ASSOCIATED(VMO_DG)) THEN
      SCALAR_INT_BUF(69) = 1
      CALL fstarpu_block_data_register(handles(69), 0, C_LOC(VMO_DG(LBOUND(VMO_DG,1),LBOUND(VMO_DG,2),LBOUND(VMO_DG,3))), SIZE(VMO_DG,1), SIZE(VMO_DG,1)*SIZE(VMO_DG,2), SIZE(VMO_DG,1), SIZE(VMO_DG,2), SIZE(VMO_DG,3), C_SIZEOF(VMO_DG(LBOUND(VMO_DG,1),LBOUND(VMO_DG,2),LBOUND(VMO_DG,3))))
    ELSE
      SCALAR_INT_BUF(69) = 0
      CALL fstarpu_block_data_register(handles(69), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(VMO_DG)
    IF (ASSOCIATED(XLEN)) THEN
      SCALAR_INT_BUF(70) = 1
      CALL fstarpu_vector_data_register(handles(70), 0, C_LOC(XLEN(LBOUND(XLEN,1))), SIZE(XLEN,1), C_SIZEOF(XLEN(LBOUND(XLEN,1))))
    ELSE
      SCALAR_INT_BUF(70) = 0
      CALL fstarpu_vector_data_register(handles(70), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(XLEN)
    IF (ASSOCIATED(HB)) THEN
      SCALAR_INT_BUF(71) = 1
      CALL fstarpu_block_data_register(handles(71), 0, C_LOC(HB(LBOUND(HB,1),LBOUND(HB,2),LBOUND(HB,3))), SIZE(HB,1), SIZE(HB,1)*SIZE(HB,2), SIZE(HB,1), SIZE(HB,2), SIZE(HB,3), C_SIZEOF(HB(LBOUND(HB,1),LBOUND(HB,2),LBOUND(HB,3))))
    ELSE
      SCALAR_INT_BUF(71) = 0
      CALL fstarpu_block_data_register(handles(71), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(HB)
    IF (ASSOCIATED(MANN)) THEN
      SCALAR_INT_BUF(72) = 1
      CALL fstarpu_matrix_data_register(handles(72), 0, C_LOC(MANN(LBOUND(MANN,1),LBOUND(MANN,2))), SIZE(MANN,1), SIZE(MANN,1), SIZE(MANN,2), C_SIZEOF(MANN(LBOUND(MANN,1),LBOUND(MANN,2))))
    ELSE
      SCALAR_INT_BUF(72) = 0
      CALL fstarpu_matrix_data_register(handles(72), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(MANN)
    IF (ASSOCIATED(IBHT)) THEN
      SCALAR_INT_BUF(73) = 1
      CALL fstarpu_vector_data_register(handles(73), 0, C_LOC(IBHT(LBOUND(IBHT,1))), SIZE(IBHT,1), C_SIZEOF(IBHT(LBOUND(IBHT,1))))
    ELSE
      SCALAR_INT_BUF(73) = 0
      CALL fstarpu_vector_data_register(handles(73), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(IBHT)
    IF (ASSOCIATED(EBHT)) THEN
      SCALAR_INT_BUF(74) = 1
      CALL fstarpu_vector_data_register(handles(74), 0, C_LOC(EBHT(LBOUND(EBHT,1))), SIZE(EBHT,1), C_SIZEOF(EBHT(LBOUND(EBHT,1))))
    ELSE
      SCALAR_INT_BUF(74) = 0
      CALL fstarpu_vector_data_register(handles(74), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(EBHT)
    IF (ASSOCIATED(EBCFSP)) THEN
      SCALAR_INT_BUF(75) = 1
      CALL fstarpu_vector_data_register(handles(75), 0, C_LOC(EBCFSP(LBOUND(EBCFSP,1))), SIZE(EBCFSP,1), C_SIZEOF(EBCFSP(LBOUND(EBCFSP,1))))
    ELSE
      SCALAR_INT_BUF(75) = 0
      CALL fstarpu_vector_data_register(handles(75), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(EBCFSP)
    IF (ASSOCIATED(IBCFSP)) THEN
      SCALAR_INT_BUF(76) = 1
      CALL fstarpu_vector_data_register(handles(76), 0, C_LOC(IBCFSP(LBOUND(IBCFSP,1))), SIZE(IBCFSP,1), C_SIZEOF(IBCFSP(LBOUND(IBCFSP,1))))
    ELSE
      SCALAR_INT_BUF(76) = 0
      CALL fstarpu_vector_data_register(handles(76), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(IBCFSP)
    IF (ASSOCIATED(IBCFSB)) THEN
      SCALAR_INT_BUF(77) = 1
      CALL fstarpu_vector_data_register(handles(77), 0, C_LOC(IBCFSB(LBOUND(IBCFSB,1))), SIZE(IBCFSB,1), C_SIZEOF(IBCFSB(LBOUND(IBCFSB,1))))
    ELSE
      SCALAR_INT_BUF(77) = 0
      CALL fstarpu_vector_data_register(handles(77), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(IBCFSB)
    IF (ASSOCIATED(JACOBI)) THEN
      SCALAR_INT_BUF(78) = 1
      CALL fstarpu_tensor_data_register(handles(78), 0, C_LOC(JACOBI(LBOUND(JACOBI,1),LBOUND(JACOBI,2),LBOUND(JACOBI,3),LBOUND(JACOBI,4))), SIZE(JACOBI,1), SIZE(JACOBI,1)*SIZE(JACOBI,2), SIZE(JACOBI,1)*SIZE(JACOBI,2)*SIZE(JACOBI,3), SIZE(JACOBI,1), SIZE(JACOBI,2), SIZE(JACOBI,3), SIZE(JACOBI,4), C_SIZEOF(JACOBI(LBOUND(JACOBI,1),LBOUND(JACOBI,2),LBOUND(JACOBI,3),LBOUND(JACOBI,4))))
    ELSE
      SCALAR_INT_BUF(78) = 0
      CALL fstarpu_tensor_data_register(handles(78), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(JACOBI)
    IF (ASSOCIATED(M_INV)) THEN
      SCALAR_INT_BUF(79) = 1
      CALL fstarpu_matrix_data_register(handles(79), 0, C_LOC(M_INV(LBOUND(M_INV,1),LBOUND(M_INV,2))), SIZE(M_INV,1), SIZE(M_INV,1), SIZE(M_INV,2), C_SIZEOF(M_INV(LBOUND(M_INV,1),LBOUND(M_INV,2))))
    ELSE
      SCALAR_INT_BUF(79) = 0
      CALL fstarpu_matrix_data_register(handles(79), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(M_INV)
    IF (ASSOCIATED(phi_edge_fixed)) THEN
      SCALAR_INT_BUF(80) = 1
      CALL fstarpu_block_data_register(handles(80), 0, C_LOC(phi_edge_fixed(LBOUND(phi_edge_fixed,1),LBOUND(phi_edge_fixed,2),LBOUND(phi_edge_fixed,3))), SIZE(phi_edge_fixed,1), SIZE(phi_edge_fixed,1)*SIZE(phi_edge_fixed,2), SIZE(phi_edge_fixed,1), SIZE(phi_edge_fixed,2), SIZE(phi_edge_fixed,3), C_SIZEOF(phi_edge_fixed(LBOUND(phi_edge_fixed,1),LBOUND(phi_edge_fixed,2),LBOUND(phi_edge_fixed,3))))
    ELSE
      SCALAR_INT_BUF(80) = 0
      CALL fstarpu_block_data_register(handles(80), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(phi_edge_fixed)
    IF (ASSOCIATED(PHI_AREA)) THEN
      SCALAR_INT_BUF(81) = 1
      CALL fstarpu_block_data_register(handles(81), 0, C_LOC(PHI_AREA(LBOUND(PHI_AREA,1),LBOUND(PHI_AREA,2),LBOUND(PHI_AREA,3))), SIZE(PHI_AREA,1), SIZE(PHI_AREA,1)*SIZE(PHI_AREA,2), SIZE(PHI_AREA,1), SIZE(PHI_AREA,2), SIZE(PHI_AREA,3), C_SIZEOF(PHI_AREA(LBOUND(PHI_AREA,1),LBOUND(PHI_AREA,2),LBOUND(PHI_AREA,3))))
    ELSE
      SCALAR_INT_BUF(81) = 0
      CALL fstarpu_block_data_register(handles(81), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(PHI_AREA)
    IF (ASSOCIATED(PHI_EDGE)) THEN
      SCALAR_INT_BUF(82) = 1
      CALL fstarpu_tensor_data_register(handles(82), 0, C_LOC(PHI_EDGE(LBOUND(PHI_EDGE,1),LBOUND(PHI_EDGE,2),LBOUND(PHI_EDGE,3),LBOUND(PHI_EDGE,4))), SIZE(PHI_EDGE,1), SIZE(PHI_EDGE,1)*SIZE(PHI_EDGE,2), SIZE(PHI_EDGE,1)*SIZE(PHI_EDGE,2)*SIZE(PHI_EDGE,3), SIZE(PHI_EDGE,1), SIZE(PHI_EDGE,2), SIZE(PHI_EDGE,3), SIZE(PHI_EDGE,4), C_SIZEOF(PHI_EDGE(LBOUND(PHI_EDGE,1),LBOUND(PHI_EDGE,2),LBOUND(PHI_EDGE,3),LBOUND(PHI_EDGE,4))))
    ELSE
      SCALAR_INT_BUF(82) = 0
      CALL fstarpu_tensor_data_register(handles(82), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(PHI_EDGE)
    IF (ASSOCIATED(PHI_CENTER)) THEN
      SCALAR_INT_BUF(83) = 1
      CALL fstarpu_matrix_data_register(handles(83), 0, C_LOC(PHI_CENTER(LBOUND(PHI_CENTER,1),LBOUND(PHI_CENTER,2))), SIZE(PHI_CENTER,1), SIZE(PHI_CENTER,1), SIZE(PHI_CENTER,2), C_SIZEOF(PHI_CENTER(LBOUND(PHI_CENTER,1),LBOUND(PHI_CENTER,2))))
    ELSE
      SCALAR_INT_BUF(83) = 0
      CALL fstarpu_matrix_data_register(handles(83), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(PHI_CENTER)
    IF (ASSOCIATED(PHI_CORNER)) THEN
      SCALAR_INT_BUF(84) = 1
      CALL fstarpu_block_data_register(handles(84), 0, C_LOC(PHI_CORNER(LBOUND(PHI_CORNER,1),LBOUND(PHI_CORNER,2),LBOUND(PHI_CORNER,3))), SIZE(PHI_CORNER,1), SIZE(PHI_CORNER,1)*SIZE(PHI_CORNER,2), SIZE(PHI_CORNER,1), SIZE(PHI_CORNER,2), SIZE(PHI_CORNER,3), C_SIZEOF(PHI_CORNER(LBOUND(PHI_CORNER,1),LBOUND(PHI_CORNER,2),LBOUND(PHI_CORNER,3))))
    ELSE
      SCALAR_INT_BUF(84) = 0
      CALL fstarpu_block_data_register(handles(84), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(PHI_CORNER)
    IF (ASSOCIATED(PHI_CHECK)) THEN
      SCALAR_INT_BUF(85) = 1
      CALL fstarpu_block_data_register(handles(85), 0, C_LOC(PHI_CHECK(LBOUND(PHI_CHECK,1),LBOUND(PHI_CHECK,2),LBOUND(PHI_CHECK,3))), SIZE(PHI_CHECK,1), SIZE(PHI_CHECK,1)*SIZE(PHI_CHECK,2), SIZE(PHI_CHECK,1), SIZE(PHI_CHECK,2), SIZE(PHI_CHECK,3), C_SIZEOF(PHI_CHECK(LBOUND(PHI_CHECK,1),LBOUND(PHI_CHECK,2),LBOUND(PHI_CHECK,3))))
    ELSE
      SCALAR_INT_BUF(85) = 0
      CALL fstarpu_block_data_register(handles(85), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(PHI_CHECK)
    IF (ASSOCIATED(PHI_CORNER1)) THEN
      SCALAR_INT_BUF(86) = 1
      CALL fstarpu_tensor_data_register(handles(86), 0, C_LOC(PHI_CORNER1(LBOUND(PHI_CORNER1,1),LBOUND(PHI_CORNER1,2),LBOUND(PHI_CORNER1,3),LBOUND(PHI_CORNER1,4))), SIZE(PHI_CORNER1,1), SIZE(PHI_CORNER1,1)*SIZE(PHI_CORNER1,2), SIZE(PHI_CORNER1,1)*SIZE(PHI_CORNER1,2)*SIZE(PHI_CORNER1,3), SIZE(PHI_CORNER1,1), SIZE(PHI_CORNER1,2), SIZE(PHI_CORNER1,3), SIZE(PHI_CORNER1,4), C_SIZEOF(PHI_CORNER1(LBOUND(PHI_CORNER1,1),LBOUND(PHI_CORNER1,2),LBOUND(PHI_CORNER1,3),LBOUND(PHI_CORNER1,4))))
    ELSE
      SCALAR_INT_BUF(86) = 0
      CALL fstarpu_tensor_data_register(handles(86), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(PHI_CORNER1)
    IF (ASSOCIATED(PHI_MID)) THEN
      SCALAR_INT_BUF(87) = 1
      CALL fstarpu_block_data_register(handles(87), 0, C_LOC(PHI_MID(LBOUND(PHI_MID,1),LBOUND(PHI_MID,2),LBOUND(PHI_MID,3))), SIZE(PHI_MID,1), SIZE(PHI_MID,1)*SIZE(PHI_MID,2), SIZE(PHI_MID,1), SIZE(PHI_MID,2), SIZE(PHI_MID,3), C_SIZEOF(PHI_MID(LBOUND(PHI_MID,1),LBOUND(PHI_MID,2),LBOUND(PHI_MID,3))))
    ELSE
      SCALAR_INT_BUF(87) = 0
      CALL fstarpu_block_data_register(handles(87), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(PHI_MID)
    IF (ASSOCIATED(PHI_INTEGRATED)) THEN
      SCALAR_INT_BUF(88) = 1
      CALL fstarpu_matrix_data_register(handles(88), 0, C_LOC(PHI_INTEGRATED(LBOUND(PHI_INTEGRATED,1),LBOUND(PHI_INTEGRATED,2))), SIZE(PHI_INTEGRATED,1), SIZE(PHI_INTEGRATED,1), SIZE(PHI_INTEGRATED,2), C_SIZEOF(PHI_INTEGRATED(LBOUND(PHI_INTEGRATED,1),LBOUND(PHI_INTEGRATED,2))))
    ELSE
      SCALAR_INT_BUF(88) = 0
      CALL fstarpu_matrix_data_register(handles(88), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(PHI_INTEGRATED)
    IF (ASSOCIATED(PSI_CHECK)) THEN
      SCALAR_INT_BUF(89) = 1
      CALL fstarpu_matrix_data_register(handles(89), 0, C_LOC(PSI_CHECK(LBOUND(PSI_CHECK,1),LBOUND(PSI_CHECK,2))), SIZE(PSI_CHECK,1), SIZE(PSI_CHECK,1), SIZE(PSI_CHECK,2), C_SIZEOF(PSI_CHECK(LBOUND(PSI_CHECK,1),LBOUND(PSI_CHECK,2))))
    ELSE
      SCALAR_INT_BUF(89) = 0
      CALL fstarpu_matrix_data_register(handles(89), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(PSI_CHECK)
    IF (ASSOCIATED(PSI1)) THEN
      SCALAR_INT_BUF(90) = 1
      CALL fstarpu_matrix_data_register(handles(90), 0, C_LOC(PSI1(LBOUND(PSI1,1),LBOUND(PSI1,2))), SIZE(PSI1,1), SIZE(PSI1,1), SIZE(PSI1,2), C_SIZEOF(PSI1(LBOUND(PSI1,1),LBOUND(PSI1,2))))
    ELSE
      SCALAR_INT_BUF(90) = 0
      CALL fstarpu_matrix_data_register(handles(90), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(PSI1)
    IF (ASSOCIATED(PSI2)) THEN
      SCALAR_INT_BUF(91) = 1
      CALL fstarpu_matrix_data_register(handles(91), 0, C_LOC(PSI2(LBOUND(PSI2,1),LBOUND(PSI2,2))), SIZE(PSI2,1), SIZE(PSI2,1), SIZE(PSI2,2), C_SIZEOF(PSI2(LBOUND(PSI2,1),LBOUND(PSI2,2))))
    ELSE
      SCALAR_INT_BUF(91) = 0
      CALL fstarpu_matrix_data_register(handles(91), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(PSI2)
    IF (ASSOCIATED(PSI3)) THEN
      SCALAR_INT_BUF(92) = 1
      CALL fstarpu_matrix_data_register(handles(92), 0, C_LOC(PSI3(LBOUND(PSI3,1),LBOUND(PSI3,2))), SIZE(PSI3,1), SIZE(PSI3,1), SIZE(PSI3,2), C_SIZEOF(PSI3(LBOUND(PSI3,1),LBOUND(PSI3,2))))
    ELSE
      SCALAR_INT_BUF(92) = 0
      CALL fstarpu_matrix_data_register(handles(92), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(PSI3)
    IF (ASSOCIATED(Q_HAT)) THEN
      SCALAR_INT_BUF(93) = 1
      CALL fstarpu_vector_data_register(handles(93), 0, C_LOC(Q_HAT(LBOUND(Q_HAT,1))), SIZE(Q_HAT,1), C_SIZEOF(Q_HAT(LBOUND(Q_HAT,1))))
    ELSE
      SCALAR_INT_BUF(93) = 0
      CALL fstarpu_vector_data_register(handles(93), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(Q_HAT)
    IF (ASSOCIATED(QIB)) THEN
      SCALAR_INT_BUF(94) = 1
      CALL fstarpu_vector_data_register(handles(94), 0, C_LOC(QIB(LBOUND(QIB,1))), SIZE(QIB,1), C_SIZEOF(QIB(LBOUND(QIB,1))))
    ELSE
      SCALAR_INT_BUF(94) = 0
      CALL fstarpu_vector_data_register(handles(94), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(QIB)
    IF (ASSOCIATED(QX)) THEN
      SCALAR_INT_BUF(95) = 1
      CALL fstarpu_block_data_register(handles(95), 0, C_LOC(QX(LBOUND(QX,1),LBOUND(QX,2),LBOUND(QX,3))), SIZE(QX,1), SIZE(QX,1)*SIZE(QX,2), SIZE(QX,1), SIZE(QX,2), SIZE(QX,3), C_SIZEOF(QX(LBOUND(QX,1),LBOUND(QX,2),LBOUND(QX,3))))
    ELSE
      SCALAR_INT_BUF(95) = 0
      CALL fstarpu_block_data_register(handles(95), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(QX)
    IF (ASSOCIATED(QY)) THEN
      SCALAR_INT_BUF(96) = 1
      CALL fstarpu_block_data_register(handles(96), 0, C_LOC(QY(LBOUND(QY,1),LBOUND(QY,2),LBOUND(QY,3))), SIZE(QY,1), SIZE(QY,1)*SIZE(QY,2), SIZE(QY,1), SIZE(QY,2), SIZE(QY,3), C_SIZEOF(QY(LBOUND(QY,1),LBOUND(QY,2),LBOUND(QY,3))))
    ELSE
      SCALAR_INT_BUF(96) = 0
      CALL fstarpu_block_data_register(handles(96), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(QY)
    IF (ASSOCIATED(ZE)) THEN
      SCALAR_INT_BUF(97) = 1
      CALL fstarpu_block_data_register(handles(97), 0, C_LOC(ZE(LBOUND(ZE,1),LBOUND(ZE,2),LBOUND(ZE,3))), SIZE(ZE,1), SIZE(ZE,1)*SIZE(ZE,2), SIZE(ZE,1), SIZE(ZE,2), SIZE(ZE,3), C_SIZEOF(ZE(LBOUND(ZE,1),LBOUND(ZE,2),LBOUND(ZE,3))))
    ELSE
      SCALAR_INT_BUF(97) = 0
      CALL fstarpu_block_data_register(handles(97), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(ZE)
    IF (ASSOCIATED(ze_edge)) THEN
      SCALAR_INT_BUF(98) = 1
      CALL fstarpu_block_data_register(handles(98), 0, C_LOC(ze_edge(LBOUND(ze_edge,1),LBOUND(ze_edge,2),LBOUND(ze_edge,3))), SIZE(ze_edge,1), SIZE(ze_edge,1)*SIZE(ze_edge,2), SIZE(ze_edge,1), SIZE(ze_edge,2), SIZE(ze_edge,3), C_SIZEOF(ze_edge(LBOUND(ze_edge,1),LBOUND(ze_edge,2),LBOUND(ze_edge,3))))
    ELSE
      SCALAR_INT_BUF(98) = 0
      CALL fstarpu_block_data_register(handles(98), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(ze_edge)
    IF (ASSOCIATED(qx_edge)) THEN
      SCALAR_INT_BUF(99) = 1
      CALL fstarpu_block_data_register(handles(99), 0, C_LOC(qx_edge(LBOUND(qx_edge,1),LBOUND(qx_edge,2),LBOUND(qx_edge,3))), SIZE(qx_edge,1), SIZE(qx_edge,1)*SIZE(qx_edge,2), SIZE(qx_edge,1), SIZE(qx_edge,2), SIZE(qx_edge,3), C_SIZEOF(qx_edge(LBOUND(qx_edge,1),LBOUND(qx_edge,2),LBOUND(qx_edge,3))))
    ELSE
      SCALAR_INT_BUF(99) = 0
      CALL fstarpu_block_data_register(handles(99), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(qx_edge)
    IF (ASSOCIATED(qy_edge)) THEN
      SCALAR_INT_BUF(100) = 1
      CALL fstarpu_block_data_register(handles(100), 0, C_LOC(qy_edge(LBOUND(qy_edge,1),LBOUND(qy_edge,2),LBOUND(qy_edge,3))), SIZE(qy_edge,1), SIZE(qy_edge,1)*SIZE(qy_edge,2), SIZE(qy_edge,1), SIZE(qy_edge,2), SIZE(qy_edge,3), C_SIZEOF(qy_edge(LBOUND(qy_edge,1),LBOUND(qy_edge,2),LBOUND(qy_edge,3))))
    ELSE
      SCALAR_INT_BUF(100) = 0
      CALL fstarpu_block_data_register(handles(100), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(qy_edge)
    IF (ASSOCIATED(elem_edge)) THEN
      SCALAR_INT_BUF(101) = 1
      CALL fstarpu_matrix_data_register(handles(101), 0, C_LOC(elem_edge(LBOUND(elem_edge,1),LBOUND(elem_edge,2))), SIZE(elem_edge,1), SIZE(elem_edge,1), SIZE(elem_edge,2), C_SIZEOF(elem_edge(LBOUND(elem_edge,1),LBOUND(elem_edge,2))))
    ELSE
      SCALAR_INT_BUF(101) = 0
      CALL fstarpu_matrix_data_register(handles(101), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(elem_edge)
    IF (ASSOCIATED(nieds_count)) THEN
      SCALAR_INT_BUF(102) = 1
      CALL fstarpu_vector_data_register(handles(102), 0, C_LOC(nieds_count(LBOUND(nieds_count,1))), SIZE(nieds_count,1), C_SIZEOF(nieds_count(LBOUND(nieds_count,1))))
    ELSE
      SCALAR_INT_BUF(102) = 0
      CALL fstarpu_vector_data_register(handles(102), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(nieds_count)
    IF (ASSOCIATED(bed)) THEN
      SCALAR_INT_BUF(103) = 1
      CALL fstarpu_tensor_data_register(handles(103), 0, C_LOC(bed(LBOUND(bed,1),LBOUND(bed,2),LBOUND(bed,3),LBOUND(bed,4))), SIZE(bed,1), SIZE(bed,1)*SIZE(bed,2), SIZE(bed,1)*SIZE(bed,2)*SIZE(bed,3), SIZE(bed,1), SIZE(bed,2), SIZE(bed,3), SIZE(bed,4), C_SIZEOF(bed(LBOUND(bed,1),LBOUND(bed,2),LBOUND(bed,3),LBOUND(bed,4))))
    ELSE
      SCALAR_INT_BUF(103) = 0
      CALL fstarpu_tensor_data_register(handles(103), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(bed)
    IF (ASSOCIATED(dynP)) THEN
      SCALAR_INT_BUF(104) = 1
      CALL fstarpu_block_data_register(handles(104), 0, C_LOC(dynP(LBOUND(dynP,1),LBOUND(dynP,2),LBOUND(dynP,3))), SIZE(dynP,1), SIZE(dynP,1)*SIZE(dynP,2), SIZE(dynP,1), SIZE(dynP,2), SIZE(dynP,3), C_SIZEOF(dynP(LBOUND(dynP,1),LBOUND(dynP,2),LBOUND(dynP,3))))
    ELSE
      SCALAR_INT_BUF(104) = 0
      CALL fstarpu_block_data_register(handles(104), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(dynP)
    IF (ASSOCIATED(dynP_MAX)) THEN
      SCALAR_INT_BUF(105) = 1
      CALL fstarpu_vector_data_register(handles(105), 0, C_LOC(dynP_MAX(LBOUND(dynP_MAX,1))), SIZE(dynP_MAX,1), C_SIZEOF(dynP_MAX(LBOUND(dynP_MAX,1))))
    ELSE
      SCALAR_INT_BUF(105) = 0
      CALL fstarpu_vector_data_register(handles(105), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(dynP_MAX)
    IF (ASSOCIATED(dynP_MIN)) THEN
      SCALAR_INT_BUF(106) = 1
      CALL fstarpu_vector_data_register(handles(106), 0, C_LOC(dynP_MIN(LBOUND(dynP_MIN,1))), SIZE(dynP_MIN,1), C_SIZEOF(dynP_MIN(LBOUND(dynP_MIN,1))))
    ELSE
      SCALAR_INT_BUF(106) = 0
      CALL fstarpu_vector_data_register(handles(106), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(dynP_MIN)
    IF (ASSOCIATED(iota)) THEN
      SCALAR_INT_BUF(107) = 1
      CALL fstarpu_block_data_register(handles(107), 0, C_LOC(iota(LBOUND(iota,1),LBOUND(iota,2),LBOUND(iota,3))), SIZE(iota,1), SIZE(iota,1)*SIZE(iota,2), SIZE(iota,1), SIZE(iota,2), SIZE(iota,3), C_SIZEOF(iota(LBOUND(iota,1),LBOUND(iota,2),LBOUND(iota,3))))
    ELSE
      SCALAR_INT_BUF(107) = 0
      CALL fstarpu_block_data_register(handles(107), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(iota)
    IF (ASSOCIATED(iotaa)) THEN
      SCALAR_INT_BUF(108) = 1
      CALL fstarpu_block_data_register(handles(108), 0, C_LOC(iotaa(LBOUND(iotaa,1),LBOUND(iotaa,2),LBOUND(iotaa,3))), SIZE(iotaa,1), SIZE(iotaa,1)*SIZE(iotaa,2), SIZE(iotaa,1), SIZE(iotaa,2), SIZE(iotaa,3), C_SIZEOF(iotaa(LBOUND(iotaa,1),LBOUND(iotaa,2),LBOUND(iotaa,3))))
    ELSE
      SCALAR_INT_BUF(108) = 0
      CALL fstarpu_block_data_register(handles(108), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(iotaa)
    IF (ASSOCIATED(iota2)) THEN
      SCALAR_INT_BUF(109) = 1
      CALL fstarpu_block_data_register(handles(109), 0, C_LOC(iota2(LBOUND(iota2,1),LBOUND(iota2,2),LBOUND(iota2,3))), SIZE(iota2,1), SIZE(iota2,1)*SIZE(iota2,2), SIZE(iota2,1), SIZE(iota2,2), SIZE(iota2,3), C_SIZEOF(iota2(LBOUND(iota2,1),LBOUND(iota2,2),LBOUND(iota2,3))))
    ELSE
      SCALAR_INT_BUF(109) = 0
      CALL fstarpu_block_data_register(handles(109), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(iota2)
    IF (ASSOCIATED(iota_MAX)) THEN
      SCALAR_INT_BUF(110) = 1
      CALL fstarpu_vector_data_register(handles(110), 0, C_LOC(iota_MAX(LBOUND(iota_MAX,1))), SIZE(iota_MAX,1), C_SIZEOF(iota_MAX(LBOUND(iota_MAX,1))))
    ELSE
      SCALAR_INT_BUF(110) = 0
      CALL fstarpu_vector_data_register(handles(110), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(iota_MAX)
    IF (ASSOCIATED(iota_MIN)) THEN
      SCALAR_INT_BUF(111) = 1
      CALL fstarpu_vector_data_register(handles(111), 0, C_LOC(iota_MIN(LBOUND(iota_MIN,1))), SIZE(iota_MIN,1), C_SIZEOF(iota_MIN(LBOUND(iota_MIN,1))))
    ELSE
      SCALAR_INT_BUF(111) = 0
      CALL fstarpu_vector_data_register(handles(111), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(iota_MIN)
    IF (ASSOCIATED(iotaa2)) THEN
      SCALAR_INT_BUF(112) = 1
      CALL fstarpu_block_data_register(handles(112), 0, C_LOC(iotaa2(LBOUND(iotaa2,1),LBOUND(iotaa2,2),LBOUND(iotaa2,3))), SIZE(iotaa2,1), SIZE(iotaa2,1)*SIZE(iotaa2,2), SIZE(iotaa2,1), SIZE(iotaa2,2), SIZE(iotaa2,3), C_SIZEOF(iotaa2(LBOUND(iotaa2,1),LBOUND(iotaa2,2),LBOUND(iotaa2,3))))
    ELSE
      SCALAR_INT_BUF(112) = 0
      CALL fstarpu_block_data_register(handles(112), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(iotaa2)
    IF (ASSOCIATED(iotaa3)) THEN
      SCALAR_INT_BUF(113) = 1
      CALL fstarpu_block_data_register(handles(113), 0, C_LOC(iotaa3(LBOUND(iotaa3,1),LBOUND(iotaa3,2),LBOUND(iotaa3,3))), SIZE(iotaa3,1), SIZE(iotaa3,1)*SIZE(iotaa3,2), SIZE(iotaa3,1), SIZE(iotaa3,2), SIZE(iotaa3,3), C_SIZEOF(iotaa3(LBOUND(iotaa3,1),LBOUND(iotaa3,2),LBOUND(iotaa3,3))))
    ELSE
      SCALAR_INT_BUF(113) = 0
      CALL fstarpu_block_data_register(handles(113), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(iotaa3)
    IF (ASSOCIATED(iota2_MAX)) THEN
      SCALAR_INT_BUF(114) = 1
      CALL fstarpu_vector_data_register(handles(114), 0, C_LOC(iota2_MAX(LBOUND(iota2_MAX,1))), SIZE(iota2_MAX,1), C_SIZEOF(iota2_MAX(LBOUND(iota2_MAX,1))))
    ELSE
      SCALAR_INT_BUF(114) = 0
      CALL fstarpu_vector_data_register(handles(114), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(iota2_MAX)
    IF (ASSOCIATED(iota2_MIN)) THEN
      SCALAR_INT_BUF(115) = 1
      CALL fstarpu_vector_data_register(handles(115), 0, C_LOC(iota2_MIN(LBOUND(iota2_MIN,1))), SIZE(iota2_MIN,1), C_SIZEOF(iota2_MIN(LBOUND(iota2_MIN,1))))
    ELSE
      SCALAR_INT_BUF(115) = 0
      CALL fstarpu_vector_data_register(handles(115), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(iota2_MIN)
    IF (ASSOCIATED(arrayfix)) THEN
      SCALAR_INT_BUF(116) = 1
      CALL fstarpu_block_data_register(handles(116), 0, C_LOC(arrayfix(LBOUND(arrayfix,1),LBOUND(arrayfix,2),LBOUND(arrayfix,3))), SIZE(arrayfix,1), SIZE(arrayfix,1)*SIZE(arrayfix,2), SIZE(arrayfix,1), SIZE(arrayfix,2), SIZE(arrayfix,3), C_SIZEOF(arrayfix(LBOUND(arrayfix,1),LBOUND(arrayfix,2),LBOUND(arrayfix,3))))
    ELSE
      SCALAR_INT_BUF(116) = 0
      CALL fstarpu_block_data_register(handles(116), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(arrayfix)
    IF (ASSOCIATED(CORI_EL)) THEN
      SCALAR_INT_BUF(117) = 1
      CALL fstarpu_vector_data_register(handles(117), 0, C_LOC(CORI_EL(LBOUND(CORI_EL,1))), SIZE(CORI_EL,1), C_SIZEOF(CORI_EL(LBOUND(CORI_EL,1))))
    ELSE
      SCALAR_INT_BUF(117) = 0
      CALL fstarpu_vector_data_register(handles(117), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(CORI_EL)
    IF (ASSOCIATED(FRIC_EL)) THEN
      SCALAR_INT_BUF(118) = 1
      CALL fstarpu_vector_data_register(handles(118), 0, C_LOC(FRIC_EL(LBOUND(FRIC_EL,1))), SIZE(FRIC_EL,1), C_SIZEOF(FRIC_EL(LBOUND(FRIC_EL,1))))
    ELSE
      SCALAR_INT_BUF(118) = 0
      CALL fstarpu_vector_data_register(handles(118), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(FRIC_EL)
    IF (ASSOCIATED(ZE_MAX)) THEN
      SCALAR_INT_BUF(119) = 1
      CALL fstarpu_vector_data_register(handles(119), 0, C_LOC(ZE_MAX(LBOUND(ZE_MAX,1))), SIZE(ZE_MAX,1), C_SIZEOF(ZE_MAX(LBOUND(ZE_MAX,1))))
    ELSE
      SCALAR_INT_BUF(119) = 0
      CALL fstarpu_vector_data_register(handles(119), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(ZE_MAX)
    IF (ASSOCIATED(ZE_MIN)) THEN
      SCALAR_INT_BUF(120) = 1
      CALL fstarpu_vector_data_register(handles(120), 0, C_LOC(ZE_MIN(LBOUND(ZE_MIN,1))), SIZE(ZE_MIN,1), C_SIZEOF(ZE_MIN(LBOUND(ZE_MIN,1))))
    ELSE
      SCALAR_INT_BUF(120) = 0
      CALL fstarpu_vector_data_register(handles(120), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(ZE_MIN)
    IF (ASSOCIATED(DPE_MIN)) THEN
      SCALAR_INT_BUF(121) = 1
      CALL fstarpu_vector_data_register(handles(121), 0, C_LOC(DPE_MIN(LBOUND(DPE_MIN,1))), SIZE(DPE_MIN,1), C_SIZEOF(DPE_MIN(LBOUND(DPE_MIN,1))))
    ELSE
      SCALAR_INT_BUF(121) = 0
      CALL fstarpu_vector_data_register(handles(121), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(DPE_MIN)
    IF (ASSOCIATED(WATER_DEPTH_OLD)) THEN
      SCALAR_INT_BUF(122) = 1
      CALL fstarpu_matrix_data_register(handles(122), 0, C_LOC(WATER_DEPTH_OLD(LBOUND(WATER_DEPTH_OLD,1),LBOUND(WATER_DEPTH_OLD,2))), SIZE(WATER_DEPTH_OLD,1), SIZE(WATER_DEPTH_OLD,1), SIZE(WATER_DEPTH_OLD,2), C_SIZEOF(WATER_DEPTH_OLD(LBOUND(WATER_DEPTH_OLD,1),LBOUND(WATER_DEPTH_OLD,2))))
    ELSE
      SCALAR_INT_BUF(122) = 0
      CALL fstarpu_matrix_data_register(handles(122), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(WATER_DEPTH_OLD)
    IF (ASSOCIATED(WATER_DEPTH)) THEN
      SCALAR_INT_BUF(123) = 1
      CALL fstarpu_matrix_data_register(handles(123), 0, C_LOC(WATER_DEPTH(LBOUND(WATER_DEPTH,1),LBOUND(WATER_DEPTH,2))), SIZE(WATER_DEPTH,1), SIZE(WATER_DEPTH,1), SIZE(WATER_DEPTH,2), C_SIZEOF(WATER_DEPTH(LBOUND(WATER_DEPTH,1),LBOUND(WATER_DEPTH,2))))
    ELSE
      SCALAR_INT_BUF(123) = 0
      CALL fstarpu_matrix_data_register(handles(123), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(WATER_DEPTH)
    IF (ASSOCIATED(ADVECTQX)) THEN
      SCALAR_INT_BUF(124) = 1
      CALL fstarpu_vector_data_register(handles(124), 0, C_LOC(ADVECTQX(LBOUND(ADVECTQX,1))), SIZE(ADVECTQX,1), C_SIZEOF(ADVECTQX(LBOUND(ADVECTQX,1))))
    ELSE
      SCALAR_INT_BUF(124) = 0
      CALL fstarpu_vector_data_register(handles(124), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(ADVECTQX)
    IF (ASSOCIATED(ADVECTQY)) THEN
      SCALAR_INT_BUF(125) = 1
      CALL fstarpu_vector_data_register(handles(125), 0, C_LOC(ADVECTQY(LBOUND(ADVECTQY,1))), SIZE(ADVECTQY,1), C_SIZEOF(ADVECTQY(LBOUND(ADVECTQY,1))))
    ELSE
      SCALAR_INT_BUF(125) = 0
      CALL fstarpu_vector_data_register(handles(125), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(ADVECTQY)
    IF (ASSOCIATED(SOURCEQX)) THEN
      SCALAR_INT_BUF(126) = 1
      CALL fstarpu_vector_data_register(handles(126), 0, C_LOC(SOURCEQX(LBOUND(SOURCEQX,1))), SIZE(SOURCEQX,1), C_SIZEOF(SOURCEQX(LBOUND(SOURCEQX,1))))
    ELSE
      SCALAR_INT_BUF(126) = 0
      CALL fstarpu_vector_data_register(handles(126), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(SOURCEQX)
    IF (ASSOCIATED(SOURCEQY)) THEN
      SCALAR_INT_BUF(127) = 1
      CALL fstarpu_vector_data_register(handles(127), 0, C_LOC(SOURCEQY(LBOUND(SOURCEQY,1))), SIZE(SOURCEQY,1), C_SIZEOF(SOURCEQY(LBOUND(SOURCEQY,1))))
    ELSE
      SCALAR_INT_BUF(127) = 0
      CALL fstarpu_vector_data_register(handles(127), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(SOURCEQY)
    IF (ASSOCIATED(LZ)) THEN
      SCALAR_INT_BUF(128) = 1
      CALL fstarpu_tensor_data_register(handles(128), 0, C_LOC(LZ(LBOUND(LZ,1),LBOUND(LZ,2),LBOUND(LZ,3),LBOUND(LZ,4))), SIZE(LZ,1), SIZE(LZ,1)*SIZE(LZ,2), SIZE(LZ,1)*SIZE(LZ,2)*SIZE(LZ,3), SIZE(LZ,1), SIZE(LZ,2), SIZE(LZ,3), SIZE(LZ,4), C_SIZEOF(LZ(LBOUND(LZ,1),LBOUND(LZ,2),LBOUND(LZ,3),LBOUND(LZ,4))))
    ELSE
      SCALAR_INT_BUF(128) = 0
      CALL fstarpu_tensor_data_register(handles(128), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(LZ)
    IF (ASSOCIATED(MZ)) THEN
      SCALAR_INT_BUF(129) = 1
      CALL fstarpu_tensor_data_register(handles(129), 0, C_LOC(MZ(LBOUND(MZ,1),LBOUND(MZ,2),LBOUND(MZ,3),LBOUND(MZ,4))), SIZE(MZ,1), SIZE(MZ,1)*SIZE(MZ,2), SIZE(MZ,1)*SIZE(MZ,2)*SIZE(MZ,3), SIZE(MZ,1), SIZE(MZ,2), SIZE(MZ,3), SIZE(MZ,4), C_SIZEOF(MZ(LBOUND(MZ,1),LBOUND(MZ,2),LBOUND(MZ,3),LBOUND(MZ,4))))
    ELSE
      SCALAR_INT_BUF(129) = 0
      CALL fstarpu_tensor_data_register(handles(129), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(MZ)
    IF (ASSOCIATED(HZ)) THEN
      SCALAR_INT_BUF(130) = 1
      CALL fstarpu_tensor_data_register(handles(130), 0, C_LOC(HZ(LBOUND(HZ,1),LBOUND(HZ,2),LBOUND(HZ,3),LBOUND(HZ,4))), SIZE(HZ,1), SIZE(HZ,1)*SIZE(HZ,2), SIZE(HZ,1)*SIZE(HZ,2)*SIZE(HZ,3), SIZE(HZ,1), SIZE(HZ,2), SIZE(HZ,3), SIZE(HZ,4), C_SIZEOF(HZ(LBOUND(HZ,1),LBOUND(HZ,2),LBOUND(HZ,3),LBOUND(HZ,4))))
    ELSE
      SCALAR_INT_BUF(130) = 0
      CALL fstarpu_tensor_data_register(handles(130), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(HZ)
    IF (ASSOCIATED(TZ)) THEN
      SCALAR_INT_BUF(131) = 1
      CALL fstarpu_tensor_data_register(handles(131), 0, C_LOC(TZ(LBOUND(TZ,1),LBOUND(TZ,2),LBOUND(TZ,3),LBOUND(TZ,4))), SIZE(TZ,1), SIZE(TZ,1)*SIZE(TZ,2), SIZE(TZ,1)*SIZE(TZ,2)*SIZE(TZ,3), SIZE(TZ,1), SIZE(TZ,2), SIZE(TZ,3), SIZE(TZ,4), C_SIZEOF(TZ(LBOUND(TZ,1),LBOUND(TZ,2),LBOUND(TZ,3),LBOUND(TZ,4))))
    ELSE
      SCALAR_INT_BUF(131) = 0
      CALL fstarpu_tensor_data_register(handles(131), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(TZ)
    IF (ASSOCIATED(QNAM_DG)) THEN
      SCALAR_INT_BUF(132) = 1
      CALL fstarpu_block_data_register(handles(132), 0, C_LOC(QNAM_DG(LBOUND(QNAM_DG,1),LBOUND(QNAM_DG,2),LBOUND(QNAM_DG,3))), SIZE(QNAM_DG,1), SIZE(QNAM_DG,1)*SIZE(QNAM_DG,2), SIZE(QNAM_DG,1), SIZE(QNAM_DG,2), SIZE(QNAM_DG,3), C_SIZEOF(QNAM_DG(LBOUND(QNAM_DG,1),LBOUND(QNAM_DG,2),LBOUND(QNAM_DG,3))))
    ELSE
      SCALAR_INT_BUF(132) = 0
      CALL fstarpu_block_data_register(handles(132), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(QNAM_DG)
    IF (ASSOCIATED(QNPH_DG)) THEN
      SCALAR_INT_BUF(133) = 1
      CALL fstarpu_block_data_register(handles(133), 0, C_LOC(QNPH_DG(LBOUND(QNPH_DG,1),LBOUND(QNPH_DG,2),LBOUND(QNPH_DG,3))), SIZE(QNPH_DG,1), SIZE(QNPH_DG,1)*SIZE(QNPH_DG,2), SIZE(QNPH_DG,1), SIZE(QNPH_DG,2), SIZE(QNPH_DG,3), C_SIZEOF(QNPH_DG(LBOUND(QNPH_DG,1),LBOUND(QNPH_DG,2),LBOUND(QNPH_DG,3))))
    ELSE
      SCALAR_INT_BUF(133) = 0
      CALL fstarpu_block_data_register(handles(133), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(QNPH_DG)
    IF (ASSOCIATED(RHS_ZE)) THEN
      SCALAR_INT_BUF(134) = 1
      CALL fstarpu_block_data_register(handles(134), 0, C_LOC(RHS_ZE(LBOUND(RHS_ZE,1),LBOUND(RHS_ZE,2),LBOUND(RHS_ZE,3))), SIZE(RHS_ZE,1), SIZE(RHS_ZE,1)*SIZE(RHS_ZE,2), SIZE(RHS_ZE,1), SIZE(RHS_ZE,2), SIZE(RHS_ZE,3), C_SIZEOF(RHS_ZE(LBOUND(RHS_ZE,1),LBOUND(RHS_ZE,2),LBOUND(RHS_ZE,3))))
    ELSE
      SCALAR_INT_BUF(134) = 0
      CALL fstarpu_block_data_register(handles(134), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(RHS_ZE)
    IF (ASSOCIATED(RHS_bed)) THEN
      SCALAR_INT_BUF(135) = 1
      CALL fstarpu_tensor_data_register(handles(135), 0, C_LOC(RHS_bed(LBOUND(RHS_bed,1),LBOUND(RHS_bed,2),LBOUND(RHS_bed,3),LBOUND(RHS_bed,4))), SIZE(RHS_bed,1), SIZE(RHS_bed,1)*SIZE(RHS_bed,2), SIZE(RHS_bed,1)*SIZE(RHS_bed,2)*SIZE(RHS_bed,3), SIZE(RHS_bed,1), SIZE(RHS_bed,2), SIZE(RHS_bed,3), SIZE(RHS_bed,4), C_SIZEOF(RHS_bed(LBOUND(RHS_bed,1),LBOUND(RHS_bed,2),LBOUND(RHS_bed,3),LBOUND(RHS_bed,4))))
    ELSE
      SCALAR_INT_BUF(135) = 0
      CALL fstarpu_tensor_data_register(handles(135), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(RHS_bed)
    IF (ASSOCIATED(RHS_QX)) THEN
      SCALAR_INT_BUF(136) = 1
      CALL fstarpu_block_data_register(handles(136), 0, C_LOC(RHS_QX(LBOUND(RHS_QX,1),LBOUND(RHS_QX,2),LBOUND(RHS_QX,3))), SIZE(RHS_QX,1), SIZE(RHS_QX,1)*SIZE(RHS_QX,2), SIZE(RHS_QX,1), SIZE(RHS_QX,2), SIZE(RHS_QX,3), C_SIZEOF(RHS_QX(LBOUND(RHS_QX,1),LBOUND(RHS_QX,2),LBOUND(RHS_QX,3))))
    ELSE
      SCALAR_INT_BUF(136) = 0
      CALL fstarpu_block_data_register(handles(136), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(RHS_QX)
    IF (ASSOCIATED(RHS_QY)) THEN
      SCALAR_INT_BUF(137) = 1
      CALL fstarpu_block_data_register(handles(137), 0, C_LOC(RHS_QY(LBOUND(RHS_QY,1),LBOUND(RHS_QY,2),LBOUND(RHS_QY,3))), SIZE(RHS_QY,1), SIZE(RHS_QY,1)*SIZE(RHS_QY,2), SIZE(RHS_QY,1), SIZE(RHS_QY,2), SIZE(RHS_QY,3), C_SIZEOF(RHS_QY(LBOUND(RHS_QY,1),LBOUND(RHS_QY,2),LBOUND(RHS_QY,3))))
    ELSE
      SCALAR_INT_BUF(137) = 0
      CALL fstarpu_block_data_register(handles(137), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(RHS_QY)
    IF (ASSOCIATED(RHS_iota)) THEN
      SCALAR_INT_BUF(138) = 1
      CALL fstarpu_block_data_register(handles(138), 0, C_LOC(RHS_iota(LBOUND(RHS_iota,1),LBOUND(RHS_iota,2),LBOUND(RHS_iota,3))), SIZE(RHS_iota,1), SIZE(RHS_iota,1)*SIZE(RHS_iota,2), SIZE(RHS_iota,1), SIZE(RHS_iota,2), SIZE(RHS_iota,3), C_SIZEOF(RHS_iota(LBOUND(RHS_iota,1),LBOUND(RHS_iota,2),LBOUND(RHS_iota,3))))
    ELSE
      SCALAR_INT_BUF(138) = 0
      CALL fstarpu_block_data_register(handles(138), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(RHS_iota)
    IF (ASSOCIATED(RHS_iota2)) THEN
      SCALAR_INT_BUF(139) = 1
      CALL fstarpu_block_data_register(handles(139), 0, C_LOC(RHS_iota2(LBOUND(RHS_iota2,1),LBOUND(RHS_iota2,2),LBOUND(RHS_iota2,3))), SIZE(RHS_iota2,1), SIZE(RHS_iota2,1)*SIZE(RHS_iota2,2), SIZE(RHS_iota2,1), SIZE(RHS_iota2,2), SIZE(RHS_iota2,3), C_SIZEOF(RHS_iota2(LBOUND(RHS_iota2,1),LBOUND(RHS_iota2,2),LBOUND(RHS_iota2,3))))
    ELSE
      SCALAR_INT_BUF(139) = 0
      CALL fstarpu_block_data_register(handles(139), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(RHS_iota2)
    IF (ASSOCIATED(RHS_dynP)) THEN
      SCALAR_INT_BUF(140) = 1
      CALL fstarpu_block_data_register(handles(140), 0, C_LOC(RHS_dynP(LBOUND(RHS_dynP,1),LBOUND(RHS_dynP,2),LBOUND(RHS_dynP,3))), SIZE(RHS_dynP,1), SIZE(RHS_dynP,1)*SIZE(RHS_dynP,2), SIZE(RHS_dynP,1), SIZE(RHS_dynP,2), SIZE(RHS_dynP,3), C_SIZEOF(RHS_dynP(LBOUND(RHS_dynP,1),LBOUND(RHS_dynP,2),LBOUND(RHS_dynP,3))))
    ELSE
      SCALAR_INT_BUF(140) = 0
      CALL fstarpu_block_data_register(handles(140), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(RHS_dynP)
    IF (ASSOCIATED(RHS_bed_IN)) THEN
      SCALAR_INT_BUF(141) = 1
      CALL fstarpu_matrix_data_register(handles(141), 0, C_LOC(RHS_bed_IN(LBOUND(RHS_bed_IN,1),LBOUND(RHS_bed_IN,2))), SIZE(RHS_bed_IN,1), SIZE(RHS_bed_IN,1), SIZE(RHS_bed_IN,2), C_SIZEOF(RHS_bed_IN(LBOUND(RHS_bed_IN,1),LBOUND(RHS_bed_IN,2))))
    ELSE
      SCALAR_INT_BUF(141) = 0
      CALL fstarpu_matrix_data_register(handles(141), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(RHS_bed_IN)
    IF (ASSOCIATED(RHS_bed_EX)) THEN
      SCALAR_INT_BUF(142) = 1
      CALL fstarpu_matrix_data_register(handles(142), 0, C_LOC(RHS_bed_EX(LBOUND(RHS_bed_EX,1),LBOUND(RHS_bed_EX,2))), SIZE(RHS_bed_EX,1), SIZE(RHS_bed_EX,1), SIZE(RHS_bed_EX,2), C_SIZEOF(RHS_bed_EX(LBOUND(RHS_bed_EX,1),LBOUND(RHS_bed_EX,2))))
    ELSE
      SCALAR_INT_BUF(142) = 0
      CALL fstarpu_matrix_data_register(handles(142), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(RHS_bed_EX)
    IF (ASSOCIATED(bed_HAT_O)) THEN
      SCALAR_INT_BUF(143) = 1
      CALL fstarpu_vector_data_register(handles(143), 0, C_LOC(bed_HAT_O(LBOUND(bed_HAT_O,1))), SIZE(bed_HAT_O,1), C_SIZEOF(bed_HAT_O(LBOUND(bed_HAT_O,1))))
    ELSE
      SCALAR_INT_BUF(143) = 0
      CALL fstarpu_vector_data_register(handles(143), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(bed_HAT_O)
    IF (ASSOCIATED(XAGP)) THEN
      SCALAR_INT_BUF(144) = 1
      CALL fstarpu_matrix_data_register(handles(144), 0, C_LOC(XAGP(LBOUND(XAGP,1),LBOUND(XAGP,2))), SIZE(XAGP,1), SIZE(XAGP,1), SIZE(XAGP,2), C_SIZEOF(XAGP(LBOUND(XAGP,1),LBOUND(XAGP,2))))
    ELSE
      SCALAR_INT_BUF(144) = 0
      CALL fstarpu_matrix_data_register(handles(144), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(XAGP)
    IF (ASSOCIATED(YAGP)) THEN
      SCALAR_INT_BUF(145) = 1
      CALL fstarpu_matrix_data_register(handles(145), 0, C_LOC(YAGP(LBOUND(YAGP,1),LBOUND(YAGP,2))), SIZE(YAGP,1), SIZE(YAGP,1), SIZE(YAGP,2), C_SIZEOF(YAGP(LBOUND(YAGP,1),LBOUND(YAGP,2))))
    ELSE
      SCALAR_INT_BUF(145) = 0
      CALL fstarpu_matrix_data_register(handles(145), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(YAGP)
    IF (ASSOCIATED(WAGP)) THEN
      SCALAR_INT_BUF(146) = 1
      CALL fstarpu_matrix_data_register(handles(146), 0, C_LOC(WAGP(LBOUND(WAGP,1),LBOUND(WAGP,2))), SIZE(WAGP,1), SIZE(WAGP,1), SIZE(WAGP,2), C_SIZEOF(WAGP(LBOUND(WAGP,1),LBOUND(WAGP,2))))
    ELSE
      SCALAR_INT_BUF(146) = 0
      CALL fstarpu_matrix_data_register(handles(146), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(WAGP)
    IF (ASSOCIATED(XEGP)) THEN
      SCALAR_INT_BUF(147) = 1
      CALL fstarpu_matrix_data_register(handles(147), 0, C_LOC(XEGP(LBOUND(XEGP,1),LBOUND(XEGP,2))), SIZE(XEGP,1), SIZE(XEGP,1), SIZE(XEGP,2), C_SIZEOF(XEGP(LBOUND(XEGP,1),LBOUND(XEGP,2))))
    ELSE
      SCALAR_INT_BUF(147) = 0
      CALL fstarpu_matrix_data_register(handles(147), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(XEGP)
    IF (ASSOCIATED(YEGP)) THEN
      SCALAR_INT_BUF(148) = 1
      CALL fstarpu_matrix_data_register(handles(148), 0, C_LOC(YEGP(LBOUND(YEGP,1),LBOUND(YEGP,2))), SIZE(YEGP,1), SIZE(YEGP,1), SIZE(YEGP,2), C_SIZEOF(YEGP(LBOUND(YEGP,1),LBOUND(YEGP,2))))
    ELSE
      SCALAR_INT_BUF(148) = 0
      CALL fstarpu_matrix_data_register(handles(148), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(YEGP)
    IF (ASSOCIATED(WEGP)) THEN
      SCALAR_INT_BUF(149) = 1
      CALL fstarpu_matrix_data_register(handles(149), 0, C_LOC(WEGP(LBOUND(WEGP,1),LBOUND(WEGP,2))), SIZE(WEGP,1), SIZE(WEGP,1), SIZE(WEGP,2), C_SIZEOF(WEGP(LBOUND(WEGP,1),LBOUND(WEGP,2))))
    ELSE
      SCALAR_INT_BUF(149) = 0
      CALL fstarpu_matrix_data_register(handles(149), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(WEGP)
    IF (ASSOCIATED(SL3)) THEN
      SCALAR_INT_BUF(150) = 1
      CALL fstarpu_matrix_data_register(handles(150), 0, C_LOC(SL3(LBOUND(SL3,1),LBOUND(SL3,2))), SIZE(SL3,1), SIZE(SL3,1), SIZE(SL3,2), C_SIZEOF(SL3(LBOUND(SL3,1),LBOUND(SL3,2))))
    ELSE
      SCALAR_INT_BUF(150) = 0
      CALL fstarpu_matrix_data_register(handles(150), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(SL3)
    IF (ASSOCIATED(XBC)) THEN
      SCALAR_INT_BUF(151) = 1
      CALL fstarpu_vector_data_register(handles(151), 0, C_LOC(XBC(LBOUND(XBC,1))), SIZE(XBC,1), C_SIZEOF(XBC(LBOUND(XBC,1))))
    ELSE
      SCALAR_INT_BUF(151) = 0
      CALL fstarpu_vector_data_register(handles(151), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(XBC)
    IF (ASSOCIATED(YBC)) THEN
      SCALAR_INT_BUF(152) = 1
      CALL fstarpu_vector_data_register(handles(152), 0, C_LOC(YBC(LBOUND(YBC,1))), SIZE(YBC,1), C_SIZEOF(YBC(LBOUND(YBC,1))))
    ELSE
      SCALAR_INT_BUF(152) = 0
      CALL fstarpu_vector_data_register(handles(152), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(YBC)
    IF (ASSOCIATED(XFAC)) THEN
      SCALAR_INT_BUF(153) = 1
      CALL fstarpu_tensor_data_register(handles(153), 0, C_LOC(XFAC(LBOUND(XFAC,1),LBOUND(XFAC,2),LBOUND(XFAC,3),LBOUND(XFAC,4))), SIZE(XFAC,1), SIZE(XFAC,1)*SIZE(XFAC,2), SIZE(XFAC,1)*SIZE(XFAC,2)*SIZE(XFAC,3), SIZE(XFAC,1), SIZE(XFAC,2), SIZE(XFAC,3), SIZE(XFAC,4), C_SIZEOF(XFAC(LBOUND(XFAC,1),LBOUND(XFAC,2),LBOUND(XFAC,3),LBOUND(XFAC,4))))
    ELSE
      SCALAR_INT_BUF(153) = 0
      CALL fstarpu_tensor_data_register(handles(153), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(XFAC)
    IF (ASSOCIATED(YFAC)) THEN
      SCALAR_INT_BUF(154) = 1
      CALL fstarpu_tensor_data_register(handles(154), 0, C_LOC(YFAC(LBOUND(YFAC,1),LBOUND(YFAC,2),LBOUND(YFAC,3),LBOUND(YFAC,4))), SIZE(YFAC,1), SIZE(YFAC,1)*SIZE(YFAC,2), SIZE(YFAC,1)*SIZE(YFAC,2)*SIZE(YFAC,3), SIZE(YFAC,1), SIZE(YFAC,2), SIZE(YFAC,3), SIZE(YFAC,4), C_SIZEOF(YFAC(LBOUND(YFAC,1),LBOUND(YFAC,2),LBOUND(YFAC,3),LBOUND(YFAC,4))))
    ELSE
      SCALAR_INT_BUF(154) = 0
      CALL fstarpu_tensor_data_register(handles(154), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(YFAC)
    IF (ASSOCIATED(SRFAC)) THEN
      SCALAR_INT_BUF(155) = 1
      CALL fstarpu_tensor_data_register(handles(155), 0, C_LOC(SRFAC(LBOUND(SRFAC,1),LBOUND(SRFAC,2),LBOUND(SRFAC,3),LBOUND(SRFAC,4))), SIZE(SRFAC,1), SIZE(SRFAC,1)*SIZE(SRFAC,2), SIZE(SRFAC,1)*SIZE(SRFAC,2)*SIZE(SRFAC,3), SIZE(SRFAC,1), SIZE(SRFAC,2), SIZE(SRFAC,3), SIZE(SRFAC,4), C_SIZEOF(SRFAC(LBOUND(SRFAC,1),LBOUND(SRFAC,2),LBOUND(SRFAC,3),LBOUND(SRFAC,4))))
    ELSE
      SCALAR_INT_BUF(155) = 0
      CALL fstarpu_tensor_data_register(handles(155), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(SRFAC)
    IF (ASSOCIATED(EDGEQ)) THEN
      SCALAR_INT_BUF(156) = 1
      CALL fstarpu_tensor_data_register(handles(156), 0, C_LOC(EDGEQ(LBOUND(EDGEQ,1),LBOUND(EDGEQ,2),LBOUND(EDGEQ,3),LBOUND(EDGEQ,4))), SIZE(EDGEQ,1), SIZE(EDGEQ,1)*SIZE(EDGEQ,2), SIZE(EDGEQ,1)*SIZE(EDGEQ,2)*SIZE(EDGEQ,3), SIZE(EDGEQ,1), SIZE(EDGEQ,2), SIZE(EDGEQ,3), SIZE(EDGEQ,4), C_SIZEOF(EDGEQ(LBOUND(EDGEQ,1),LBOUND(EDGEQ,2),LBOUND(EDGEQ,3),LBOUND(EDGEQ,4))))
    ELSE
      SCALAR_INT_BUF(156) = 0
      CALL fstarpu_tensor_data_register(handles(156), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(EDGEQ)
    IF (ASSOCIATED(PHI)) THEN
      SCALAR_INT_BUF(157) = 1
      CALL fstarpu_vector_data_register(handles(157), 0, C_LOC(PHI(LBOUND(PHI,1))), SIZE(PHI,1), C_SIZEOF(PHI(LBOUND(PHI,1))))
    ELSE
      SCALAR_INT_BUF(157) = 0
      CALL fstarpu_vector_data_register(handles(157), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(PHI)
    IF (ASSOCIATED(DPHIDZ1)) THEN
      SCALAR_INT_BUF(158) = 1
      CALL fstarpu_vector_data_register(handles(158), 0, C_LOC(DPHIDZ1(LBOUND(DPHIDZ1,1))), SIZE(DPHIDZ1,1), C_SIZEOF(DPHIDZ1(LBOUND(DPHIDZ1,1))))
    ELSE
      SCALAR_INT_BUF(158) = 0
      CALL fstarpu_vector_data_register(handles(158), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(DPHIDZ1)
    IF (ASSOCIATED(DPHIDZ2)) THEN
      SCALAR_INT_BUF(159) = 1
      CALL fstarpu_vector_data_register(handles(159), 0, C_LOC(DPHIDZ2(LBOUND(DPHIDZ2,1))), SIZE(DPHIDZ2,1), C_SIZEOF(DPHIDZ2(LBOUND(DPHIDZ2,1))))
    ELSE
      SCALAR_INT_BUF(159) = 0
      CALL fstarpu_vector_data_register(handles(159), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(DPHIDZ2)
    IF (ASSOCIATED(PHI_STAE)) THEN
      SCALAR_INT_BUF(160) = 1
      CALL fstarpu_matrix_data_register(handles(160), 0, C_LOC(PHI_STAE(LBOUND(PHI_STAE,1),LBOUND(PHI_STAE,2))), SIZE(PHI_STAE,1), SIZE(PHI_STAE,1), SIZE(PHI_STAE,2), C_SIZEOF(PHI_STAE(LBOUND(PHI_STAE,1),LBOUND(PHI_STAE,2))))
    ELSE
      SCALAR_INT_BUF(160) = 0
      CALL fstarpu_matrix_data_register(handles(160), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(PHI_STAE)
    IF (ASSOCIATED(PHI_STAV)) THEN
      SCALAR_INT_BUF(161) = 1
      CALL fstarpu_matrix_data_register(handles(161), 0, C_LOC(PHI_STAV(LBOUND(PHI_STAV,1),LBOUND(PHI_STAV,2))), SIZE(PHI_STAV,1), SIZE(PHI_STAV,1), SIZE(PHI_STAV,2), C_SIZEOF(PHI_STAV(LBOUND(PHI_STAV,1),LBOUND(PHI_STAV,2))))
    ELSE
      SCALAR_INT_BUF(161) = 0
      CALL fstarpu_matrix_data_register(handles(161), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(PHI_STAV)
    IF (ASSOCIATED(bed_IN)) THEN
      SCALAR_INT_BUF(162) = 1
      CALL fstarpu_vector_data_register(handles(162), 0, C_LOC(bed_IN(LBOUND(bed_IN,1))), SIZE(bed_IN,1), C_SIZEOF(bed_IN(LBOUND(bed_IN,1))))
    ELSE
      SCALAR_INT_BUF(162) = 0
      CALL fstarpu_vector_data_register(handles(162), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(bed_IN)
    IF (ASSOCIATED(bed_EX)) THEN
      SCALAR_INT_BUF(163) = 1
      CALL fstarpu_vector_data_register(handles(163), 0, C_LOC(bed_EX(LBOUND(bed_EX,1))), SIZE(bed_EX,1), C_SIZEOF(bed_EX(LBOUND(bed_EX,1))))
    ELSE
      SCALAR_INT_BUF(163) = 0
      CALL fstarpu_vector_data_register(handles(163), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(bed_EX)
    IF (ASSOCIATED(bed_HAT)) THEN
      SCALAR_INT_BUF(164) = 1
      CALL fstarpu_vector_data_register(handles(164), 0, C_LOC(bed_HAT(LBOUND(bed_HAT,1))), SIZE(bed_HAT,1), C_SIZEOF(bed_HAT(LBOUND(bed_HAT,1))))
    ELSE
      SCALAR_INT_BUF(164) = 0
      CALL fstarpu_vector_data_register(handles(164), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(bed_HAT)
    IF (ASSOCIATED(fact)) THEN
      SCALAR_INT_BUF(165) = 1
      CALL fstarpu_vector_data_register(handles(165), 0, C_LOC(fact(LBOUND(fact,1))), SIZE(fact,1), C_SIZEOF(fact(LBOUND(fact,1))))
    ELSE
      SCALAR_INT_BUF(165) = 0
      CALL fstarpu_vector_data_register(handles(165), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(fact)
    IF (ASSOCIATED(focal_neigh)) THEN
      SCALAR_INT_BUF(166) = 1
      CALL fstarpu_matrix_data_register(handles(166), 0, C_LOC(focal_neigh(LBOUND(focal_neigh,1),LBOUND(focal_neigh,2))), SIZE(focal_neigh,1), SIZE(focal_neigh,1), SIZE(focal_neigh,2), C_SIZEOF(focal_neigh(LBOUND(focal_neigh,1),LBOUND(focal_neigh,2))))
    ELSE
      SCALAR_INT_BUF(166) = 0
      CALL fstarpu_matrix_data_register(handles(166), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(focal_neigh)
    IF (ASSOCIATED(focal_up)) THEN
      SCALAR_INT_BUF(167) = 1
      CALL fstarpu_vector_data_register(handles(167), 0, C_LOC(focal_up(LBOUND(focal_up,1))), SIZE(focal_up,1), C_SIZEOF(focal_up(LBOUND(focal_up,1))))
    ELSE
      SCALAR_INT_BUF(167) = 0
      CALL fstarpu_vector_data_register(handles(167), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(focal_up)
    IF (ASSOCIATED(bi)) THEN
      SCALAR_INT_BUF(168) = 1
      CALL fstarpu_vector_data_register(handles(168), 0, C_LOC(bi(LBOUND(bi,1))), SIZE(bi,1), C_SIZEOF(bi(LBOUND(bi,1))))
    ELSE
      SCALAR_INT_BUF(168) = 0
      CALL fstarpu_vector_data_register(handles(168), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(bi)
    IF (ASSOCIATED(bj)) THEN
      SCALAR_INT_BUF(169) = 1
      CALL fstarpu_vector_data_register(handles(169), 0, C_LOC(bj(LBOUND(bj,1))), SIZE(bj,1), C_SIZEOF(bj(LBOUND(bj,1))))
    ELSE
      SCALAR_INT_BUF(169) = 0
      CALL fstarpu_vector_data_register(handles(169), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(bj)
    IF (ASSOCIATED(XBCb)) THEN
      SCALAR_INT_BUF(170) = 1
      CALL fstarpu_vector_data_register(handles(170), 0, C_LOC(XBCb(LBOUND(XBCb,1))), SIZE(XBCb,1), C_SIZEOF(XBCb(LBOUND(XBCb,1))))
    ELSE
      SCALAR_INT_BUF(170) = 0
      CALL fstarpu_vector_data_register(handles(170), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(XBCb)
    IF (ASSOCIATED(YBCb)) THEN
      SCALAR_INT_BUF(171) = 1
      CALL fstarpu_vector_data_register(handles(171), 0, C_LOC(YBCb(LBOUND(YBCb,1))), SIZE(YBCb,1), C_SIZEOF(YBCb(LBOUND(YBCb,1))))
    ELSE
      SCALAR_INT_BUF(171) = 0
      CALL fstarpu_vector_data_register(handles(171), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(YBCb)
    IF (ASSOCIATED(xi1)) THEN
      SCALAR_INT_BUF(172) = 1
      CALL fstarpu_matrix_data_register(handles(172), 0, C_LOC(xi1(LBOUND(xi1,1),LBOUND(xi1,2))), SIZE(xi1,1), SIZE(xi1,1), SIZE(xi1,2), C_SIZEOF(xi1(LBOUND(xi1,1),LBOUND(xi1,2))))
    ELSE
      SCALAR_INT_BUF(172) = 0
      CALL fstarpu_matrix_data_register(handles(172), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(xi1)
    IF (ASSOCIATED(xi2)) THEN
      SCALAR_INT_BUF(173) = 1
      CALL fstarpu_matrix_data_register(handles(173), 0, C_LOC(xi2(LBOUND(xi2,1),LBOUND(xi2,2))), SIZE(xi2,1), SIZE(xi2,1), SIZE(xi2,2), C_SIZEOF(xi2(LBOUND(xi2,1),LBOUND(xi2,2))))
    ELSE
      SCALAR_INT_BUF(173) = 0
      CALL fstarpu_matrix_data_register(handles(173), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(xi2)
    IF (ASSOCIATED(xtransform)) THEN
      SCALAR_INT_BUF(174) = 1
      CALL fstarpu_matrix_data_register(handles(174), 0, C_LOC(xtransform(LBOUND(xtransform,1),LBOUND(xtransform,2))), SIZE(xtransform,1), SIZE(xtransform,1), SIZE(xtransform,2), C_SIZEOF(xtransform(LBOUND(xtransform,1),LBOUND(xtransform,2))))
    ELSE
      SCALAR_INT_BUF(174) = 0
      CALL fstarpu_matrix_data_register(handles(174), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(xtransform)
    IF (ASSOCIATED(ytransform)) THEN
      SCALAR_INT_BUF(175) = 1
      CALL fstarpu_matrix_data_register(handles(175), 0, C_LOC(ytransform(LBOUND(ytransform,1),LBOUND(ytransform,2))), SIZE(ytransform,1), SIZE(ytransform,1), SIZE(ytransform,2), C_SIZEOF(ytransform(LBOUND(ytransform,1),LBOUND(ytransform,2))))
    ELSE
      SCALAR_INT_BUF(175) = 0
      CALL fstarpu_matrix_data_register(handles(175), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(ytransform)
    IF (ASSOCIATED(xi1BCb)) THEN
      SCALAR_INT_BUF(176) = 1
      CALL fstarpu_vector_data_register(handles(176), 0, C_LOC(xi1BCb(LBOUND(xi1BCb,1))), SIZE(xi1BCb,1), C_SIZEOF(xi1BCb(LBOUND(xi1BCb,1))))
    ELSE
      SCALAR_INT_BUF(176) = 0
      CALL fstarpu_vector_data_register(handles(176), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(xi1BCb)
    IF (ASSOCIATED(xi2BCb)) THEN
      SCALAR_INT_BUF(177) = 1
      CALL fstarpu_vector_data_register(handles(177), 0, C_LOC(xi2BCb(LBOUND(xi2BCb,1))), SIZE(xi2BCb,1), C_SIZEOF(xi2BCb(LBOUND(xi2BCb,1))))
    ELSE
      SCALAR_INT_BUF(177) = 0
      CALL fstarpu_vector_data_register(handles(177), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(xi2BCb)
    IF (ASSOCIATED(xi1vert)) THEN
      SCALAR_INT_BUF(178) = 1
      CALL fstarpu_matrix_data_register(handles(178), 0, C_LOC(xi1vert(LBOUND(xi1vert,1),LBOUND(xi1vert,2))), SIZE(xi1vert,1), SIZE(xi1vert,1), SIZE(xi1vert,2), C_SIZEOF(xi1vert(LBOUND(xi1vert,1),LBOUND(xi1vert,2))))
    ELSE
      SCALAR_INT_BUF(178) = 0
      CALL fstarpu_matrix_data_register(handles(178), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(xi1vert)
    IF (ASSOCIATED(xi2vert)) THEN
      SCALAR_INT_BUF(179) = 1
      CALL fstarpu_matrix_data_register(handles(179), 0, C_LOC(xi2vert(LBOUND(xi2vert,1),LBOUND(xi2vert,2))), SIZE(xi2vert,1), SIZE(xi2vert,1), SIZE(xi2vert,2), C_SIZEOF(xi2vert(LBOUND(xi2vert,1),LBOUND(xi2vert,2))))
    ELSE
      SCALAR_INT_BUF(179) = 0
      CALL fstarpu_matrix_data_register(handles(179), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(xi2vert)
    IF (ASSOCIATED(xtransformv)) THEN
      SCALAR_INT_BUF(180) = 1
      CALL fstarpu_matrix_data_register(handles(180), 0, C_LOC(xtransformv(LBOUND(xtransformv,1),LBOUND(xtransformv,2))), SIZE(xtransformv,1), SIZE(xtransformv,1), SIZE(xtransformv,2), C_SIZEOF(xtransformv(LBOUND(xtransformv,1),LBOUND(xtransformv,2))))
    ELSE
      SCALAR_INT_BUF(180) = 0
      CALL fstarpu_matrix_data_register(handles(180), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(xtransformv)
    IF (ASSOCIATED(ytransformv)) THEN
      SCALAR_INT_BUF(181) = 1
      CALL fstarpu_matrix_data_register(handles(181), 0, C_LOC(ytransformv(LBOUND(ytransformv,1),LBOUND(ytransformv,2))), SIZE(ytransformv,1), SIZE(ytransformv,1), SIZE(ytransformv,2), C_SIZEOF(ytransformv(LBOUND(ytransformv,1),LBOUND(ytransformv,2))))
    ELSE
      SCALAR_INT_BUF(181) = 0
      CALL fstarpu_matrix_data_register(handles(181), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(ytransformv)
    IF (ASSOCIATED(XBCv)) THEN
      SCALAR_INT_BUF(182) = 1
      CALL fstarpu_matrix_data_register(handles(182), 0, C_LOC(XBCv(LBOUND(XBCv,1),LBOUND(XBCv,2))), SIZE(XBCv,1), SIZE(XBCv,1), SIZE(XBCv,2), C_SIZEOF(XBCv(LBOUND(XBCv,1),LBOUND(XBCv,2))))
    ELSE
      SCALAR_INT_BUF(182) = 0
      CALL fstarpu_matrix_data_register(handles(182), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(XBCv)
    IF (ASSOCIATED(YBCv)) THEN
      SCALAR_INT_BUF(183) = 1
      CALL fstarpu_matrix_data_register(handles(183), 0, C_LOC(YBCv(LBOUND(YBCv,1),LBOUND(YBCv,2))), SIZE(YBCv,1), SIZE(YBCv,1), SIZE(YBCv,2), C_SIZEOF(YBCv(LBOUND(YBCv,1),LBOUND(YBCv,2))))
    ELSE
      SCALAR_INT_BUF(183) = 0
      CALL fstarpu_matrix_data_register(handles(183), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(YBCv)
    IF (ASSOCIATED(xi1BCv)) THEN
      SCALAR_INT_BUF(184) = 1
      CALL fstarpu_matrix_data_register(handles(184), 0, C_LOC(xi1BCv(LBOUND(xi1BCv,1),LBOUND(xi1BCv,2))), SIZE(xi1BCv,1), SIZE(xi1BCv,1), SIZE(xi1BCv,2), C_SIZEOF(xi1BCv(LBOUND(xi1BCv,1),LBOUND(xi1BCv,2))))
    ELSE
      SCALAR_INT_BUF(184) = 0
      CALL fstarpu_matrix_data_register(handles(184), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(xi1BCv)
    IF (ASSOCIATED(xi2BCv)) THEN
      SCALAR_INT_BUF(185) = 1
      CALL fstarpu_matrix_data_register(handles(185), 0, C_LOC(xi2BCv(LBOUND(xi2BCv,1),LBOUND(xi2BCv,2))), SIZE(xi2BCv,1), SIZE(xi2BCv,1), SIZE(xi2BCv,2), C_SIZEOF(xi2BCv(LBOUND(xi2BCv,1),LBOUND(xi2BCv,2))))
    ELSE
      SCALAR_INT_BUF(185) = 0
      CALL fstarpu_matrix_data_register(handles(185), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(xi2BCv)
    IF (ASSOCIATED(Area_integral)) THEN
      SCALAR_INT_BUF(186) = 1
      CALL fstarpu_block_data_register(handles(186), 0, C_LOC(Area_integral(LBOUND(Area_integral,1),LBOUND(Area_integral,2),LBOUND(Area_integral,3))), SIZE(Area_integral,1), SIZE(Area_integral,1)*SIZE(Area_integral,2), SIZE(Area_integral,1), SIZE(Area_integral,2), SIZE(Area_integral,3), C_SIZEOF(Area_integral(LBOUND(Area_integral,1),LBOUND(Area_integral,2),LBOUND(Area_integral,3))))
    ELSE
      SCALAR_INT_BUF(186) = 0
      CALL fstarpu_block_data_register(handles(186), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(Area_integral)
    IF (ASSOCIATED(f)) THEN
      SCALAR_INT_BUF(187) = 1
      CALL fstarpu_tensor_data_register(handles(187), 0, C_LOC(f(LBOUND(f,1),LBOUND(f,2),LBOUND(f,3),LBOUND(f,4))), SIZE(f,1), SIZE(f,1)*SIZE(f,2), SIZE(f,1)*SIZE(f,2)*SIZE(f,3), SIZE(f,1), SIZE(f,2), SIZE(f,3), SIZE(f,4), C_SIZEOF(f(LBOUND(f,1),LBOUND(f,2),LBOUND(f,3),LBOUND(f,4))))
    ELSE
      SCALAR_INT_BUF(187) = 0
      CALL fstarpu_tensor_data_register(handles(187), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(f)
    IF (ASSOCIATED(g0)) THEN
      SCALAR_INT_BUF(188) = 1
      CALL fstarpu_tensor_data_register(handles(188), 0, C_LOC(g0(LBOUND(g0,1),LBOUND(g0,2),LBOUND(g0,3),LBOUND(g0,4))), SIZE(g0,1), SIZE(g0,1)*SIZE(g0,2), SIZE(g0,1)*SIZE(g0,2)*SIZE(g0,3), SIZE(g0,1), SIZE(g0,2), SIZE(g0,3), SIZE(g0,4), C_SIZEOF(g0(LBOUND(g0,1),LBOUND(g0,2),LBOUND(g0,3),LBOUND(g0,4))))
    ELSE
      SCALAR_INT_BUF(188) = 0
      CALL fstarpu_tensor_data_register(handles(188), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(g0)
    IF (ASSOCIATED(varsigma0)) THEN
      SCALAR_INT_BUF(189) = 1
      CALL fstarpu_tensor_data_register(handles(189), 0, C_LOC(varsigma0(LBOUND(varsigma0,1),LBOUND(varsigma0,2),LBOUND(varsigma0,3),LBOUND(varsigma0,4))), SIZE(varsigma0,1), SIZE(varsigma0,1)*SIZE(varsigma0,2), SIZE(varsigma0,1)*SIZE(varsigma0,2)*SIZE(varsigma0,3), SIZE(varsigma0,1), SIZE(varsigma0,2), SIZE(varsigma0,3), SIZE(varsigma0,4), C_SIZEOF(varsigma0(LBOUND(varsigma0,1),LBOUND(varsigma0,2),LBOUND(varsigma0,3),LBOUND(varsigma0,4))))
    ELSE
      SCALAR_INT_BUF(189) = 0
      CALL fstarpu_tensor_data_register(handles(189), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(varsigma0)
    IF (ASSOCIATED(fv)) THEN
      SCALAR_INT_BUF(190) = 1
      CALL fstarpu_tensor_data_register(handles(190), 0, C_LOC(fv(LBOUND(fv,1),LBOUND(fv,2),LBOUND(fv,3),LBOUND(fv,4))), SIZE(fv,1), SIZE(fv,1)*SIZE(fv,2), SIZE(fv,1)*SIZE(fv,2)*SIZE(fv,3), SIZE(fv,1), SIZE(fv,2), SIZE(fv,3), SIZE(fv,4), C_SIZEOF(fv(LBOUND(fv,1),LBOUND(fv,2),LBOUND(fv,3),LBOUND(fv,4))))
    ELSE
      SCALAR_INT_BUF(190) = 0
      CALL fstarpu_tensor_data_register(handles(190), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(fv)
    IF (ASSOCIATED(g0v)) THEN
      SCALAR_INT_BUF(191) = 1
      CALL fstarpu_tensor_data_register(handles(191), 0, C_LOC(g0v(LBOUND(g0v,1),LBOUND(g0v,2),LBOUND(g0v,3),LBOUND(g0v,4))), SIZE(g0v,1), SIZE(g0v,1)*SIZE(g0v,2), SIZE(g0v,1)*SIZE(g0v,2)*SIZE(g0v,3), SIZE(g0v,1), SIZE(g0v,2), SIZE(g0v,3), SIZE(g0v,4), C_SIZEOF(g0v(LBOUND(g0v,1),LBOUND(g0v,2),LBOUND(g0v,3),LBOUND(g0v,4))))
    ELSE
      SCALAR_INT_BUF(191) = 0
      CALL fstarpu_tensor_data_register(handles(191), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(g0v)
    IF (ASSOCIATED(varsigma0v)) THEN
      SCALAR_INT_BUF(192) = 1
      CALL fstarpu_tensor_data_register(handles(192), 0, C_LOC(varsigma0v(LBOUND(varsigma0v,1),LBOUND(varsigma0v,2),LBOUND(varsigma0v,3),LBOUND(varsigma0v,4))), SIZE(varsigma0v,1), SIZE(varsigma0v,1)*SIZE(varsigma0v,2), SIZE(varsigma0v,1)*SIZE(varsigma0v,2)*SIZE(varsigma0v,3), SIZE(varsigma0v,1), SIZE(varsigma0v,2), SIZE(varsigma0v,3), SIZE(varsigma0v,4), C_SIZEOF(varsigma0v(LBOUND(varsigma0v,1),LBOUND(varsigma0v,2),LBOUND(varsigma0v,3),LBOUND(varsigma0v,4))))
    ELSE
      SCALAR_INT_BUF(192) = 0
      CALL fstarpu_tensor_data_register(handles(192), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(varsigma0v)
    IF (ASSOCIATED(var2sigmag)) THEN
      SCALAR_INT_BUF(193) = 1
      CALL fstarpu_block_data_register(handles(193), 0, C_LOC(var2sigmag(LBOUND(var2sigmag,1),LBOUND(var2sigmag,2),LBOUND(var2sigmag,3))), SIZE(var2sigmag,1), SIZE(var2sigmag,1)*SIZE(var2sigmag,2), SIZE(var2sigmag,1), SIZE(var2sigmag,2), SIZE(var2sigmag,3), C_SIZEOF(var2sigmag(LBOUND(var2sigmag,1),LBOUND(var2sigmag,2),LBOUND(var2sigmag,3))))
    ELSE
      SCALAR_INT_BUF(193) = 0
      CALL fstarpu_block_data_register(handles(193), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(var2sigmag)
    IF (ASSOCIATED(var2sigmav)) THEN
      SCALAR_INT_BUF(194) = 1
      CALL fstarpu_block_data_register(handles(194), 0, C_LOC(var2sigmav(LBOUND(var2sigmav,1),LBOUND(var2sigmav,2),LBOUND(var2sigmav,3))), SIZE(var2sigmav,1), SIZE(var2sigmav,1)*SIZE(var2sigmav,2), SIZE(var2sigmav,1), SIZE(var2sigmav,2), SIZE(var2sigmav,3), C_SIZEOF(var2sigmav(LBOUND(var2sigmav,1),LBOUND(var2sigmav,2),LBOUND(var2sigmav,3))))
    ELSE
      SCALAR_INT_BUF(194) = 0
      CALL fstarpu_block_data_register(handles(194), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(var2sigmav)
    IF (ASSOCIATED(Nmatrix)) THEN
      SCALAR_INT_BUF(195) = 1
      CALL fstarpu_tensor_data_register(handles(195), 0, C_LOC(Nmatrix(LBOUND(Nmatrix,1),LBOUND(Nmatrix,2),LBOUND(Nmatrix,3),LBOUND(Nmatrix,4))), SIZE(Nmatrix,1), SIZE(Nmatrix,1)*SIZE(Nmatrix,2), SIZE(Nmatrix,1)*SIZE(Nmatrix,2)*SIZE(Nmatrix,3), SIZE(Nmatrix,1), SIZE(Nmatrix,2), SIZE(Nmatrix,3), SIZE(Nmatrix,4), C_SIZEOF(Nmatrix(LBOUND(Nmatrix,1),LBOUND(Nmatrix,2),LBOUND(Nmatrix,3),LBOUND(Nmatrix,4))))
    ELSE
      SCALAR_INT_BUF(195) = 0
      CALL fstarpu_tensor_data_register(handles(195), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(Nmatrix)
    IF (ASSOCIATED(NmatrixInv)) THEN
      SCALAR_INT_BUF(196) = 1
      CALL fstarpu_tensor_data_register(handles(196), 0, C_LOC(NmatrixInv(LBOUND(NmatrixInv,1),LBOUND(NmatrixInv,2),LBOUND(NmatrixInv,3),LBOUND(NmatrixInv,4))), SIZE(NmatrixInv,1), SIZE(NmatrixInv,1)*SIZE(NmatrixInv,2), SIZE(NmatrixInv,1)*SIZE(NmatrixInv,2)*SIZE(NmatrixInv,3), SIZE(NmatrixInv,1), SIZE(NmatrixInv,2), SIZE(NmatrixInv,3), SIZE(NmatrixInv,4), C_SIZEOF(NmatrixInv(LBOUND(NmatrixInv,1),LBOUND(NmatrixInv,2),LBOUND(NmatrixInv,3),LBOUND(NmatrixInv,4))))
    ELSE
      SCALAR_INT_BUF(196) = 0
      CALL fstarpu_tensor_data_register(handles(196), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(NmatrixInv)
    IF (ASSOCIATED(deltx)) THEN
      SCALAR_INT_BUF(197) = 1
      CALL fstarpu_vector_data_register(handles(197), 0, C_LOC(deltx(LBOUND(deltx,1))), SIZE(deltx,1), C_SIZEOF(deltx(LBOUND(deltx,1))))
    ELSE
      SCALAR_INT_BUF(197) = 0
      CALL fstarpu_vector_data_register(handles(197), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(deltx)
    IF (ASSOCIATED(delty)) THEN
      SCALAR_INT_BUF(198) = 1
      CALL fstarpu_vector_data_register(handles(198), 0, C_LOC(delty(LBOUND(delty,1))), SIZE(delty,1), C_SIZEOF(delty(LBOUND(delty,1))))
    ELSE
      SCALAR_INT_BUF(198) = 0
      CALL fstarpu_vector_data_register(handles(198), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(delty)
    IF (ASSOCIATED(pmatrix)) THEN
      SCALAR_INT_BUF(199) = 1
      CALL fstarpu_block_data_register(handles(199), 0, C_LOC(pmatrix(LBOUND(pmatrix,1),LBOUND(pmatrix,2),LBOUND(pmatrix,3))), SIZE(pmatrix,1), SIZE(pmatrix,1)*SIZE(pmatrix,2), SIZE(pmatrix,1), SIZE(pmatrix,2), SIZE(pmatrix,3), C_SIZEOF(pmatrix(LBOUND(pmatrix,1),LBOUND(pmatrix,2),LBOUND(pmatrix,3))))
    ELSE
      SCALAR_INT_BUF(199) = 0
      CALL fstarpu_block_data_register(handles(199), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(pmatrix)
    IF (ASSOCIATED(ZEmin)) THEN
      SCALAR_INT_BUF(200) = 1
      CALL fstarpu_matrix_data_register(handles(200), 0, C_LOC(ZEmin(LBOUND(ZEmin,1),LBOUND(ZEmin,2))), SIZE(ZEmin,1), SIZE(ZEmin,1), SIZE(ZEmin,2), C_SIZEOF(ZEmin(LBOUND(ZEmin,1),LBOUND(ZEmin,2))))
    ELSE
      SCALAR_INT_BUF(200) = 0
      CALL fstarpu_matrix_data_register(handles(200), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(ZEmin)
    IF (ASSOCIATED(ZEmax)) THEN
      SCALAR_INT_BUF(201) = 1
      CALL fstarpu_matrix_data_register(handles(201), 0, C_LOC(ZEmax(LBOUND(ZEmax,1),LBOUND(ZEmax,2))), SIZE(ZEmax,1), SIZE(ZEmax,1), SIZE(ZEmax,2), C_SIZEOF(ZEmax(LBOUND(ZEmax,1),LBOUND(ZEmax,2))))
    ELSE
      SCALAR_INT_BUF(201) = 0
      CALL fstarpu_matrix_data_register(handles(201), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(ZEmax)
    IF (ASSOCIATED(QXmin)) THEN
      SCALAR_INT_BUF(202) = 1
      CALL fstarpu_matrix_data_register(handles(202), 0, C_LOC(QXmin(LBOUND(QXmin,1),LBOUND(QXmin,2))), SIZE(QXmin,1), SIZE(QXmin,1), SIZE(QXmin,2), C_SIZEOF(QXmin(LBOUND(QXmin,1),LBOUND(QXmin,2))))
    ELSE
      SCALAR_INT_BUF(202) = 0
      CALL fstarpu_matrix_data_register(handles(202), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(QXmin)
    IF (ASSOCIATED(QXmax)) THEN
      SCALAR_INT_BUF(203) = 1
      CALL fstarpu_matrix_data_register(handles(203), 0, C_LOC(QXmax(LBOUND(QXmax,1),LBOUND(QXmax,2))), SIZE(QXmax,1), SIZE(QXmax,1), SIZE(QXmax,2), C_SIZEOF(QXmax(LBOUND(QXmax,1),LBOUND(QXmax,2))))
    ELSE
      SCALAR_INT_BUF(203) = 0
      CALL fstarpu_matrix_data_register(handles(203), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(QXmax)
    IF (ASSOCIATED(QYmin)) THEN
      SCALAR_INT_BUF(204) = 1
      CALL fstarpu_matrix_data_register(handles(204), 0, C_LOC(QYmin(LBOUND(QYmin,1),LBOUND(QYmin,2))), SIZE(QYmin,1), SIZE(QYmin,1), SIZE(QYmin,2), C_SIZEOF(QYmin(LBOUND(QYmin,1),LBOUND(QYmin,2))))
    ELSE
      SCALAR_INT_BUF(204) = 0
      CALL fstarpu_matrix_data_register(handles(204), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(QYmin)
    IF (ASSOCIATED(QYmax)) THEN
      SCALAR_INT_BUF(205) = 1
      CALL fstarpu_matrix_data_register(handles(205), 0, C_LOC(QYmax(LBOUND(QYmax,1),LBOUND(QYmax,2))), SIZE(QYmax,1), SIZE(QYmax,1), SIZE(QYmax,2), C_SIZEOF(QYmax(LBOUND(QYmax,1),LBOUND(QYmax,2))))
    ELSE
      SCALAR_INT_BUF(205) = 0
      CALL fstarpu_matrix_data_register(handles(205), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(QYmax)
    IF (ASSOCIATED(iotamin)) THEN
      SCALAR_INT_BUF(206) = 1
      CALL fstarpu_matrix_data_register(handles(206), 0, C_LOC(iotamin(LBOUND(iotamin,1),LBOUND(iotamin,2))), SIZE(iotamin,1), SIZE(iotamin,1), SIZE(iotamin,2), C_SIZEOF(iotamin(LBOUND(iotamin,1),LBOUND(iotamin,2))))
    ELSE
      SCALAR_INT_BUF(206) = 0
      CALL fstarpu_matrix_data_register(handles(206), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(iotamin)
    IF (ASSOCIATED(iotamax)) THEN
      SCALAR_INT_BUF(207) = 1
      CALL fstarpu_matrix_data_register(handles(207), 0, C_LOC(iotamax(LBOUND(iotamax,1),LBOUND(iotamax,2))), SIZE(iotamax,1), SIZE(iotamax,1), SIZE(iotamax,2), C_SIZEOF(iotamax(LBOUND(iotamax,1),LBOUND(iotamax,2))))
    ELSE
      SCALAR_INT_BUF(207) = 0
      CALL fstarpu_matrix_data_register(handles(207), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(iotamax)
    IF (ASSOCIATED(iota2min)) THEN
      SCALAR_INT_BUF(208) = 1
      CALL fstarpu_matrix_data_register(handles(208), 0, C_LOC(iota2min(LBOUND(iota2min,1),LBOUND(iota2min,2))), SIZE(iota2min,1), SIZE(iota2min,1), SIZE(iota2min,2), C_SIZEOF(iota2min(LBOUND(iota2min,1),LBOUND(iota2min,2))))
    ELSE
      SCALAR_INT_BUF(208) = 0
      CALL fstarpu_matrix_data_register(handles(208), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(iota2min)
    IF (ASSOCIATED(iota2max)) THEN
      SCALAR_INT_BUF(209) = 1
      CALL fstarpu_matrix_data_register(handles(209), 0, C_LOC(iota2max(LBOUND(iota2max,1),LBOUND(iota2max,2))), SIZE(iota2max,1), SIZE(iota2max,1), SIZE(iota2max,2), C_SIZEOF(iota2max(LBOUND(iota2max,1),LBOUND(iota2max,2))))
    ELSE
      SCALAR_INT_BUF(209) = 0
      CALL fstarpu_matrix_data_register(handles(209), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(iota2max)
    IF (ASSOCIATED(ZEtaylor)) THEN
      SCALAR_INT_BUF(210) = 1
      CALL fstarpu_block_data_register(handles(210), 0, C_LOC(ZEtaylor(LBOUND(ZEtaylor,1),LBOUND(ZEtaylor,2),LBOUND(ZEtaylor,3))), SIZE(ZEtaylor,1), SIZE(ZEtaylor,1)*SIZE(ZEtaylor,2), SIZE(ZEtaylor,1), SIZE(ZEtaylor,2), SIZE(ZEtaylor,3), C_SIZEOF(ZEtaylor(LBOUND(ZEtaylor,1),LBOUND(ZEtaylor,2),LBOUND(ZEtaylor,3))))
    ELSE
      SCALAR_INT_BUF(210) = 0
      CALL fstarpu_block_data_register(handles(210), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(ZEtaylor)
    IF (ASSOCIATED(QXtaylor)) THEN
      SCALAR_INT_BUF(211) = 1
      CALL fstarpu_block_data_register(handles(211), 0, C_LOC(QXtaylor(LBOUND(QXtaylor,1),LBOUND(QXtaylor,2),LBOUND(QXtaylor,3))), SIZE(QXtaylor,1), SIZE(QXtaylor,1)*SIZE(QXtaylor,2), SIZE(QXtaylor,1), SIZE(QXtaylor,2), SIZE(QXtaylor,3), C_SIZEOF(QXtaylor(LBOUND(QXtaylor,1),LBOUND(QXtaylor,2),LBOUND(QXtaylor,3))))
    ELSE
      SCALAR_INT_BUF(211) = 0
      CALL fstarpu_block_data_register(handles(211), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(QXtaylor)
    IF (ASSOCIATED(QYtaylor)) THEN
      SCALAR_INT_BUF(212) = 1
      CALL fstarpu_block_data_register(handles(212), 0, C_LOC(QYtaylor(LBOUND(QYtaylor,1),LBOUND(QYtaylor,2),LBOUND(QYtaylor,3))), SIZE(QYtaylor,1), SIZE(QYtaylor,1)*SIZE(QYtaylor,2), SIZE(QYtaylor,1), SIZE(QYtaylor,2), SIZE(QYtaylor,3), C_SIZEOF(QYtaylor(LBOUND(QYtaylor,1),LBOUND(QYtaylor,2),LBOUND(QYtaylor,3))))
    ELSE
      SCALAR_INT_BUF(212) = 0
      CALL fstarpu_block_data_register(handles(212), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(QYtaylor)
    IF (ASSOCIATED(iotataylor)) THEN
      SCALAR_INT_BUF(213) = 1
      CALL fstarpu_block_data_register(handles(213), 0, C_LOC(iotataylor(LBOUND(iotataylor,1),LBOUND(iotataylor,2),LBOUND(iotataylor,3))), SIZE(iotataylor,1), SIZE(iotataylor,1)*SIZE(iotataylor,2), SIZE(iotataylor,1), SIZE(iotataylor,2), SIZE(iotataylor,3), C_SIZEOF(iotataylor(LBOUND(iotataylor,1),LBOUND(iotataylor,2),LBOUND(iotataylor,3))))
    ELSE
      SCALAR_INT_BUF(213) = 0
      CALL fstarpu_block_data_register(handles(213), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(iotataylor)
    IF (ASSOCIATED(iota2taylor)) THEN
      SCALAR_INT_BUF(214) = 1
      CALL fstarpu_block_data_register(handles(214), 0, C_LOC(iota2taylor(LBOUND(iota2taylor,1),LBOUND(iota2taylor,2),LBOUND(iota2taylor,3))), SIZE(iota2taylor,1), SIZE(iota2taylor,1)*SIZE(iota2taylor,2), SIZE(iota2taylor,1), SIZE(iota2taylor,2), SIZE(iota2taylor,3), C_SIZEOF(iota2taylor(LBOUND(iota2taylor,1),LBOUND(iota2taylor,2),LBOUND(iota2taylor,3))))
    ELSE
      SCALAR_INT_BUF(214) = 0
      CALL fstarpu_block_data_register(handles(214), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(iota2taylor)
    IF (ASSOCIATED(ZEtaylorvert)) THEN
      SCALAR_INT_BUF(215) = 1
      CALL fstarpu_block_data_register(handles(215), 0, C_LOC(ZEtaylorvert(LBOUND(ZEtaylorvert,1),LBOUND(ZEtaylorvert,2),LBOUND(ZEtaylorvert,3))), SIZE(ZEtaylorvert,1), SIZE(ZEtaylorvert,1)*SIZE(ZEtaylorvert,2), SIZE(ZEtaylorvert,1), SIZE(ZEtaylorvert,2), SIZE(ZEtaylorvert,3), C_SIZEOF(ZEtaylorvert(LBOUND(ZEtaylorvert,1),LBOUND(ZEtaylorvert,2),LBOUND(ZEtaylorvert,3))))
    ELSE
      SCALAR_INT_BUF(215) = 0
      CALL fstarpu_block_data_register(handles(215), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(ZEtaylorvert)
    IF (ASSOCIATED(QXtaylorvert)) THEN
      SCALAR_INT_BUF(216) = 1
      CALL fstarpu_block_data_register(handles(216), 0, C_LOC(QXtaylorvert(LBOUND(QXtaylorvert,1),LBOUND(QXtaylorvert,2),LBOUND(QXtaylorvert,3))), SIZE(QXtaylorvert,1), SIZE(QXtaylorvert,1)*SIZE(QXtaylorvert,2), SIZE(QXtaylorvert,1), SIZE(QXtaylorvert,2), SIZE(QXtaylorvert,3), C_SIZEOF(QXtaylorvert(LBOUND(QXtaylorvert,1),LBOUND(QXtaylorvert,2),LBOUND(QXtaylorvert,3))))
    ELSE
      SCALAR_INT_BUF(216) = 0
      CALL fstarpu_block_data_register(handles(216), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(QXtaylorvert)
    IF (ASSOCIATED(QYtaylorvert)) THEN
      SCALAR_INT_BUF(217) = 1
      CALL fstarpu_block_data_register(handles(217), 0, C_LOC(QYtaylorvert(LBOUND(QYtaylorvert,1),LBOUND(QYtaylorvert,2),LBOUND(QYtaylorvert,3))), SIZE(QYtaylorvert,1), SIZE(QYtaylorvert,1)*SIZE(QYtaylorvert,2), SIZE(QYtaylorvert,1), SIZE(QYtaylorvert,2), SIZE(QYtaylorvert,3), C_SIZEOF(QYtaylorvert(LBOUND(QYtaylorvert,1),LBOUND(QYtaylorvert,2),LBOUND(QYtaylorvert,3))))
    ELSE
      SCALAR_INT_BUF(217) = 0
      CALL fstarpu_block_data_register(handles(217), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(QYtaylorvert)
    IF (ASSOCIATED(iotataylorvert)) THEN
      SCALAR_INT_BUF(218) = 1
      CALL fstarpu_block_data_register(handles(218), 0, C_LOC(iotataylorvert(LBOUND(iotataylorvert,1),LBOUND(iotataylorvert,2),LBOUND(iotataylorvert,3))), SIZE(iotataylorvert,1), SIZE(iotataylorvert,1)*SIZE(iotataylorvert,2), SIZE(iotataylorvert,1), SIZE(iotataylorvert,2), SIZE(iotataylorvert,3), C_SIZEOF(iotataylorvert(LBOUND(iotataylorvert,1),LBOUND(iotataylorvert,2),LBOUND(iotataylorvert,3))))
    ELSE
      SCALAR_INT_BUF(218) = 0
      CALL fstarpu_block_data_register(handles(218), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(iotataylorvert)
    IF (ASSOCIATED(iota2taylorvert)) THEN
      SCALAR_INT_BUF(219) = 1
      CALL fstarpu_block_data_register(handles(219), 0, C_LOC(iota2taylorvert(LBOUND(iota2taylorvert,1),LBOUND(iota2taylorvert,2),LBOUND(iota2taylorvert,3))), SIZE(iota2taylorvert,1), SIZE(iota2taylorvert,1)*SIZE(iota2taylorvert,2), SIZE(iota2taylorvert,1), SIZE(iota2taylorvert,2), SIZE(iota2taylorvert,3), C_SIZEOF(iota2taylorvert(LBOUND(iota2taylorvert,1),LBOUND(iota2taylorvert,2),LBOUND(iota2taylorvert,3))))
    ELSE
      SCALAR_INT_BUF(219) = 0
      CALL fstarpu_block_data_register(handles(219), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(iota2taylorvert)
    IF (ASSOCIATED(alphaZE0)) THEN
      SCALAR_INT_BUF(220) = 1
      CALL fstarpu_block_data_register(handles(220), 0, C_LOC(alphaZE0(LBOUND(alphaZE0,1),LBOUND(alphaZE0,2),LBOUND(alphaZE0,3))), SIZE(alphaZE0,1), SIZE(alphaZE0,1)*SIZE(alphaZE0,2), SIZE(alphaZE0,1), SIZE(alphaZE0,2), SIZE(alphaZE0,3), C_SIZEOF(alphaZE0(LBOUND(alphaZE0,1),LBOUND(alphaZE0,2),LBOUND(alphaZE0,3))))
    ELSE
      SCALAR_INT_BUF(220) = 0
      CALL fstarpu_block_data_register(handles(220), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(alphaZE0)
    IF (ASSOCIATED(alphaQX0)) THEN
      SCALAR_INT_BUF(221) = 1
      CALL fstarpu_block_data_register(handles(221), 0, C_LOC(alphaQX0(LBOUND(alphaQX0,1),LBOUND(alphaQX0,2),LBOUND(alphaQX0,3))), SIZE(alphaQX0,1), SIZE(alphaQX0,1)*SIZE(alphaQX0,2), SIZE(alphaQX0,1), SIZE(alphaQX0,2), SIZE(alphaQX0,3), C_SIZEOF(alphaQX0(LBOUND(alphaQX0,1),LBOUND(alphaQX0,2),LBOUND(alphaQX0,3))))
    ELSE
      SCALAR_INT_BUF(221) = 0
      CALL fstarpu_block_data_register(handles(221), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(alphaQX0)
    IF (ASSOCIATED(alphaQY0)) THEN
      SCALAR_INT_BUF(222) = 1
      CALL fstarpu_block_data_register(handles(222), 0, C_LOC(alphaQY0(LBOUND(alphaQY0,1),LBOUND(alphaQY0,2),LBOUND(alphaQY0,3))), SIZE(alphaQY0,1), SIZE(alphaQY0,1)*SIZE(alphaQY0,2), SIZE(alphaQY0,1), SIZE(alphaQY0,2), SIZE(alphaQY0,3), C_SIZEOF(alphaQY0(LBOUND(alphaQY0,1),LBOUND(alphaQY0,2),LBOUND(alphaQY0,3))))
    ELSE
      SCALAR_INT_BUF(222) = 0
      CALL fstarpu_block_data_register(handles(222), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(alphaQY0)
    IF (ASSOCIATED(alphaiota0)) THEN
      SCALAR_INT_BUF(223) = 1
      CALL fstarpu_block_data_register(handles(223), 0, C_LOC(alphaiota0(LBOUND(alphaiota0,1),LBOUND(alphaiota0,2),LBOUND(alphaiota0,3))), SIZE(alphaiota0,1), SIZE(alphaiota0,1)*SIZE(alphaiota0,2), SIZE(alphaiota0,1), SIZE(alphaiota0,2), SIZE(alphaiota0,3), C_SIZEOF(alphaiota0(LBOUND(alphaiota0,1),LBOUND(alphaiota0,2),LBOUND(alphaiota0,3))))
    ELSE
      SCALAR_INT_BUF(223) = 0
      CALL fstarpu_block_data_register(handles(223), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(alphaiota0)
    IF (ASSOCIATED(alphaiota20)) THEN
      SCALAR_INT_BUF(224) = 1
      CALL fstarpu_block_data_register(handles(224), 0, C_LOC(alphaiota20(LBOUND(alphaiota20,1),LBOUND(alphaiota20,2),LBOUND(alphaiota20,3))), SIZE(alphaiota20,1), SIZE(alphaiota20,1)*SIZE(alphaiota20,2), SIZE(alphaiota20,1), SIZE(alphaiota20,2), SIZE(alphaiota20,3), C_SIZEOF(alphaiota20(LBOUND(alphaiota20,1),LBOUND(alphaiota20,2),LBOUND(alphaiota20,3))))
    ELSE
      SCALAR_INT_BUF(224) = 0
      CALL fstarpu_block_data_register(handles(224), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(alphaiota20)
    IF (ASSOCIATED(alphaZE)) THEN
      SCALAR_INT_BUF(225) = 1
      CALL fstarpu_matrix_data_register(handles(225), 0, C_LOC(alphaZE(LBOUND(alphaZE,1),LBOUND(alphaZE,2))), SIZE(alphaZE,1), SIZE(alphaZE,1), SIZE(alphaZE,2), C_SIZEOF(alphaZE(LBOUND(alphaZE,1),LBOUND(alphaZE,2))))
    ELSE
      SCALAR_INT_BUF(225) = 0
      CALL fstarpu_matrix_data_register(handles(225), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(alphaZE)
    IF (ASSOCIATED(alphaQX)) THEN
      SCALAR_INT_BUF(226) = 1
      CALL fstarpu_matrix_data_register(handles(226), 0, C_LOC(alphaQX(LBOUND(alphaQX,1),LBOUND(alphaQX,2))), SIZE(alphaQX,1), SIZE(alphaQX,1), SIZE(alphaQX,2), C_SIZEOF(alphaQX(LBOUND(alphaQX,1),LBOUND(alphaQX,2))))
    ELSE
      SCALAR_INT_BUF(226) = 0
      CALL fstarpu_matrix_data_register(handles(226), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(alphaQX)
    IF (ASSOCIATED(alphaQY)) THEN
      SCALAR_INT_BUF(227) = 1
      CALL fstarpu_matrix_data_register(handles(227), 0, C_LOC(alphaQY(LBOUND(alphaQY,1),LBOUND(alphaQY,2))), SIZE(alphaQY,1), SIZE(alphaQY,1), SIZE(alphaQY,2), C_SIZEOF(alphaQY(LBOUND(alphaQY,1),LBOUND(alphaQY,2))))
    ELSE
      SCALAR_INT_BUF(227) = 0
      CALL fstarpu_matrix_data_register(handles(227), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(alphaQY)
    IF (ASSOCIATED(alphaiota)) THEN
      SCALAR_INT_BUF(228) = 1
      CALL fstarpu_matrix_data_register(handles(228), 0, C_LOC(alphaiota(LBOUND(alphaiota,1),LBOUND(alphaiota,2))), SIZE(alphaiota,1), SIZE(alphaiota,1), SIZE(alphaiota,2), C_SIZEOF(alphaiota(LBOUND(alphaiota,1),LBOUND(alphaiota,2))))
    ELSE
      SCALAR_INT_BUF(228) = 0
      CALL fstarpu_matrix_data_register(handles(228), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(alphaiota)
    IF (ASSOCIATED(alphaiota2)) THEN
      SCALAR_INT_BUF(229) = 1
      CALL fstarpu_matrix_data_register(handles(229), 0, C_LOC(alphaiota2(LBOUND(alphaiota2,1),LBOUND(alphaiota2,2))), SIZE(alphaiota2,1), SIZE(alphaiota2,1), SIZE(alphaiota2,2), C_SIZEOF(alphaiota2(LBOUND(alphaiota2,1),LBOUND(alphaiota2,2))))
    ELSE
      SCALAR_INT_BUF(229) = 0
      CALL fstarpu_matrix_data_register(handles(229), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(alphaiota2)
    IF (ASSOCIATED(alphaZEm)) THEN
      SCALAR_INT_BUF(230) = 1
      CALL fstarpu_matrix_data_register(handles(230), 0, C_LOC(alphaZEm(LBOUND(alphaZEm,1),LBOUND(alphaZEm,2))), SIZE(alphaZEm,1), SIZE(alphaZEm,1), SIZE(alphaZEm,2), C_SIZEOF(alphaZEm(LBOUND(alphaZEm,1),LBOUND(alphaZEm,2))))
    ELSE
      SCALAR_INT_BUF(230) = 0
      CALL fstarpu_matrix_data_register(handles(230), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(alphaZEm)
    IF (ASSOCIATED(alphaQXm)) THEN
      SCALAR_INT_BUF(231) = 1
      CALL fstarpu_matrix_data_register(handles(231), 0, C_LOC(alphaQXm(LBOUND(alphaQXm,1),LBOUND(alphaQXm,2))), SIZE(alphaQXm,1), SIZE(alphaQXm,1), SIZE(alphaQXm,2), C_SIZEOF(alphaQXm(LBOUND(alphaQXm,1),LBOUND(alphaQXm,2))))
    ELSE
      SCALAR_INT_BUF(231) = 0
      CALL fstarpu_matrix_data_register(handles(231), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(alphaQXm)
    IF (ASSOCIATED(alphaQYm)) THEN
      SCALAR_INT_BUF(232) = 1
      CALL fstarpu_matrix_data_register(handles(232), 0, C_LOC(alphaQYm(LBOUND(alphaQYm,1),LBOUND(alphaQYm,2))), SIZE(alphaQYm,1), SIZE(alphaQYm,1), SIZE(alphaQYm,2), C_SIZEOF(alphaQYm(LBOUND(alphaQYm,1),LBOUND(alphaQYm,2))))
    ELSE
      SCALAR_INT_BUF(232) = 0
      CALL fstarpu_matrix_data_register(handles(232), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(alphaQYm)
    IF (ASSOCIATED(alphaiotam)) THEN
      SCALAR_INT_BUF(233) = 1
      CALL fstarpu_matrix_data_register(handles(233), 0, C_LOC(alphaiotam(LBOUND(alphaiotam,1),LBOUND(alphaiotam,2))), SIZE(alphaiotam,1), SIZE(alphaiotam,1), SIZE(alphaiotam,2), C_SIZEOF(alphaiotam(LBOUND(alphaiotam,1),LBOUND(alphaiotam,2))))
    ELSE
      SCALAR_INT_BUF(233) = 0
      CALL fstarpu_matrix_data_register(handles(233), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(alphaiotam)
    IF (ASSOCIATED(alphaiota2m)) THEN
      SCALAR_INT_BUF(234) = 1
      CALL fstarpu_matrix_data_register(handles(234), 0, C_LOC(alphaiota2m(LBOUND(alphaiota2m,1),LBOUND(alphaiota2m,2))), SIZE(alphaiota2m,1), SIZE(alphaiota2m,1), SIZE(alphaiota2m,2), C_SIZEOF(alphaiota2m(LBOUND(alphaiota2m,1),LBOUND(alphaiota2m,2))))
    ELSE
      SCALAR_INT_BUF(234) = 0
      CALL fstarpu_matrix_data_register(handles(234), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(alphaiota2m)
    IF (ASSOCIATED(alphaZE_max)) THEN
      SCALAR_INT_BUF(235) = 1
      CALL fstarpu_matrix_data_register(handles(235), 0, C_LOC(alphaZE_max(LBOUND(alphaZE_max,1),LBOUND(alphaZE_max,2))), SIZE(alphaZE_max,1), SIZE(alphaZE_max,1), SIZE(alphaZE_max,2), C_SIZEOF(alphaZE_max(LBOUND(alphaZE_max,1),LBOUND(alphaZE_max,2))))
    ELSE
      SCALAR_INT_BUF(235) = 0
      CALL fstarpu_matrix_data_register(handles(235), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(alphaZE_max)
    IF (ASSOCIATED(alphaQX_max)) THEN
      SCALAR_INT_BUF(236) = 1
      CALL fstarpu_matrix_data_register(handles(236), 0, C_LOC(alphaQX_max(LBOUND(alphaQX_max,1),LBOUND(alphaQX_max,2))), SIZE(alphaQX_max,1), SIZE(alphaQX_max,1), SIZE(alphaQX_max,2), C_SIZEOF(alphaQX_max(LBOUND(alphaQX_max,1),LBOUND(alphaQX_max,2))))
    ELSE
      SCALAR_INT_BUF(236) = 0
      CALL fstarpu_matrix_data_register(handles(236), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(alphaQX_max)
    IF (ASSOCIATED(alphaQY_max)) THEN
      SCALAR_INT_BUF(237) = 1
      CALL fstarpu_matrix_data_register(handles(237), 0, C_LOC(alphaQY_max(LBOUND(alphaQY_max,1),LBOUND(alphaQY_max,2))), SIZE(alphaQY_max,1), SIZE(alphaQY_max,1), SIZE(alphaQY_max,2), C_SIZEOF(alphaQY_max(LBOUND(alphaQY_max,1),LBOUND(alphaQY_max,2))))
    ELSE
      SCALAR_INT_BUF(237) = 0
      CALL fstarpu_matrix_data_register(handles(237), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(alphaQY_max)
    IF (ASSOCIATED(alphaiota_max)) THEN
      SCALAR_INT_BUF(238) = 1
      CALL fstarpu_matrix_data_register(handles(238), 0, C_LOC(alphaiota_max(LBOUND(alphaiota_max,1),LBOUND(alphaiota_max,2))), SIZE(alphaiota_max,1), SIZE(alphaiota_max,1), SIZE(alphaiota_max,2), C_SIZEOF(alphaiota_max(LBOUND(alphaiota_max,1),LBOUND(alphaiota_max,2))))
    ELSE
      SCALAR_INT_BUF(238) = 0
      CALL fstarpu_matrix_data_register(handles(238), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(alphaiota_max)
    IF (ASSOCIATED(alphaiota2_max)) THEN
      SCALAR_INT_BUF(239) = 1
      CALL fstarpu_matrix_data_register(handles(239), 0, C_LOC(alphaiota2_max(LBOUND(alphaiota2_max,1),LBOUND(alphaiota2_max,2))), SIZE(alphaiota2_max,1), SIZE(alphaiota2_max,1), SIZE(alphaiota2_max,2), C_SIZEOF(alphaiota2_max(LBOUND(alphaiota2_max,1),LBOUND(alphaiota2_max,2))))
    ELSE
      SCALAR_INT_BUF(239) = 0
      CALL fstarpu_matrix_data_register(handles(239), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(alphaiota2_max)
    IF (ASSOCIATED(limitZE)) THEN
      SCALAR_INT_BUF(240) = 1
      CALL fstarpu_matrix_data_register(handles(240), 0, C_LOC(limitZE(LBOUND(limitZE,1),LBOUND(limitZE,2))), SIZE(limitZE,1), SIZE(limitZE,1), SIZE(limitZE,2), C_SIZEOF(limitZE(LBOUND(limitZE,1),LBOUND(limitZE,2))))
    ELSE
      SCALAR_INT_BUF(240) = 0
      CALL fstarpu_matrix_data_register(handles(240), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(limitZE)
    IF (ASSOCIATED(limitQX)) THEN
      SCALAR_INT_BUF(241) = 1
      CALL fstarpu_matrix_data_register(handles(241), 0, C_LOC(limitQX(LBOUND(limitQX,1),LBOUND(limitQX,2))), SIZE(limitQX,1), SIZE(limitQX,1), SIZE(limitQX,2), C_SIZEOF(limitQX(LBOUND(limitQX,1),LBOUND(limitQX,2))))
    ELSE
      SCALAR_INT_BUF(241) = 0
      CALL fstarpu_matrix_data_register(handles(241), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(limitQX)
    IF (ASSOCIATED(limitQY)) THEN
      SCALAR_INT_BUF(242) = 1
      CALL fstarpu_matrix_data_register(handles(242), 0, C_LOC(limitQY(LBOUND(limitQY,1),LBOUND(limitQY,2))), SIZE(limitQY,1), SIZE(limitQY,1), SIZE(limitQY,2), C_SIZEOF(limitQY(LBOUND(limitQY,1),LBOUND(limitQY,2))))
    ELSE
      SCALAR_INT_BUF(242) = 0
      CALL fstarpu_matrix_data_register(handles(242), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(limitQY)
    IF (ASSOCIATED(limitiota)) THEN
      SCALAR_INT_BUF(243) = 1
      CALL fstarpu_matrix_data_register(handles(243), 0, C_LOC(limitiota(LBOUND(limitiota,1),LBOUND(limitiota,2))), SIZE(limitiota,1), SIZE(limitiota,1), SIZE(limitiota,2), C_SIZEOF(limitiota(LBOUND(limitiota,1),LBOUND(limitiota,2))))
    ELSE
      SCALAR_INT_BUF(243) = 0
      CALL fstarpu_matrix_data_register(handles(243), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(limitiota)
    IF (ASSOCIATED(limitiota2)) THEN
      SCALAR_INT_BUF(244) = 1
      CALL fstarpu_matrix_data_register(handles(244), 0, C_LOC(limitiota2(LBOUND(limitiota2,1),LBOUND(limitiota2,2))), SIZE(limitiota2,1), SIZE(limitiota2,1), SIZE(limitiota2,2), C_SIZEOF(limitiota2(LBOUND(limitiota2,1),LBOUND(limitiota2,2))))
    ELSE
      SCALAR_INT_BUF(244) = 0
      CALL fstarpu_matrix_data_register(handles(244), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(limitiota2)
    IF (ASSOCIATED(ZEconst)) THEN
      SCALAR_INT_BUF(245) = 1
      CALL fstarpu_matrix_data_register(handles(245), 0, C_LOC(ZEconst(LBOUND(ZEconst,1),LBOUND(ZEconst,2))), SIZE(ZEconst,1), SIZE(ZEconst,1), SIZE(ZEconst,2), C_SIZEOF(ZEconst(LBOUND(ZEconst,1),LBOUND(ZEconst,2))))
    ELSE
      SCALAR_INT_BUF(245) = 0
      CALL fstarpu_matrix_data_register(handles(245), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(ZEconst)
    IF (ASSOCIATED(QXconst)) THEN
      SCALAR_INT_BUF(246) = 1
      CALL fstarpu_matrix_data_register(handles(246), 0, C_LOC(QXconst(LBOUND(QXconst,1),LBOUND(QXconst,2))), SIZE(QXconst,1), SIZE(QXconst,1), SIZE(QXconst,2), C_SIZEOF(QXconst(LBOUND(QXconst,1),LBOUND(QXconst,2))))
    ELSE
      SCALAR_INT_BUF(246) = 0
      CALL fstarpu_matrix_data_register(handles(246), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(QXconst)
    IF (ASSOCIATED(QYconst)) THEN
      SCALAR_INT_BUF(247) = 1
      CALL fstarpu_matrix_data_register(handles(247), 0, C_LOC(QYconst(LBOUND(QYconst,1),LBOUND(QYconst,2))), SIZE(QYconst,1), SIZE(QYconst,1), SIZE(QYconst,2), C_SIZEOF(QYconst(LBOUND(QYconst,1),LBOUND(QYconst,2))))
    ELSE
      SCALAR_INT_BUF(247) = 0
      CALL fstarpu_matrix_data_register(handles(247), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(QYconst)
    IF (ASSOCIATED(iotaconst)) THEN
      SCALAR_INT_BUF(248) = 1
      CALL fstarpu_matrix_data_register(handles(248), 0, C_LOC(iotaconst(LBOUND(iotaconst,1),LBOUND(iotaconst,2))), SIZE(iotaconst,1), SIZE(iotaconst,1), SIZE(iotaconst,2), C_SIZEOF(iotaconst(LBOUND(iotaconst,1),LBOUND(iotaconst,2))))
    ELSE
      SCALAR_INT_BUF(248) = 0
      CALL fstarpu_matrix_data_register(handles(248), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(iotaconst)
    IF (ASSOCIATED(iota2const)) THEN
      SCALAR_INT_BUF(249) = 1
      CALL fstarpu_matrix_data_register(handles(249), 0, C_LOC(iota2const(LBOUND(iota2const,1),LBOUND(iota2const,2))), SIZE(iota2const,1), SIZE(iota2const,1), SIZE(iota2const,2), C_SIZEOF(iota2const(LBOUND(iota2const,1),LBOUND(iota2const,2))))
    ELSE
      SCALAR_INT_BUF(249) = 0
      CALL fstarpu_matrix_data_register(handles(249), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(iota2const)
    IF (ASSOCIATED(NODES_LG)) THEN
      SCALAR_INT_BUF(250) = 1
      CALL fstarpu_vector_data_register(handles(250), 0, C_LOC(NODES_LG(LBOUND(NODES_LG,1))), SIZE(NODES_LG,1), C_SIZEOF(NODES_LG(LBOUND(NODES_LG,1))))
    ELSE
      SCALAR_INT_BUF(250) = 0
      CALL fstarpu_vector_data_register(handles(250), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(NODES_LG)
    IF (ASSOCIATED(ANGTAB)) THEN
      SCALAR_INT_BUF(251) = 1
      CALL fstarpu_matrix_data_register(handles(251), 0, C_LOC(ANGTAB(LBOUND(ANGTAB,1),LBOUND(ANGTAB,2))), SIZE(ANGTAB,1), SIZE(ANGTAB,1), SIZE(ANGTAB,2), C_SIZEOF(ANGTAB(LBOUND(ANGTAB,1),LBOUND(ANGTAB,2))))
    ELSE
      SCALAR_INT_BUF(251) = 0
      CALL fstarpu_matrix_data_register(handles(251), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(ANGTAB)
    IF (ASSOCIATED(CENTAB)) THEN
      SCALAR_INT_BUF(252) = 1
      CALL fstarpu_matrix_data_register(handles(252), 0, C_LOC(CENTAB(LBOUND(CENTAB,1),LBOUND(CENTAB,2))), SIZE(CENTAB,1), SIZE(CENTAB,1), SIZE(CENTAB,2), C_SIZEOF(CENTAB(LBOUND(CENTAB,1),LBOUND(CENTAB,2))))
    ELSE
      SCALAR_INT_BUF(252) = 0
      CALL fstarpu_matrix_data_register(handles(252), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(CENTAB)
    IF (ASSOCIATED(ELETAB)) THEN
      SCALAR_INT_BUF(253) = 1
      CALL fstarpu_matrix_data_register(handles(253), 0, C_LOC(ELETAB(LBOUND(ELETAB,1),LBOUND(ELETAB,2))), SIZE(ELETAB,1), SIZE(ELETAB,1), SIZE(ELETAB,2), C_SIZEOF(ELETAB(LBOUND(ELETAB,1),LBOUND(ELETAB,2))))
    ELSE
      SCALAR_INT_BUF(253) = 0
      CALL fstarpu_matrix_data_register(handles(253), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(ELETAB)
    IF (ASSOCIATED(DG_ANG)) THEN
      SCALAR_INT_BUF(254) = 1
      CALL fstarpu_vector_data_register(handles(254), 0, C_LOC(DG_ANG(LBOUND(DG_ANG,1))), SIZE(DG_ANG,1), C_SIZEOF(DG_ANG(LBOUND(DG_ANG,1))))
    ELSE
      SCALAR_INT_BUF(254) = 0
      CALL fstarpu_vector_data_register(handles(254), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(DG_ANG)
    IF (ASSOCIATED(DP_DG)) THEN
      SCALAR_INT_BUF(255) = 1
      CALL fstarpu_vector_data_register(handles(255), 0, C_LOC(DP_DG(LBOUND(DP_DG,1))), SIZE(DP_DG,1), C_SIZEOF(DP_DG(LBOUND(DP_DG,1))))
    ELSE
      SCALAR_INT_BUF(255) = 0
      CALL fstarpu_vector_data_register(handles(255), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(DP_DG)
    IF (ASSOCIATED(EL_COUNT)) THEN
      SCALAR_INT_BUF(256) = 1
      CALL fstarpu_vector_data_register(handles(256), 0, C_LOC(EL_COUNT(LBOUND(EL_COUNT,1))), SIZE(EL_COUNT,1), C_SIZEOF(EL_COUNT(LBOUND(EL_COUNT,1))))
    ELSE
      SCALAR_INT_BUF(256) = 0
      CALL fstarpu_vector_data_register(handles(256), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(EL_COUNT)
    IF (ASSOCIATED(NNOEL)) THEN
      SCALAR_INT_BUF(257) = 1
      CALL fstarpu_matrix_data_register(handles(257), 0, C_LOC(NNOEL(LBOUND(NNOEL,1),LBOUND(NNOEL,2))), SIZE(NNOEL,1), SIZE(NNOEL,1), SIZE(NNOEL,2), C_SIZEOF(NNOEL(LBOUND(NNOEL,1),LBOUND(NNOEL,2))))
    ELSE
      SCALAR_INT_BUF(257) = 0
      CALL fstarpu_matrix_data_register(handles(257), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(NNOEL)
    IF (ASSOCIATED(NNDEL)) THEN
      SCALAR_INT_BUF(258) = 1
      CALL fstarpu_vector_data_register(handles(258), 0, C_LOC(NNDEL(LBOUND(NNDEL,1))), SIZE(NNDEL,1), C_SIZEOF(NNDEL(LBOUND(NNDEL,1))))
    ELSE
      SCALAR_INT_BUF(258) = 0
      CALL fstarpu_vector_data_register(handles(258), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(NNDEL)
    IF (ASSOCIATED(NDEL)) THEN
      SCALAR_INT_BUF(259) = 1
      CALL fstarpu_matrix_data_register(handles(259), 0, C_LOC(NDEL(LBOUND(NDEL,1),LBOUND(NDEL,2))), SIZE(NDEL,1), SIZE(NDEL,1), SIZE(NDEL,2), C_SIZEOF(NDEL(LBOUND(NDEL,1),LBOUND(NDEL,2))))
    ELSE
      SCALAR_INT_BUF(259) = 0
      CALL fstarpu_matrix_data_register(handles(259), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(NDEL)
    IF (ASSOCIATED(FX_MID)) THEN
      SCALAR_INT_BUF(260) = 1
      CALL fstarpu_matrix_data_register(handles(260), 0, C_LOC(FX_MID(LBOUND(FX_MID,1),LBOUND(FX_MID,2))), SIZE(FX_MID,1), SIZE(FX_MID,1), SIZE(FX_MID,2), C_SIZEOF(FX_MID(LBOUND(FX_MID,1),LBOUND(FX_MID,2))))
    ELSE
      SCALAR_INT_BUF(260) = 0
      CALL fstarpu_matrix_data_register(handles(260), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(FX_MID)
    IF (ASSOCIATED(GX_MID)) THEN
      SCALAR_INT_BUF(261) = 1
      CALL fstarpu_matrix_data_register(handles(261), 0, C_LOC(GX_MID(LBOUND(GX_MID,1),LBOUND(GX_MID,2))), SIZE(GX_MID,1), SIZE(GX_MID,1), SIZE(GX_MID,2), C_SIZEOF(GX_MID(LBOUND(GX_MID,1),LBOUND(GX_MID,2))))
    ELSE
      SCALAR_INT_BUF(261) = 0
      CALL fstarpu_matrix_data_register(handles(261), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(GX_MID)
    IF (ASSOCIATED(HX_MID)) THEN
      SCALAR_INT_BUF(262) = 1
      CALL fstarpu_matrix_data_register(handles(262), 0, C_LOC(HX_MID(LBOUND(HX_MID,1),LBOUND(HX_MID,2))), SIZE(HX_MID,1), SIZE(HX_MID,1), SIZE(HX_MID,2), C_SIZEOF(HX_MID(LBOUND(HX_MID,1),LBOUND(HX_MID,2))))
    ELSE
      SCALAR_INT_BUF(262) = 0
      CALL fstarpu_matrix_data_register(handles(262), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(HX_MID)
    IF (ASSOCIATED(FY_MID)) THEN
      SCALAR_INT_BUF(263) = 1
      CALL fstarpu_matrix_data_register(handles(263), 0, C_LOC(FY_MID(LBOUND(FY_MID,1),LBOUND(FY_MID,2))), SIZE(FY_MID,1), SIZE(FY_MID,1), SIZE(FY_MID,2), C_SIZEOF(FY_MID(LBOUND(FY_MID,1),LBOUND(FY_MID,2))))
    ELSE
      SCALAR_INT_BUF(263) = 0
      CALL fstarpu_matrix_data_register(handles(263), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(FY_MID)
    IF (ASSOCIATED(GY_MID)) THEN
      SCALAR_INT_BUF(264) = 1
      CALL fstarpu_matrix_data_register(handles(264), 0, C_LOC(GY_MID(LBOUND(GY_MID,1),LBOUND(GY_MID,2))), SIZE(GY_MID,1), SIZE(GY_MID,1), SIZE(GY_MID,2), C_SIZEOF(GY_MID(LBOUND(GY_MID,1),LBOUND(GY_MID,2))))
    ELSE
      SCALAR_INT_BUF(264) = 0
      CALL fstarpu_matrix_data_register(handles(264), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(GY_MID)
    IF (ASSOCIATED(HY_MID)) THEN
      SCALAR_INT_BUF(265) = 1
      CALL fstarpu_matrix_data_register(handles(265), 0, C_LOC(HY_MID(LBOUND(HY_MID,1),LBOUND(HY_MID,2))), SIZE(HY_MID,1), SIZE(HY_MID,1), SIZE(HY_MID,2), C_SIZEOF(HY_MID(LBOUND(HY_MID,1),LBOUND(HY_MID,2))))
    ELSE
      SCALAR_INT_BUF(265) = 0
      CALL fstarpu_matrix_data_register(handles(265), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(HY_MID)
    IF (ASSOCIATED(ZE_C)) THEN
      SCALAR_INT_BUF(266) = 1
      CALL fstarpu_vector_data_register(handles(266), 0, C_LOC(ZE_C(LBOUND(ZE_C,1))), SIZE(ZE_C,1), C_SIZEOF(ZE_C(LBOUND(ZE_C,1))))
    ELSE
      SCALAR_INT_BUF(266) = 0
      CALL fstarpu_vector_data_register(handles(266), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(ZE_C)
    IF (ASSOCIATED(QX_C)) THEN
      SCALAR_INT_BUF(267) = 1
      CALL fstarpu_vector_data_register(handles(267), 0, C_LOC(QX_C(LBOUND(QX_C,1))), SIZE(QX_C,1), C_SIZEOF(QX_C(LBOUND(QX_C,1))))
    ELSE
      SCALAR_INT_BUF(267) = 0
      CALL fstarpu_vector_data_register(handles(267), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(QX_C)
    IF (ASSOCIATED(QY_C)) THEN
      SCALAR_INT_BUF(268) = 1
      CALL fstarpu_vector_data_register(handles(268), 0, C_LOC(QY_C(LBOUND(QY_C,1))), SIZE(QY_C,1), C_SIZEOF(QY_C(LBOUND(QY_C,1))))
    ELSE
      SCALAR_INT_BUF(268) = 0
      CALL fstarpu_vector_data_register(handles(268), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(QY_C)
    IF (ASSOCIATED(dynP_DG)) THEN
      SCALAR_INT_BUF(269) = 1
      CALL fstarpu_vector_data_register(handles(269), 0, C_LOC(dynP_DG(LBOUND(dynP_DG,1))), SIZE(dynP_DG,1), C_SIZEOF(dynP_DG(LBOUND(dynP_DG,1))))
    ELSE
      SCALAR_INT_BUF(269) = 0
      CALL fstarpu_vector_data_register(handles(269), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(dynP_DG)
    IF (ASSOCIATED(iota2_DG)) THEN
      SCALAR_INT_BUF(270) = 1
      CALL fstarpu_vector_data_register(handles(270), 0, C_LOC(iota2_DG(LBOUND(iota2_DG,1))), SIZE(iota2_DG,1), C_SIZEOF(iota2_DG(LBOUND(iota2_DG,1))))
    ELSE
      SCALAR_INT_BUF(270) = 0
      CALL fstarpu_vector_data_register(handles(270), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(iota2_DG)
    IF (ASSOCIATED(iota_DG)) THEN
      SCALAR_INT_BUF(271) = 1
      CALL fstarpu_vector_data_register(handles(271), 0, C_LOC(iota_DG(LBOUND(iota_DG,1))), SIZE(iota_DG,1), C_SIZEOF(iota_DG(LBOUND(iota_DG,1))))
    ELSE
      SCALAR_INT_BUF(271) = 0
      CALL fstarpu_vector_data_register(handles(271), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(iota_DG)
    IF (ASSOCIATED(iotaa_DG)) THEN
      SCALAR_INT_BUF(272) = 1
      CALL fstarpu_vector_data_register(handles(272), 0, C_LOC(iotaa_DG(LBOUND(iotaa_DG,1))), SIZE(iotaa_DG,1), C_SIZEOF(iotaa_DG(LBOUND(iotaa_DG,1))))
    ELSE
      SCALAR_INT_BUF(272) = 0
      CALL fstarpu_vector_data_register(handles(272), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(iotaa_DG)
    IF (ASSOCIATED(bed_DG)) THEN
      SCALAR_INT_BUF(273) = 1
      CALL fstarpu_matrix_data_register(handles(273), 0, C_LOC(bed_DG(LBOUND(bed_DG,1),LBOUND(bed_DG,2))), SIZE(bed_DG,1), SIZE(bed_DG,1), SIZE(bed_DG,2), C_SIZEOF(bed_DG(LBOUND(bed_DG,1),LBOUND(bed_DG,2))))
    ELSE
      SCALAR_INT_BUF(273) = 0
      CALL fstarpu_matrix_data_register(handles(273), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(bed_DG)
    IF (ASSOCIATED(bed_N_int)) THEN
      SCALAR_INT_BUF(274) = 1
      CALL fstarpu_vector_data_register(handles(274), 0, C_LOC(bed_N_int(LBOUND(bed_N_int,1))), SIZE(bed_N_int,1), C_SIZEOF(bed_N_int(LBOUND(bed_N_int,1))))
    ELSE
      SCALAR_INT_BUF(274) = 0
      CALL fstarpu_vector_data_register(handles(274), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(bed_N_int)
    IF (ASSOCIATED(bed_N_ext)) THEN
      SCALAR_INT_BUF(275) = 1
      CALL fstarpu_vector_data_register(handles(275), 0, C_LOC(bed_N_ext(LBOUND(bed_N_ext,1))), SIZE(bed_N_ext,1), C_SIZEOF(bed_N_ext(LBOUND(bed_N_ext,1))))
    ELSE
      SCALAR_INT_BUF(275) = 0
      CALL fstarpu_vector_data_register(handles(275), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(bed_N_ext)
    IF (ASSOCIATED(IDUMY)) THEN
      SCALAR_INT_BUF(276) = 1
      CALL fstarpu_vector_data_register(handles(276), 0, C_LOC(IDUMY(LBOUND(IDUMY,1))), SIZE(IDUMY,1), C_SIZEOF(IDUMY(LBOUND(IDUMY,1))))
    ELSE
      SCALAR_INT_BUF(276) = 0
      CALL fstarpu_vector_data_register(handles(276), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(IDUMY)
    IF (ASSOCIATED(DUMY1)) THEN
      SCALAR_INT_BUF(277) = 1
      CALL fstarpu_vector_data_register(handles(277), 0, C_LOC(DUMY1(LBOUND(DUMY1,1))), SIZE(DUMY1,1), C_SIZEOF(DUMY1(LBOUND(DUMY1,1))))
    ELSE
      SCALAR_INT_BUF(277) = 0
      CALL fstarpu_vector_data_register(handles(277), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(DUMY1)
    IF (ASSOCIATED(DUMY2)) THEN
      SCALAR_INT_BUF(278) = 1
      CALL fstarpu_vector_data_register(handles(278), 0, C_LOC(DUMY2(LBOUND(DUMY2,1))), SIZE(DUMY2,1), C_SIZEOF(DUMY2(LBOUND(DUMY2,1))))
    ELSE
      SCALAR_INT_BUF(278) = 0
      CALL fstarpu_vector_data_register(handles(278), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(DUMY2)
    IF (ASSOCIATED(DGDUMY1)) THEN
      SCALAR_INT_BUF(279) = 1
      CALL fstarpu_block_data_register(handles(279), 0, C_LOC(DGDUMY1(LBOUND(DGDUMY1,1),LBOUND(DGDUMY1,2),LBOUND(DGDUMY1,3))), SIZE(DGDUMY1,1), SIZE(DGDUMY1,1)*SIZE(DGDUMY1,2), SIZE(DGDUMY1,1), SIZE(DGDUMY1,2), SIZE(DGDUMY1,3), C_SIZEOF(DGDUMY1(LBOUND(DGDUMY1,1),LBOUND(DGDUMY1,2),LBOUND(DGDUMY1,3))))
    ELSE
      SCALAR_INT_BUF(279) = 0
      CALL fstarpu_block_data_register(handles(279), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(DGDUMY1)
    IF (ASSOCIATED(DGDUMY2)) THEN
      SCALAR_INT_BUF(280) = 1
      CALL fstarpu_block_data_register(handles(280), 0, C_LOC(DGDUMY2(LBOUND(DGDUMY2,1),LBOUND(DGDUMY2,2),LBOUND(DGDUMY2,3))), SIZE(DGDUMY2,1), SIZE(DGDUMY2,1)*SIZE(DGDUMY2,2), SIZE(DGDUMY2,1), SIZE(DGDUMY2,2), SIZE(DGDUMY2,3), C_SIZEOF(DGDUMY2(LBOUND(DGDUMY2,1),LBOUND(DGDUMY2,2),LBOUND(DGDUMY2,3))))
    ELSE
      SCALAR_INT_BUF(280) = 0
      CALL fstarpu_block_data_register(handles(280), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(DGDUMY2)
    IF (ASSOCIATED(pdg_el)) THEN
      SCALAR_INT_BUF(281) = 1
      CALL fstarpu_vector_data_register(handles(281), 0, C_LOC(pdg_el(LBOUND(pdg_el,1))), SIZE(pdg_el,1), C_SIZEOF(pdg_el(LBOUND(pdg_el,1))))
    ELSE
      SCALAR_INT_BUF(281) = 0
      CALL fstarpu_vector_data_register(handles(281), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(pdg_el)
    IF (ASSOCIATED(ETAS)) THEN
      SCALAR_INT_BUF(282) = 1
      CALL fstarpu_vector_data_register(handles(282), 0, C_LOC(ETAS(LBOUND(ETAS,1))), SIZE(ETAS,1), C_SIZEOF(ETAS(LBOUND(ETAS,1))))
    ELSE
      SCALAR_INT_BUF(282) = 0
      CALL fstarpu_vector_data_register(handles(282), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(ETAS)
    IF (ASSOCIATED(ETA1)) THEN
      SCALAR_INT_BUF(283) = 1
      CALL fstarpu_vector_data_register(handles(283), 0, C_LOC(ETA1(LBOUND(ETA1,1))), SIZE(ETA1,1), C_SIZEOF(ETA1(LBOUND(ETA1,1))))
    ELSE
      SCALAR_INT_BUF(283) = 0
      CALL fstarpu_vector_data_register(handles(283), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(ETA1)
    IF (ASSOCIATED(ETA2)) THEN
      SCALAR_INT_BUF(284) = 1
      CALL fstarpu_vector_data_register(handles(284), 0, C_LOC(ETA2(LBOUND(ETA2,1))), SIZE(ETA2,1), C_SIZEOF(ETA2(LBOUND(ETA2,1))))
    ELSE
      SCALAR_INT_BUF(284) = 0
      CALL fstarpu_vector_data_register(handles(284), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(ETA2)
    IF (ASSOCIATED(ETAMAX)) THEN
      SCALAR_INT_BUF(285) = 1
      CALL fstarpu_vector_data_register(handles(285), 0, C_LOC(ETAMAX(LBOUND(ETAMAX,1))), SIZE(ETAMAX,1), C_SIZEOF(ETAMAX(LBOUND(ETAMAX,1))))
    ELSE
      SCALAR_INT_BUF(285) = 0
      CALL fstarpu_vector_data_register(handles(285), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(ETAMAX)
    IF (ASSOCIATED(entrop)) THEN
      SCALAR_INT_BUF(286) = 1
      CALL fstarpu_matrix_data_register(handles(286), 0, C_LOC(entrop(LBOUND(entrop,1),LBOUND(entrop,2))), SIZE(entrop,1), SIZE(entrop,1), SIZE(entrop,2), C_SIZEOF(entrop(LBOUND(entrop,1),LBOUND(entrop,2))))
    ELSE
      SCALAR_INT_BUF(286) = 0
      CALL fstarpu_matrix_data_register(handles(286), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(entrop)
    IF (ASSOCIATED(tracer)) THEN
      SCALAR_INT_BUF(287) = 1
      CALL fstarpu_vector_data_register(handles(287), 0, C_LOC(tracer(LBOUND(tracer,1))), SIZE(tracer,1), C_SIZEOF(tracer(LBOUND(tracer,1))))
    ELSE
      SCALAR_INT_BUF(287) = 0
      CALL fstarpu_vector_data_register(handles(287), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(tracer)
    IF (ASSOCIATED(tracer2)) THEN
      SCALAR_INT_BUF(288) = 1
      CALL fstarpu_vector_data_register(handles(288), 0, C_LOC(tracer2(LBOUND(tracer2,1))), SIZE(tracer2,1), C_SIZEOF(tracer2(LBOUND(tracer2,1))))
    ELSE
      SCALAR_INT_BUF(288) = 0
      CALL fstarpu_vector_data_register(handles(288), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(tracer2)
    IF (ASSOCIATED(MassMax)) THEN
      SCALAR_INT_BUF(289) = 1
      CALL fstarpu_vector_data_register(handles(289), 0, C_LOC(MassMax(LBOUND(MassMax,1))), SIZE(MassMax,1), C_SIZEOF(MassMax(LBOUND(MassMax,1))))
    ELSE
      SCALAR_INT_BUF(289) = 0
      CALL fstarpu_vector_data_register(handles(289), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(MassMax)
    IF (ASSOCIATED(bed_int)) THEN
      SCALAR_INT_BUF(290) = 1
      CALL fstarpu_matrix_data_register(handles(290), 0, C_LOC(bed_int(LBOUND(bed_int,1),LBOUND(bed_int,2))), SIZE(bed_int,1), SIZE(bed_int,1), SIZE(bed_int,2), C_SIZEOF(bed_int(LBOUND(bed_int,1),LBOUND(bed_int,2))))
    ELSE
      SCALAR_INT_BUF(290) = 0
      CALL fstarpu_matrix_data_register(handles(290), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(bed_int)
    IF (ASSOCIATED(UU1)) THEN
      SCALAR_INT_BUF(291) = 1
      CALL fstarpu_vector_data_register(handles(291), 0, C_LOC(UU1(LBOUND(UU1,1))), SIZE(UU1,1), C_SIZEOF(UU1(LBOUND(UU1,1))))
    ELSE
      SCALAR_INT_BUF(291) = 0
      CALL fstarpu_vector_data_register(handles(291), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(UU1)
    IF (ASSOCIATED(UU2)) THEN
      SCALAR_INT_BUF(292) = 1
      CALL fstarpu_vector_data_register(handles(292), 0, C_LOC(UU2(LBOUND(UU2,1))), SIZE(UU2,1), C_SIZEOF(UU2(LBOUND(UU2,1))))
    ELSE
      SCALAR_INT_BUF(292) = 0
      CALL fstarpu_vector_data_register(handles(292), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(UU2)
    IF (ASSOCIATED(VV1)) THEN
      SCALAR_INT_BUF(293) = 1
      CALL fstarpu_vector_data_register(handles(293), 0, C_LOC(VV1(LBOUND(VV1,1))), SIZE(VV1,1), C_SIZEOF(VV1(LBOUND(VV1,1))))
    ELSE
      SCALAR_INT_BUF(293) = 0
      CALL fstarpu_vector_data_register(handles(293), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(VV1)
    IF (ASSOCIATED(VV2)) THEN
      SCALAR_INT_BUF(294) = 1
      CALL fstarpu_vector_data_register(handles(294), 0, C_LOC(VV2(LBOUND(VV2,1))), SIZE(VV2,1), C_SIZEOF(VV2(LBOUND(VV2,1))))
    ELSE
      SCALAR_INT_BUF(294) = 0
      CALL fstarpu_vector_data_register(handles(294), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(VV2)
    IF (ASSOCIATED(dyn_P)) THEN
      SCALAR_INT_BUF(295) = 1
      CALL fstarpu_vector_data_register(handles(295), 0, C_LOC(dyn_P(LBOUND(dyn_P,1))), SIZE(dyn_P,1), C_SIZEOF(dyn_P(LBOUND(dyn_P,1))))
    ELSE
      SCALAR_INT_BUF(295) = 0
      CALL fstarpu_vector_data_register(handles(295), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(dyn_P)
    IF (ASSOCIATED(DP)) THEN
      SCALAR_INT_BUF(296) = 1
      CALL fstarpu_vector_data_register(handles(296), 0, C_LOC(DP(LBOUND(DP,1))), SIZE(DP,1), C_SIZEOF(DP(LBOUND(DP,1))))
    ELSE
      SCALAR_INT_BUF(296) = 0
      CALL fstarpu_vector_data_register(handles(296), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(DP)
    IF (ASSOCIATED(DP0)) THEN
      SCALAR_INT_BUF(297) = 1
      CALL fstarpu_vector_data_register(handles(297), 0, C_LOC(DP0(LBOUND(DP0,1))), SIZE(DP0,1), C_SIZEOF(DP0(LBOUND(DP0,1))))
    ELSE
      SCALAR_INT_BUF(297) = 0
      CALL fstarpu_vector_data_register(handles(297), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(DP0)
    IF (ASSOCIATED(DPe)) THEN
      SCALAR_INT_BUF(298) = 1
      CALL fstarpu_vector_data_register(handles(298), 0, C_LOC(DPe(LBOUND(DPe,1))), SIZE(DPe,1), C_SIZEOF(DPe(LBOUND(DPe,1))))
    ELSE
      SCALAR_INT_BUF(298) = 0
      CALL fstarpu_vector_data_register(handles(298), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(DPe)
    IF (ASSOCIATED(SFAC)) THEN
      SCALAR_INT_BUF(299) = 1
      CALL fstarpu_vector_data_register(handles(299), 0, C_LOC(SFAC(LBOUND(SFAC,1))), SIZE(SFAC,1), C_SIZEOF(SFAC(LBOUND(SFAC,1))))
    ELSE
      SCALAR_INT_BUF(299) = 0
      CALL fstarpu_vector_data_register(handles(299), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(SFAC)
    IF (ASSOCIATED(QU)) THEN
      SCALAR_INT_BUF(300) = 1
      CALL fstarpu_vector_data_register(handles(300), 0, C_LOC(QU(LBOUND(QU,1))), SIZE(QU,1), C_SIZEOF(QU(LBOUND(QU,1))))
    ELSE
      SCALAR_INT_BUF(300) = 0
      CALL fstarpu_vector_data_register(handles(300), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(QU)
    IF (ASSOCIATED(QV)) THEN
      SCALAR_INT_BUF(301) = 1
      CALL fstarpu_vector_data_register(handles(301), 0, C_LOC(QV(LBOUND(QV,1))), SIZE(QV,1), C_SIZEOF(QV(LBOUND(QV,1))))
    ELSE
      SCALAR_INT_BUF(301) = 0
      CALL fstarpu_vector_data_register(handles(301), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(QV)
    IF (ASSOCIATED(QW)) THEN
      SCALAR_INT_BUF(302) = 1
      CALL fstarpu_vector_data_register(handles(302), 0, C_LOC(QW(LBOUND(QW,1))), SIZE(QW,1), C_SIZEOF(QW(LBOUND(QW,1))))
    ELSE
      SCALAR_INT_BUF(302) = 0
      CALL fstarpu_vector_data_register(handles(302), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(QW)
    IF (ASSOCIATED(CORIF)) THEN
      SCALAR_INT_BUF(303) = 1
      CALL fstarpu_vector_data_register(handles(303), 0, C_LOC(CORIF(LBOUND(CORIF,1))), SIZE(CORIF,1), C_SIZEOF(CORIF(LBOUND(CORIF,1))))
    ELSE
      SCALAR_INT_BUF(303) = 0
      CALL fstarpu_vector_data_register(handles(303), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(CORIF)
    IF (ASSOCIATED(TPK)) THEN
      SCALAR_INT_BUF(304) = 1
      CALL fstarpu_vector_data_register(handles(304), 0, C_LOC(TPK(LBOUND(TPK,1))), SIZE(TPK,1), C_SIZEOF(TPK(LBOUND(TPK,1))))
    ELSE
      SCALAR_INT_BUF(304) = 0
      CALL fstarpu_vector_data_register(handles(304), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(TPK)
    IF (ASSOCIATED(FFT)) THEN
      SCALAR_INT_BUF(305) = 1
      CALL fstarpu_vector_data_register(handles(305), 0, C_LOC(FFT(LBOUND(FFT,1))), SIZE(FFT,1), C_SIZEOF(FFT(LBOUND(FFT,1))))
    ELSE
      SCALAR_INT_BUF(305) = 0
      CALL fstarpu_vector_data_register(handles(305), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(FFT)
    IF (ASSOCIATED(FACET)) THEN
      SCALAR_INT_BUF(306) = 1
      CALL fstarpu_vector_data_register(handles(306), 0, C_LOC(FACET(LBOUND(FACET,1))), SIZE(FACET,1), C_SIZEOF(FACET(LBOUND(FACET,1))))
    ELSE
      SCALAR_INT_BUF(306) = 0
      CALL fstarpu_vector_data_register(handles(306), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(FACET)
    IF (ASSOCIATED(ETRF)) THEN
      SCALAR_INT_BUF(307) = 1
      CALL fstarpu_vector_data_register(handles(307), 0, C_LOC(ETRF(LBOUND(ETRF,1))), SIZE(ETRF,1), C_SIZEOF(ETRF(LBOUND(ETRF,1))))
    ELSE
      SCALAR_INT_BUF(307) = 0
      CALL fstarpu_vector_data_register(handles(307), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(ETRF)
    IF (ASSOCIATED(ESBIN1)) THEN
      SCALAR_INT_BUF(308) = 1
      CALL fstarpu_vector_data_register(handles(308), 0, C_LOC(ESBIN1(LBOUND(ESBIN1,1))), SIZE(ESBIN1,1), C_SIZEOF(ESBIN1(LBOUND(ESBIN1,1))))
    ELSE
      SCALAR_INT_BUF(308) = 0
      CALL fstarpu_vector_data_register(handles(308), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(ESBIN1)
    IF (ASSOCIATED(ESBIN2)) THEN
      SCALAR_INT_BUF(309) = 1
      CALL fstarpu_vector_data_register(handles(309), 0, C_LOC(ESBIN2(LBOUND(ESBIN2,1))), SIZE(ESBIN2,1), C_SIZEOF(ESBIN2(LBOUND(ESBIN2,1))))
    ELSE
      SCALAR_INT_BUF(309) = 0
      CALL fstarpu_vector_data_register(handles(309), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(ESBIN2)
    IF (ASSOCIATED(QTEMA)) THEN
      SCALAR_INT_BUF(310) = 1
      CALL fstarpu_matrix_data_register(handles(310), 0, C_LOC(QTEMA(LBOUND(QTEMA,1),LBOUND(QTEMA,2))), SIZE(QTEMA,1), SIZE(QTEMA,1), SIZE(QTEMA,2), C_SIZEOF(QTEMA(LBOUND(QTEMA,1),LBOUND(QTEMA,2))))
    ELSE
      SCALAR_INT_BUF(310) = 0
      CALL fstarpu_matrix_data_register(handles(310), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(QTEMA)
    IF (ASSOCIATED(QTEMB)) THEN
      SCALAR_INT_BUF(311) = 1
      CALL fstarpu_matrix_data_register(handles(311), 0, C_LOC(QTEMB(LBOUND(QTEMB,1),LBOUND(QTEMB,2))), SIZE(QTEMB,1), SIZE(QTEMB,1), SIZE(QTEMB,2), C_SIZEOF(QTEMB(LBOUND(QTEMB,1),LBOUND(QTEMB,2))))
    ELSE
      SCALAR_INT_BUF(311) = 0
      CALL fstarpu_matrix_data_register(handles(311), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(QTEMB)
    IF (ASSOCIATED(QN2)) THEN
      SCALAR_INT_BUF(312) = 1
      CALL fstarpu_vector_data_register(handles(312), 0, C_LOC(QN2(LBOUND(QN2,1))), SIZE(QN2,1), C_SIZEOF(QN2(LBOUND(QN2,1))))
    ELSE
      SCALAR_INT_BUF(312) = 0
      CALL fstarpu_vector_data_register(handles(312), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(QN2)
    IF (ASSOCIATED(BNDLEN2O3)) THEN
      SCALAR_INT_BUF(313) = 1
      CALL fstarpu_vector_data_register(handles(313), 0, C_LOC(BNDLEN2O3(LBOUND(BNDLEN2O3,1))), SIZE(BNDLEN2O3,1), C_SIZEOF(BNDLEN2O3(LBOUND(BNDLEN2O3,1))))
    ELSE
      SCALAR_INT_BUF(313) = 0
      CALL fstarpu_vector_data_register(handles(313), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(BNDLEN2O3)
    IF (ASSOCIATED(CSII)) THEN
      SCALAR_INT_BUF(314) = 1
      CALL fstarpu_vector_data_register(handles(314), 0, C_LOC(CSII(LBOUND(CSII,1))), SIZE(CSII,1), C_SIZEOF(CSII(LBOUND(CSII,1))))
    ELSE
      SCALAR_INT_BUF(314) = 0
      CALL fstarpu_vector_data_register(handles(314), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(CSII)
    IF (ASSOCIATED(SIII)) THEN
      SCALAR_INT_BUF(315) = 1
      CALL fstarpu_vector_data_register(handles(315), 0, C_LOC(SIII(LBOUND(SIII,1))), SIZE(SIII,1), C_SIZEOF(SIII(LBOUND(SIII,1))))
    ELSE
      SCALAR_INT_BUF(315) = 0
      CALL fstarpu_vector_data_register(handles(315), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(SIII)
    IF (ASSOCIATED(QNAM)) THEN
      SCALAR_INT_BUF(316) = 1
      CALL fstarpu_matrix_data_register(handles(316), 0, C_LOC(QNAM(LBOUND(QNAM,1),LBOUND(QNAM,2))), SIZE(QNAM,1), SIZE(QNAM,1), SIZE(QNAM,2), C_SIZEOF(QNAM(LBOUND(QNAM,1),LBOUND(QNAM,2))))
    ELSE
      SCALAR_INT_BUF(316) = 0
      CALL fstarpu_matrix_data_register(handles(316), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(QNAM)
    IF (ASSOCIATED(QNPH)) THEN
      SCALAR_INT_BUF(317) = 1
      CALL fstarpu_matrix_data_register(handles(317), 0, C_LOC(QNPH(LBOUND(QNPH,1),LBOUND(QNPH,2))), SIZE(QNPH,1), SIZE(QNPH,1), SIZE(QNPH,2), C_SIZEOF(QNPH(LBOUND(QNPH,1),LBOUND(QNPH,2))))
    ELSE
      SCALAR_INT_BUF(317) = 0
      CALL fstarpu_matrix_data_register(handles(317), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(QNPH)
    IF (ASSOCIATED(QNIN1)) THEN
      SCALAR_INT_BUF(318) = 1
      CALL fstarpu_vector_data_register(handles(318), 0, C_LOC(QNIN1(LBOUND(QNIN1,1))), SIZE(QNIN1,1), C_SIZEOF(QNIN1(LBOUND(QNIN1,1))))
    ELSE
      SCALAR_INT_BUF(318) = 0
      CALL fstarpu_vector_data_register(handles(318), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(QNIN1)
    IF (ASSOCIATED(QNIN2)) THEN
      SCALAR_INT_BUF(319) = 1
      CALL fstarpu_vector_data_register(handles(319), 0, C_LOC(QNIN2(LBOUND(QNIN2,1))), SIZE(QNIN2,1), C_SIZEOF(QNIN2(LBOUND(QNIN2,1))))
    ELSE
      SCALAR_INT_BUF(319) = 0
      CALL fstarpu_vector_data_register(handles(319), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(QNIN2)
    IF (ASSOCIATED(CSI)) THEN
      SCALAR_INT_BUF(320) = 1
      CALL fstarpu_vector_data_register(handles(320), 0, C_LOC(CSI(LBOUND(CSI,1))), SIZE(CSI,1), C_SIZEOF(CSI(LBOUND(CSI,1))))
    ELSE
      SCALAR_INT_BUF(320) = 0
      CALL fstarpu_vector_data_register(handles(320), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(CSI)
    IF (ASSOCIATED(SII)) THEN
      SCALAR_INT_BUF(321) = 1
      CALL fstarpu_vector_data_register(handles(321), 0, C_LOC(SII(LBOUND(SII,1))), SIZE(SII,1), C_SIZEOF(SII(LBOUND(SII,1))))
    ELSE
      SCALAR_INT_BUF(321) = 0
      CALL fstarpu_vector_data_register(handles(321), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(SII)
    IF (ASSOCIATED(ET00)) THEN
      SCALAR_INT_BUF(322) = 1
      CALL fstarpu_vector_data_register(handles(322), 0, C_LOC(ET00(LBOUND(ET00,1))), SIZE(ET00,1), C_SIZEOF(ET00(LBOUND(ET00,1))))
    ELSE
      SCALAR_INT_BUF(322) = 0
      CALL fstarpu_vector_data_register(handles(322), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(ET00)
    IF (ASSOCIATED(BT00)) THEN
      SCALAR_INT_BUF(323) = 1
      CALL fstarpu_vector_data_register(handles(323), 0, C_LOC(BT00(LBOUND(BT00,1))), SIZE(BT00,1), C_SIZEOF(BT00(LBOUND(BT00,1))))
    ELSE
      SCALAR_INT_BUF(323) = 0
      CALL fstarpu_vector_data_register(handles(323), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(BT00)
    IF (ASSOCIATED(STAIE1)) THEN
      SCALAR_INT_BUF(324) = 1
      CALL fstarpu_vector_data_register(handles(324), 0, C_LOC(STAIE1(LBOUND(STAIE1,1))), SIZE(STAIE1,1), C_SIZEOF(STAIE1(LBOUND(STAIE1,1))))
    ELSE
      SCALAR_INT_BUF(324) = 0
      CALL fstarpu_vector_data_register(handles(324), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(STAIE1)
    IF (ASSOCIATED(STAIE2)) THEN
      SCALAR_INT_BUF(325) = 1
      CALL fstarpu_vector_data_register(handles(325), 0, C_LOC(STAIE2(LBOUND(STAIE2,1))), SIZE(STAIE2,1), C_SIZEOF(STAIE2(LBOUND(STAIE2,1))))
    ELSE
      SCALAR_INT_BUF(325) = 0
      CALL fstarpu_vector_data_register(handles(325), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(STAIE2)
    IF (ASSOCIATED(STAIE3)) THEN
      SCALAR_INT_BUF(326) = 1
      CALL fstarpu_vector_data_register(handles(326), 0, C_LOC(STAIE3(LBOUND(STAIE3,1))), SIZE(STAIE3,1), C_SIZEOF(STAIE3(LBOUND(STAIE3,1))))
    ELSE
      SCALAR_INT_BUF(326) = 0
      CALL fstarpu_vector_data_register(handles(326), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(STAIE3)
    IF (ASSOCIATED(XEV)) THEN
      SCALAR_INT_BUF(327) = 1
      CALL fstarpu_vector_data_register(handles(327), 0, C_LOC(XEV(LBOUND(XEV,1))), SIZE(XEV,1), C_SIZEOF(XEV(LBOUND(XEV,1))))
    ELSE
      SCALAR_INT_BUF(327) = 0
      CALL fstarpu_vector_data_register(handles(327), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(XEV)
    IF (ASSOCIATED(YEV)) THEN
      SCALAR_INT_BUF(328) = 1
      CALL fstarpu_vector_data_register(handles(328), 0, C_LOC(YEV(LBOUND(YEV,1))), SIZE(YEV,1), C_SIZEOF(YEV(LBOUND(YEV,1))))
    ELSE
      SCALAR_INT_BUF(328) = 0
      CALL fstarpu_vector_data_register(handles(328), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(YEV)
    IF (ASSOCIATED(SLEV)) THEN
      SCALAR_INT_BUF(329) = 1
      CALL fstarpu_vector_data_register(handles(329), 0, C_LOC(SLEV(LBOUND(SLEV,1))), SIZE(SLEV,1), C_SIZEOF(SLEV(LBOUND(SLEV,1))))
    ELSE
      SCALAR_INT_BUF(329) = 0
      CALL fstarpu_vector_data_register(handles(329), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(SLEV)
    IF (ASSOCIATED(SFEV)) THEN
      SCALAR_INT_BUF(330) = 1
      CALL fstarpu_vector_data_register(handles(330), 0, C_LOC(SFEV(LBOUND(SFEV,1))), SIZE(SFEV,1), C_SIZEOF(SFEV(LBOUND(SFEV,1))))
    ELSE
      SCALAR_INT_BUF(330) = 0
      CALL fstarpu_vector_data_register(handles(330), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(SFEV)
    IF (ASSOCIATED(UU00)) THEN
      SCALAR_INT_BUF(331) = 1
      CALL fstarpu_vector_data_register(handles(331), 0, C_LOC(UU00(LBOUND(UU00,1))), SIZE(UU00,1), C_SIZEOF(UU00(LBOUND(UU00,1))))
    ELSE
      SCALAR_INT_BUF(331) = 0
      CALL fstarpu_vector_data_register(handles(331), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(UU00)
    IF (ASSOCIATED(VV00)) THEN
      SCALAR_INT_BUF(332) = 1
      CALL fstarpu_vector_data_register(handles(332), 0, C_LOC(VV00(LBOUND(VV00,1))), SIZE(VV00,1), C_SIZEOF(VV00(LBOUND(VV00,1))))
    ELSE
      SCALAR_INT_BUF(332) = 0
      CALL fstarpu_vector_data_register(handles(332), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(VV00)
    IF (ASSOCIATED(STAIV1)) THEN
      SCALAR_INT_BUF(333) = 1
      CALL fstarpu_vector_data_register(handles(333), 0, C_LOC(STAIV1(LBOUND(STAIV1,1))), SIZE(STAIV1,1), C_SIZEOF(STAIV1(LBOUND(STAIV1,1))))
    ELSE
      SCALAR_INT_BUF(333) = 0
      CALL fstarpu_vector_data_register(handles(333), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(STAIV1)
    IF (ASSOCIATED(STAIV2)) THEN
      SCALAR_INT_BUF(334) = 1
      CALL fstarpu_vector_data_register(handles(334), 0, C_LOC(STAIV2(LBOUND(STAIV2,1))), SIZE(STAIV2,1), C_SIZEOF(STAIV2(LBOUND(STAIV2,1))))
    ELSE
      SCALAR_INT_BUF(334) = 0
      CALL fstarpu_vector_data_register(handles(334), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(STAIV2)
    IF (ASSOCIATED(STAIV3)) THEN
      SCALAR_INT_BUF(335) = 1
      CALL fstarpu_vector_data_register(handles(335), 0, C_LOC(STAIV3(LBOUND(STAIV3,1))), SIZE(STAIV3,1), C_SIZEOF(STAIV3(LBOUND(STAIV3,1))))
    ELSE
      SCALAR_INT_BUF(335) = 0
      CALL fstarpu_vector_data_register(handles(335), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(STAIV3)
    IF (ASSOCIATED(XEC)) THEN
      SCALAR_INT_BUF(336) = 1
      CALL fstarpu_vector_data_register(handles(336), 0, C_LOC(XEC(LBOUND(XEC,1))), SIZE(XEC,1), C_SIZEOF(XEC(LBOUND(XEC,1))))
    ELSE
      SCALAR_INT_BUF(336) = 0
      CALL fstarpu_vector_data_register(handles(336), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(XEC)
    IF (ASSOCIATED(YEC)) THEN
      SCALAR_INT_BUF(337) = 1
      CALL fstarpu_vector_data_register(handles(337), 0, C_LOC(YEC(LBOUND(YEC,1))), SIZE(YEC,1), C_SIZEOF(YEC(LBOUND(YEC,1))))
    ELSE
      SCALAR_INT_BUF(337) = 0
      CALL fstarpu_vector_data_register(handles(337), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(YEC)
    IF (ASSOCIATED(SLEC)) THEN
      SCALAR_INT_BUF(338) = 1
      CALL fstarpu_vector_data_register(handles(338), 0, C_LOC(SLEC(LBOUND(SLEC,1))), SIZE(SLEC,1), C_SIZEOF(SLEC(LBOUND(SLEC,1))))
    ELSE
      SCALAR_INT_BUF(338) = 0
      CALL fstarpu_vector_data_register(handles(338), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(SLEC)
    IF (ASSOCIATED(SFEC)) THEN
      SCALAR_INT_BUF(339) = 1
      CALL fstarpu_vector_data_register(handles(339), 0, C_LOC(SFEC(LBOUND(SFEC,1))), SIZE(SFEC,1), C_SIZEOF(SFEC(LBOUND(SFEC,1))))
    ELSE
      SCALAR_INT_BUF(339) = 0
      CALL fstarpu_vector_data_register(handles(339), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(SFEC)
    IF (ASSOCIATED(CC00)) THEN
      SCALAR_INT_BUF(340) = 1
      CALL fstarpu_vector_data_register(handles(340), 0, C_LOC(CC00(LBOUND(CC00,1))), SIZE(CC00,1), C_SIZEOF(CC00(LBOUND(CC00,1))))
    ELSE
      SCALAR_INT_BUF(340) = 0
      CALL fstarpu_vector_data_register(handles(340), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(CC00)
    IF (ASSOCIATED(STAIC1)) THEN
      SCALAR_INT_BUF(341) = 1
      CALL fstarpu_vector_data_register(handles(341), 0, C_LOC(STAIC1(LBOUND(STAIC1,1))), SIZE(STAIC1,1), C_SIZEOF(STAIC1(LBOUND(STAIC1,1))))
    ELSE
      SCALAR_INT_BUF(341) = 0
      CALL fstarpu_vector_data_register(handles(341), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(STAIC1)
    IF (ASSOCIATED(STAIC2)) THEN
      SCALAR_INT_BUF(342) = 1
      CALL fstarpu_vector_data_register(handles(342), 0, C_LOC(STAIC2(LBOUND(STAIC2,1))), SIZE(STAIC2,1), C_SIZEOF(STAIC2(LBOUND(STAIC2,1))))
    ELSE
      SCALAR_INT_BUF(342) = 0
      CALL fstarpu_vector_data_register(handles(342), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(STAIC2)
    IF (ASSOCIATED(STAIC3)) THEN
      SCALAR_INT_BUF(343) = 1
      CALL fstarpu_vector_data_register(handles(343), 0, C_LOC(STAIC3(LBOUND(STAIC3,1))), SIZE(STAIC3,1), C_SIZEOF(STAIC3(LBOUND(STAIC3,1))))
    ELSE
      SCALAR_INT_BUF(343) = 0
      CALL fstarpu_vector_data_register(handles(343), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(STAIC3)
    IF (ASSOCIATED(XEM)) THEN
      SCALAR_INT_BUF(344) = 1
      CALL fstarpu_vector_data_register(handles(344), 0, C_LOC(XEM(LBOUND(XEM,1))), SIZE(XEM,1), C_SIZEOF(XEM(LBOUND(XEM,1))))
    ELSE
      SCALAR_INT_BUF(344) = 0
      CALL fstarpu_vector_data_register(handles(344), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(XEM)
    IF (ASSOCIATED(YEM)) THEN
      SCALAR_INT_BUF(345) = 1
      CALL fstarpu_vector_data_register(handles(345), 0, C_LOC(YEM(LBOUND(YEM,1))), SIZE(YEM,1), C_SIZEOF(YEM(LBOUND(YEM,1))))
    ELSE
      SCALAR_INT_BUF(345) = 0
      CALL fstarpu_vector_data_register(handles(345), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(YEM)
    IF (ASSOCIATED(SLEM)) THEN
      SCALAR_INT_BUF(346) = 1
      CALL fstarpu_vector_data_register(handles(346), 0, C_LOC(SLEM(LBOUND(SLEM,1))), SIZE(SLEM,1), C_SIZEOF(SLEM(LBOUND(SLEM,1))))
    ELSE
      SCALAR_INT_BUF(346) = 0
      CALL fstarpu_vector_data_register(handles(346), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(SLEM)
    IF (ASSOCIATED(SFEM)) THEN
      SCALAR_INT_BUF(347) = 1
      CALL fstarpu_vector_data_register(handles(347), 0, C_LOC(SFEM(LBOUND(SFEM,1))), SIZE(SFEM,1), C_SIZEOF(SFEM(LBOUND(SFEM,1))))
    ELSE
      SCALAR_INT_BUF(347) = 0
      CALL fstarpu_vector_data_register(handles(347), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(SFEM)
    IF (ASSOCIATED(RMU00)) THEN
      SCALAR_INT_BUF(348) = 1
      CALL fstarpu_vector_data_register(handles(348), 0, C_LOC(RMU00(LBOUND(RMU00,1))), SIZE(RMU00,1), C_SIZEOF(RMU00(LBOUND(RMU00,1))))
    ELSE
      SCALAR_INT_BUF(348) = 0
      CALL fstarpu_vector_data_register(handles(348), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(RMU00)
    IF (ASSOCIATED(RMV00)) THEN
      SCALAR_INT_BUF(349) = 1
      CALL fstarpu_vector_data_register(handles(349), 0, C_LOC(RMV00(LBOUND(RMV00,1))), SIZE(RMV00,1), C_SIZEOF(RMV00(LBOUND(RMV00,1))))
    ELSE
      SCALAR_INT_BUF(349) = 0
      CALL fstarpu_vector_data_register(handles(349), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(RMV00)
    IF (ASSOCIATED(RMP00)) THEN
      SCALAR_INT_BUF(350) = 1
      CALL fstarpu_vector_data_register(handles(350), 0, C_LOC(RMP00(LBOUND(RMP00,1))), SIZE(RMP00,1), C_SIZEOF(RMP00(LBOUND(RMP00,1))))
    ELSE
      SCALAR_INT_BUF(350) = 0
      CALL fstarpu_vector_data_register(handles(350), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(RMP00)
    IF (ASSOCIATED(STAIM1)) THEN
      SCALAR_INT_BUF(351) = 1
      CALL fstarpu_vector_data_register(handles(351), 0, C_LOC(STAIM1(LBOUND(STAIM1,1))), SIZE(STAIM1,1), C_SIZEOF(STAIM1(LBOUND(STAIM1,1))))
    ELSE
      SCALAR_INT_BUF(351) = 0
      CALL fstarpu_vector_data_register(handles(351), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(STAIM1)
    IF (ASSOCIATED(STAIM2)) THEN
      SCALAR_INT_BUF(352) = 1
      CALL fstarpu_vector_data_register(handles(352), 0, C_LOC(STAIM2(LBOUND(STAIM2,1))), SIZE(STAIM2,1), C_SIZEOF(STAIM2(LBOUND(STAIM2,1))))
    ELSE
      SCALAR_INT_BUF(352) = 0
      CALL fstarpu_vector_data_register(handles(352), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(STAIM2)
    IF (ASSOCIATED(STAIM3)) THEN
      SCALAR_INT_BUF(353) = 1
      CALL fstarpu_vector_data_register(handles(353), 0, C_LOC(STAIM3(LBOUND(STAIM3,1))), SIZE(STAIM3,1), C_SIZEOF(STAIM3(LBOUND(STAIM3,1))))
    ELSE
      SCALAR_INT_BUF(353) = 0
      CALL fstarpu_vector_data_register(handles(353), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(STAIM3)
    IF (ASSOCIATED(CH1)) THEN
      SCALAR_INT_BUF(354) = 1
      CALL fstarpu_vector_data_register(handles(354), 0, C_LOC(CH1(LBOUND(CH1,1))), SIZE(CH1,1), C_SIZEOF(CH1(LBOUND(CH1,1))))
    ELSE
      SCALAR_INT_BUF(354) = 0
      CALL fstarpu_vector_data_register(handles(354), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(CH1)
    IF (ASSOCIATED(QB)) THEN
      SCALAR_INT_BUF(355) = 1
      CALL fstarpu_vector_data_register(handles(355), 0, C_LOC(QB(LBOUND(QB,1))), SIZE(QB,1), C_SIZEOF(QB(LBOUND(QB,1))))
    ELSE
      SCALAR_INT_BUF(355) = 0
      CALL fstarpu_vector_data_register(handles(355), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(QB)
    IF (ASSOCIATED(QA)) THEN
      SCALAR_INT_BUF(356) = 1
      CALL fstarpu_vector_data_register(handles(356), 0, C_LOC(QA(LBOUND(QA,1))), SIZE(QA,1), C_SIZEOF(QA(LBOUND(QA,1))))
    ELSE
      SCALAR_INT_BUF(356) = 0
      CALL fstarpu_vector_data_register(handles(356), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(QA)
    IF (ASSOCIATED(SOURSIN)) THEN
      SCALAR_INT_BUF(357) = 1
      CALL fstarpu_vector_data_register(handles(357), 0, C_LOC(SOURSIN(LBOUND(SOURSIN,1))), SIZE(SOURSIN,1), C_SIZEOF(SOURSIN(LBOUND(SOURSIN,1))))
    ELSE
      SCALAR_INT_BUF(357) = 0
      CALL fstarpu_vector_data_register(handles(357), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(SOURSIN)
    IF (ASSOCIATED(WSX1)) THEN
      SCALAR_INT_BUF(358) = 1
      CALL fstarpu_vector_data_register(handles(358), 0, C_LOC(WSX1(LBOUND(WSX1,1))), SIZE(WSX1,1), C_SIZEOF(WSX1(LBOUND(WSX1,1))))
    ELSE
      SCALAR_INT_BUF(358) = 0
      CALL fstarpu_vector_data_register(handles(358), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(WSX1)
    IF (ASSOCIATED(WSY1)) THEN
      SCALAR_INT_BUF(359) = 1
      CALL fstarpu_vector_data_register(handles(359), 0, C_LOC(WSY1(LBOUND(WSY1,1))), SIZE(WSY1,1), C_SIZEOF(WSY1(LBOUND(WSY1,1))))
    ELSE
      SCALAR_INT_BUF(359) = 0
      CALL fstarpu_vector_data_register(handles(359), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(WSY1)
    IF (ASSOCIATED(PR1)) THEN
      SCALAR_INT_BUF(360) = 1
      CALL fstarpu_vector_data_register(handles(360), 0, C_LOC(PR1(LBOUND(PR1,1))), SIZE(PR1,1), C_SIZEOF(PR1(LBOUND(PR1,1))))
    ELSE
      SCALAR_INT_BUF(360) = 0
      CALL fstarpu_vector_data_register(handles(360), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(PR1)
    IF (ASSOCIATED(WSX2)) THEN
      SCALAR_INT_BUF(361) = 1
      CALL fstarpu_vector_data_register(handles(361), 0, C_LOC(WSX2(LBOUND(WSX2,1))), SIZE(WSX2,1), C_SIZEOF(WSX2(LBOUND(WSX2,1))))
    ELSE
      SCALAR_INT_BUF(361) = 0
      CALL fstarpu_vector_data_register(handles(361), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(WSX2)
    IF (ASSOCIATED(WSY2)) THEN
      SCALAR_INT_BUF(362) = 1
      CALL fstarpu_vector_data_register(handles(362), 0, C_LOC(WSY2(LBOUND(WSY2,1))), SIZE(WSY2,1), C_SIZEOF(WSY2(LBOUND(WSY2,1))))
    ELSE
      SCALAR_INT_BUF(362) = 0
      CALL fstarpu_vector_data_register(handles(362), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(WSY2)
    IF (ASSOCIATED(PR2)) THEN
      SCALAR_INT_BUF(363) = 1
      CALL fstarpu_vector_data_register(handles(363), 0, C_LOC(PR2(LBOUND(PR2,1))), SIZE(PR2,1), C_SIZEOF(PR2(LBOUND(PR2,1))))
    ELSE
      SCALAR_INT_BUF(363) = 0
      CALL fstarpu_vector_data_register(handles(363), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(PR2)
    IF (ASSOCIATED(WVNX1)) THEN
      SCALAR_INT_BUF(364) = 1
      CALL fstarpu_vector_data_register(handles(364), 0, C_LOC(WVNX1(LBOUND(WVNX1,1))), SIZE(WVNX1,1), C_SIZEOF(WVNX1(LBOUND(WVNX1,1))))
    ELSE
      SCALAR_INT_BUF(364) = 0
      CALL fstarpu_vector_data_register(handles(364), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(WVNX1)
    IF (ASSOCIATED(WVNY1)) THEN
      SCALAR_INT_BUF(365) = 1
      CALL fstarpu_vector_data_register(handles(365), 0, C_LOC(WVNY1(LBOUND(WVNY1,1))), SIZE(WVNY1,1), C_SIZEOF(WVNY1(LBOUND(WVNY1,1))))
    ELSE
      SCALAR_INT_BUF(365) = 0
      CALL fstarpu_vector_data_register(handles(365), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(WVNY1)
    IF (ASSOCIATED(PRN1)) THEN
      SCALAR_INT_BUF(366) = 1
      CALL fstarpu_vector_data_register(handles(366), 0, C_LOC(PRN1(LBOUND(PRN1,1))), SIZE(PRN1,1), C_SIZEOF(PRN1(LBOUND(PRN1,1))))
    ELSE
      SCALAR_INT_BUF(366) = 0
      CALL fstarpu_vector_data_register(handles(366), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(PRN1)
    IF (ASSOCIATED(WVNX2)) THEN
      SCALAR_INT_BUF(367) = 1
      CALL fstarpu_vector_data_register(handles(367), 0, C_LOC(WVNX2(LBOUND(WVNX2,1))), SIZE(WVNX2,1), C_SIZEOF(WVNX2(LBOUND(WVNX2,1))))
    ELSE
      SCALAR_INT_BUF(367) = 0
      CALL fstarpu_vector_data_register(handles(367), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(WVNX2)
    IF (ASSOCIATED(WVNY2)) THEN
      SCALAR_INT_BUF(368) = 1
      CALL fstarpu_vector_data_register(handles(368), 0, C_LOC(WVNY2(LBOUND(WVNY2,1))), SIZE(WVNY2,1), C_SIZEOF(WVNY2(LBOUND(WVNY2,1))))
    ELSE
      SCALAR_INT_BUF(368) = 0
      CALL fstarpu_vector_data_register(handles(368), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(WVNY2)
    IF (ASSOCIATED(PRN2)) THEN
      SCALAR_INT_BUF(369) = 1
      CALL fstarpu_vector_data_register(handles(369), 0, C_LOC(PRN2(LBOUND(PRN2,1))), SIZE(PRN2,1), C_SIZEOF(PRN2(LBOUND(PRN2,1))))
    ELSE
      SCALAR_INT_BUF(369) = 0
      CALL fstarpu_vector_data_register(handles(369), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(PRN2)
    IF (ASSOCIATED(RSNX1)) THEN
      SCALAR_INT_BUF(370) = 1
      CALL fstarpu_vector_data_register(handles(370), 0, C_LOC(RSNX1(LBOUND(RSNX1,1))), SIZE(RSNX1,1), C_SIZEOF(RSNX1(LBOUND(RSNX1,1))))
    ELSE
      SCALAR_INT_BUF(370) = 0
      CALL fstarpu_vector_data_register(handles(370), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(RSNX1)
    IF (ASSOCIATED(RSNY1)) THEN
      SCALAR_INT_BUF(371) = 1
      CALL fstarpu_vector_data_register(handles(371), 0, C_LOC(RSNY1(LBOUND(RSNY1,1))), SIZE(RSNY1,1), C_SIZEOF(RSNY1(LBOUND(RSNY1,1))))
    ELSE
      SCALAR_INT_BUF(371) = 0
      CALL fstarpu_vector_data_register(handles(371), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(RSNY1)
    IF (ASSOCIATED(RSNX2)) THEN
      SCALAR_INT_BUF(372) = 1
      CALL fstarpu_vector_data_register(handles(372), 0, C_LOC(RSNX2(LBOUND(RSNX2,1))), SIZE(RSNX2,1), C_SIZEOF(RSNX2(LBOUND(RSNX2,1))))
    ELSE
      SCALAR_INT_BUF(372) = 0
      CALL fstarpu_vector_data_register(handles(372), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(RSNX2)
    IF (ASSOCIATED(RSNY2)) THEN
      SCALAR_INT_BUF(373) = 1
      CALL fstarpu_vector_data_register(handles(373), 0, C_LOC(RSNY2(LBOUND(RSNY2,1))), SIZE(RSNY2,1), C_SIZEOF(RSNY2(LBOUND(RSNY2,1))))
    ELSE
      SCALAR_INT_BUF(373) = 0
      CALL fstarpu_vector_data_register(handles(373), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(RSNY2)
    IF (ASSOCIATED(RSNXOUT)) THEN
      SCALAR_INT_BUF(374) = 1
      CALL fstarpu_vector_data_register(handles(374), 0, C_LOC(RSNXOUT(LBOUND(RSNXOUT,1))), SIZE(RSNXOUT,1), C_SIZEOF(RSNXOUT(LBOUND(RSNXOUT,1))))
    ELSE
      SCALAR_INT_BUF(374) = 0
      CALL fstarpu_vector_data_register(handles(374), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(RSNXOUT)
    IF (ASSOCIATED(RSNYOUT)) THEN
      SCALAR_INT_BUF(375) = 1
      CALL fstarpu_vector_data_register(handles(375), 0, C_LOC(RSNYOUT(LBOUND(RSNYOUT,1))), SIZE(RSNYOUT,1), C_SIZEOF(RSNYOUT(LBOUND(RSNYOUT,1))))
    ELSE
      SCALAR_INT_BUF(375) = 0
      CALL fstarpu_vector_data_register(handles(375), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(RSNYOUT)
    IF (ASSOCIATED(WAVE_T1)) THEN
      SCALAR_INT_BUF(376) = 1
      CALL fstarpu_vector_data_register(handles(376), 0, C_LOC(WAVE_T1(LBOUND(WAVE_T1,1))), SIZE(WAVE_T1,1), C_SIZEOF(WAVE_T1(LBOUND(WAVE_T1,1))))
    ELSE
      SCALAR_INT_BUF(376) = 0
      CALL fstarpu_vector_data_register(handles(376), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(WAVE_T1)
    IF (ASSOCIATED(WAVE_H1)) THEN
      SCALAR_INT_BUF(377) = 1
      CALL fstarpu_vector_data_register(handles(377), 0, C_LOC(WAVE_H1(LBOUND(WAVE_H1,1))), SIZE(WAVE_H1,1), C_SIZEOF(WAVE_H1(LBOUND(WAVE_H1,1))))
    ELSE
      SCALAR_INT_BUF(377) = 0
      CALL fstarpu_vector_data_register(handles(377), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(WAVE_H1)
    IF (ASSOCIATED(WAVE_A1)) THEN
      SCALAR_INT_BUF(378) = 1
      CALL fstarpu_vector_data_register(handles(378), 0, C_LOC(WAVE_A1(LBOUND(WAVE_A1,1))), SIZE(WAVE_A1,1), C_SIZEOF(WAVE_A1(LBOUND(WAVE_A1,1))))
    ELSE
      SCALAR_INT_BUF(378) = 0
      CALL fstarpu_vector_data_register(handles(378), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(WAVE_A1)
    IF (ASSOCIATED(WAVE_D1)) THEN
      SCALAR_INT_BUF(379) = 1
      CALL fstarpu_vector_data_register(handles(379), 0, C_LOC(WAVE_D1(LBOUND(WAVE_D1,1))), SIZE(WAVE_D1,1), C_SIZEOF(WAVE_D1(LBOUND(WAVE_D1,1))))
    ELSE
      SCALAR_INT_BUF(379) = 0
      CALL fstarpu_vector_data_register(handles(379), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(WAVE_D1)
    IF (ASSOCIATED(WAVE_T2)) THEN
      SCALAR_INT_BUF(380) = 1
      CALL fstarpu_vector_data_register(handles(380), 0, C_LOC(WAVE_T2(LBOUND(WAVE_T2,1))), SIZE(WAVE_T2,1), C_SIZEOF(WAVE_T2(LBOUND(WAVE_T2,1))))
    ELSE
      SCALAR_INT_BUF(380) = 0
      CALL fstarpu_vector_data_register(handles(380), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(WAVE_T2)
    IF (ASSOCIATED(WAVE_H2)) THEN
      SCALAR_INT_BUF(381) = 1
      CALL fstarpu_vector_data_register(handles(381), 0, C_LOC(WAVE_H2(LBOUND(WAVE_H2,1))), SIZE(WAVE_H2,1), C_SIZEOF(WAVE_H2(LBOUND(WAVE_H2,1))))
    ELSE
      SCALAR_INT_BUF(381) = 0
      CALL fstarpu_vector_data_register(handles(381), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(WAVE_H2)
    IF (ASSOCIATED(WAVE_A2)) THEN
      SCALAR_INT_BUF(382) = 1
      CALL fstarpu_vector_data_register(handles(382), 0, C_LOC(WAVE_A2(LBOUND(WAVE_A2,1))), SIZE(WAVE_A2,1), C_SIZEOF(WAVE_A2(LBOUND(WAVE_A2,1))))
    ELSE
      SCALAR_INT_BUF(382) = 0
      CALL fstarpu_vector_data_register(handles(382), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(WAVE_A2)
    IF (ASSOCIATED(WAVE_D2)) THEN
      SCALAR_INT_BUF(383) = 1
      CALL fstarpu_vector_data_register(handles(383), 0, C_LOC(WAVE_D2(LBOUND(WAVE_D2,1))), SIZE(WAVE_D2,1), C_SIZEOF(WAVE_D2(LBOUND(WAVE_D2,1))))
    ELSE
      SCALAR_INT_BUF(383) = 0
      CALL fstarpu_vector_data_register(handles(383), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(WAVE_D2)
    IF (ASSOCIATED(WAVE_T)) THEN
      SCALAR_INT_BUF(384) = 1
      CALL fstarpu_vector_data_register(handles(384), 0, C_LOC(WAVE_T(LBOUND(WAVE_T,1))), SIZE(WAVE_T,1), C_SIZEOF(WAVE_T(LBOUND(WAVE_T,1))))
    ELSE
      SCALAR_INT_BUF(384) = 0
      CALL fstarpu_vector_data_register(handles(384), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(WAVE_T)
    IF (ASSOCIATED(WAVE_H)) THEN
      SCALAR_INT_BUF(385) = 1
      CALL fstarpu_vector_data_register(handles(385), 0, C_LOC(WAVE_H(LBOUND(WAVE_H,1))), SIZE(WAVE_H,1), C_SIZEOF(WAVE_H(LBOUND(WAVE_H,1))))
    ELSE
      SCALAR_INT_BUF(385) = 0
      CALL fstarpu_vector_data_register(handles(385), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(WAVE_H)
    IF (ASSOCIATED(WAVE_A)) THEN
      SCALAR_INT_BUF(386) = 1
      CALL fstarpu_vector_data_register(handles(386), 0, C_LOC(WAVE_A(LBOUND(WAVE_A,1))), SIZE(WAVE_A,1), C_SIZEOF(WAVE_A(LBOUND(WAVE_A,1))))
    ELSE
      SCALAR_INT_BUF(386) = 0
      CALL fstarpu_vector_data_register(handles(386), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(WAVE_A)
    IF (ASSOCIATED(WAVE_D)) THEN
      SCALAR_INT_BUF(387) = 1
      CALL fstarpu_vector_data_register(handles(387), 0, C_LOC(WAVE_D(LBOUND(WAVE_D,1))), SIZE(WAVE_D,1), C_SIZEOF(WAVE_D(LBOUND(WAVE_D,1))))
    ELSE
      SCALAR_INT_BUF(387) = 0
      CALL fstarpu_vector_data_register(handles(387), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(WAVE_D)
    IF (ASSOCIATED(WB)) THEN
      SCALAR_INT_BUF(388) = 1
      CALL fstarpu_vector_data_register(handles(388), 0, C_LOC(WB(LBOUND(WB,1))), SIZE(WB,1), C_SIZEOF(WB(LBOUND(WB,1))))
    ELSE
      SCALAR_INT_BUF(388) = 0
      CALL fstarpu_vector_data_register(handles(388), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(WB)
    IF (ASSOCIATED(WVNXOUT)) THEN
      SCALAR_INT_BUF(389) = 1
      CALL fstarpu_vector_data_register(handles(389), 0, C_LOC(WVNXOUT(LBOUND(WVNXOUT,1))), SIZE(WVNXOUT,1), C_SIZEOF(WVNXOUT(LBOUND(WVNXOUT,1))))
    ELSE
      SCALAR_INT_BUF(389) = 0
      CALL fstarpu_vector_data_register(handles(389), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(WVNXOUT)
    IF (ASSOCIATED(WVNYOUT)) THEN
      SCALAR_INT_BUF(390) = 1
      CALL fstarpu_vector_data_register(handles(390), 0, C_LOC(WVNYOUT(LBOUND(WVNYOUT,1))), SIZE(WVNYOUT,1), C_SIZEOF(WVNYOUT(LBOUND(WVNYOUT,1))))
    ELSE
      SCALAR_INT_BUF(390) = 0
      CALL fstarpu_vector_data_register(handles(390), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(WVNYOUT)
    IF (ASSOCIATED(TKXX)) THEN
      SCALAR_INT_BUF(391) = 1
      CALL fstarpu_vector_data_register(handles(391), 0, C_LOC(TKXX(LBOUND(TKXX,1))), SIZE(TKXX,1), C_SIZEOF(TKXX(LBOUND(TKXX,1))))
    ELSE
      SCALAR_INT_BUF(391) = 0
      CALL fstarpu_vector_data_register(handles(391), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(TKXX)
    IF (ASSOCIATED(TKYY)) THEN
      SCALAR_INT_BUF(392) = 1
      CALL fstarpu_vector_data_register(handles(392), 0, C_LOC(TKYY(LBOUND(TKYY,1))), SIZE(TKYY,1), C_SIZEOF(TKYY(LBOUND(TKYY,1))))
    ELSE
      SCALAR_INT_BUF(392) = 0
      CALL fstarpu_vector_data_register(handles(392), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(TKYY)
    IF (ASSOCIATED(TKXY)) THEN
      SCALAR_INT_BUF(393) = 1
      CALL fstarpu_vector_data_register(handles(393), 0, C_LOC(TKXY(LBOUND(TKXY,1))), SIZE(TKXY,1), C_SIZEOF(TKXY(LBOUND(TKXY,1))))
    ELSE
      SCALAR_INT_BUF(393) = 0
      CALL fstarpu_vector_data_register(handles(393), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(TKXY)
    IF (ASSOCIATED(EMO)) THEN
      SCALAR_INT_BUF(394) = 1
      CALL fstarpu_matrix_data_register(handles(394), 0, C_LOC(EMO(LBOUND(EMO,1),LBOUND(EMO,2))), SIZE(EMO,1), SIZE(EMO,1), SIZE(EMO,2), C_SIZEOF(EMO(LBOUND(EMO,1),LBOUND(EMO,2))))
    ELSE
      SCALAR_INT_BUF(394) = 0
      CALL fstarpu_matrix_data_register(handles(394), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(EMO)
    IF (ASSOCIATED(EFA)) THEN
      SCALAR_INT_BUF(395) = 1
      CALL fstarpu_matrix_data_register(handles(395), 0, C_LOC(EFA(LBOUND(EFA,1),LBOUND(EFA,2))), SIZE(EFA,1), SIZE(EFA,1), SIZE(EFA,2), C_SIZEOF(EFA(LBOUND(EFA,1),LBOUND(EFA,2))))
    ELSE
      SCALAR_INT_BUF(395) = 0
      CALL fstarpu_matrix_data_register(handles(395), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(EFA)
    IF (ASSOCIATED(UMO)) THEN
      SCALAR_INT_BUF(396) = 1
      CALL fstarpu_matrix_data_register(handles(396), 0, C_LOC(UMO(LBOUND(UMO,1),LBOUND(UMO,2))), SIZE(UMO,1), SIZE(UMO,1), SIZE(UMO,2), C_SIZEOF(UMO(LBOUND(UMO,1),LBOUND(UMO,2))))
    ELSE
      SCALAR_INT_BUF(396) = 0
      CALL fstarpu_matrix_data_register(handles(396), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(UMO)
    IF (ASSOCIATED(UFA)) THEN
      SCALAR_INT_BUF(397) = 1
      CALL fstarpu_matrix_data_register(handles(397), 0, C_LOC(UFA(LBOUND(UFA,1),LBOUND(UFA,2))), SIZE(UFA,1), SIZE(UFA,1), SIZE(UFA,2), C_SIZEOF(UFA(LBOUND(UFA,1),LBOUND(UFA,2))))
    ELSE
      SCALAR_INT_BUF(397) = 0
      CALL fstarpu_matrix_data_register(handles(397), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(UFA)
    IF (ASSOCIATED(VMO)) THEN
      SCALAR_INT_BUF(398) = 1
      CALL fstarpu_matrix_data_register(handles(398), 0, C_LOC(VMO(LBOUND(VMO,1),LBOUND(VMO,2))), SIZE(VMO,1), SIZE(VMO,1), SIZE(VMO,2), C_SIZEOF(VMO(LBOUND(VMO,1),LBOUND(VMO,2))))
    ELSE
      SCALAR_INT_BUF(398) = 0
      CALL fstarpu_matrix_data_register(handles(398), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(VMO)
    IF (ASSOCIATED(VFA)) THEN
      SCALAR_INT_BUF(399) = 1
      CALL fstarpu_matrix_data_register(handles(399), 0, C_LOC(VFA(LBOUND(VFA,1),LBOUND(VFA,2))), SIZE(VFA,1), SIZE(VFA,1), SIZE(VFA,2), C_SIZEOF(VFA(LBOUND(VFA,1),LBOUND(VFA,2))))
    ELSE
      SCALAR_INT_BUF(399) = 0
      CALL fstarpu_matrix_data_register(handles(399), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(VFA)
    IF (ASSOCIATED(XEL)) THEN
      SCALAR_INT_BUF(400) = 1
      CALL fstarpu_vector_data_register(handles(400), 0, C_LOC(XEL(LBOUND(XEL,1))), SIZE(XEL,1), C_SIZEOF(XEL(LBOUND(XEL,1))))
    ELSE
      SCALAR_INT_BUF(400) = 0
      CALL fstarpu_vector_data_register(handles(400), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(XEL)
    IF (ASSOCIATED(YEL)) THEN
      SCALAR_INT_BUF(401) = 1
      CALL fstarpu_vector_data_register(handles(401), 0, C_LOC(YEL(LBOUND(YEL,1))), SIZE(YEL,1), C_SIZEOF(YEL(LBOUND(YEL,1))))
    ELSE
      SCALAR_INT_BUF(401) = 0
      CALL fstarpu_vector_data_register(handles(401), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(YEL)
    IF (ASSOCIATED(SLEL)) THEN
      SCALAR_INT_BUF(402) = 1
      CALL fstarpu_vector_data_register(handles(402), 0, C_LOC(SLEL(LBOUND(SLEL,1))), SIZE(SLEL,1), C_SIZEOF(SLEL(LBOUND(SLEL,1))))
    ELSE
      SCALAR_INT_BUF(402) = 0
      CALL fstarpu_vector_data_register(handles(402), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(SLEL)
    IF (ASSOCIATED(SFEL)) THEN
      SCALAR_INT_BUF(403) = 1
      CALL fstarpu_vector_data_register(handles(403), 0, C_LOC(SFEL(LBOUND(SFEL,1))), SIZE(SFEL,1), C_SIZEOF(SFEL(LBOUND(SFEL,1))))
    ELSE
      SCALAR_INT_BUF(403) = 0
      CALL fstarpu_vector_data_register(handles(403), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(SFEL)
    IF (ASSOCIATED(AREAS)) THEN
      SCALAR_INT_BUF(404) = 1
      CALL fstarpu_vector_data_register(handles(404), 0, C_LOC(AREAS(LBOUND(AREAS,1))), SIZE(AREAS,1), C_SIZEOF(AREAS(LBOUND(AREAS,1))))
    ELSE
      SCALAR_INT_BUF(404) = 0
      CALL fstarpu_vector_data_register(handles(404), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(AREAS)
    IF (ASSOCIATED(SFACDUB)) THEN
      SCALAR_INT_BUF(405) = 1
      CALL fstarpu_matrix_data_register(handles(405), 0, C_LOC(SFACDUB(LBOUND(SFACDUB,1),LBOUND(SFACDUB,2))), SIZE(SFACDUB,1), SIZE(SFACDUB,1), SIZE(SFACDUB,2), C_SIZEOF(SFACDUB(LBOUND(SFACDUB,1),LBOUND(SFACDUB,2))))
    ELSE
      SCALAR_INT_BUF(405) = 0
      CALL fstarpu_matrix_data_register(handles(405), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(SFACDUB)
    IF (ASSOCIATED(YDUB)) THEN
      SCALAR_INT_BUF(406) = 1
      CALL fstarpu_block_data_register(handles(406), 0, C_LOC(YDUB(LBOUND(YDUB,1),LBOUND(YDUB,2),LBOUND(YDUB,3))), SIZE(YDUB,1), SIZE(YDUB,1)*SIZE(YDUB,2), SIZE(YDUB,1), SIZE(YDUB,2), SIZE(YDUB,3), C_SIZEOF(YDUB(LBOUND(YDUB,1),LBOUND(YDUB,2),LBOUND(YDUB,3))))
    ELSE
      SCALAR_INT_BUF(406) = 0
      CALL fstarpu_block_data_register(handles(406), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(YDUB)
    IF (ASSOCIATED(RTEMP2)) THEN
      SCALAR_INT_BUF(407) = 1
      CALL fstarpu_vector_data_register(handles(407), 0, C_LOC(RTEMP2(LBOUND(RTEMP2,1))), SIZE(RTEMP2,1), C_SIZEOF(RTEMP2(LBOUND(RTEMP2,1))))
    ELSE
      SCALAR_INT_BUF(407) = 0
      CALL fstarpu_vector_data_register(handles(407), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(RTEMP2)
    IF (ASSOCIATED(AUV11)) THEN
      SCALAR_INT_BUF(408) = 1
      CALL fstarpu_vector_data_register(handles(408), 0, C_LOC(AUV11(LBOUND(AUV11,1))), SIZE(AUV11,1), C_SIZEOF(AUV11(LBOUND(AUV11,1))))
    ELSE
      SCALAR_INT_BUF(408) = 0
      CALL fstarpu_vector_data_register(handles(408), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(AUV11)
    IF (ASSOCIATED(AUV12)) THEN
      SCALAR_INT_BUF(409) = 1
      CALL fstarpu_vector_data_register(handles(409), 0, C_LOC(AUV12(LBOUND(AUV12,1))), SIZE(AUV12,1), C_SIZEOF(AUV12(LBOUND(AUV12,1))))
    ELSE
      SCALAR_INT_BUF(409) = 0
      CALL fstarpu_vector_data_register(handles(409), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(AUV12)
    IF (ASSOCIATED(AUV13)) THEN
      SCALAR_INT_BUF(410) = 1
      CALL fstarpu_vector_data_register(handles(410), 0, C_LOC(AUV13(LBOUND(AUV13,1))), SIZE(AUV13,1), C_SIZEOF(AUV13(LBOUND(AUV13,1))))
    ELSE
      SCALAR_INT_BUF(410) = 0
      CALL fstarpu_vector_data_register(handles(410), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(AUV13)
    IF (ASSOCIATED(AUV14)) THEN
      SCALAR_INT_BUF(411) = 1
      CALL fstarpu_vector_data_register(handles(411), 0, C_LOC(AUV14(LBOUND(AUV14,1))), SIZE(AUV14,1), C_SIZEOF(AUV14(LBOUND(AUV14,1))))
    ELSE
      SCALAR_INT_BUF(411) = 0
      CALL fstarpu_vector_data_register(handles(411), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(AUV14)
    IF (ASSOCIATED(AUVXX)) THEN
      SCALAR_INT_BUF(412) = 1
      CALL fstarpu_vector_data_register(handles(412), 0, C_LOC(AUVXX(LBOUND(AUVXX,1))), SIZE(AUVXX,1), C_SIZEOF(AUVXX(LBOUND(AUVXX,1))))
    ELSE
      SCALAR_INT_BUF(412) = 0
      CALL fstarpu_vector_data_register(handles(412), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(AUVXX)
    IF (ASSOCIATED(AUVYY)) THEN
      SCALAR_INT_BUF(413) = 1
      CALL fstarpu_vector_data_register(handles(413), 0, C_LOC(AUVYY(LBOUND(AUVYY,1))), SIZE(AUVYY,1), C_SIZEOF(AUVYY(LBOUND(AUVYY,1))))
    ELSE
      SCALAR_INT_BUF(413) = 0
      CALL fstarpu_vector_data_register(handles(413), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(AUVYY)
    IF (ASSOCIATED(AUVXY)) THEN
      SCALAR_INT_BUF(414) = 1
      CALL fstarpu_vector_data_register(handles(414), 0, C_LOC(AUVXY(LBOUND(AUVXY,1))), SIZE(AUVXY,1), C_SIZEOF(AUVXY(LBOUND(AUVXY,1))))
    ELSE
      SCALAR_INT_BUF(414) = 0
      CALL fstarpu_vector_data_register(handles(414), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(AUVXY)
    IF (ASSOCIATED(AUVYX)) THEN
      SCALAR_INT_BUF(415) = 1
      CALL fstarpu_vector_data_register(handles(415), 0, C_LOC(AUVYX(LBOUND(AUVYX,1))), SIZE(AUVYX,1), C_SIZEOF(AUVYX(LBOUND(AUVYX,1))))
    ELSE
      SCALAR_INT_BUF(415) = 0
      CALL fstarpu_vector_data_register(handles(415), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(AUVYX)
    IF (ASSOCIATED(DUU1)) THEN
      SCALAR_INT_BUF(416) = 1
      CALL fstarpu_vector_data_register(handles(416), 0, C_LOC(DUU1(LBOUND(DUU1,1))), SIZE(DUU1,1), C_SIZEOF(DUU1(LBOUND(DUU1,1))))
    ELSE
      SCALAR_INT_BUF(416) = 0
      CALL fstarpu_vector_data_register(handles(416), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(DUU1)
    IF (ASSOCIATED(DUV1)) THEN
      SCALAR_INT_BUF(417) = 1
      CALL fstarpu_vector_data_register(handles(417), 0, C_LOC(DUV1(LBOUND(DUV1,1))), SIZE(DUV1,1), C_SIZEOF(DUV1(LBOUND(DUV1,1))))
    ELSE
      SCALAR_INT_BUF(417) = 0
      CALL fstarpu_vector_data_register(handles(417), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(DUV1)
    IF (ASSOCIATED(DVV1)) THEN
      SCALAR_INT_BUF(418) = 1
      CALL fstarpu_vector_data_register(handles(418), 0, C_LOC(DVV1(LBOUND(DVV1,1))), SIZE(DVV1,1), C_SIZEOF(DVV1(LBOUND(DVV1,1))))
    ELSE
      SCALAR_INT_BUF(418) = 0
      CALL fstarpu_vector_data_register(handles(418), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(DVV1)
    IF (ASSOCIATED(BSX1)) THEN
      SCALAR_INT_BUF(419) = 1
      CALL fstarpu_vector_data_register(handles(419), 0, C_LOC(BSX1(LBOUND(BSX1,1))), SIZE(BSX1,1), C_SIZEOF(BSX1(LBOUND(BSX1,1))))
    ELSE
      SCALAR_INT_BUF(419) = 0
      CALL fstarpu_vector_data_register(handles(419), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(BSX1)
    IF (ASSOCIATED(BSY1)) THEN
      SCALAR_INT_BUF(420) = 1
      CALL fstarpu_vector_data_register(handles(420), 0, C_LOC(BSY1(LBOUND(BSY1,1))), SIZE(BSY1,1), C_SIZEOF(BSY1(LBOUND(BSY1,1))))
    ELSE
      SCALAR_INT_BUF(420) = 0
      CALL fstarpu_vector_data_register(handles(420), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(BSY1)
    IF (ASSOCIATED(TIP1)) THEN
      SCALAR_INT_BUF(421) = 1
      CALL fstarpu_vector_data_register(handles(421), 0, C_LOC(TIP1(LBOUND(TIP1,1))), SIZE(TIP1,1), C_SIZEOF(TIP1(LBOUND(TIP1,1))))
    ELSE
      SCALAR_INT_BUF(421) = 0
      CALL fstarpu_vector_data_register(handles(421), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(TIP1)
    IF (ASSOCIATED(TIP2)) THEN
      SCALAR_INT_BUF(422) = 1
      CALL fstarpu_vector_data_register(handles(422), 0, C_LOC(TIP2(LBOUND(TIP2,1))), SIZE(TIP2,1), C_SIZEOF(TIP2(LBOUND(TIP2,1))))
    ELSE
      SCALAR_INT_BUF(422) = 0
      CALL fstarpu_vector_data_register(handles(422), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(TIP2)
    IF (ASSOCIATED(SALTAMP)) THEN
      SCALAR_INT_BUF(423) = 1
      CALL fstarpu_matrix_data_register(handles(423), 0, C_LOC(SALTAMP(LBOUND(SALTAMP,1),LBOUND(SALTAMP,2))), SIZE(SALTAMP,1), SIZE(SALTAMP,1), SIZE(SALTAMP,2), C_SIZEOF(SALTAMP(LBOUND(SALTAMP,1),LBOUND(SALTAMP,2))))
    ELSE
      SCALAR_INT_BUF(423) = 0
      CALL fstarpu_matrix_data_register(handles(423), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(SALTAMP)
    IF (ASSOCIATED(SALTPHA)) THEN
      SCALAR_INT_BUF(424) = 1
      CALL fstarpu_matrix_data_register(handles(424), 0, C_LOC(SALTPHA(LBOUND(SALTPHA,1),LBOUND(SALTPHA,2))), SIZE(SALTPHA,1), SIZE(SALTPHA,1), SIZE(SALTPHA,2), C_SIZEOF(SALTPHA(LBOUND(SALTPHA,1),LBOUND(SALTPHA,2))))
    ELSE
      SCALAR_INT_BUF(424) = 0
      CALL fstarpu_matrix_data_register(handles(424), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(SALTPHA)
    IF (ASSOCIATED(OBCCOEF)) THEN
      SCALAR_INT_BUF(425) = 1
      CALL fstarpu_matrix_data_register(handles(425), 0, C_LOC(OBCCOEF(LBOUND(OBCCOEF,1),LBOUND(OBCCOEF,2))), SIZE(OBCCOEF,1), SIZE(OBCCOEF,1), SIZE(OBCCOEF,2), C_SIZEOF(OBCCOEF(LBOUND(OBCCOEF,1),LBOUND(OBCCOEF,2))))
    ELSE
      SCALAR_INT_BUF(425) = 0
      CALL fstarpu_matrix_data_register(handles(425), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(OBCCOEF)
    IF (ASSOCIATED(COEF)) THEN
      SCALAR_INT_BUF(426) = 1
      CALL fstarpu_matrix_data_register(handles(426), 0, C_LOC(COEF(LBOUND(COEF,1),LBOUND(COEF,2))), SIZE(COEF,1), SIZE(COEF,1), SIZE(COEF,2), C_SIZEOF(COEF(LBOUND(COEF,1),LBOUND(COEF,2))))
    ELSE
      SCALAR_INT_BUF(426) = 0
      CALL fstarpu_matrix_data_register(handles(426), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(COEF)
    IF (ASSOCIATED(WKSP)) THEN
      SCALAR_INT_BUF(427) = 1
      CALL fstarpu_vector_data_register(handles(427), 0, C_LOC(WKSP(LBOUND(WKSP,1))), SIZE(WKSP,1), C_SIZEOF(WKSP(LBOUND(WKSP,1))))
    ELSE
      SCALAR_INT_BUF(427) = 0
      CALL fstarpu_vector_data_register(handles(427), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(WKSP)
    IF (ASSOCIATED(RPARM)) THEN
      SCALAR_INT_BUF(428) = 1
      CALL fstarpu_vector_data_register(handles(428), 0, C_LOC(RPARM(LBOUND(RPARM,1))), SIZE(RPARM,1), C_SIZEOF(RPARM(LBOUND(RPARM,1))))
    ELSE
      SCALAR_INT_BUF(428) = 0
      CALL fstarpu_vector_data_register(handles(428), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(RPARM)
    IF (ASSOCIATED(ABD)) THEN
      SCALAR_INT_BUF(429) = 1
      CALL fstarpu_matrix_data_register(handles(429), 0, C_LOC(ABD(LBOUND(ABD,1),LBOUND(ABD,2))), SIZE(ABD,1), SIZE(ABD,1), SIZE(ABD,2), C_SIZEOF(ABD(LBOUND(ABD,1),LBOUND(ABD,2))))
    ELSE
      SCALAR_INT_BUF(429) = 0
      CALL fstarpu_matrix_data_register(handles(429), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(ABD)
    IF (ASSOCIATED(ZX)) THEN
      SCALAR_INT_BUF(430) = 1
      CALL fstarpu_vector_data_register(handles(430), 0, C_LOC(ZX(LBOUND(ZX,1))), SIZE(ZX,1), C_SIZEOF(ZX(LBOUND(ZX,1))))
    ELSE
      SCALAR_INT_BUF(430) = 0
      CALL fstarpu_vector_data_register(handles(430), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(ZX)
    IF (ASSOCIATED(GRAVX)) THEN
      SCALAR_INT_BUF(431) = 1
      CALL fstarpu_vector_data_register(handles(431), 0, C_LOC(GRAVX(LBOUND(GRAVX,1))), SIZE(GRAVX,1), C_SIZEOF(GRAVX(LBOUND(GRAVX,1))))
    ELSE
      SCALAR_INT_BUF(431) = 0
      CALL fstarpu_vector_data_register(handles(431), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(GRAVX)
    IF (ASSOCIATED(GRAVY)) THEN
      SCALAR_INT_BUF(432) = 1
      CALL fstarpu_vector_data_register(handles(432), 0, C_LOC(GRAVY(LBOUND(GRAVY,1))), SIZE(GRAVY,1), C_SIZEOF(GRAVY(LBOUND(GRAVY,1))))
    ELSE
      SCALAR_INT_BUF(432) = 0
      CALL fstarpu_vector_data_register(handles(432), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(GRAVY)
    IF (ASSOCIATED(ME2GW)) THEN
      SCALAR_INT_BUF(433) = 1
      CALL fstarpu_vector_data_register(handles(433), 0, C_LOC(ME2GW(LBOUND(ME2GW,1))), SIZE(ME2GW,1), C_SIZEOF(ME2GW(LBOUND(ME2GW,1))))
    ELSE
      SCALAR_INT_BUF(433) = 0
      CALL fstarpu_vector_data_register(handles(433), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(ME2GW)
    IF (ASSOCIATED(NBV)) THEN
      SCALAR_INT_BUF(434) = 1
      CALL fstarpu_vector_data_register(handles(434), 0, C_LOC(NBV(LBOUND(NBV,1))), SIZE(NBV,1), C_SIZEOF(NBV(LBOUND(NBV,1))))
    ELSE
      SCALAR_INT_BUF(434) = 0
      CALL fstarpu_vector_data_register(handles(434), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(NBV)
    IF (ASSOCIATED(LBCODEI)) THEN
      SCALAR_INT_BUF(435) = 1
      CALL fstarpu_vector_data_register(handles(435), 0, C_LOC(LBCODEI(LBOUND(LBCODEI,1))), SIZE(LBCODEI,1), C_SIZEOF(LBCODEI(LBOUND(LBCODEI,1))))
    ELSE
      SCALAR_INT_BUF(435) = 0
      CALL fstarpu_vector_data_register(handles(435), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(LBCODEI)
    IF (ASSOCIATED(NNODECODE)) THEN
      SCALAR_INT_BUF(436) = 1
      CALL fstarpu_vector_data_register(handles(436), 0, C_LOC(NNODECODE(LBOUND(NNODECODE,1))), SIZE(NNODECODE,1), C_SIZEOF(NNODECODE(LBOUND(NNODECODE,1))))
    ELSE
      SCALAR_INT_BUF(436) = 0
      CALL fstarpu_vector_data_register(handles(436), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(NNODECODE)
    IF (ASSOCIATED(NODECODE)) THEN
      SCALAR_INT_BUF(437) = 1
      CALL fstarpu_vector_data_register(handles(437), 0, C_LOC(NODECODE(LBOUND(NODECODE,1))), SIZE(NODECODE,1), C_SIZEOF(NODECODE(LBOUND(NODECODE,1))))
    ELSE
      SCALAR_INT_BUF(437) = 0
      CALL fstarpu_vector_data_register(handles(437), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(NODECODE)
    IF (ASSOCIATED(NODEREP)) THEN
      SCALAR_INT_BUF(438) = 1
      CALL fstarpu_vector_data_register(handles(438), 0, C_LOC(NODEREP(LBOUND(NODEREP,1))), SIZE(NODEREP,1), C_SIZEOF(NODEREP(LBOUND(NODEREP,1))))
    ELSE
      SCALAR_INT_BUF(438) = 0
      CALL fstarpu_vector_data_register(handles(438), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(NODEREP)
    IF (ASSOCIATED(NIBCNT)) THEN
      SCALAR_INT_BUF(439) = 1
      CALL fstarpu_vector_data_register(handles(439), 0, C_LOC(NIBCNT(LBOUND(NIBCNT,1))), SIZE(NIBCNT,1), C_SIZEOF(NIBCNT(LBOUND(NIBCNT,1))))
    ELSE
      SCALAR_INT_BUF(439) = 0
      CALL fstarpu_vector_data_register(handles(439), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(NIBCNT)
    IF (ASSOCIATED(NM)) THEN
      SCALAR_INT_BUF(440) = 1
      CALL fstarpu_matrix_data_register(handles(440), 0, C_LOC(NM(LBOUND(NM,1),LBOUND(NM,2))), SIZE(NM,1), SIZE(NM,1), SIZE(NM,2), C_SIZEOF(NM(LBOUND(NM,1),LBOUND(NM,2))))
    ELSE
      SCALAR_INT_BUF(440) = 0
      CALL fstarpu_matrix_data_register(handles(440), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(NM)
    IF (ASSOCIATED(NNEIGH)) THEN
      SCALAR_INT_BUF(441) = 1
      CALL fstarpu_vector_data_register(handles(441), 0, C_LOC(NNEIGH(LBOUND(NNEIGH,1))), SIZE(NNEIGH,1), C_SIZEOF(NNEIGH(LBOUND(NNEIGH,1))))
    ELSE
      SCALAR_INT_BUF(441) = 0
      CALL fstarpu_vector_data_register(handles(441), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(NNEIGH)
    IF (ASSOCIATED(MJU)) THEN
      SCALAR_INT_BUF(442) = 1
      CALL fstarpu_vector_data_register(handles(442), 0, C_LOC(MJU(LBOUND(MJU,1))), SIZE(MJU,1), C_SIZEOF(MJU(LBOUND(MJU,1))))
    ELSE
      SCALAR_INT_BUF(442) = 0
      CALL fstarpu_vector_data_register(handles(442), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(MJU)
    IF (ASSOCIATED(NODELE)) THEN
      SCALAR_INT_BUF(443) = 1
      CALL fstarpu_vector_data_register(handles(443), 0, C_LOC(NODELE(LBOUND(NODELE,1))), SIZE(NODELE,1), C_SIZEOF(NODELE(LBOUND(NODELE,1))))
    ELSE
      SCALAR_INT_BUF(443) = 0
      CALL fstarpu_vector_data_register(handles(443), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(NODELE)
    IF (ASSOCIATED(NEITAB)) THEN
      SCALAR_INT_BUF(444) = 1
      CALL fstarpu_matrix_data_register(handles(444), 0, C_LOC(NEITAB(LBOUND(NEITAB,1),LBOUND(NEITAB,2))), SIZE(NEITAB,1), SIZE(NEITAB,1), SIZE(NEITAB,2), C_SIZEOF(NEITAB(LBOUND(NEITAB,1),LBOUND(NEITAB,2))))
    ELSE
      SCALAR_INT_BUF(444) = 0
      CALL fstarpu_matrix_data_register(handles(444), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(NEITAB)
    IF (ASSOCIATED(NNEIGH_ELEM)) THEN
      SCALAR_INT_BUF(445) = 1
      CALL fstarpu_vector_data_register(handles(445), 0, C_LOC(NNEIGH_ELEM(LBOUND(NNEIGH_ELEM,1))), SIZE(NNEIGH_ELEM,1), C_SIZEOF(NNEIGH_ELEM(LBOUND(NNEIGH_ELEM,1))))
    ELSE
      SCALAR_INT_BUF(445) = 0
      CALL fstarpu_vector_data_register(handles(445), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(NNEIGH_ELEM)
    IF (ASSOCIATED(NIBNODECODE)) THEN
      SCALAR_INT_BUF(446) = 1
      CALL fstarpu_vector_data_register(handles(446), 0, C_LOC(NIBNODECODE(LBOUND(NIBNODECODE,1))), SIZE(NIBNODECODE,1), C_SIZEOF(NIBNODECODE(LBOUND(NIBNODECODE,1))))
    ELSE
      SCALAR_INT_BUF(446) = 0
      CALL fstarpu_vector_data_register(handles(446), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(NIBNODECODE)
    IF (ASSOCIATED(NEIGH_ELEM)) THEN
      SCALAR_INT_BUF(447) = 1
      CALL fstarpu_matrix_data_register(handles(447), 0, C_LOC(NEIGH_ELEM(LBOUND(NEIGH_ELEM,1),LBOUND(NEIGH_ELEM,2))), SIZE(NEIGH_ELEM,1), SIZE(NEIGH_ELEM,1), SIZE(NEIGH_ELEM,2), C_SIZEOF(NEIGH_ELEM(LBOUND(NEIGH_ELEM,1),LBOUND(NEIGH_ELEM,2))))
    ELSE
      SCALAR_INT_BUF(447) = 0
      CALL fstarpu_matrix_data_register(handles(447), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(NEIGH_ELEM)
    IF (ASSOCIATED(LBCODE)) THEN
      SCALAR_INT_BUF(448) = 1
      CALL fstarpu_vector_data_register(handles(448), 0, C_LOC(LBCODE(LBOUND(LBCODE,1))), SIZE(LBCODE,1), C_SIZEOF(LBCODE(LBOUND(LBCODE,1))))
    ELSE
      SCALAR_INT_BUF(448) = 0
      CALL fstarpu_vector_data_register(handles(448), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(LBCODE)
    IF (ASSOCIATED(NNC)) THEN
      SCALAR_INT_BUF(449) = 1
      CALL fstarpu_vector_data_register(handles(449), 0, C_LOC(NNC(LBOUND(NNC,1))), SIZE(NNC,1), C_SIZEOF(NNC(LBOUND(NNC,1))))
    ELSE
      SCALAR_INT_BUF(449) = 0
      CALL fstarpu_vector_data_register(handles(449), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(NNC)
    IF (ASSOCIATED(NNE)) THEN
      SCALAR_INT_BUF(450) = 1
      CALL fstarpu_vector_data_register(handles(450), 0, C_LOC(NNE(LBOUND(NNE,1))), SIZE(NNE,1), C_SIZEOF(NNE(LBOUND(NNE,1))))
    ELSE
      SCALAR_INT_BUF(450) = 0
      CALL fstarpu_vector_data_register(handles(450), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(NNE)
    IF (ASSOCIATED(NNV)) THEN
      SCALAR_INT_BUF(451) = 1
      CALL fstarpu_vector_data_register(handles(451), 0, C_LOC(NNV(LBOUND(NNV,1))), SIZE(NNV,1), C_SIZEOF(NNV(LBOUND(NNV,1))))
    ELSE
      SCALAR_INT_BUF(451) = 0
      CALL fstarpu_vector_data_register(handles(451), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(NNV)
    IF (ASSOCIATED(NNM)) THEN
      SCALAR_INT_BUF(452) = 1
      CALL fstarpu_vector_data_register(handles(452), 0, C_LOC(NNM(LBOUND(NNM,1))), SIZE(NNM,1), C_SIZEOF(NNM(LBOUND(NNM,1))))
    ELSE
      SCALAR_INT_BUF(452) = 0
      CALL fstarpu_vector_data_register(handles(452), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(NNM)
    IF (ASSOCIATED(IWKSP)) THEN
      SCALAR_INT_BUF(453) = 1
      CALL fstarpu_vector_data_register(handles(453), 0, C_LOC(IWKSP(LBOUND(IWKSP,1))), SIZE(IWKSP,1), C_SIZEOF(IWKSP(LBOUND(IWKSP,1))))
    ELSE
      SCALAR_INT_BUF(453) = 0
      CALL fstarpu_vector_data_register(handles(453), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(IWKSP)
    IF (ASSOCIATED(IPARM)) THEN
      SCALAR_INT_BUF(454) = 1
      CALL fstarpu_vector_data_register(handles(454), 0, C_LOC(IPARM(LBOUND(IPARM,1))), SIZE(IPARM,1), C_SIZEOF(IPARM(LBOUND(IPARM,1))))
    ELSE
      SCALAR_INT_BUF(454) = 0
      CALL fstarpu_vector_data_register(handles(454), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(IPARM)
    IF (ASSOCIATED(IPV)) THEN
      SCALAR_INT_BUF(455) = 1
      CALL fstarpu_vector_data_register(handles(455), 0, C_LOC(IPV(LBOUND(IPV,1))), SIZE(IPV,1), C_SIZEOF(IPV(LBOUND(IPV,1))))
    ELSE
      SCALAR_INT_BUF(455) = 0
      CALL fstarpu_vector_data_register(handles(455), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(IPV)
    IF (ASSOCIATED(NVDLL)) THEN
      SCALAR_INT_BUF(456) = 1
      CALL fstarpu_vector_data_register(handles(456), 0, C_LOC(NVDLL(LBOUND(NVDLL,1))), SIZE(NVDLL,1), C_SIZEOF(NVDLL(LBOUND(NVDLL,1))))
    ELSE
      SCALAR_INT_BUF(456) = 0
      CALL fstarpu_vector_data_register(handles(456), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(NVDLL)
    IF (ASSOCIATED(NBD)) THEN
      SCALAR_INT_BUF(457) = 1
      CALL fstarpu_vector_data_register(handles(457), 0, C_LOC(NBD(LBOUND(NBD,1))), SIZE(NBD,1), C_SIZEOF(NBD(LBOUND(NBD,1))))
    ELSE
      SCALAR_INT_BUF(457) = 0
      CALL fstarpu_vector_data_register(handles(457), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(NBD)
    IF (ASSOCIATED(NBDV)) THEN
      SCALAR_INT_BUF(458) = 1
      CALL fstarpu_matrix_data_register(handles(458), 0, C_LOC(NBDV(LBOUND(NBDV,1),LBOUND(NBDV,2))), SIZE(NBDV,1), SIZE(NBDV,1), SIZE(NBDV,2), C_SIZEOF(NBDV(LBOUND(NBDV,1),LBOUND(NBDV,2))))
    ELSE
      SCALAR_INT_BUF(458) = 0
      CALL fstarpu_matrix_data_register(handles(458), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(NBDV)
    IF (ASSOCIATED(NVELL)) THEN
      SCALAR_INT_BUF(459) = 1
      CALL fstarpu_vector_data_register(handles(459), 0, C_LOC(NVELL(LBOUND(NVELL,1))), SIZE(NVELL,1), C_SIZEOF(NVELL(LBOUND(NVELL,1))))
    ELSE
      SCALAR_INT_BUF(459) = 0
      CALL fstarpu_vector_data_register(handles(459), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(NVELL)
    IF (ASSOCIATED(NBVV)) THEN
      SCALAR_INT_BUF(460) = 1
      CALL fstarpu_matrix_data_register(handles(460), 0, C_LOC(NBVV(LBOUND(NBVV,1),LBOUND(NBVV,2))), SIZE(NBVV,1), SIZE(NBVV,1), SIZE(NBVV,2), C_SIZEOF(NBVV(LBOUND(NBVV,1),LBOUND(NBVV,2))))
    ELSE
      SCALAR_INT_BUF(460) = 0
      CALL fstarpu_matrix_data_register(handles(460), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(NBVV)
    IF (ASSOCIATED(NELED)) THEN
      SCALAR_INT_BUF(461) = 1
      CALL fstarpu_matrix_data_register(handles(461), 0, C_LOC(NELED(LBOUND(NELED,1),LBOUND(NELED,2))), SIZE(NELED,1), SIZE(NELED,1), SIZE(NELED,2), C_SIZEOF(NELED(LBOUND(NELED,1),LBOUND(NELED,2))))
    ELSE
      SCALAR_INT_BUF(461) = 0
      CALL fstarpu_matrix_data_register(handles(461), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(NELED)
    IF (ASSOCIATED(SEGTYPE)) THEN
      SCALAR_INT_BUF(462) = 1
      CALL fstarpu_vector_data_register(handles(462), 0, C_LOC(SEGTYPE(LBOUND(SEGTYPE,1))), SIZE(SEGTYPE,1), C_SIZEOF(SEGTYPE(LBOUND(SEGTYPE,1))))
    ELSE
      SCALAR_INT_BUF(462) = 0
      CALL fstarpu_vector_data_register(handles(462), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(SEGTYPE)
    IF (ASSOCIATED(NOT_AN_EDGE)) THEN
      SCALAR_INT_BUF(463) = 1
      CALL fstarpu_vector_data_register(handles(463), 0, C_LOC(NOT_AN_EDGE(LBOUND(NOT_AN_EDGE,1))), SIZE(NOT_AN_EDGE,1), C_SIZEOF(NOT_AN_EDGE(LBOUND(NOT_AN_EDGE,1))))
    ELSE
      SCALAR_INT_BUF(463) = 0
      CALL fstarpu_vector_data_register(handles(463), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(NOT_AN_EDGE)
    IF (ASSOCIATED(WEIR_BUDDY_NODE)) THEN
      SCALAR_INT_BUF(464) = 1
      CALL fstarpu_matrix_data_register(handles(464), 0, C_LOC(WEIR_BUDDY_NODE(LBOUND(WEIR_BUDDY_NODE,1),LBOUND(WEIR_BUDDY_NODE,2))), SIZE(WEIR_BUDDY_NODE,1), SIZE(WEIR_BUDDY_NODE,1), SIZE(WEIR_BUDDY_NODE,2), C_SIZEOF(WEIR_BUDDY_NODE(LBOUND(WEIR_BUDDY_NODE,1),LBOUND(WEIR_BUDDY_NODE,2))))
    ELSE
      SCALAR_INT_BUF(464) = 0
      CALL fstarpu_matrix_data_register(handles(464), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(WEIR_BUDDY_NODE)
    IF (ASSOCIATED(ONE_OR_TWO)) THEN
      SCALAR_INT_BUF(465) = 1
      CALL fstarpu_vector_data_register(handles(465), 0, C_LOC(ONE_OR_TWO(LBOUND(ONE_OR_TWO,1))), SIZE(ONE_OR_TWO,1), C_SIZEOF(ONE_OR_TWO(LBOUND(ONE_OR_TWO,1))))
    ELSE
      SCALAR_INT_BUF(465) = 0
      CALL fstarpu_vector_data_register(handles(465), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(ONE_OR_TWO)
    IF (ASSOCIATED(EDFLG)) THEN
      SCALAR_INT_BUF(466) = 1
      CALL fstarpu_matrix_data_register(handles(466), 0, C_LOC(EDFLG(LBOUND(EDFLG,1),LBOUND(EDFLG,2))), SIZE(EDFLG,1), SIZE(EDFLG,1), SIZE(EDFLG,2), C_SIZEOF(EDFLG(LBOUND(EDFLG,1),LBOUND(EDFLG,2))))
    ELSE
      SCALAR_INT_BUF(466) = 0
      CALL fstarpu_matrix_data_register(handles(466), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(EDFLG)
    IF (ASSOCIATED(BARLANHTR)) THEN
      SCALAR_INT_BUF(467) = 1
      CALL fstarpu_vector_data_register(handles(467), 0, C_LOC(BARLANHTR(LBOUND(BARLANHTR,1))), SIZE(BARLANHTR,1), C_SIZEOF(BARLANHTR(LBOUND(BARLANHTR,1))))
    ELSE
      SCALAR_INT_BUF(467) = 0
      CALL fstarpu_vector_data_register(handles(467), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(BARLANHTR)
    IF (ASSOCIATED(BARLANCFSPR)) THEN
      SCALAR_INT_BUF(468) = 1
      CALL fstarpu_vector_data_register(handles(468), 0, C_LOC(BARLANCFSPR(LBOUND(BARLANCFSPR,1))), SIZE(BARLANCFSPR,1), C_SIZEOF(BARLANCFSPR(LBOUND(BARLANCFSPR,1))))
    ELSE
      SCALAR_INT_BUF(468) = 0
      CALL fstarpu_vector_data_register(handles(468), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(BARLANCFSPR)
    IF (ASSOCIATED(BARINHTR)) THEN
      SCALAR_INT_BUF(469) = 1
      CALL fstarpu_vector_data_register(handles(469), 0, C_LOC(BARINHTR(LBOUND(BARINHTR,1))), SIZE(BARINHTR,1), C_SIZEOF(BARINHTR(LBOUND(BARINHTR,1))))
    ELSE
      SCALAR_INT_BUF(469) = 0
      CALL fstarpu_vector_data_register(handles(469), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(BARINHTR)
    IF (ASSOCIATED(BARINCFSBR)) THEN
      SCALAR_INT_BUF(470) = 1
      CALL fstarpu_vector_data_register(handles(470), 0, C_LOC(BARINCFSBR(LBOUND(BARINCFSBR,1))), SIZE(BARINCFSBR,1), C_SIZEOF(BARINCFSBR(LBOUND(BARINCFSBR,1))))
    ELSE
      SCALAR_INT_BUF(470) = 0
      CALL fstarpu_vector_data_register(handles(470), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(BARINCFSBR)
    IF (ASSOCIATED(BARINCFSPR)) THEN
      SCALAR_INT_BUF(471) = 1
      CALL fstarpu_vector_data_register(handles(471), 0, C_LOC(BARINCFSPR(LBOUND(BARINCFSPR,1))), SIZE(BARINCFSPR,1), C_SIZEOF(BARINCFSPR(LBOUND(BARINCFSPR,1))))
    ELSE
      SCALAR_INT_BUF(471) = 0
      CALL fstarpu_vector_data_register(handles(471), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(BARINCFSPR)
    IF (ASSOCIATED(PIPEHTR)) THEN
      SCALAR_INT_BUF(472) = 1
      CALL fstarpu_vector_data_register(handles(472), 0, C_LOC(PIPEHTR(LBOUND(PIPEHTR,1))), SIZE(PIPEHTR,1), C_SIZEOF(PIPEHTR(LBOUND(PIPEHTR,1))))
    ELSE
      SCALAR_INT_BUF(472) = 0
      CALL fstarpu_vector_data_register(handles(472), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(PIPEHTR)
    IF (ASSOCIATED(PIPECOEFR)) THEN
      SCALAR_INT_BUF(473) = 1
      CALL fstarpu_vector_data_register(handles(473), 0, C_LOC(PIPECOEFR(LBOUND(PIPECOEFR,1))), SIZE(PIPECOEFR,1), C_SIZEOF(PIPECOEFR(LBOUND(PIPECOEFR,1))))
    ELSE
      SCALAR_INT_BUF(473) = 0
      CALL fstarpu_vector_data_register(handles(473), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(PIPECOEFR)
    IF (ASSOCIATED(PIPEDIAMR)) THEN
      SCALAR_INT_BUF(474) = 1
      CALL fstarpu_vector_data_register(handles(474), 0, C_LOC(PIPEDIAMR(LBOUND(PIPEDIAMR,1))), SIZE(PIPEDIAMR,1), C_SIZEOF(PIPEDIAMR(LBOUND(PIPEDIAMR,1))))
    ELSE
      SCALAR_INT_BUF(474) = 0
      CALL fstarpu_vector_data_register(handles(474), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(PIPEDIAMR)
    IF (ASSOCIATED(BARLANHT)) THEN
      SCALAR_INT_BUF(475) = 1
      CALL fstarpu_vector_data_register(handles(475), 0, C_LOC(BARLANHT(LBOUND(BARLANHT,1))), SIZE(BARLANHT,1), C_SIZEOF(BARLANHT(LBOUND(BARLANHT,1))))
    ELSE
      SCALAR_INT_BUF(475) = 0
      CALL fstarpu_vector_data_register(handles(475), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(BARLANHT)
    IF (ASSOCIATED(BARLANCFSP)) THEN
      SCALAR_INT_BUF(476) = 1
      CALL fstarpu_vector_data_register(handles(476), 0, C_LOC(BARLANCFSP(LBOUND(BARLANCFSP,1))), SIZE(BARLANCFSP,1), C_SIZEOF(BARLANCFSP(LBOUND(BARLANCFSP,1))))
    ELSE
      SCALAR_INT_BUF(476) = 0
      CALL fstarpu_vector_data_register(handles(476), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(BARLANCFSP)
    IF (ASSOCIATED(FFF)) THEN
      SCALAR_INT_BUF(477) = 1
      CALL fstarpu_vector_data_register(handles(477), 0, C_LOC(FFF(LBOUND(FFF,1))), SIZE(FFF,1), C_SIZEOF(FFF(LBOUND(FFF,1))))
    ELSE
      SCALAR_INT_BUF(477) = 0
      CALL fstarpu_vector_data_register(handles(477), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(FFF)
    IF (ASSOCIATED(FFACE)) THEN
      SCALAR_INT_BUF(478) = 1
      CALL fstarpu_vector_data_register(handles(478), 0, C_LOC(FFACE(LBOUND(FFACE,1))), SIZE(FFACE,1), C_SIZEOF(FFACE(LBOUND(FFACE,1))))
    ELSE
      SCALAR_INT_BUF(478) = 0
      CALL fstarpu_vector_data_register(handles(478), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(FFACE)
    IF (ASSOCIATED(BTRAN3)) THEN
      SCALAR_INT_BUF(479) = 1
      CALL fstarpu_vector_data_register(handles(479), 0, C_LOC(BTRAN3(LBOUND(BTRAN3,1))), SIZE(BTRAN3,1), C_SIZEOF(BTRAN3(LBOUND(BTRAN3,1))))
    ELSE
      SCALAR_INT_BUF(479) = 0
      CALL fstarpu_vector_data_register(handles(479), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(BTRAN3)
    IF (ASSOCIATED(BTRAN4)) THEN
      SCALAR_INT_BUF(480) = 1
      CALL fstarpu_vector_data_register(handles(480), 0, C_LOC(BTRAN4(LBOUND(BTRAN4,1))), SIZE(BTRAN4,1), C_SIZEOF(BTRAN4(LBOUND(BTRAN4,1))))
    ELSE
      SCALAR_INT_BUF(480) = 0
      CALL fstarpu_vector_data_register(handles(480), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(BTRAN4)
    IF (ASSOCIATED(BTRAN5)) THEN
      SCALAR_INT_BUF(481) = 1
      CALL fstarpu_vector_data_register(handles(481), 0, C_LOC(BTRAN5(LBOUND(BTRAN5,1))), SIZE(BTRAN5,1), C_SIZEOF(BTRAN5(LBOUND(BTRAN5,1))))
    ELSE
      SCALAR_INT_BUF(481) = 0
      CALL fstarpu_vector_data_register(handles(481), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(BTRAN5)
    IF (ASSOCIATED(BTRAN6)) THEN
      SCALAR_INT_BUF(482) = 1
      CALL fstarpu_vector_data_register(handles(482), 0, C_LOC(BTRAN6(LBOUND(BTRAN6,1))), SIZE(BTRAN6,1), C_SIZEOF(BTRAN6(LBOUND(BTRAN6,1))))
    ELSE
      SCALAR_INT_BUF(482) = 0
      CALL fstarpu_vector_data_register(handles(482), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(BTRAN6)
    IF (ASSOCIATED(BTRAN7)) THEN
      SCALAR_INT_BUF(483) = 1
      CALL fstarpu_vector_data_register(handles(483), 0, C_LOC(BTRAN7(LBOUND(BTRAN7,1))), SIZE(BTRAN7,1), C_SIZEOF(BTRAN7(LBOUND(BTRAN7,1))))
    ELSE
      SCALAR_INT_BUF(483) = 0
      CALL fstarpu_vector_data_register(handles(483), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(BTRAN7)
    IF (ASSOCIATED(BTRAN8)) THEN
      SCALAR_INT_BUF(484) = 1
      CALL fstarpu_vector_data_register(handles(484), 0, C_LOC(BTRAN8(LBOUND(BTRAN8,1))), SIZE(BTRAN8,1), C_SIZEOF(BTRAN8(LBOUND(BTRAN8,1))))
    ELSE
      SCALAR_INT_BUF(484) = 0
      CALL fstarpu_vector_data_register(handles(484), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(BTRAN8)
    IF (ASSOCIATED(BARINHT)) THEN
      SCALAR_INT_BUF(485) = 1
      CALL fstarpu_vector_data_register(handles(485), 0, C_LOC(BARINHT(LBOUND(BARINHT,1))), SIZE(BARINHT,1), C_SIZEOF(BARINHT(LBOUND(BARINHT,1))))
    ELSE
      SCALAR_INT_BUF(485) = 0
      CALL fstarpu_vector_data_register(handles(485), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(BARINHT)
    IF (ASSOCIATED(BARINCFSB)) THEN
      SCALAR_INT_BUF(486) = 1
      CALL fstarpu_vector_data_register(handles(486), 0, C_LOC(BARINCFSB(LBOUND(BARINCFSB,1))), SIZE(BARINCFSB,1), C_SIZEOF(BARINCFSB(LBOUND(BARINCFSB,1))))
    ELSE
      SCALAR_INT_BUF(486) = 0
      CALL fstarpu_vector_data_register(handles(486), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(BARINCFSB)
    IF (ASSOCIATED(BARINCFSP)) THEN
      SCALAR_INT_BUF(487) = 1
      CALL fstarpu_vector_data_register(handles(487), 0, C_LOC(BARINCFSP(LBOUND(BARINCFSP,1))), SIZE(BARINCFSP,1), C_SIZEOF(BARINCFSP(LBOUND(BARINCFSP,1))))
    ELSE
      SCALAR_INT_BUF(487) = 0
      CALL fstarpu_vector_data_register(handles(487), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(BARINCFSP)
    IF (ASSOCIATED(PIPEHT)) THEN
      SCALAR_INT_BUF(488) = 1
      CALL fstarpu_vector_data_register(handles(488), 0, C_LOC(PIPEHT(LBOUND(PIPEHT,1))), SIZE(PIPEHT,1), C_SIZEOF(PIPEHT(LBOUND(PIPEHT,1))))
    ELSE
      SCALAR_INT_BUF(488) = 0
      CALL fstarpu_vector_data_register(handles(488), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(PIPEHT)
    IF (ASSOCIATED(PIPECOEF)) THEN
      SCALAR_INT_BUF(489) = 1
      CALL fstarpu_vector_data_register(handles(489), 0, C_LOC(PIPECOEF(LBOUND(PIPECOEF,1))), SIZE(PIPECOEF,1), C_SIZEOF(PIPECOEF(LBOUND(PIPECOEF,1))))
    ELSE
      SCALAR_INT_BUF(489) = 0
      CALL fstarpu_vector_data_register(handles(489), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(PIPECOEF)
    IF (ASSOCIATED(PIPEDIAM)) THEN
      SCALAR_INT_BUF(490) = 1
      CALL fstarpu_vector_data_register(handles(490), 0, C_LOC(PIPEDIAM(LBOUND(PIPEDIAM,1))), SIZE(PIPEDIAM,1), C_SIZEOF(PIPEDIAM(LBOUND(PIPEDIAM,1))))
    ELSE
      SCALAR_INT_BUF(490) = 0
      CALL fstarpu_vector_data_register(handles(490), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(PIPEDIAM)
    IF (ASSOCIATED(RBARWL1AVG)) THEN
      SCALAR_INT_BUF(491) = 1
      CALL fstarpu_vector_data_register(handles(491), 0, C_LOC(RBARWL1AVG(LBOUND(RBARWL1AVG,1))), SIZE(RBARWL1AVG,1), C_SIZEOF(RBARWL1AVG(LBOUND(RBARWL1AVG,1))))
    ELSE
      SCALAR_INT_BUF(491) = 0
      CALL fstarpu_vector_data_register(handles(491), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(RBARWL1AVG)
    IF (ASSOCIATED(RBARWL2AVG)) THEN
      SCALAR_INT_BUF(492) = 1
      CALL fstarpu_vector_data_register(handles(492), 0, C_LOC(RBARWL2AVG(LBOUND(RBARWL2AVG,1))), SIZE(RBARWL2AVG,1), C_SIZEOF(RBARWL2AVG(LBOUND(RBARWL2AVG,1))))
    ELSE
      SCALAR_INT_BUF(492) = 0
      CALL fstarpu_vector_data_register(handles(492), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(RBARWL2AVG)
    IF (ASSOCIATED(ELEXLEN)) THEN
      SCALAR_INT_BUF(493) = 1
      CALL fstarpu_matrix_data_register(handles(493), 0, C_LOC(ELEXLEN(LBOUND(ELEXLEN,1),LBOUND(ELEXLEN,2))), SIZE(ELEXLEN,1), SIZE(ELEXLEN,1), SIZE(ELEXLEN,2), C_SIZEOF(ELEXLEN(LBOUND(ELEXLEN,1),LBOUND(ELEXLEN,2))))
    ELSE
      SCALAR_INT_BUF(493) = 0
      CALL fstarpu_matrix_data_register(handles(493), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(ELEXLEN)
    IF (ASSOCIATED(IBCONN)) THEN
      SCALAR_INT_BUF(494) = 1
      CALL fstarpu_vector_data_register(handles(494), 0, C_LOC(IBCONN(LBOUND(IBCONN,1))), SIZE(IBCONN,1), C_SIZEOF(IBCONN(LBOUND(IBCONN,1))))
    ELSE
      SCALAR_INT_BUF(494) = 0
      CALL fstarpu_vector_data_register(handles(494), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(IBCONN)
    IF (ASSOCIATED(IBCONNR)) THEN
      SCALAR_INT_BUF(495) = 1
      CALL fstarpu_vector_data_register(handles(495), 0, C_LOC(IBCONNR(LBOUND(IBCONNR,1))), SIZE(IBCONNR,1), C_SIZEOF(IBCONNR(LBOUND(IBCONNR,1))))
    ELSE
      SCALAR_INT_BUF(495) = 0
      CALL fstarpu_vector_data_register(handles(495), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(IBCONNR)
    IF (ASSOCIATED(NTRAN1)) THEN
      SCALAR_INT_BUF(496) = 1
      CALL fstarpu_vector_data_register(handles(496), 0, C_LOC(NTRAN1(LBOUND(NTRAN1,1))), SIZE(NTRAN1,1), C_SIZEOF(NTRAN1(LBOUND(NTRAN1,1))))
    ELSE
      SCALAR_INT_BUF(496) = 0
      CALL fstarpu_vector_data_register(handles(496), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(NTRAN1)
    IF (ASSOCIATED(NTRAN2)) THEN
      SCALAR_INT_BUF(497) = 1
      CALL fstarpu_vector_data_register(handles(497), 0, C_LOC(NTRAN2(LBOUND(NTRAN2,1))), SIZE(NTRAN2,1), C_SIZEOF(NTRAN2(LBOUND(NTRAN2,1))))
    ELSE
      SCALAR_INT_BUF(497) = 0
      CALL fstarpu_vector_data_register(handles(497), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(NTRAN2)
    IF (ASSOCIATED(BK)) THEN
      SCALAR_INT_BUF(498) = 1
      CALL fstarpu_vector_data_register(handles(498), 0, C_LOC(BK(LBOUND(BK,1))), SIZE(BK,1), C_SIZEOF(BK(LBOUND(BK,1))))
    ELSE
      SCALAR_INT_BUF(498) = 0
      CALL fstarpu_vector_data_register(handles(498), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(BK)
    IF (ASSOCIATED(BALPHA)) THEN
      SCALAR_INT_BUF(499) = 1
      CALL fstarpu_vector_data_register(handles(499), 0, C_LOC(BALPHA(LBOUND(BALPHA,1))), SIZE(BALPHA,1), C_SIZEOF(BALPHA(LBOUND(BALPHA,1))))
    ELSE
      SCALAR_INT_BUF(499) = 0
      CALL fstarpu_vector_data_register(handles(499), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(BALPHA)
    IF (ASSOCIATED(BDELX)) THEN
      SCALAR_INT_BUF(500) = 1
      CALL fstarpu_vector_data_register(handles(500), 0, C_LOC(BDELX(LBOUND(BDELX,1))), SIZE(BDELX,1), C_SIZEOF(BDELX(LBOUND(BDELX,1))))
    ELSE
      SCALAR_INT_BUF(500) = 0
      CALL fstarpu_vector_data_register(handles(500), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(BDELX)
    IF (ASSOCIATED(NBNNUM)) THEN
      SCALAR_INT_BUF(501) = 1
      CALL fstarpu_vector_data_register(handles(501), 0, C_LOC(NBNNUM(LBOUND(NBNNUM,1))), SIZE(NBNNUM,1), C_SIZEOF(NBNNUM(LBOUND(NBNNUM,1))))
    ELSE
      SCALAR_INT_BUF(501) = 0
      CALL fstarpu_vector_data_register(handles(501), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(NBNNUM)
    IF (ASSOCIATED(AMIG)) THEN
      SCALAR_INT_BUF(502) = 1
      CALL fstarpu_vector_data_register(handles(502), 0, C_LOC(AMIG(LBOUND(AMIG,1))), SIZE(AMIG,1), C_SIZEOF(AMIG(LBOUND(AMIG,1))))
    ELSE
      SCALAR_INT_BUF(502) = 0
      CALL fstarpu_vector_data_register(handles(502), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(AMIG)
    IF (ASSOCIATED(AMIGT)) THEN
      SCALAR_INT_BUF(503) = 1
      CALL fstarpu_vector_data_register(handles(503), 0, C_LOC(AMIGT(LBOUND(AMIGT,1))), SIZE(AMIGT,1), C_SIZEOF(AMIGT(LBOUND(AMIGT,1))))
    ELSE
      SCALAR_INT_BUF(503) = 0
      CALL fstarpu_vector_data_register(handles(503), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(AMIGT)
    IF (ASSOCIATED(FAMIG)) THEN
      SCALAR_INT_BUF(504) = 1
      CALL fstarpu_vector_data_register(handles(504), 0, C_LOC(FAMIG(LBOUND(FAMIG,1))), SIZE(FAMIG,1), C_SIZEOF(FAMIG(LBOUND(FAMIG,1))))
    ELSE
      SCALAR_INT_BUF(504) = 0
      CALL fstarpu_vector_data_register(handles(504), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(FAMIG)
    IF (ASSOCIATED(PER)) THEN
      SCALAR_INT_BUF(505) = 1
      CALL fstarpu_vector_data_register(handles(505), 0, C_LOC(PER(LBOUND(PER,1))), SIZE(PER,1), C_SIZEOF(PER(LBOUND(PER,1))))
    ELSE
      SCALAR_INT_BUF(505) = 0
      CALL fstarpu_vector_data_register(handles(505), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(PER)
    IF (ASSOCIATED(PERT)) THEN
      SCALAR_INT_BUF(506) = 1
      CALL fstarpu_vector_data_register(handles(506), 0, C_LOC(PERT(LBOUND(PERT,1))), SIZE(PERT,1), C_SIZEOF(PERT(LBOUND(PERT,1))))
    ELSE
      SCALAR_INT_BUF(506) = 0
      CALL fstarpu_vector_data_register(handles(506), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(PERT)
    IF (ASSOCIATED(FPER)) THEN
      SCALAR_INT_BUF(507) = 1
      CALL fstarpu_vector_data_register(handles(507), 0, C_LOC(FPER(LBOUND(FPER,1))), SIZE(FPER,1), C_SIZEOF(FPER(LBOUND(FPER,1))))
    ELSE
      SCALAR_INT_BUF(507) = 0
      CALL fstarpu_vector_data_register(handles(507), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(FPER)
    IF (ASSOCIATED(FREQ)) THEN
      SCALAR_INT_BUF(508) = 1
      CALL fstarpu_vector_data_register(handles(508), 0, C_LOC(FREQ(LBOUND(FREQ,1))), SIZE(FREQ,1), C_SIZEOF(FREQ(LBOUND(FREQ,1))))
    ELSE
      SCALAR_INT_BUF(508) = 0
      CALL fstarpu_vector_data_register(handles(508), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(FREQ)
    IF (ASSOCIATED(FF)) THEN
      SCALAR_INT_BUF(509) = 1
      CALL fstarpu_vector_data_register(handles(509), 0, C_LOC(FF(LBOUND(FF,1))), SIZE(FF,1), C_SIZEOF(FF(LBOUND(FF,1))))
    ELSE
      SCALAR_INT_BUF(509) = 0
      CALL fstarpu_vector_data_register(handles(509), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(FF)
    IF (ASSOCIATED(FACE)) THEN
      SCALAR_INT_BUF(510) = 1
      CALL fstarpu_vector_data_register(handles(510), 0, C_LOC(FACE(LBOUND(FACE,1))), SIZE(FACE,1), C_SIZEOF(FACE(LBOUND(FACE,1))))
    ELSE
      SCALAR_INT_BUF(510) = 0
      CALL fstarpu_vector_data_register(handles(510), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(FACE)
    IF (ASSOCIATED(SLAM)) THEN
      SCALAR_INT_BUF(511) = 1
      CALL fstarpu_vector_data_register(handles(511), 0, C_LOC(SLAM(LBOUND(SLAM,1))), SIZE(SLAM,1), C_SIZEOF(SLAM(LBOUND(SLAM,1))))
    ELSE
      SCALAR_INT_BUF(511) = 0
      CALL fstarpu_vector_data_register(handles(511), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(SLAM)
    IF (ASSOCIATED(SFEA)) THEN
      SCALAR_INT_BUF(512) = 1
      CALL fstarpu_vector_data_register(handles(512), 0, C_LOC(SFEA(LBOUND(SFEA,1))), SIZE(SFEA,1), C_SIZEOF(SFEA(LBOUND(SFEA,1))))
    ELSE
      SCALAR_INT_BUF(512) = 0
      CALL fstarpu_vector_data_register(handles(512), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(SFEA)
    IF (ASSOCIATED(X)) THEN
      SCALAR_INT_BUF(513) = 1
      CALL fstarpu_vector_data_register(handles(513), 0, C_LOC(X(LBOUND(X,1))), SIZE(X,1), C_SIZEOF(X(LBOUND(X,1))))
    ELSE
      SCALAR_INT_BUF(513) = 0
      CALL fstarpu_vector_data_register(handles(513), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(X)
    IF (ASSOCIATED(Y)) THEN
      SCALAR_INT_BUF(514) = 1
      CALL fstarpu_vector_data_register(handles(514), 0, C_LOC(Y(LBOUND(Y,1))), SIZE(Y,1), C_SIZEOF(Y(LBOUND(Y,1))))
    ELSE
      SCALAR_INT_BUF(514) = 0
      CALL fstarpu_vector_data_register(handles(514), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))
    END IF
    NULLIFY(Y)
    SCALAR_INT_BUF(515) = MNPROC
    SCALAR_INT_BUF(516) = MNE
    SCALAR_INT_BUF(517) = MNP
    SCALAR_INT_BUF(518) = MNEI
    SCALAR_INT_BUF(519) = MNOPE
    SCALAR_INT_BUF(520) = MNETA
    SCALAR_INT_BUF(521) = MNBOU
    SCALAR_INT_BUF(522) = MNVEL
    SCALAR_INT_BUF(523) = MNTIF
    SCALAR_INT_BUF(524) = MNBFR
    SCALAR_INT_BUF(525) = MNFFR
    SCALAR_INT_BUF(526) = MNSTAE
    SCALAR_INT_BUF(527) = MNSTAV
    SCALAR_INT_BUF(528) = MNSTAC
    SCALAR_INT_BUF(529) = MNSTAM
    SCALAR_INT_BUF(530) = MNHARF
    SCALAR_INT_BUF(531) = layers
    SCALAR_INT_BUF(532) = MNNDEL
    SCALAR_INT_BUF(533) = MYPROC
    SCALAR_INT_BUF(534) = LNAME
    SCALAR_INT_BUF(535) = rainfall
    SCALAR_INT_BUF(536) = DGFLAG
    SCALAR_INT_BUF(537) = DGHOT
    SCALAR_INT_BUF(538) = DGHOTSPOOL
    SCALAR_INT_BUF(539) = DOF
    SCALAR_INT_BUF(540) = dofh
    SCALAR_INT_BUF(541) = dofl
    SCALAR_INT_BUF(542) = dofx
    SCALAR_INT_BUF(543) = EL
    SCALAR_INT_BUF(544) = MNES
    SCALAR_INT_BUF(545) = artdif
    SCALAR_INT_BUF(546) = tune_by_hand
    SCALAR_INT_BUF(547) = IRK
    SCALAR_INT_BUF(548) = J1
    SCALAR_INT_BUF(549) = J2
    SCALAR_INT_BUF(550) = J3
    SCALAR_INT_BUF(551) = negp_fixed
    SCALAR_INT_BUF(552) = nagp_fixed
    SCALAR_INT_BUF(553) = NAGP(1)
    SCALAR_INT_BUF(554) = NAGP(2)
    SCALAR_INT_BUF(555) = NAGP(3)
    SCALAR_INT_BUF(556) = NAGP(4)
    SCALAR_INT_BUF(557) = NAGP(5)
    SCALAR_INT_BUF(558) = NAGP(6)
    SCALAR_INT_BUF(559) = NAGP(7)
    SCALAR_INT_BUF(560) = NAGP(8)
    SCALAR_INT_BUF(561) = NCHECK(1)
    SCALAR_INT_BUF(562) = NCHECK(2)
    SCALAR_INT_BUF(563) = NCHECK(3)
    SCALAR_INT_BUF(564) = NCHECK(4)
    SCALAR_INT_BUF(565) = NCHECK(5)
    SCALAR_INT_BUF(566) = NCHECK(6)
    SCALAR_INT_BUF(567) = NCHECK(7)
    SCALAR_INT_BUF(568) = NCHECK(8)
    SCALAR_INT_BUF(569) = NEGP(1)
    SCALAR_INT_BUF(570) = NEGP(2)
    SCALAR_INT_BUF(571) = NEGP(3)
    SCALAR_INT_BUF(572) = NEGP(4)
    SCALAR_INT_BUF(573) = NEGP(5)
    SCALAR_INT_BUF(574) = NEGP(6)
    SCALAR_INT_BUF(575) = NEGP(7)
    SCALAR_INT_BUF(576) = NEGP(8)
    SCALAR_INT_BUF(577) = NEDGES
    SCALAR_INT_BUF(578) = NRK
    SCALAR_INT_BUF(579) = NIEDS
    SCALAR_INT_BUF(580) = NLEDS
    SCALAR_INT_BUF(581) = NEEDS
    SCALAR_INT_BUF(582) = NFEDS
    SCALAR_INT_BUF(583) = NREDS
    SCALAR_INT_BUF(584) = NEBEDS
    SCALAR_INT_BUF(585) = NIBEDS
    SCALAR_INT_BUF(586) = NIBSEG
    SCALAR_INT_BUF(587) = NEBSEG
    SCALAR_INT_BUF(588) = MNED
    SCALAR_INT_BUF(589) = MNLED
    SCALAR_INT_BUF(590) = MNSED
    SCALAR_INT_BUF(591) = MNRAED
    SCALAR_INT_BUF(592) = MNRIED
    SCALAR_INT_BUF(593) = MODAL_IC
    SCALAR_INT_BUF(594) = P_READ
    SCALAR_INT_BUF(595) = P_READ2
    SCALAR_INT_BUF(596) = SLOPEFLAG
    SCALAR_INT_BUF(597) = test_el
    SCALAR_INT_BUF(598) = FLUXTYPE
    SCALAR_INT_BUF(599) = RK_STAGE
    SCALAR_INT_BUF(600) = RK_ORDER
    SCALAR_INT_BUF(601) = padapt
    SCALAR_INT_BUF(602) = pflag
    SCALAR_INT_BUF(603) = pl
    SCALAR_INT_BUF(604) = ph
    SCALAR_INT_BUF(605) = px
    SCALAR_INT_BUF(606) = lebesgueP
    SCALAR_INT_BUF(607) = gflag
    SCALAR_INT_BUF(608) = pa
    SCALAR_INT_BUF(609) = iwrite
    SCALAR_INT_BUF(610) = lim_count
    SCALAR_INT_BUF(611) = lim_count_roll
    SCALAR_INT_BUF(612) = SEDFLAG
    SCALAR_INT_BUF(613) = MAXEL
    SCALAR_INT_BUF(614) = ELEM_ED
    SCALAR_INT_BUF(615) = NBOR_ED
    SCALAR_INT_BUF(616) = NBOR_EL
    SCALAR_INT_BUF(617) = ITDG
    SCALAR_INT_BUF(618) = ModetoNode
    SCALAR_INT_BUF(619) = tracer_flag
    SCALAR_INT_BUF(620) = chem_flag
    SCALAR_INT_BUF(621) = N1
    SCALAR_INT_BUF(622) = N2
    SCALAR_INT_BUF(623) = NO_NBORS
    SCALAR_INT_BUF(624) = NBOR
    SCALAR_INT_BUF(625) = SEDFLAG_W
    SCALAR_INT_BUF(626) = OPEN_INDEX
    SCALAR_INT_BUF(627) = DG_TO_CG
    SCALAR_INT_BUF(628) = NSCREEN_INC
    SCALAR_INT_BUF(629) = ScreenUnit
    SCALAR_INT_BUF(630) = FluxSettlingIT
    SCALAR_INT_BUF(631) = DGSWE
    SCALAR_INT_BUF(632) = EL_IN
    SCALAR_INT_BUF(633) = EL_EX
    SCALAR_INT_BUF(634) = SD_IN
    SCALAR_INT_BUF(635) = SD_EX
    SCALAR_INT_BUF(636) = EDGE(1)
    SCALAR_INT_BUF(637) = EDGE(2)
    SCALAR_INT_BUF(638) = EDGE(3)
    SCALAR_INT_BUF(639) = SIDE(1)
    SCALAR_INT_BUF(640) = SIDE(2)
    SCALAR_INT_BUF(641) = TESTPROBLEM
    SCALAR_INT_BUF(642) = NBPNODES
    SCALAR_INT_BUF(643) = NP
    SCALAR_INT_BUF(644) = NOLICA
    SCALAR_INT_BUF(645) = NOLIFA
    SCALAR_INT_BUF(646) = NSCREEN
    SCALAR_INT_BUF(647) = IHOT
    SCALAR_INT_BUF(648) = ICS
    SCALAR_INT_BUF(649) = FRW
    SCALAR_INT_BUF(650) = NODEDRYMIN
    SCALAR_INT_BUF(651) = NODEWETMIN
    SCALAR_INT_BUF(652) = IBTYPE
    SCALAR_INT_BUF(653) = ICK
    SCALAR_INT_BUF(654) = IDR
    SCALAR_INT_BUF(655) = IM
    SCALAR_INT_BUF(656) = IPRBI
    SCALAR_INT_BUF(657) = JGW
    SCALAR_INT_BUF(658) = JKI
    SCALAR_INT_BUF(659) = JME
    SCALAR_INT_BUF(660) = JNMM
    SCALAR_INT_BUF(661) = KMIN
    SCALAR_INT_BUF(662) = N3
    SCALAR_INT_BUF(663) = NABOUT
    SCALAR_INT_BUF(664) = NBFR
    SCALAR_INT_BUF(665) = NBOU
    SCALAR_INT_BUF(666) = NBVI
    SCALAR_INT_BUF(667) = NBVJ
    SCALAR_INT_BUF(668) = NCOR
    SCALAR_INT_BUF(669) = NE
    SCALAR_INT_BUF(670) = NE2
    SCALAR_INT_BUF(671) = NP2
    SCALAR_INT_BUF(672) = NEIMIN
    SCALAR_INT_BUF(673) = NEIMAX
    SCALAR_INT_BUF(674) = NETA
    SCALAR_INT_BUF(675) = NFFR
    SCALAR_INT_BUF(676) = NFLUXB
    SCALAR_INT_BUF(677) = NFLUXF
    SCALAR_INT_BUF(678) = NFLUXIB
    SCALAR_INT_BUF(679) = NFLUXRBC
    SCALAR_INT_BUF(680) = NFLUXIBP
    SCALAR_INT_BUF(681) = NPIPE
    SCALAR_INT_BUF(682) = NFOVER
    SCALAR_INT_BUF(683) = NHG
    SCALAR_INT_BUF(684) = NHY
    SCALAR_INT_BUF(685) = NOLICAT
    SCALAR_INT_BUF(686) = NOPE
    SCALAR_INT_BUF(687) = NOUTC
    SCALAR_INT_BUF(688) = NOUTE
    SCALAR_INT_BUF(689) = NSPOOLE
    SCALAR_INT_BUF(690) = NOUTV
    SCALAR_INT_BUF(691) = NSPOOLV
    SCALAR_INT_BUF(692) = NPRBI
    SCALAR_INT_BUF(693) = NRAMP
    SCALAR_INT_BUF(694) = NRS
    SCALAR_INT_BUF(695) = NSTAE
    SCALAR_INT_BUF(696) = NSTARTDRY
    SCALAR_INT_BUF(697) = NSTAV
    SCALAR_INT_BUF(698) = NT
    SCALAR_INT_BUF(699) = NTCYFE
    SCALAR_INT_BUF(700) = NTCYFV
    SCALAR_INT_BUF(701) = NTCYSE
    SCALAR_INT_BUF(702) = NTCYSV
    SCALAR_INT_BUF(703) = NTIF
    SCALAR_INT_BUF(704) = NTIP
    SCALAR_INT_BUF(705) = NTRSPE
    SCALAR_INT_BUF(706) = NTRSPV
    SCALAR_INT_BUF(707) = NVEL
    SCALAR_INT_BUF(708) = NVELEXT
    SCALAR_INT_BUF(709) = NVELME
    SCALAR_INT_BUF(710) = NWLAT
    SCALAR_INT_BUF(711) = NWLON
    SCALAR_INT_BUF(712) = NWS
    SCALAR_INT_BUF(713) = IBSTART
    SCALAR_INT_BUF(714) = ICHA
    SCALAR_INT_BUF(715) = ICSTP
    SCALAR_INT_BUF(716) = IDSETFLG
    SCALAR_INT_BUF(717) = IE
    SCALAR_INT_BUF(718) = IER
    SCALAR_INT_BUF(719) = IESTP
    SCALAR_INT_BUF(720) = IFNLCAT
    SCALAR_INT_BUF(721) = IFNLCT
    SCALAR_INT_BUF(722) = IFNLFA
    SCALAR_INT_BUF(723) = IFWIND
    SCALAR_INT_BUF(724) = IGCP
    SCALAR_INT_BUF(725) = IGEP
    SCALAR_INT_BUF(726) = IGPP
    SCALAR_INT_BUF(727) = IGVP
    SCALAR_INT_BUF(728) = IGWP
    SCALAR_INT_BUF(729) = IHABEG
    SCALAR_INT_BUF(730) = IGRadS
    SCALAR_INT_BUF(731) = IHOTSTP
    SCALAR_INT_BUF(732) = IHSFIL
    SCALAR_INT_BUF(733) = IJ
    SCALAR_INT_BUF(734) = ILUMP
    SCALAR_INT_BUF(735) = IMHS
    SCALAR_INT_BUF(736) = IPSTP
    SCALAR_INT_BUF(737) = IREFYR
    SCALAR_INT_BUF(738) = IREFMO
    SCALAR_INT_BUF(739) = IREFDAY
    SCALAR_INT_BUF(740) = IREFHR
    SCALAR_INT_BUF(741) = IREFMIN
    SCALAR_INT_BUF(742) = ISLDIA
    SCALAR_INT_BUF(743) = ITIME_A
    SCALAR_INT_BUF(744) = ITEMPSTP
    SCALAR_INT_BUF(745) = ITEST
    SCALAR_INT_BUF(746) = ITHS
    SCALAR_INT_BUF(747) = ITITER
    SCALAR_INT_BUF(748) = ITMAX
    SCALAR_INT_BUF(749) = IVSTP
    SCALAR_INT_BUF(750) = IWSTP
    SCALAR_INT_BUF(751) = IWTIME
    SCALAR_INT_BUF(752) = IWTIMEP
    SCALAR_INT_BUF(753) = IWYR
    SCALAR_INT_BUF(754) = J12
    SCALAR_INT_BUF(755) = J13
    SCALAR_INT_BUF(756) = J21
    SCALAR_INT_BUF(757) = J23
    SCALAR_INT_BUF(758) = J31
    SCALAR_INT_BUF(759) = J32
    SCALAR_INT_BUF(760) = JN
    SCALAR_INT_BUF(761) = KEMAX
    SCALAR_INT_BUF(762) = KVMAX
    SCALAR_INT_BUF(763) = LRC
    SCALAR_INT_BUF(764) = LUMPT
    SCALAR_INT_BUF(765) = MMAX
    SCALAR_INT_BUF(766) = MBW
    SCALAR_INT_BUF(767) = MDF
    SCALAR_INT_BUF(768) = MMIN
    SCALAR_INT_BUF(769) = NA
    SCALAR_INT_BUF(770) = NBDI
    SCALAR_INT_BUF(771) = NBDJ
    SCALAR_INT_BUF(772) = NBNCTOT
    SCALAR_INT_BUF(773) = NBW
    SCALAR_INT_BUF(774) = NC1
    SCALAR_INT_BUF(775) = NC2
    SCALAR_INT_BUF(776) = NC3
    SCALAR_INT_BUF(777) = NCBND
    SCALAR_INT_BUF(778) = NCELE
    SCALAR_INT_BUF(779) = NCI
    SCALAR_INT_BUF(780) = NCJ
    SCALAR_INT_BUF(781) = NCTOT
    SCALAR_INT_BUF(782) = NCYC
    SCALAR_INT_BUF(783) = NDRY
    SCALAR_INT_BUF(784) = NDSETSC
    SCALAR_INT_BUF(785) = NDSETSE
    SCALAR_INT_BUF(786) = NDSETSV
    SCALAR_INT_BUF(787) = NDSETSW
    SCALAR_INT_BUF(788) = NHSINC
    SCALAR_INT_BUF(789) = NHSTAR
    SCALAR_INT_BUF(790) = NM1
    SCALAR_INT_BUF(791) = NM123
    SCALAR_INT_BUF(792) = NM2
    SCALAR_INT_BUF(793) = NM3
    SCALAR_INT_BUF(794) = NMI1
    SCALAR_INT_BUF(795) = NMI2
    SCALAR_INT_BUF(796) = NMI3
    SCALAR_INT_BUF(797) = NMJ1
    SCALAR_INT_BUF(798) = NMJ2
    SCALAR_INT_BUF(799) = NMJ3
    SCALAR_INT_BUF(800) = NNBB
    SCALAR_INT_BUF(801) = NNBB1
    SCALAR_INT_BUF(802) = NNBB2
    SCALAR_INT_BUF(803) = NOUTGC
    SCALAR_INT_BUF(804) = NOUTGE
    SCALAR_INT_BUF(805) = NOUTGV
    SCALAR_INT_BUF(806) = NOUTGW
    SCALAR_INT_BUF(807) = NOUTM
    SCALAR_INT_BUF(808) = NSCOUC
    SCALAR_INT_BUF(809) = NSCOUE
    SCALAR_INT_BUF(810) = NSCOUGC
    SCALAR_INT_BUF(811) = NSCOUGE
    SCALAR_INT_BUF(812) = NSCOUGV
    SCALAR_INT_BUF(813) = NSCOUGW
    SCALAR_INT_BUF(814) = NSCOUM
    SCALAR_INT_BUF(815) = NSCOUV
    SCALAR_INT_BUF(816) = NSPOOLC
    SCALAR_INT_BUF(817) = NSPOOLGC
    SCALAR_INT_BUF(818) = NSPOOLGE
    SCALAR_INT_BUF(819) = NSPOOLGV
    SCALAR_INT_BUF(820) = NSPOOLGW
    SCALAR_INT_BUF(821) = NSPOOLM
    SCALAR_INT_BUF(822) = NSTAC
    SCALAR_INT_BUF(823) = NSTAM
    SCALAR_INT_BUF(824) = NTCYFC
    SCALAR_INT_BUF(825) = NTCYFGC
    SCALAR_INT_BUF(826) = NTCYFGE
    SCALAR_INT_BUF(827) = NTCYFGV
    SCALAR_INT_BUF(828) = NTCYFGW
    SCALAR_INT_BUF(829) = NTCYFM
    SCALAR_INT_BUF(830) = NTCYSC
    SCALAR_INT_BUF(831) = NTCYSGC
    SCALAR_INT_BUF(832) = NTCYSGE
    SCALAR_INT_BUF(833) = NTCYSGV
    SCALAR_INT_BUF(834) = NTCYSGW
    SCALAR_INT_BUF(835) = NTCYSM
    SCALAR_INT_BUF(836) = NTRSPC
    SCALAR_INT_BUF(837) = NTRSPM
    SCALAR_INT_BUF(838) = NUMITR
    SCALAR_INT_BUF(839) = NW
    SCALAR_INT_BUF(840) = NWET
    SCALAR_INT_BUF(841) = NWSEGWI
    SCALAR_INT_BUF(842) = NWSGGWI
    SCALAR_INT_BUF(843) = NCCHANGE
    SCALAR_INT_BUF(844) = IRAMPING
    IF (vertexslope) THEN
      SCALAR_INT_BUF(847) = 1
    ELSE
      SCALAR_INT_BUF(847) = 0
    END IF
    CALL fstarpu_vector_data_register(handles(515), 0, C_LOC(SCALAR_INT_BUF(1)), NUM_INT_SCALARS, C_SIZEOF(SCALAR_INT_BUF(1)))
    NULLIFY(SCALAR_INT_BUF)
    SCALAR_REAL_BUF(1) = C13
    SCALAR_REAL_BUF(2) = C16
    SCALAR_REAL_BUF(3) = diorism
    SCALAR_REAL_BUF(4) = porosity
    SCALAR_REAL_BUF(5) = SEVDM
    SCALAR_REAL_BUF(6) = DOT
    SCALAR_REAL_BUF(7) = DHB_X
    SCALAR_REAL_BUF(8) = DHB_Y
    SCALAR_REAL_BUF(9) = DPHIDX
    SCALAR_REAL_BUF(10) = DPHIDY
    SCALAR_REAL_BUF(11) = slimit
    SCALAR_REAL_BUF(12) = plimit
    SCALAR_REAL_BUF(13) = pflag2con1
    SCALAR_REAL_BUF(14) = pflag2con2
    SCALAR_REAL_BUF(15) = EFA_GP
    SCALAR_REAL_BUF(16) = EMO_GP
    SCALAR_REAL_BUF(17) = slimit1
    SCALAR_REAL_BUF(18) = slimit2
    SCALAR_REAL_BUF(19) = slimit3
    SCALAR_REAL_BUF(20) = EL_ANG
    SCALAR_REAL_BUF(21) = slimit4
    SCALAR_REAL_BUF(22) = bg_dif
    SCALAR_REAL_BUF(23) = trc_dif
    SCALAR_REAL_BUF(24) = slimit5
    SCALAR_REAL_BUF(25) = FG_L
    SCALAR_REAL_BUF(26) = l2er_global
    SCALAR_REAL_BUF(27) = temperg
    SCALAR_REAL_BUF(28) = slope_weight
    SCALAR_REAL_BUF(29) = H_TRI
    SCALAR_REAL_BUF(30) = MAG1
    SCALAR_REAL_BUF(31) = MAG2
    SCALAR_REAL_BUF(32) = kappa
    SCALAR_REAL_BUF(33) = s0
    SCALAR_REAL_BUF(34) = uniform_dif
    SCALAR_REAL_BUF(35) = RAMPDG
    SCALAR_REAL_BUF(36) = S1
    SCALAR_REAL_BUF(37) = S2
    SCALAR_REAL_BUF(38) = SAV
    SCALAR_REAL_BUF(39) = SOURCE_X
    SCALAR_REAL_BUF(40) = SOURCE_Y
    SCALAR_REAL_BUF(41) = TIMEDG
    SCALAR_REAL_BUF(42) = TIMEH_DG
    SCALAR_REAL_BUF(43) = TK
    SCALAR_REAL_BUF(44) = QNAM_GP
    SCALAR_REAL_BUF(45) = QNPH_GP
    SCALAR_REAL_BUF(46) = SL2_M
    SCALAR_REAL_BUF(47) = SL2_NYU
    SCALAR_REAL_BUF(48) = SL3_MD
    SCALAR_REAL_BUF(49) = EVMAvg
    SCALAR_REAL_BUF(50) = SEVDMAvg
    SCALAR_REAL_BUF(51) = UMAG
    SCALAR_REAL_BUF(52) = WSX_GP
    SCALAR_REAL_BUF(53) = WSY_GP
    SCALAR_REAL_BUF(54) = QMag_IN
    SCALAR_REAL_BUF(55) = QMag_EX
    SCALAR_REAL_BUF(56) = subphi_IN
    SCALAR_REAL_BUF(57) = subphi_EX
    SCALAR_REAL_BUF(58) = iota_EX
    SCALAR_REAL_BUF(59) = iota_IN
    SCALAR_REAL_BUF(60) = iota2_EX
    SCALAR_REAL_BUF(61) = iota2_IN
    SCALAR_REAL_BUF(62) = valx(1)
    SCALAR_REAL_BUF(63) = valx(2)
    SCALAR_REAL_BUF(64) = valx(3)
    SCALAR_REAL_BUF(65) = valx(4)
    SCALAR_REAL_BUF(66) = valy(1)
    SCALAR_REAL_BUF(67) = valy(2)
    SCALAR_REAL_BUF(68) = valy(3)
    SCALAR_REAL_BUF(69) = valy(4)
    SCALAR_REAL_BUF(70) = DRPSI(1)
    SCALAR_REAL_BUF(71) = DRPSI(2)
    SCALAR_REAL_BUF(72) = DRPSI(3)
    SCALAR_REAL_BUF(73) = DSPSI(1)
    SCALAR_REAL_BUF(74) = DSPSI(2)
    SCALAR_REAL_BUF(75) = DSPSI(3)
    SCALAR_REAL_BUF(76) = VEC1(1)
    SCALAR_REAL_BUF(77) = VEC1(2)
    SCALAR_REAL_BUF(78) = VEC2(1)
    SCALAR_REAL_BUF(79) = VEC2(2)
    SCALAR_REAL_BUF(80) = rhoAir
    SCALAR_REAL_BUF(81) = windReduction
    SCALAR_REAL_BUF(82) = one2ten
    SCALAR_REAL_BUF(83) = ten2one
    SCALAR_REAL_BUF(84) = WaveWindMultiplier
    SCALAR_REAL_BUF(85) = DEPAVG
    SCALAR_REAL_BUF(86) = DEPMAX
    SCALAR_REAL_BUF(87) = DEPMIN
    SCALAR_REAL_BUF(88) = AREA_SUM
    SCALAR_REAL_BUF(89) = CEN_SUM
    SCALAR_REAL_BUF(90) = NLEQ
    SCALAR_REAL_BUF(91) = LEQ
    SCALAR_REAL_BUF(92) = NLEQG
    SCALAR_REAL_BUF(93) = reaction_rate
    SCALAR_REAL_BUF(94) = FluxSettlingTime
    SCALAR_REAL_BUF(95) = RampExtFlux
    SCALAR_REAL_BUF(96) = DRampExtFlux
    SCALAR_REAL_BUF(97) = RampIntFlux
    SCALAR_REAL_BUF(98) = DRampIntFlux
    SCALAR_REAL_BUF(99) = RampElev
    SCALAR_REAL_BUF(100) = DRampElev
    SCALAR_REAL_BUF(101) = RampTip
    SCALAR_REAL_BUF(102) = DRampTip
    SCALAR_REAL_BUF(103) = RampMete
    SCALAR_REAL_BUF(104) = DRampMete
    SCALAR_REAL_BUF(105) = RampWRad
    SCALAR_REAL_BUF(106) = DRampWRad
    SCALAR_REAL_BUF(107) = FX_IN
    SCALAR_REAL_BUF(108) = FY_IN
    SCALAR_REAL_BUF(109) = GX_IN
    SCALAR_REAL_BUF(110) = GY_IN
    SCALAR_REAL_BUF(111) = HX_IN
    SCALAR_REAL_BUF(112) = HY_IN
    SCALAR_REAL_BUF(113) = FX_EX
    SCALAR_REAL_BUF(114) = FY_EX
    SCALAR_REAL_BUF(115) = GX_EX
    SCALAR_REAL_BUF(116) = GY_EX
    SCALAR_REAL_BUF(117) = HX_EX
    SCALAR_REAL_BUF(118) = HY_EX
    SCALAR_REAL_BUF(119) = F_AVG
    SCALAR_REAL_BUF(120) = G_AVG
    SCALAR_REAL_BUF(121) = H_AVG
    SCALAR_REAL_BUF(122) = JUMP(1)
    SCALAR_REAL_BUF(123) = JUMP(2)
    SCALAR_REAL_BUF(124) = JUMP(3)
    SCALAR_REAL_BUF(125) = JUMP(4)
    SCALAR_REAL_BUF(126) = HT_IN
    SCALAR_REAL_BUF(127) = HT_EX
    SCALAR_REAL_BUF(128) = UMag_IN
    SCALAR_REAL_BUF(129) = UMag_EX
    SCALAR_REAL_BUF(130) = ZE_ROE
    SCALAR_REAL_BUF(131) = QX_ROE
    SCALAR_REAL_BUF(132) = QY_ROE
    SCALAR_REAL_BUF(133) = bed_ROE
    SCALAR_REAL_BUF(134) = Q_N
    SCALAR_REAL_BUF(135) = Q_T
    SCALAR_REAL_BUF(136) = U_N
    SCALAR_REAL_BUF(137) = U_T
    SCALAR_REAL_BUF(138) = U_IN
    SCALAR_REAL_BUF(139) = U_EX
    SCALAR_REAL_BUF(140) = V_IN
    SCALAR_REAL_BUF(141) = V_EX
    SCALAR_REAL_BUF(142) = ZE_SUM
    SCALAR_REAL_BUF(143) = QX_SUM
    SCALAR_REAL_BUF(144) = QY_SUM
    SCALAR_REAL_BUF(145) = DG_MAX
    SCALAR_REAL_BUF(146) = DG_MIN
    SCALAR_REAL_BUF(147) = U_N_EXT
    SCALAR_REAL_BUF(148) = U_T_EXT
    SCALAR_REAL_BUF(149) = Q_N_INT
    SCALAR_REAL_BUF(150) = Q_T_INT
    SCALAR_REAL_BUF(151) = U_N_INT
    SCALAR_REAL_BUF(152) = U_T_INT
    SCALAR_REAL_BUF(153) = Q_N_EXT
    SCALAR_REAL_BUF(154) = Q_T_EXT
    SCALAR_REAL_BUF(155) = BX_INT
    SCALAR_REAL_BUF(156) = BY_INT
    SCALAR_REAL_BUF(157) = SOURCE_1
    SCALAR_REAL_BUF(158) = SOURCE_2
    SCALAR_REAL_BUF(159) = SOURCE_SUM
    SCALAR_REAL_BUF(160) = k_hat
    SCALAR_REAL_BUF(161) = FRIC_AVG
    SCALAR_REAL_BUF(162) = DP_MID
    SCALAR_REAL_BUF(163) = i_hat
    SCALAR_REAL_BUF(164) = j_hat
    SCALAR_REAL_BUF(165) = INFLOW_ZE
    SCALAR_REAL_BUF(166) = INFLOW_QX
    SCALAR_REAL_BUF(167) = INFLOW_QY
    SCALAR_REAL_BUF(168) = H_LEN
    SCALAR_REAL_BUF(169) = INFLOW_LEN
    SCALAR_REAL_BUF(170) = ZE_NORM
    SCALAR_REAL_BUF(171) = QX_NORM
    SCALAR_REAL_BUF(172) = QY_NORM
    SCALAR_REAL_BUF(173) = ZE_DECT
    SCALAR_REAL_BUF(174) = QX_DECT
    SCALAR_REAL_BUF(175) = QY_DECT
    SCALAR_REAL_BUF(176) = POAN
    SCALAR_REAL_BUF(177) = Fr
    SCALAR_REAL_BUF(178) = FRICBP
    SCALAR_REAL_BUF(179) = STATIM
    SCALAR_REAL_BUF(180) = REFTIM
    SCALAR_REAL_BUF(181) = TIME_A
    SCALAR_REAL_BUF(182) = DT
    SCALAR_REAL_BUF(183) = DTDP
    SCALAR_REAL_BUF(184) = TIMEH
    SCALAR_REAL_BUF(185) = vdtdp
    SCALAR_REAL_BUF(186) = cfl_max
    SCALAR_REAL_BUF(187) = AVGXY
    SCALAR_REAL_BUF(188) = DIF1R
    SCALAR_REAL_BUF(189) = DIF2R
    SCALAR_REAL_BUF(190) = DIF3R
    SCALAR_REAL_BUF(191) = AEMIN
    SCALAR_REAL_BUF(192) = AE
    SCALAR_REAL_BUF(193) = AA
    SCALAR_REAL_BUF(194) = A1
    SCALAR_REAL_BUF(195) = A2
    SCALAR_REAL_BUF(196) = A3
    SCALAR_REAL_BUF(197) = X1
    SCALAR_REAL_BUF(198) = X2
    SCALAR_REAL_BUF(199) = X3
    SCALAR_REAL_BUF(200) = X4
    SCALAR_REAL_BUF(201) = Y1
    SCALAR_REAL_BUF(202) = Y2
    SCALAR_REAL_BUF(203) = Y3
    SCALAR_REAL_BUF(204) = Y4
    SCALAR_REAL_BUF(205) = FDX1
    SCALAR_REAL_BUF(206) = FDX2
    SCALAR_REAL_BUF(207) = FDX3
    SCALAR_REAL_BUF(208) = FDY1
    SCALAR_REAL_BUF(209) = FDY2
    SCALAR_REAL_BUF(210) = FDY3
    SCALAR_REAL_BUF(211) = FDX1OA
    SCALAR_REAL_BUF(212) = FDX2OA
    SCALAR_REAL_BUF(213) = FDX3OA
    SCALAR_REAL_BUF(214) = FDY1OA
    SCALAR_REAL_BUF(215) = FDY2OA
    SCALAR_REAL_BUF(216) = FDY3OA
    SCALAR_REAL_BUF(217) = AREAIE
    SCALAR_REAL_BUF(218) = DDX1
    SCALAR_REAL_BUF(219) = DDX2
    SCALAR_REAL_BUF(220) = DDX3
    SCALAR_REAL_BUF(221) = DDY1
    SCALAR_REAL_BUF(222) = DDY2
    SCALAR_REAL_BUF(223) = DDY3
    SCALAR_REAL_BUF(224) = DXX11
    SCALAR_REAL_BUF(225) = DXX12
    SCALAR_REAL_BUF(226) = DXX13
    SCALAR_REAL_BUF(227) = DXX21
    SCALAR_REAL_BUF(228) = DXX22
    SCALAR_REAL_BUF(229) = DXX23
    SCALAR_REAL_BUF(230) = DXX31
    SCALAR_REAL_BUF(231) = DXX32
    SCALAR_REAL_BUF(232) = DXX33
    SCALAR_REAL_BUF(233) = DYY11
    SCALAR_REAL_BUF(234) = DYY12
    SCALAR_REAL_BUF(235) = DYY13
    SCALAR_REAL_BUF(236) = DYY21
    SCALAR_REAL_BUF(237) = DYY22
    SCALAR_REAL_BUF(238) = DYY23
    SCALAR_REAL_BUF(239) = DYY31
    SCALAR_REAL_BUF(240) = DYY32
    SCALAR_REAL_BUF(241) = DYY33
    SCALAR_REAL_BUF(242) = DXY11
    SCALAR_REAL_BUF(243) = DXY12
    SCALAR_REAL_BUF(244) = DXY13
    SCALAR_REAL_BUF(245) = DXY21
    SCALAR_REAL_BUF(246) = DXY22
    SCALAR_REAL_BUF(247) = DXY23
    SCALAR_REAL_BUF(248) = DXY31
    SCALAR_REAL_BUF(249) = DXY32
    SCALAR_REAL_BUF(250) = DXY33
    SCALAR_REAL_BUF(251) = XL0
    SCALAR_REAL_BUF(252) = XL1
    SCALAR_REAL_BUF(253) = XL2
    SCALAR_REAL_BUF(254) = YL0
    SCALAR_REAL_BUF(255) = YL1
    SCALAR_REAL_BUF(256) = YL2
    SCALAR_REAL_BUF(257) = SLAM0
    SCALAR_REAL_BUF(258) = SFEA0
    SCALAR_REAL_BUF(259) = WREFTIM
    SCALAR_REAL_BUF(260) = WTIMED
    SCALAR_REAL_BUF(261) = WTIME2
    SCALAR_REAL_BUF(262) = WTIME1
    SCALAR_REAL_BUF(263) = WTIMINC
    SCALAR_REAL_BUF(264) = QTIME1
    SCALAR_REAL_BUF(265) = QTIME2
    SCALAR_REAL_BUF(266) = FTIMINC
    SCALAR_REAL_BUF(267) = ETIMINC
    SCALAR_REAL_BUF(268) = RSTIME1
    SCALAR_REAL_BUF(269) = RSTIME2
    SCALAR_REAL_BUF(270) = RSTIMINC
    SCALAR_REAL_BUF(271) = DELX
    SCALAR_REAL_BUF(272) = DELY
    SCALAR_REAL_BUF(273) = DIST
    SCALAR_REAL_BUF(274) = DELDIST
    SCALAR_REAL_BUF(275) = DELETA
    SCALAR_REAL_BUF(276) = ADVECX
    SCALAR_REAL_BUF(277) = ADVECY
    SCALAR_REAL_BUF(278) = AGIRD
    SCALAR_REAL_BUF(279) = AH
    SCALAR_REAL_BUF(280) = AO12
    SCALAR_REAL_BUF(281) = AO6
    SCALAR_REAL_BUF(282) = ARG
    SCALAR_REAL_BUF(283) = ARG1
    SCALAR_REAL_BUF(284) = ARG2
    SCALAR_REAL_BUF(285) = ARGJ
    SCALAR_REAL_BUF(286) = ARGJ1
    SCALAR_REAL_BUF(287) = ARGJ2
    SCALAR_REAL_BUF(288) = ARGSALT
    SCALAR_REAL_BUF(289) = ARGT
    SCALAR_REAL_BUF(290) = ARGTP
    SCALAR_REAL_BUF(291) = AUV21
    SCALAR_REAL_BUF(292) = AUV22
    SCALAR_REAL_BUF(293) = BARAVGWT
    SCALAR_REAL_BUF(294) = BEDSTR
    SCALAR_REAL_BUF(295) = BNDLEN2O3NC
    SCALAR_REAL_BUF(296) = BSXN1
    SCALAR_REAL_BUF(297) = BSXN2
    SCALAR_REAL_BUF(298) = BSXN3
    SCALAR_REAL_BUF(299) = BSXPP3
    SCALAR_REAL_BUF(300) = BSYN1
    SCALAR_REAL_BUF(301) = BSYN2
    SCALAR_REAL_BUF(302) = BSYN3
    SCALAR_REAL_BUF(303) = BSYPP3
    SCALAR_REAL_BUF(304) = C1
    SCALAR_REAL_BUF(305) = C2
    SCALAR_REAL_BUF(306) = C3
    SCALAR_REAL_BUF(307) = CBEDSTRD
    SCALAR_REAL_BUF(308) = CBEDSTRE
    SCALAR_REAL_BUF(309) = CCRITD
    SCALAR_REAL_BUF(310) = CCSFEA
    SCALAR_REAL_BUF(311) = CELERITY
    SCALAR_REAL_BUF(312) = CH1N1
    SCALAR_REAL_BUF(313) = CH1N2
    SCALAR_REAL_BUF(314) = CH1N3
    SCALAR_REAL_BUF(315) = CHSUM
    SCALAR_REAL_BUF(316) = COND
    SCALAR_REAL_BUF(317) = CONVCR
    SCALAR_REAL_BUF(318) = CORIFPP
    SCALAR_REAL_BUF(319) = DDU
    SCALAR_REAL_BUF(320) = DHDX
    SCALAR_REAL_BUF(321) = DHDY
    SCALAR_REAL_BUF(322) = DISPERX
    SCALAR_REAL_BUF(323) = DISPERY
    SCALAR_REAL_BUF(324) = DT2
    SCALAR_REAL_BUF(325) = DTO2
    SCALAR_REAL_BUF(326) = DTOHPP
    SCALAR_REAL_BUF(327) = DUU1N1
    SCALAR_REAL_BUF(328) = DUU1N2
    SCALAR_REAL_BUF(329) = DUU1N3
    SCALAR_REAL_BUF(330) = DUV1N1
    SCALAR_REAL_BUF(331) = DUV1N2
    SCALAR_REAL_BUF(332) = DUV1N3
    SCALAR_REAL_BUF(333) = DVV1N1
    SCALAR_REAL_BUF(334) = DVV1N2
    SCALAR_REAL_BUF(335) = DVV1N3
    SCALAR_REAL_BUF(336) = DXXYY11
    SCALAR_REAL_BUF(337) = DXXYY12
    SCALAR_REAL_BUF(338) = DXXYY13
    SCALAR_REAL_BUF(339) = DXXYY21
    SCALAR_REAL_BUF(340) = DXXYY22
    SCALAR_REAL_BUF(341) = DXXYY23
    SCALAR_REAL_BUF(342) = DXXYY31
    SCALAR_REAL_BUF(343) = DXXYY32
    SCALAR_REAL_BUF(344) = DXXYY33
    SCALAR_REAL_BUF(345) = DXYH11
    SCALAR_REAL_BUF(346) = DXYH12
    SCALAR_REAL_BUF(347) = DXYH13
    SCALAR_REAL_BUF(348) = DXYH21
    SCALAR_REAL_BUF(349) = DXYH22
    SCALAR_REAL_BUF(350) = DXYH23
    SCALAR_REAL_BUF(351) = DXYH31
    SCALAR_REAL_BUF(352) = DXYH32
    SCALAR_REAL_BUF(353) = DXYH33
    SCALAR_REAL_BUF(354) = E0N1
    SCALAR_REAL_BUF(355) = E0N2
    SCALAR_REAL_BUF(356) = E0N3
    SCALAR_REAL_BUF(357) = E1N1
    SCALAR_REAL_BUF(358) = E1N1SQ
    SCALAR_REAL_BUF(359) = E1N2
    SCALAR_REAL_BUF(360) = E1N2SQ
    SCALAR_REAL_BUF(361) = E1N3
    SCALAR_REAL_BUF(362) = E1N3SQ
    SCALAR_REAL_BUF(363) = ECONST
    SCALAR_REAL_BUF(364) = EE1
    SCALAR_REAL_BUF(365) = EE2
    SCALAR_REAL_BUF(366) = EE3
    SCALAR_REAL_BUF(367) = ELMAX
    SCALAR_REAL_BUF(368) = EP
    SCALAR_REAL_BUF(369) = ESN1
    SCALAR_REAL_BUF(370) = ESN2
    SCALAR_REAL_BUF(371) = BE1
    SCALAR_REAL_BUF(372) = BE2
    SCALAR_REAL_BUF(373) = BE3
    SCALAR_REAL_BUF(374) = ESN3
    SCALAR_REAL_BUF(375) = ETIME1
    SCALAR_REAL_BUF(376) = ETIME2
    SCALAR_REAL_BUF(377) = ETRATIO
    SCALAR_REAL_BUF(378) = EVC1
    SCALAR_REAL_BUF(379) = EVC2
    SCALAR_REAL_BUF(380) = EVC3
    SCALAR_REAL_BUF(381) = EVCEA
    SCALAR_REAL_BUF(382) = EVMPPODT
    SCALAR_REAL_BUF(383) = EVMPPDT
    SCALAR_REAL_BUF(384) = FDDD
    SCALAR_REAL_BUF(385) = FDDDODT
    SCALAR_REAL_BUF(386) = FDDOD
    SCALAR_REAL_BUF(387) = FDDODODT
    SCALAR_REAL_BUF(388) = FIIN
    SCALAR_REAL_BUF(389) = G
    SCALAR_REAL_BUF(390) = GA00
    SCALAR_REAL_BUF(391) = GB00A00
    SCALAR_REAL_BUF(392) = GC00
    SCALAR_REAL_BUF(393) = GDTO2
    SCALAR_REAL_BUF(394) = GFAO2
    SCALAR_REAL_BUF(395) = GHPP
    SCALAR_REAL_BUF(396) = GO3
    SCALAR_REAL_BUF(397) = HABSMIN
    SCALAR_REAL_BUF(398) = HEA
    SCALAR_REAL_BUF(399) = HH1
    SCALAR_REAL_BUF(400) = HH1N1
    SCALAR_REAL_BUF(401) = HH1N2
    SCALAR_REAL_BUF(402) = HH1N3
    SCALAR_REAL_BUF(403) = HH2
    SCALAR_REAL_BUF(404) = HH2N1
    SCALAR_REAL_BUF(405) = HH2N2
    SCALAR_REAL_BUF(406) = HH2N3
    SCALAR_REAL_BUF(407) = HHU1N1
    SCALAR_REAL_BUF(408) = HHU1N2
    SCALAR_REAL_BUF(409) = HHU1N3
    SCALAR_REAL_BUF(410) = HHV1N1
    SCALAR_REAL_BUF(411) = HHV1N2
    SCALAR_REAL_BUF(412) = HHV1N3
    SCALAR_REAL_BUF(413) = HPP
    SCALAR_REAL_BUF(414) = HSD
    SCALAR_REAL_BUF(415) = HSE
    SCALAR_REAL_BUF(416) = HTOT
    SCALAR_REAL_BUF(417) = G2ROOT
    SCALAR_REAL_BUF(418) = P11
    SCALAR_REAL_BUF(419) = P22
    SCALAR_REAL_BUF(420) = P33
    SCALAR_REAL_BUF(421) = PR1N1
    SCALAR_REAL_BUF(422) = PR1N2
    SCALAR_REAL_BUF(423) = PR1N3
    SCALAR_REAL_BUF(424) = QFORCEI
    SCALAR_REAL_BUF(425) = QFORCEJ
    SCALAR_REAL_BUF(426) = QTEMA1
    SCALAR_REAL_BUF(427) = QTEMA2
    SCALAR_REAL_BUF(428) = QTEMA3
    SCALAR_REAL_BUF(429) = QTEMB1
    SCALAR_REAL_BUF(430) = QTEMB2
    SCALAR_REAL_BUF(431) = QTEMB3
    SCALAR_REAL_BUF(432) = QTRATIO
    SCALAR_REAL_BUF(433) = QUNORM
    SCALAR_REAL_BUF(434) = QVNORM
    SCALAR_REAL_BUF(435) = RAMP1
    SCALAR_REAL_BUF(436) = RAMP2
    SCALAR_REAL_BUF(437) = RBARWL
    SCALAR_REAL_BUF(438) = RBARWL1
    SCALAR_REAL_BUF(439) = RBARWL1F
    SCALAR_REAL_BUF(440) = RBARWL2
    SCALAR_REAL_BUF(441) = RBARWL2F
    SCALAR_REAL_BUF(442) = RFF
    SCALAR_REAL_BUF(443) = RFF1
    SCALAR_REAL_BUF(444) = RFF2
    SCALAR_REAL_BUF(445) = RHO0
    SCALAR_REAL_BUF(446) = RSTRATIO
    SCALAR_REAL_BUF(447) = RSX
    SCALAR_REAL_BUF(448) = RSY
    SCALAR_REAL_BUF(449) = S2SFEA
    SCALAR_REAL_BUF(450) = SADVDTO3
    SCALAR_REAL_BUF(451) = SALTMUL
    SCALAR_REAL_BUF(452) = SFACPP
    SCALAR_REAL_BUF(453) = SS1N1
    SCALAR_REAL_BUF(454) = SS1N2
    SCALAR_REAL_BUF(455) = SS1N3
    SCALAR_REAL_BUF(456) = T0N1
    SCALAR_REAL_BUF(457) = T0N2
    SCALAR_REAL_BUF(458) = T0N3
    SCALAR_REAL_BUF(459) = T0XN1
    SCALAR_REAL_BUF(460) = T0XN2
    SCALAR_REAL_BUF(461) = T0XN3
    SCALAR_REAL_BUF(462) = T0XPP3
    SCALAR_REAL_BUF(463) = T0YN1
    SCALAR_REAL_BUF(464) = T0YN2
    SCALAR_REAL_BUF(465) = T0YN3
    SCALAR_REAL_BUF(466) = T0YPP3
    SCALAR_REAL_BUF(467) = TADVODT
    SCALAR_REAL_BUF(468) = TAU0AVG
    SCALAR_REAL_BUF(469) = THENALLDSSSTUP
    SCALAR_REAL_BUF(470) = TIMEIT
    SCALAR_REAL_BUF(471) = TIPN1
    SCALAR_REAL_BUF(472) = TOUTFC
    SCALAR_REAL_BUF(473) = TIPN2
    SCALAR_REAL_BUF(474) = TIPN3
    SCALAR_REAL_BUF(475) = TKWET
    SCALAR_REAL_BUF(476) = TOUTFGC
    SCALAR_REAL_BUF(477) = TOUTFGE
    SCALAR_REAL_BUF(478) = TOUTFGV
    SCALAR_REAL_BUF(479) = TOUTFGW
    SCALAR_REAL_BUF(480) = TOUTFM
    SCALAR_REAL_BUF(481) = TOUTSGC
    SCALAR_REAL_BUF(482) = TOUTSGE
    SCALAR_REAL_BUF(483) = TOUTSGV
    SCALAR_REAL_BUF(484) = TOUTSGW
    SCALAR_REAL_BUF(485) = TOUTSM
    SCALAR_REAL_BUF(486) = TPMUL
    SCALAR_REAL_BUF(487) = TT0L
    SCALAR_REAL_BUF(488) = TT0R
    SCALAR_REAL_BUF(489) = U11
    SCALAR_REAL_BUF(490) = U1N1
    SCALAR_REAL_BUF(491) = U1N2
    SCALAR_REAL_BUF(492) = U1N3
    SCALAR_REAL_BUF(493) = U22
    SCALAR_REAL_BUF(494) = U33
    SCALAR_REAL_BUF(495) = UEA
    SCALAR_REAL_BUF(496) = UHPP
    SCALAR_REAL_BUF(497) = UHPP3
    SCALAR_REAL_BUF(498) = UN1
    SCALAR_REAL_BUF(499) = UPEA
    SCALAR_REAL_BUF(500) = UPP
    SCALAR_REAL_BUF(501) = UPPDT
    SCALAR_REAL_BUF(502) = UPPDTDDX1
    SCALAR_REAL_BUF(503) = UPPDTDDX2
    SCALAR_REAL_BUF(504) = UPPDTDDX3
    SCALAR_REAL_BUF(505) = UV1
    SCALAR_REAL_BUF(506) = V11
    SCALAR_REAL_BUF(507) = V1N1
    SCALAR_REAL_BUF(508) = UVW1
    SCALAR_REAL_BUF(509) = UVW2
    SCALAR_REAL_BUF(510) = COSA1
    SCALAR_REAL_BUF(511) = SINA1
    SCALAR_REAL_BUF(512) = WD1
    SCALAR_REAL_BUF(513) = WD
    SCALAR_REAL_BUF(514) = WDXX
    SCALAR_REAL_BUF(515) = WDYY
    SCALAR_REAL_BUF(516) = WDXY
    SCALAR_REAL_BUF(517) = CHYBR
    SCALAR_REAL_BUF(518) = VCOEFXX
    SCALAR_REAL_BUF(519) = VCOEFYY
    SCALAR_REAL_BUF(520) = VCOEFXY
    SCALAR_REAL_BUF(521) = VCOEFYX
    SCALAR_REAL_BUF(522) = VCOEF2
    SCALAR_REAL_BUF(523) = V1N2
    SCALAR_REAL_BUF(524) = V1N3
    SCALAR_REAL_BUF(525) = V22
    SCALAR_REAL_BUF(526) = V33
    SCALAR_REAL_BUF(527) = VCOEF3N1
    SCALAR_REAL_BUF(528) = VCOEF3N2
    SCALAR_REAL_BUF(529) = VCOEF3N3
    SCALAR_REAL_BUF(530) = VCOEF3X
    SCALAR_REAL_BUF(531) = VCOEF3Y
    SCALAR_REAL_BUF(532) = VEA
    SCALAR_REAL_BUF(533) = VEL
    SCALAR_REAL_BUF(534) = VELABS
    SCALAR_REAL_BUF(535) = VELMAX
    SCALAR_REAL_BUF(536) = VELNORM
    SCALAR_REAL_BUF(537) = VELTAN
    SCALAR_REAL_BUF(538) = VHPP
    SCALAR_REAL_BUF(539) = VHPP3
    SCALAR_REAL_BUF(540) = VPEA
    SCALAR_REAL_BUF(541) = VPP
    SCALAR_REAL_BUF(542) = VPPDT
    SCALAR_REAL_BUF(543) = VPPDTDDY1
    SCALAR_REAL_BUF(544) = VPPDTDDY2
    SCALAR_REAL_BUF(545) = VPPDTDDY3
    SCALAR_REAL_BUF(546) = WDRAGCO
    SCALAR_REAL_BUF(547) = WINDMAG
    SCALAR_REAL_BUF(548) = WINDX
    SCALAR_REAL_BUF(549) = WINDY
    SCALAR_REAL_BUF(550) = WS
    SCALAR_REAL_BUF(551) = WSMOD
    SCALAR_REAL_BUF(552) = WSX
    SCALAR_REAL_BUF(553) = WSXN1
    SCALAR_REAL_BUF(554) = WSXN2
    SCALAR_REAL_BUF(555) = WSXN3
    SCALAR_REAL_BUF(556) = WSY
    SCALAR_REAL_BUF(557) = WSYN1
    SCALAR_REAL_BUF(558) = WSYN2
    SCALAR_REAL_BUF(559) = WSYN3
    SCALAR_REAL_BUF(560) = WTRATIO
    SCALAR_REAL_BUF(561) = A00
    SCALAR_REAL_BUF(562) = B00
    SCALAR_REAL_BUF(563) = C00
    SCALAR_REAL_BUF(564) = ANGINN
    SCALAR_REAL_BUF(565) = CORI
    SCALAR_REAL_BUF(566) = COSTHETA
    SCALAR_REAL_BUF(567) = COSTHETA1
    SCALAR_REAL_BUF(568) = COSTSET
    SCALAR_REAL_BUF(569) = CROSS
    SCALAR_REAL_BUF(570) = CROSS1
    SCALAR_REAL_BUF(571) = DAY
    SCALAR_REAL_BUF(572) = DOTVEC
    SCALAR_REAL_BUF(573) = DRAMP
    SCALAR_REAL_BUF(574) = DUM1
    SCALAR_REAL_BUF(575) = DUM2
    SCALAR_REAL_BUF(576) = EVMSUM
    SCALAR_REAL_BUF(577) = H0
    SCALAR_REAL_BUF(578) = H0L
    SCALAR_REAL_BUF(579) = H0H
    SCALAR_REAL_BUF(580) = RNDAY
    SCALAR_REAL_BUF(581) = THETA
    SCALAR_REAL_BUF(582) = THETA1
    SCALAR_REAL_BUF(583) = TOUTSC
    SCALAR_REAL_BUF(584) = RAMP
    SCALAR_REAL_BUF(585) = RHOWAT0
    SCALAR_REAL_BUF(586) = TOUTSE
    SCALAR_REAL_BUF(587) = TOUTFE
    SCALAR_REAL_BUF(588) = TOUTSV
    SCALAR_REAL_BUF(589) = TOUTFV
    SCALAR_REAL_BUF(590) = XL
    SCALAR_REAL_BUF(591) = VECNORM
    SCALAR_REAL_BUF(592) = VL1X
    SCALAR_REAL_BUF(593) = VL1Y
    SCALAR_REAL_BUF(594) = VL2X
    SCALAR_REAL_BUF(595) = VL2Y
    SCALAR_REAL_BUF(596) = WLATMAX
    SCALAR_REAL_BUF(597) = WLONMIN
    SCALAR_REAL_BUF(598) = WLATINC
    SCALAR_REAL_BUF(599) = WLONINC
    SCALAR_REAL_BUF(600) = VELMIN
    SCALAR_REAL_BUF(601) = RNP_GLOBAL
    SCALAR_REAL_BUF(602) = REFSEC
    CALL fstarpu_vector_data_register(handles(516), 0, C_LOC(SCALAR_REAL_BUF(1)), NUM_REAL_SCALARS, C_SIZEOF(SCALAR_REAL_BUF(1)))
    NULLIFY(SCALAR_REAL_BUF)
  END SUBROUTINE DGSWEM_STATE_REGISTER

  SUBROUTINE DGSWEM_STATE_ACTIVATE(buffers)
    TYPE(C_PTR), VALUE, INTENT(IN) :: buffers
    TYPE(C_PTR) :: curr_ptr
    INTEGER, POINTER :: SCALAR_INT_BUF(:)
    REAL(SZ), POINTER :: SCALAR_REAL_BUF(:)
    curr_ptr = fstarpu_vector_get_ptr(buffers, 514)
    CALL c_f_pointer(curr_ptr, SCALAR_INT_BUF, shape=[NUM_INT_SCALARS])
    MNPROC = SCALAR_INT_BUF(515)
    MNE = SCALAR_INT_BUF(516)
    MNP = SCALAR_INT_BUF(517)
    MNEI = SCALAR_INT_BUF(518)
    MNOPE = SCALAR_INT_BUF(519)
    MNETA = SCALAR_INT_BUF(520)
    MNBOU = SCALAR_INT_BUF(521)
    MNVEL = SCALAR_INT_BUF(522)
    MNTIF = SCALAR_INT_BUF(523)
    MNBFR = SCALAR_INT_BUF(524)
    MNFFR = SCALAR_INT_BUF(525)
    MNSTAE = SCALAR_INT_BUF(526)
    MNSTAV = SCALAR_INT_BUF(527)
    MNSTAC = SCALAR_INT_BUF(528)
    MNSTAM = SCALAR_INT_BUF(529)
    MNHARF = SCALAR_INT_BUF(530)
    layers = SCALAR_INT_BUF(531)
    MNNDEL = SCALAR_INT_BUF(532)
    MYPROC = SCALAR_INT_BUF(533)
    LNAME = SCALAR_INT_BUF(534)
    rainfall = SCALAR_INT_BUF(535)
    DGFLAG = SCALAR_INT_BUF(536)
    DGHOT = SCALAR_INT_BUF(537)
    DGHOTSPOOL = SCALAR_INT_BUF(538)
    DOF = SCALAR_INT_BUF(539)
    dofh = SCALAR_INT_BUF(540)
    dofl = SCALAR_INT_BUF(541)
    dofx = SCALAR_INT_BUF(542)
    EL = SCALAR_INT_BUF(543)
    MNES = SCALAR_INT_BUF(544)
    artdif = SCALAR_INT_BUF(545)
    tune_by_hand = SCALAR_INT_BUF(546)
    IRK = SCALAR_INT_BUF(547)
    J1 = SCALAR_INT_BUF(548)
    J2 = SCALAR_INT_BUF(549)
    J3 = SCALAR_INT_BUF(550)
    negp_fixed = SCALAR_INT_BUF(551)
    nagp_fixed = SCALAR_INT_BUF(552)
    NAGP(1) = SCALAR_INT_BUF(553)
    NAGP(2) = SCALAR_INT_BUF(554)
    NAGP(3) = SCALAR_INT_BUF(555)
    NAGP(4) = SCALAR_INT_BUF(556)
    NAGP(5) = SCALAR_INT_BUF(557)
    NAGP(6) = SCALAR_INT_BUF(558)
    NAGP(7) = SCALAR_INT_BUF(559)
    NAGP(8) = SCALAR_INT_BUF(560)
    NCHECK(1) = SCALAR_INT_BUF(561)
    NCHECK(2) = SCALAR_INT_BUF(562)
    NCHECK(3) = SCALAR_INT_BUF(563)
    NCHECK(4) = SCALAR_INT_BUF(564)
    NCHECK(5) = SCALAR_INT_BUF(565)
    NCHECK(6) = SCALAR_INT_BUF(566)
    NCHECK(7) = SCALAR_INT_BUF(567)
    NCHECK(8) = SCALAR_INT_BUF(568)
    NEGP(1) = SCALAR_INT_BUF(569)
    NEGP(2) = SCALAR_INT_BUF(570)
    NEGP(3) = SCALAR_INT_BUF(571)
    NEGP(4) = SCALAR_INT_BUF(572)
    NEGP(5) = SCALAR_INT_BUF(573)
    NEGP(6) = SCALAR_INT_BUF(574)
    NEGP(7) = SCALAR_INT_BUF(575)
    NEGP(8) = SCALAR_INT_BUF(576)
    NEDGES = SCALAR_INT_BUF(577)
    NRK = SCALAR_INT_BUF(578)
    NIEDS = SCALAR_INT_BUF(579)
    NLEDS = SCALAR_INT_BUF(580)
    NEEDS = SCALAR_INT_BUF(581)
    NFEDS = SCALAR_INT_BUF(582)
    NREDS = SCALAR_INT_BUF(583)
    NEBEDS = SCALAR_INT_BUF(584)
    NIBEDS = SCALAR_INT_BUF(585)
    NIBSEG = SCALAR_INT_BUF(586)
    NEBSEG = SCALAR_INT_BUF(587)
    MNED = SCALAR_INT_BUF(588)
    MNLED = SCALAR_INT_BUF(589)
    MNSED = SCALAR_INT_BUF(590)
    MNRAED = SCALAR_INT_BUF(591)
    MNRIED = SCALAR_INT_BUF(592)
    MODAL_IC = SCALAR_INT_BUF(593)
    P_READ = SCALAR_INT_BUF(594)
    P_READ2 = SCALAR_INT_BUF(595)
    SLOPEFLAG = SCALAR_INT_BUF(596)
    test_el = SCALAR_INT_BUF(597)
    FLUXTYPE = SCALAR_INT_BUF(598)
    RK_STAGE = SCALAR_INT_BUF(599)
    RK_ORDER = SCALAR_INT_BUF(600)
    padapt = SCALAR_INT_BUF(601)
    pflag = SCALAR_INT_BUF(602)
    pl = SCALAR_INT_BUF(603)
    ph = SCALAR_INT_BUF(604)
    px = SCALAR_INT_BUF(605)
    lebesgueP = SCALAR_INT_BUF(606)
    gflag = SCALAR_INT_BUF(607)
    pa = SCALAR_INT_BUF(608)
    iwrite = SCALAR_INT_BUF(609)
    lim_count = SCALAR_INT_BUF(610)
    lim_count_roll = SCALAR_INT_BUF(611)
    SEDFLAG = SCALAR_INT_BUF(612)
    MAXEL = SCALAR_INT_BUF(613)
    ELEM_ED = SCALAR_INT_BUF(614)
    NBOR_ED = SCALAR_INT_BUF(615)
    NBOR_EL = SCALAR_INT_BUF(616)
    ITDG = SCALAR_INT_BUF(617)
    ModetoNode = SCALAR_INT_BUF(618)
    tracer_flag = SCALAR_INT_BUF(619)
    chem_flag = SCALAR_INT_BUF(620)
    N1 = SCALAR_INT_BUF(621)
    N2 = SCALAR_INT_BUF(622)
    NO_NBORS = SCALAR_INT_BUF(623)
    NBOR = SCALAR_INT_BUF(624)
    SEDFLAG_W = SCALAR_INT_BUF(625)
    OPEN_INDEX = SCALAR_INT_BUF(626)
    DG_TO_CG = SCALAR_INT_BUF(627)
    NSCREEN_INC = SCALAR_INT_BUF(628)
    ScreenUnit = SCALAR_INT_BUF(629)
    FluxSettlingIT = SCALAR_INT_BUF(630)
    DGSWE = SCALAR_INT_BUF(631)
    EL_IN = SCALAR_INT_BUF(632)
    EL_EX = SCALAR_INT_BUF(633)
    SD_IN = SCALAR_INT_BUF(634)
    SD_EX = SCALAR_INT_BUF(635)
    EDGE(1) = SCALAR_INT_BUF(636)
    EDGE(2) = SCALAR_INT_BUF(637)
    EDGE(3) = SCALAR_INT_BUF(638)
    SIDE(1) = SCALAR_INT_BUF(639)
    SIDE(2) = SCALAR_INT_BUF(640)
    TESTPROBLEM = SCALAR_INT_BUF(641)
    NBPNODES = SCALAR_INT_BUF(642)
    NP = SCALAR_INT_BUF(643)
    NOLICA = SCALAR_INT_BUF(644)
    NOLIFA = SCALAR_INT_BUF(645)
    NSCREEN = SCALAR_INT_BUF(646)
    IHOT = SCALAR_INT_BUF(647)
    ICS = SCALAR_INT_BUF(648)
    FRW = SCALAR_INT_BUF(649)
    NODEDRYMIN = SCALAR_INT_BUF(650)
    NODEWETMIN = SCALAR_INT_BUF(651)
    IBTYPE = SCALAR_INT_BUF(652)
    ICK = SCALAR_INT_BUF(653)
    IDR = SCALAR_INT_BUF(654)
    IM = SCALAR_INT_BUF(655)
    IPRBI = SCALAR_INT_BUF(656)
    JGW = SCALAR_INT_BUF(657)
    JKI = SCALAR_INT_BUF(658)
    JME = SCALAR_INT_BUF(659)
    JNMM = SCALAR_INT_BUF(660)
    KMIN = SCALAR_INT_BUF(661)
    N3 = SCALAR_INT_BUF(662)
    NABOUT = SCALAR_INT_BUF(663)
    NBFR = SCALAR_INT_BUF(664)
    NBOU = SCALAR_INT_BUF(665)
    NBVI = SCALAR_INT_BUF(666)
    NBVJ = SCALAR_INT_BUF(667)
    NCOR = SCALAR_INT_BUF(668)
    NE = SCALAR_INT_BUF(669)
    NE2 = SCALAR_INT_BUF(670)
    NP2 = SCALAR_INT_BUF(671)
    NEIMIN = SCALAR_INT_BUF(672)
    NEIMAX = SCALAR_INT_BUF(673)
    NETA = SCALAR_INT_BUF(674)
    NFFR = SCALAR_INT_BUF(675)
    NFLUXB = SCALAR_INT_BUF(676)
    NFLUXF = SCALAR_INT_BUF(677)
    NFLUXIB = SCALAR_INT_BUF(678)
    NFLUXRBC = SCALAR_INT_BUF(679)
    NFLUXIBP = SCALAR_INT_BUF(680)
    NPIPE = SCALAR_INT_BUF(681)
    NFOVER = SCALAR_INT_BUF(682)
    NHG = SCALAR_INT_BUF(683)
    NHY = SCALAR_INT_BUF(684)
    NOLICAT = SCALAR_INT_BUF(685)
    NOPE = SCALAR_INT_BUF(686)
    NOUTC = SCALAR_INT_BUF(687)
    NOUTE = SCALAR_INT_BUF(688)
    NSPOOLE = SCALAR_INT_BUF(689)
    NOUTV = SCALAR_INT_BUF(690)
    NSPOOLV = SCALAR_INT_BUF(691)
    NPRBI = SCALAR_INT_BUF(692)
    NRAMP = SCALAR_INT_BUF(693)
    NRS = SCALAR_INT_BUF(694)
    NSTAE = SCALAR_INT_BUF(695)
    NSTARTDRY = SCALAR_INT_BUF(696)
    NSTAV = SCALAR_INT_BUF(697)
    NT = SCALAR_INT_BUF(698)
    NTCYFE = SCALAR_INT_BUF(699)
    NTCYFV = SCALAR_INT_BUF(700)
    NTCYSE = SCALAR_INT_BUF(701)
    NTCYSV = SCALAR_INT_BUF(702)
    NTIF = SCALAR_INT_BUF(703)
    NTIP = SCALAR_INT_BUF(704)
    NTRSPE = SCALAR_INT_BUF(705)
    NTRSPV = SCALAR_INT_BUF(706)
    NVEL = SCALAR_INT_BUF(707)
    NVELEXT = SCALAR_INT_BUF(708)
    NVELME = SCALAR_INT_BUF(709)
    NWLAT = SCALAR_INT_BUF(710)
    NWLON = SCALAR_INT_BUF(711)
    NWS = SCALAR_INT_BUF(712)
    IBSTART = SCALAR_INT_BUF(713)
    ICHA = SCALAR_INT_BUF(714)
    ICSTP = SCALAR_INT_BUF(715)
    IDSETFLG = SCALAR_INT_BUF(716)
    IE = SCALAR_INT_BUF(717)
    IER = SCALAR_INT_BUF(718)
    IESTP = SCALAR_INT_BUF(719)
    IFNLCAT = SCALAR_INT_BUF(720)
    IFNLCT = SCALAR_INT_BUF(721)
    IFNLFA = SCALAR_INT_BUF(722)
    IFWIND = SCALAR_INT_BUF(723)
    IGCP = SCALAR_INT_BUF(724)
    IGEP = SCALAR_INT_BUF(725)
    IGPP = SCALAR_INT_BUF(726)
    IGVP = SCALAR_INT_BUF(727)
    IGWP = SCALAR_INT_BUF(728)
    IHABEG = SCALAR_INT_BUF(729)
    IGRadS = SCALAR_INT_BUF(730)
    IHOTSTP = SCALAR_INT_BUF(731)
    IHSFIL = SCALAR_INT_BUF(732)
    IJ = SCALAR_INT_BUF(733)
    ILUMP = SCALAR_INT_BUF(734)
    IMHS = SCALAR_INT_BUF(735)
    IPSTP = SCALAR_INT_BUF(736)
    IREFYR = SCALAR_INT_BUF(737)
    IREFMO = SCALAR_INT_BUF(738)
    IREFDAY = SCALAR_INT_BUF(739)
    IREFHR = SCALAR_INT_BUF(740)
    IREFMIN = SCALAR_INT_BUF(741)
    ISLDIA = SCALAR_INT_BUF(742)
    ITIME_A = SCALAR_INT_BUF(743)
    ITEMPSTP = SCALAR_INT_BUF(744)
    ITEST = SCALAR_INT_BUF(745)
    ITHS = SCALAR_INT_BUF(746)
    ITITER = SCALAR_INT_BUF(747)
    ITMAX = SCALAR_INT_BUF(748)
    IVSTP = SCALAR_INT_BUF(749)
    IWSTP = SCALAR_INT_BUF(750)
    IWTIME = SCALAR_INT_BUF(751)
    IWTIMEP = SCALAR_INT_BUF(752)
    IWYR = SCALAR_INT_BUF(753)
    J12 = SCALAR_INT_BUF(754)
    J13 = SCALAR_INT_BUF(755)
    J21 = SCALAR_INT_BUF(756)
    J23 = SCALAR_INT_BUF(757)
    J31 = SCALAR_INT_BUF(758)
    J32 = SCALAR_INT_BUF(759)
    JN = SCALAR_INT_BUF(760)
    KEMAX = SCALAR_INT_BUF(761)
    KVMAX = SCALAR_INT_BUF(762)
    LRC = SCALAR_INT_BUF(763)
    LUMPT = SCALAR_INT_BUF(764)
    MMAX = SCALAR_INT_BUF(765)
    MBW = SCALAR_INT_BUF(766)
    MDF = SCALAR_INT_BUF(767)
    MMIN = SCALAR_INT_BUF(768)
    NA = SCALAR_INT_BUF(769)
    NBDI = SCALAR_INT_BUF(770)
    NBDJ = SCALAR_INT_BUF(771)
    NBNCTOT = SCALAR_INT_BUF(772)
    NBW = SCALAR_INT_BUF(773)
    NC1 = SCALAR_INT_BUF(774)
    NC2 = SCALAR_INT_BUF(775)
    NC3 = SCALAR_INT_BUF(776)
    NCBND = SCALAR_INT_BUF(777)
    NCELE = SCALAR_INT_BUF(778)
    NCI = SCALAR_INT_BUF(779)
    NCJ = SCALAR_INT_BUF(780)
    NCTOT = SCALAR_INT_BUF(781)
    NCYC = SCALAR_INT_BUF(782)
    NDRY = SCALAR_INT_BUF(783)
    NDSETSC = SCALAR_INT_BUF(784)
    NDSETSE = SCALAR_INT_BUF(785)
    NDSETSV = SCALAR_INT_BUF(786)
    NDSETSW = SCALAR_INT_BUF(787)
    NHSINC = SCALAR_INT_BUF(788)
    NHSTAR = SCALAR_INT_BUF(789)
    NM1 = SCALAR_INT_BUF(790)
    NM123 = SCALAR_INT_BUF(791)
    NM2 = SCALAR_INT_BUF(792)
    NM3 = SCALAR_INT_BUF(793)
    NMI1 = SCALAR_INT_BUF(794)
    NMI2 = SCALAR_INT_BUF(795)
    NMI3 = SCALAR_INT_BUF(796)
    NMJ1 = SCALAR_INT_BUF(797)
    NMJ2 = SCALAR_INT_BUF(798)
    NMJ3 = SCALAR_INT_BUF(799)
    NNBB = SCALAR_INT_BUF(800)
    NNBB1 = SCALAR_INT_BUF(801)
    NNBB2 = SCALAR_INT_BUF(802)
    NOUTGC = SCALAR_INT_BUF(803)
    NOUTGE = SCALAR_INT_BUF(804)
    NOUTGV = SCALAR_INT_BUF(805)
    NOUTGW = SCALAR_INT_BUF(806)
    NOUTM = SCALAR_INT_BUF(807)
    NSCOUC = SCALAR_INT_BUF(808)
    NSCOUE = SCALAR_INT_BUF(809)
    NSCOUGC = SCALAR_INT_BUF(810)
    NSCOUGE = SCALAR_INT_BUF(811)
    NSCOUGV = SCALAR_INT_BUF(812)
    NSCOUGW = SCALAR_INT_BUF(813)
    NSCOUM = SCALAR_INT_BUF(814)
    NSCOUV = SCALAR_INT_BUF(815)
    NSPOOLC = SCALAR_INT_BUF(816)
    NSPOOLGC = SCALAR_INT_BUF(817)
    NSPOOLGE = SCALAR_INT_BUF(818)
    NSPOOLGV = SCALAR_INT_BUF(819)
    NSPOOLGW = SCALAR_INT_BUF(820)
    NSPOOLM = SCALAR_INT_BUF(821)
    NSTAC = SCALAR_INT_BUF(822)
    NSTAM = SCALAR_INT_BUF(823)
    NTCYFC = SCALAR_INT_BUF(824)
    NTCYFGC = SCALAR_INT_BUF(825)
    NTCYFGE = SCALAR_INT_BUF(826)
    NTCYFGV = SCALAR_INT_BUF(827)
    NTCYFGW = SCALAR_INT_BUF(828)
    NTCYFM = SCALAR_INT_BUF(829)
    NTCYSC = SCALAR_INT_BUF(830)
    NTCYSGC = SCALAR_INT_BUF(831)
    NTCYSGE = SCALAR_INT_BUF(832)
    NTCYSGV = SCALAR_INT_BUF(833)
    NTCYSGW = SCALAR_INT_BUF(834)
    NTCYSM = SCALAR_INT_BUF(835)
    NTRSPC = SCALAR_INT_BUF(836)
    NTRSPM = SCALAR_INT_BUF(837)
    NUMITR = SCALAR_INT_BUF(838)
    NW = SCALAR_INT_BUF(839)
    NWET = SCALAR_INT_BUF(840)
    NWSEGWI = SCALAR_INT_BUF(841)
    NWSGGWI = SCALAR_INT_BUF(842)
    NCCHANGE = SCALAR_INT_BUF(843)
    IRAMPING = SCALAR_INT_BUF(844)
    vertexslope = (SCALAR_INT_BUF(847) == 1)
    curr_ptr = fstarpu_vector_get_ptr(buffers, 515)
    CALL c_f_pointer(curr_ptr, SCALAR_REAL_BUF, shape=[NUM_REAL_SCALARS])
    C13 = SCALAR_REAL_BUF(1)
    C16 = SCALAR_REAL_BUF(2)
    diorism = SCALAR_REAL_BUF(3)
    porosity = SCALAR_REAL_BUF(4)
    SEVDM = SCALAR_REAL_BUF(5)
    DOT = SCALAR_REAL_BUF(6)
    DHB_X = SCALAR_REAL_BUF(7)
    DHB_Y = SCALAR_REAL_BUF(8)
    DPHIDX = SCALAR_REAL_BUF(9)
    DPHIDY = SCALAR_REAL_BUF(10)
    slimit = SCALAR_REAL_BUF(11)
    plimit = SCALAR_REAL_BUF(12)
    pflag2con1 = SCALAR_REAL_BUF(13)
    pflag2con2 = SCALAR_REAL_BUF(14)
    EFA_GP = SCALAR_REAL_BUF(15)
    EMO_GP = SCALAR_REAL_BUF(16)
    slimit1 = SCALAR_REAL_BUF(17)
    slimit2 = SCALAR_REAL_BUF(18)
    slimit3 = SCALAR_REAL_BUF(19)
    EL_ANG = SCALAR_REAL_BUF(20)
    slimit4 = SCALAR_REAL_BUF(21)
    bg_dif = SCALAR_REAL_BUF(22)
    trc_dif = SCALAR_REAL_BUF(23)
    slimit5 = SCALAR_REAL_BUF(24)
    FG_L = SCALAR_REAL_BUF(25)
    l2er_global = SCALAR_REAL_BUF(26)
    temperg = SCALAR_REAL_BUF(27)
    slope_weight = SCALAR_REAL_BUF(28)
    H_TRI = SCALAR_REAL_BUF(29)
    MAG1 = SCALAR_REAL_BUF(30)
    MAG2 = SCALAR_REAL_BUF(31)
    kappa = SCALAR_REAL_BUF(32)
    s0 = SCALAR_REAL_BUF(33)
    uniform_dif = SCALAR_REAL_BUF(34)
    RAMPDG = SCALAR_REAL_BUF(35)
    S1 = SCALAR_REAL_BUF(36)
    S2 = SCALAR_REAL_BUF(37)
    SAV = SCALAR_REAL_BUF(38)
    SOURCE_X = SCALAR_REAL_BUF(39)
    SOURCE_Y = SCALAR_REAL_BUF(40)
    TIMEDG = SCALAR_REAL_BUF(41)
    TIMEH_DG = SCALAR_REAL_BUF(42)
    TK = SCALAR_REAL_BUF(43)
    QNAM_GP = SCALAR_REAL_BUF(44)
    QNPH_GP = SCALAR_REAL_BUF(45)
    SL2_M = SCALAR_REAL_BUF(46)
    SL2_NYU = SCALAR_REAL_BUF(47)
    SL3_MD = SCALAR_REAL_BUF(48)
    EVMAvg = SCALAR_REAL_BUF(49)
    SEVDMAvg = SCALAR_REAL_BUF(50)
    UMAG = SCALAR_REAL_BUF(51)
    WSX_GP = SCALAR_REAL_BUF(52)
    WSY_GP = SCALAR_REAL_BUF(53)
    QMag_IN = SCALAR_REAL_BUF(54)
    QMag_EX = SCALAR_REAL_BUF(55)
    subphi_IN = SCALAR_REAL_BUF(56)
    subphi_EX = SCALAR_REAL_BUF(57)
    iota_EX = SCALAR_REAL_BUF(58)
    iota_IN = SCALAR_REAL_BUF(59)
    iota2_EX = SCALAR_REAL_BUF(60)
    iota2_IN = SCALAR_REAL_BUF(61)
    valx(1) = SCALAR_REAL_BUF(62)
    valx(2) = SCALAR_REAL_BUF(63)
    valx(3) = SCALAR_REAL_BUF(64)
    valx(4) = SCALAR_REAL_BUF(65)
    valy(1) = SCALAR_REAL_BUF(66)
    valy(2) = SCALAR_REAL_BUF(67)
    valy(3) = SCALAR_REAL_BUF(68)
    valy(4) = SCALAR_REAL_BUF(69)
    DRPSI(1) = SCALAR_REAL_BUF(70)
    DRPSI(2) = SCALAR_REAL_BUF(71)
    DRPSI(3) = SCALAR_REAL_BUF(72)
    DSPSI(1) = SCALAR_REAL_BUF(73)
    DSPSI(2) = SCALAR_REAL_BUF(74)
    DSPSI(3) = SCALAR_REAL_BUF(75)
    VEC1(1) = SCALAR_REAL_BUF(76)
    VEC1(2) = SCALAR_REAL_BUF(77)
    VEC2(1) = SCALAR_REAL_BUF(78)
    VEC2(2) = SCALAR_REAL_BUF(79)
    rhoAir = SCALAR_REAL_BUF(80)
    windReduction = SCALAR_REAL_BUF(81)
    one2ten = SCALAR_REAL_BUF(82)
    ten2one = SCALAR_REAL_BUF(83)
    WaveWindMultiplier = SCALAR_REAL_BUF(84)
    DEPAVG = SCALAR_REAL_BUF(85)
    DEPMAX = SCALAR_REAL_BUF(86)
    DEPMIN = SCALAR_REAL_BUF(87)
    AREA_SUM = SCALAR_REAL_BUF(88)
    CEN_SUM = SCALAR_REAL_BUF(89)
    NLEQ = SCALAR_REAL_BUF(90)
    LEQ = SCALAR_REAL_BUF(91)
    NLEQG = SCALAR_REAL_BUF(92)
    reaction_rate = SCALAR_REAL_BUF(93)
    FluxSettlingTime = SCALAR_REAL_BUF(94)
    RampExtFlux = SCALAR_REAL_BUF(95)
    DRampExtFlux = SCALAR_REAL_BUF(96)
    RampIntFlux = SCALAR_REAL_BUF(97)
    DRampIntFlux = SCALAR_REAL_BUF(98)
    RampElev = SCALAR_REAL_BUF(99)
    DRampElev = SCALAR_REAL_BUF(100)
    RampTip = SCALAR_REAL_BUF(101)
    DRampTip = SCALAR_REAL_BUF(102)
    RampMete = SCALAR_REAL_BUF(103)
    DRampMete = SCALAR_REAL_BUF(104)
    RampWRad = SCALAR_REAL_BUF(105)
    DRampWRad = SCALAR_REAL_BUF(106)
    FX_IN = SCALAR_REAL_BUF(107)
    FY_IN = SCALAR_REAL_BUF(108)
    GX_IN = SCALAR_REAL_BUF(109)
    GY_IN = SCALAR_REAL_BUF(110)
    HX_IN = SCALAR_REAL_BUF(111)
    HY_IN = SCALAR_REAL_BUF(112)
    FX_EX = SCALAR_REAL_BUF(113)
    FY_EX = SCALAR_REAL_BUF(114)
    GX_EX = SCALAR_REAL_BUF(115)
    GY_EX = SCALAR_REAL_BUF(116)
    HX_EX = SCALAR_REAL_BUF(117)
    HY_EX = SCALAR_REAL_BUF(118)
    F_AVG = SCALAR_REAL_BUF(119)
    G_AVG = SCALAR_REAL_BUF(120)
    H_AVG = SCALAR_REAL_BUF(121)
    JUMP(1) = SCALAR_REAL_BUF(122)
    JUMP(2) = SCALAR_REAL_BUF(123)
    JUMP(3) = SCALAR_REAL_BUF(124)
    JUMP(4) = SCALAR_REAL_BUF(125)
    HT_IN = SCALAR_REAL_BUF(126)
    HT_EX = SCALAR_REAL_BUF(127)
    UMag_IN = SCALAR_REAL_BUF(128)
    UMag_EX = SCALAR_REAL_BUF(129)
    ZE_ROE = SCALAR_REAL_BUF(130)
    QX_ROE = SCALAR_REAL_BUF(131)
    QY_ROE = SCALAR_REAL_BUF(132)
    bed_ROE = SCALAR_REAL_BUF(133)
    Q_N = SCALAR_REAL_BUF(134)
    Q_T = SCALAR_REAL_BUF(135)
    U_N = SCALAR_REAL_BUF(136)
    U_T = SCALAR_REAL_BUF(137)
    U_IN = SCALAR_REAL_BUF(138)
    U_EX = SCALAR_REAL_BUF(139)
    V_IN = SCALAR_REAL_BUF(140)
    V_EX = SCALAR_REAL_BUF(141)
    ZE_SUM = SCALAR_REAL_BUF(142)
    QX_SUM = SCALAR_REAL_BUF(143)
    QY_SUM = SCALAR_REAL_BUF(144)
    DG_MAX = SCALAR_REAL_BUF(145)
    DG_MIN = SCALAR_REAL_BUF(146)
    U_N_EXT = SCALAR_REAL_BUF(147)
    U_T_EXT = SCALAR_REAL_BUF(148)
    Q_N_INT = SCALAR_REAL_BUF(149)
    Q_T_INT = SCALAR_REAL_BUF(150)
    U_N_INT = SCALAR_REAL_BUF(151)
    U_T_INT = SCALAR_REAL_BUF(152)
    Q_N_EXT = SCALAR_REAL_BUF(153)
    Q_T_EXT = SCALAR_REAL_BUF(154)
    BX_INT = SCALAR_REAL_BUF(155)
    BY_INT = SCALAR_REAL_BUF(156)
    SOURCE_1 = SCALAR_REAL_BUF(157)
    SOURCE_2 = SCALAR_REAL_BUF(158)
    SOURCE_SUM = SCALAR_REAL_BUF(159)
    k_hat = SCALAR_REAL_BUF(160)
    FRIC_AVG = SCALAR_REAL_BUF(161)
    DP_MID = SCALAR_REAL_BUF(162)
    i_hat = SCALAR_REAL_BUF(163)
    j_hat = SCALAR_REAL_BUF(164)
    INFLOW_ZE = SCALAR_REAL_BUF(165)
    INFLOW_QX = SCALAR_REAL_BUF(166)
    INFLOW_QY = SCALAR_REAL_BUF(167)
    H_LEN = SCALAR_REAL_BUF(168)
    INFLOW_LEN = SCALAR_REAL_BUF(169)
    ZE_NORM = SCALAR_REAL_BUF(170)
    QX_NORM = SCALAR_REAL_BUF(171)
    QY_NORM = SCALAR_REAL_BUF(172)
    ZE_DECT = SCALAR_REAL_BUF(173)
    QX_DECT = SCALAR_REAL_BUF(174)
    QY_DECT = SCALAR_REAL_BUF(175)
    POAN = SCALAR_REAL_BUF(176)
    Fr = SCALAR_REAL_BUF(177)
    FRICBP = SCALAR_REAL_BUF(178)
    STATIM = SCALAR_REAL_BUF(179)
    REFTIM = SCALAR_REAL_BUF(180)
    TIME_A = SCALAR_REAL_BUF(181)
    DT = SCALAR_REAL_BUF(182)
    DTDP = SCALAR_REAL_BUF(183)
    TIMEH = SCALAR_REAL_BUF(184)
    vdtdp = SCALAR_REAL_BUF(185)
    cfl_max = SCALAR_REAL_BUF(186)
    AVGXY = SCALAR_REAL_BUF(187)
    DIF1R = SCALAR_REAL_BUF(188)
    DIF2R = SCALAR_REAL_BUF(189)
    DIF3R = SCALAR_REAL_BUF(190)
    AEMIN = SCALAR_REAL_BUF(191)
    AE = SCALAR_REAL_BUF(192)
    AA = SCALAR_REAL_BUF(193)
    A1 = SCALAR_REAL_BUF(194)
    A2 = SCALAR_REAL_BUF(195)
    A3 = SCALAR_REAL_BUF(196)
    X1 = SCALAR_REAL_BUF(197)
    X2 = SCALAR_REAL_BUF(198)
    X3 = SCALAR_REAL_BUF(199)
    X4 = SCALAR_REAL_BUF(200)
    Y1 = SCALAR_REAL_BUF(201)
    Y2 = SCALAR_REAL_BUF(202)
    Y3 = SCALAR_REAL_BUF(203)
    Y4 = SCALAR_REAL_BUF(204)
    FDX1 = SCALAR_REAL_BUF(205)
    FDX2 = SCALAR_REAL_BUF(206)
    FDX3 = SCALAR_REAL_BUF(207)
    FDY1 = SCALAR_REAL_BUF(208)
    FDY2 = SCALAR_REAL_BUF(209)
    FDY3 = SCALAR_REAL_BUF(210)
    FDX1OA = SCALAR_REAL_BUF(211)
    FDX2OA = SCALAR_REAL_BUF(212)
    FDX3OA = SCALAR_REAL_BUF(213)
    FDY1OA = SCALAR_REAL_BUF(214)
    FDY2OA = SCALAR_REAL_BUF(215)
    FDY3OA = SCALAR_REAL_BUF(216)
    AREAIE = SCALAR_REAL_BUF(217)
    DDX1 = SCALAR_REAL_BUF(218)
    DDX2 = SCALAR_REAL_BUF(219)
    DDX3 = SCALAR_REAL_BUF(220)
    DDY1 = SCALAR_REAL_BUF(221)
    DDY2 = SCALAR_REAL_BUF(222)
    DDY3 = SCALAR_REAL_BUF(223)
    DXX11 = SCALAR_REAL_BUF(224)
    DXX12 = SCALAR_REAL_BUF(225)
    DXX13 = SCALAR_REAL_BUF(226)
    DXX21 = SCALAR_REAL_BUF(227)
    DXX22 = SCALAR_REAL_BUF(228)
    DXX23 = SCALAR_REAL_BUF(229)
    DXX31 = SCALAR_REAL_BUF(230)
    DXX32 = SCALAR_REAL_BUF(231)
    DXX33 = SCALAR_REAL_BUF(232)
    DYY11 = SCALAR_REAL_BUF(233)
    DYY12 = SCALAR_REAL_BUF(234)
    DYY13 = SCALAR_REAL_BUF(235)
    DYY21 = SCALAR_REAL_BUF(236)
    DYY22 = SCALAR_REAL_BUF(237)
    DYY23 = SCALAR_REAL_BUF(238)
    DYY31 = SCALAR_REAL_BUF(239)
    DYY32 = SCALAR_REAL_BUF(240)
    DYY33 = SCALAR_REAL_BUF(241)
    DXY11 = SCALAR_REAL_BUF(242)
    DXY12 = SCALAR_REAL_BUF(243)
    DXY13 = SCALAR_REAL_BUF(244)
    DXY21 = SCALAR_REAL_BUF(245)
    DXY22 = SCALAR_REAL_BUF(246)
    DXY23 = SCALAR_REAL_BUF(247)
    DXY31 = SCALAR_REAL_BUF(248)
    DXY32 = SCALAR_REAL_BUF(249)
    DXY33 = SCALAR_REAL_BUF(250)
    XL0 = SCALAR_REAL_BUF(251)
    XL1 = SCALAR_REAL_BUF(252)
    XL2 = SCALAR_REAL_BUF(253)
    YL0 = SCALAR_REAL_BUF(254)
    YL1 = SCALAR_REAL_BUF(255)
    YL2 = SCALAR_REAL_BUF(256)
    SLAM0 = SCALAR_REAL_BUF(257)
    SFEA0 = SCALAR_REAL_BUF(258)
    WREFTIM = SCALAR_REAL_BUF(259)
    WTIMED = SCALAR_REAL_BUF(260)
    WTIME2 = SCALAR_REAL_BUF(261)
    WTIME1 = SCALAR_REAL_BUF(262)
    WTIMINC = SCALAR_REAL_BUF(263)
    QTIME1 = SCALAR_REAL_BUF(264)
    QTIME2 = SCALAR_REAL_BUF(265)
    FTIMINC = SCALAR_REAL_BUF(266)
    ETIMINC = SCALAR_REAL_BUF(267)
    RSTIME1 = SCALAR_REAL_BUF(268)
    RSTIME2 = SCALAR_REAL_BUF(269)
    RSTIMINC = SCALAR_REAL_BUF(270)
    DELX = SCALAR_REAL_BUF(271)
    DELY = SCALAR_REAL_BUF(272)
    DIST = SCALAR_REAL_BUF(273)
    DELDIST = SCALAR_REAL_BUF(274)
    DELETA = SCALAR_REAL_BUF(275)
    ADVECX = SCALAR_REAL_BUF(276)
    ADVECY = SCALAR_REAL_BUF(277)
    AGIRD = SCALAR_REAL_BUF(278)
    AH = SCALAR_REAL_BUF(279)
    AO12 = SCALAR_REAL_BUF(280)
    AO6 = SCALAR_REAL_BUF(281)
    ARG = SCALAR_REAL_BUF(282)
    ARG1 = SCALAR_REAL_BUF(283)
    ARG2 = SCALAR_REAL_BUF(284)
    ARGJ = SCALAR_REAL_BUF(285)
    ARGJ1 = SCALAR_REAL_BUF(286)
    ARGJ2 = SCALAR_REAL_BUF(287)
    ARGSALT = SCALAR_REAL_BUF(288)
    ARGT = SCALAR_REAL_BUF(289)
    ARGTP = SCALAR_REAL_BUF(290)
    AUV21 = SCALAR_REAL_BUF(291)
    AUV22 = SCALAR_REAL_BUF(292)
    BARAVGWT = SCALAR_REAL_BUF(293)
    BEDSTR = SCALAR_REAL_BUF(294)
    BNDLEN2O3NC = SCALAR_REAL_BUF(295)
    BSXN1 = SCALAR_REAL_BUF(296)
    BSXN2 = SCALAR_REAL_BUF(297)
    BSXN3 = SCALAR_REAL_BUF(298)
    BSXPP3 = SCALAR_REAL_BUF(299)
    BSYN1 = SCALAR_REAL_BUF(300)
    BSYN2 = SCALAR_REAL_BUF(301)
    BSYN3 = SCALAR_REAL_BUF(302)
    BSYPP3 = SCALAR_REAL_BUF(303)
    C1 = SCALAR_REAL_BUF(304)
    C2 = SCALAR_REAL_BUF(305)
    C3 = SCALAR_REAL_BUF(306)
    CBEDSTRD = SCALAR_REAL_BUF(307)
    CBEDSTRE = SCALAR_REAL_BUF(308)
    CCRITD = SCALAR_REAL_BUF(309)
    CCSFEA = SCALAR_REAL_BUF(310)
    CELERITY = SCALAR_REAL_BUF(311)
    CH1N1 = SCALAR_REAL_BUF(312)
    CH1N2 = SCALAR_REAL_BUF(313)
    CH1N3 = SCALAR_REAL_BUF(314)
    CHSUM = SCALAR_REAL_BUF(315)
    COND = SCALAR_REAL_BUF(316)
    CONVCR = SCALAR_REAL_BUF(317)
    CORIFPP = SCALAR_REAL_BUF(318)
    DDU = SCALAR_REAL_BUF(319)
    DHDX = SCALAR_REAL_BUF(320)
    DHDY = SCALAR_REAL_BUF(321)
    DISPERX = SCALAR_REAL_BUF(322)
    DISPERY = SCALAR_REAL_BUF(323)
    DT2 = SCALAR_REAL_BUF(324)
    DTO2 = SCALAR_REAL_BUF(325)
    DTOHPP = SCALAR_REAL_BUF(326)
    DUU1N1 = SCALAR_REAL_BUF(327)
    DUU1N2 = SCALAR_REAL_BUF(328)
    DUU1N3 = SCALAR_REAL_BUF(329)
    DUV1N1 = SCALAR_REAL_BUF(330)
    DUV1N2 = SCALAR_REAL_BUF(331)
    DUV1N3 = SCALAR_REAL_BUF(332)
    DVV1N1 = SCALAR_REAL_BUF(333)
    DVV1N2 = SCALAR_REAL_BUF(334)
    DVV1N3 = SCALAR_REAL_BUF(335)
    DXXYY11 = SCALAR_REAL_BUF(336)
    DXXYY12 = SCALAR_REAL_BUF(337)
    DXXYY13 = SCALAR_REAL_BUF(338)
    DXXYY21 = SCALAR_REAL_BUF(339)
    DXXYY22 = SCALAR_REAL_BUF(340)
    DXXYY23 = SCALAR_REAL_BUF(341)
    DXXYY31 = SCALAR_REAL_BUF(342)
    DXXYY32 = SCALAR_REAL_BUF(343)
    DXXYY33 = SCALAR_REAL_BUF(344)
    DXYH11 = SCALAR_REAL_BUF(345)
    DXYH12 = SCALAR_REAL_BUF(346)
    DXYH13 = SCALAR_REAL_BUF(347)
    DXYH21 = SCALAR_REAL_BUF(348)
    DXYH22 = SCALAR_REAL_BUF(349)
    DXYH23 = SCALAR_REAL_BUF(350)
    DXYH31 = SCALAR_REAL_BUF(351)
    DXYH32 = SCALAR_REAL_BUF(352)
    DXYH33 = SCALAR_REAL_BUF(353)
    E0N1 = SCALAR_REAL_BUF(354)
    E0N2 = SCALAR_REAL_BUF(355)
    E0N3 = SCALAR_REAL_BUF(356)
    E1N1 = SCALAR_REAL_BUF(357)
    E1N1SQ = SCALAR_REAL_BUF(358)
    E1N2 = SCALAR_REAL_BUF(359)
    E1N2SQ = SCALAR_REAL_BUF(360)
    E1N3 = SCALAR_REAL_BUF(361)
    E1N3SQ = SCALAR_REAL_BUF(362)
    ECONST = SCALAR_REAL_BUF(363)
    EE1 = SCALAR_REAL_BUF(364)
    EE2 = SCALAR_REAL_BUF(365)
    EE3 = SCALAR_REAL_BUF(366)
    ELMAX = SCALAR_REAL_BUF(367)
    EP = SCALAR_REAL_BUF(368)
    ESN1 = SCALAR_REAL_BUF(369)
    ESN2 = SCALAR_REAL_BUF(370)
    BE1 = SCALAR_REAL_BUF(371)
    BE2 = SCALAR_REAL_BUF(372)
    BE3 = SCALAR_REAL_BUF(373)
    ESN3 = SCALAR_REAL_BUF(374)
    ETIME1 = SCALAR_REAL_BUF(375)
    ETIME2 = SCALAR_REAL_BUF(376)
    ETRATIO = SCALAR_REAL_BUF(377)
    EVC1 = SCALAR_REAL_BUF(378)
    EVC2 = SCALAR_REAL_BUF(379)
    EVC3 = SCALAR_REAL_BUF(380)
    EVCEA = SCALAR_REAL_BUF(381)
    EVMPPODT = SCALAR_REAL_BUF(382)
    EVMPPDT = SCALAR_REAL_BUF(383)
    FDDD = SCALAR_REAL_BUF(384)
    FDDDODT = SCALAR_REAL_BUF(385)
    FDDOD = SCALAR_REAL_BUF(386)
    FDDODODT = SCALAR_REAL_BUF(387)
    FIIN = SCALAR_REAL_BUF(388)
    G = SCALAR_REAL_BUF(389)
    GA00 = SCALAR_REAL_BUF(390)
    GB00A00 = SCALAR_REAL_BUF(391)
    GC00 = SCALAR_REAL_BUF(392)
    GDTO2 = SCALAR_REAL_BUF(393)
    GFAO2 = SCALAR_REAL_BUF(394)
    GHPP = SCALAR_REAL_BUF(395)
    GO3 = SCALAR_REAL_BUF(396)
    HABSMIN = SCALAR_REAL_BUF(397)
    HEA = SCALAR_REAL_BUF(398)
    HH1 = SCALAR_REAL_BUF(399)
    HH1N1 = SCALAR_REAL_BUF(400)
    HH1N2 = SCALAR_REAL_BUF(401)
    HH1N3 = SCALAR_REAL_BUF(402)
    HH2 = SCALAR_REAL_BUF(403)
    HH2N1 = SCALAR_REAL_BUF(404)
    HH2N2 = SCALAR_REAL_BUF(405)
    HH2N3 = SCALAR_REAL_BUF(406)
    HHU1N1 = SCALAR_REAL_BUF(407)
    HHU1N2 = SCALAR_REAL_BUF(408)
    HHU1N3 = SCALAR_REAL_BUF(409)
    HHV1N1 = SCALAR_REAL_BUF(410)
    HHV1N2 = SCALAR_REAL_BUF(411)
    HHV1N3 = SCALAR_REAL_BUF(412)
    HPP = SCALAR_REAL_BUF(413)
    HSD = SCALAR_REAL_BUF(414)
    HSE = SCALAR_REAL_BUF(415)
    HTOT = SCALAR_REAL_BUF(416)
    G2ROOT = SCALAR_REAL_BUF(417)
    P11 = SCALAR_REAL_BUF(418)
    P22 = SCALAR_REAL_BUF(419)
    P33 = SCALAR_REAL_BUF(420)
    PR1N1 = SCALAR_REAL_BUF(421)
    PR1N2 = SCALAR_REAL_BUF(422)
    PR1N3 = SCALAR_REAL_BUF(423)
    QFORCEI = SCALAR_REAL_BUF(424)
    QFORCEJ = SCALAR_REAL_BUF(425)
    QTEMA1 = SCALAR_REAL_BUF(426)
    QTEMA2 = SCALAR_REAL_BUF(427)
    QTEMA3 = SCALAR_REAL_BUF(428)
    QTEMB1 = SCALAR_REAL_BUF(429)
    QTEMB2 = SCALAR_REAL_BUF(430)
    QTEMB3 = SCALAR_REAL_BUF(431)
    QTRATIO = SCALAR_REAL_BUF(432)
    QUNORM = SCALAR_REAL_BUF(433)
    QVNORM = SCALAR_REAL_BUF(434)
    RAMP1 = SCALAR_REAL_BUF(435)
    RAMP2 = SCALAR_REAL_BUF(436)
    RBARWL = SCALAR_REAL_BUF(437)
    RBARWL1 = SCALAR_REAL_BUF(438)
    RBARWL1F = SCALAR_REAL_BUF(439)
    RBARWL2 = SCALAR_REAL_BUF(440)
    RBARWL2F = SCALAR_REAL_BUF(441)
    RFF = SCALAR_REAL_BUF(442)
    RFF1 = SCALAR_REAL_BUF(443)
    RFF2 = SCALAR_REAL_BUF(444)
    RHO0 = SCALAR_REAL_BUF(445)
    RSTRATIO = SCALAR_REAL_BUF(446)
    RSX = SCALAR_REAL_BUF(447)
    RSY = SCALAR_REAL_BUF(448)
    S2SFEA = SCALAR_REAL_BUF(449)
    SADVDTO3 = SCALAR_REAL_BUF(450)
    SALTMUL = SCALAR_REAL_BUF(451)
    SFACPP = SCALAR_REAL_BUF(452)
    SS1N1 = SCALAR_REAL_BUF(453)
    SS1N2 = SCALAR_REAL_BUF(454)
    SS1N3 = SCALAR_REAL_BUF(455)
    T0N1 = SCALAR_REAL_BUF(456)
    T0N2 = SCALAR_REAL_BUF(457)
    T0N3 = SCALAR_REAL_BUF(458)
    T0XN1 = SCALAR_REAL_BUF(459)
    T0XN2 = SCALAR_REAL_BUF(460)
    T0XN3 = SCALAR_REAL_BUF(461)
    T0XPP3 = SCALAR_REAL_BUF(462)
    T0YN1 = SCALAR_REAL_BUF(463)
    T0YN2 = SCALAR_REAL_BUF(464)
    T0YN3 = SCALAR_REAL_BUF(465)
    T0YPP3 = SCALAR_REAL_BUF(466)
    TADVODT = SCALAR_REAL_BUF(467)
    TAU0AVG = SCALAR_REAL_BUF(468)
    THENALLDSSSTUP = SCALAR_REAL_BUF(469)
    TIMEIT = SCALAR_REAL_BUF(470)
    TIPN1 = SCALAR_REAL_BUF(471)
    TOUTFC = SCALAR_REAL_BUF(472)
    TIPN2 = SCALAR_REAL_BUF(473)
    TIPN3 = SCALAR_REAL_BUF(474)
    TKWET = SCALAR_REAL_BUF(475)
    TOUTFGC = SCALAR_REAL_BUF(476)
    TOUTFGE = SCALAR_REAL_BUF(477)
    TOUTFGV = SCALAR_REAL_BUF(478)
    TOUTFGW = SCALAR_REAL_BUF(479)
    TOUTFM = SCALAR_REAL_BUF(480)
    TOUTSGC = SCALAR_REAL_BUF(481)
    TOUTSGE = SCALAR_REAL_BUF(482)
    TOUTSGV = SCALAR_REAL_BUF(483)
    TOUTSGW = SCALAR_REAL_BUF(484)
    TOUTSM = SCALAR_REAL_BUF(485)
    TPMUL = SCALAR_REAL_BUF(486)
    TT0L = SCALAR_REAL_BUF(487)
    TT0R = SCALAR_REAL_BUF(488)
    U11 = SCALAR_REAL_BUF(489)
    U1N1 = SCALAR_REAL_BUF(490)
    U1N2 = SCALAR_REAL_BUF(491)
    U1N3 = SCALAR_REAL_BUF(492)
    U22 = SCALAR_REAL_BUF(493)
    U33 = SCALAR_REAL_BUF(494)
    UEA = SCALAR_REAL_BUF(495)
    UHPP = SCALAR_REAL_BUF(496)
    UHPP3 = SCALAR_REAL_BUF(497)
    UN1 = SCALAR_REAL_BUF(498)
    UPEA = SCALAR_REAL_BUF(499)
    UPP = SCALAR_REAL_BUF(500)
    UPPDT = SCALAR_REAL_BUF(501)
    UPPDTDDX1 = SCALAR_REAL_BUF(502)
    UPPDTDDX2 = SCALAR_REAL_BUF(503)
    UPPDTDDX3 = SCALAR_REAL_BUF(504)
    UV1 = SCALAR_REAL_BUF(505)
    V11 = SCALAR_REAL_BUF(506)
    V1N1 = SCALAR_REAL_BUF(507)
    UVW1 = SCALAR_REAL_BUF(508)
    UVW2 = SCALAR_REAL_BUF(509)
    COSA1 = SCALAR_REAL_BUF(510)
    SINA1 = SCALAR_REAL_BUF(511)
    WD1 = SCALAR_REAL_BUF(512)
    WD = SCALAR_REAL_BUF(513)
    WDXX = SCALAR_REAL_BUF(514)
    WDYY = SCALAR_REAL_BUF(515)
    WDXY = SCALAR_REAL_BUF(516)
    CHYBR = SCALAR_REAL_BUF(517)
    VCOEFXX = SCALAR_REAL_BUF(518)
    VCOEFYY = SCALAR_REAL_BUF(519)
    VCOEFXY = SCALAR_REAL_BUF(520)
    VCOEFYX = SCALAR_REAL_BUF(521)
    VCOEF2 = SCALAR_REAL_BUF(522)
    V1N2 = SCALAR_REAL_BUF(523)
    V1N3 = SCALAR_REAL_BUF(524)
    V22 = SCALAR_REAL_BUF(525)
    V33 = SCALAR_REAL_BUF(526)
    VCOEF3N1 = SCALAR_REAL_BUF(527)
    VCOEF3N2 = SCALAR_REAL_BUF(528)
    VCOEF3N3 = SCALAR_REAL_BUF(529)
    VCOEF3X = SCALAR_REAL_BUF(530)
    VCOEF3Y = SCALAR_REAL_BUF(531)
    VEA = SCALAR_REAL_BUF(532)
    VEL = SCALAR_REAL_BUF(533)
    VELABS = SCALAR_REAL_BUF(534)
    VELMAX = SCALAR_REAL_BUF(535)
    VELNORM = SCALAR_REAL_BUF(536)
    VELTAN = SCALAR_REAL_BUF(537)
    VHPP = SCALAR_REAL_BUF(538)
    VHPP3 = SCALAR_REAL_BUF(539)
    VPEA = SCALAR_REAL_BUF(540)
    VPP = SCALAR_REAL_BUF(541)
    VPPDT = SCALAR_REAL_BUF(542)
    VPPDTDDY1 = SCALAR_REAL_BUF(543)
    VPPDTDDY2 = SCALAR_REAL_BUF(544)
    VPPDTDDY3 = SCALAR_REAL_BUF(545)
    WDRAGCO = SCALAR_REAL_BUF(546)
    WINDMAG = SCALAR_REAL_BUF(547)
    WINDX = SCALAR_REAL_BUF(548)
    WINDY = SCALAR_REAL_BUF(549)
    WS = SCALAR_REAL_BUF(550)
    WSMOD = SCALAR_REAL_BUF(551)
    WSX = SCALAR_REAL_BUF(552)
    WSXN1 = SCALAR_REAL_BUF(553)
    WSXN2 = SCALAR_REAL_BUF(554)
    WSXN3 = SCALAR_REAL_BUF(555)
    WSY = SCALAR_REAL_BUF(556)
    WSYN1 = SCALAR_REAL_BUF(557)
    WSYN2 = SCALAR_REAL_BUF(558)
    WSYN3 = SCALAR_REAL_BUF(559)
    WTRATIO = SCALAR_REAL_BUF(560)
    A00 = SCALAR_REAL_BUF(561)
    B00 = SCALAR_REAL_BUF(562)
    C00 = SCALAR_REAL_BUF(563)
    ANGINN = SCALAR_REAL_BUF(564)
    CORI = SCALAR_REAL_BUF(565)
    COSTHETA = SCALAR_REAL_BUF(566)
    COSTHETA1 = SCALAR_REAL_BUF(567)
    COSTSET = SCALAR_REAL_BUF(568)
    CROSS = SCALAR_REAL_BUF(569)
    CROSS1 = SCALAR_REAL_BUF(570)
    DAY = SCALAR_REAL_BUF(571)
    DOTVEC = SCALAR_REAL_BUF(572)
    DRAMP = SCALAR_REAL_BUF(573)
    DUM1 = SCALAR_REAL_BUF(574)
    DUM2 = SCALAR_REAL_BUF(575)
    EVMSUM = SCALAR_REAL_BUF(576)
    H0 = SCALAR_REAL_BUF(577)
    H0L = SCALAR_REAL_BUF(578)
    H0H = SCALAR_REAL_BUF(579)
    RNDAY = SCALAR_REAL_BUF(580)
    THETA = SCALAR_REAL_BUF(581)
    THETA1 = SCALAR_REAL_BUF(582)
    TOUTSC = SCALAR_REAL_BUF(583)
    RAMP = SCALAR_REAL_BUF(584)
    RHOWAT0 = SCALAR_REAL_BUF(585)
    TOUTSE = SCALAR_REAL_BUF(586)
    TOUTFE = SCALAR_REAL_BUF(587)
    TOUTSV = SCALAR_REAL_BUF(588)
    TOUTFV = SCALAR_REAL_BUF(589)
    XL = SCALAR_REAL_BUF(590)
    VECNORM = SCALAR_REAL_BUF(591)
    VL1X = SCALAR_REAL_BUF(592)
    VL1Y = SCALAR_REAL_BUF(593)
    VL2X = SCALAR_REAL_BUF(594)
    VL2Y = SCALAR_REAL_BUF(595)
    WLATMAX = SCALAR_REAL_BUF(596)
    WLONMIN = SCALAR_REAL_BUF(597)
    WLATINC = SCALAR_REAL_BUF(598)
    WLONINC = SCALAR_REAL_BUF(599)
    VELMIN = SCALAR_REAL_BUF(600)
    RNP_GLOBAL = SCALAR_REAL_BUF(601)
    REFSEC = SCALAR_REAL_BUF(602)
    IF (SCALAR_INT_BUF(1) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 0)
      CALL c_f_pointer(curr_ptr, WDFLG, shape=[fstarpu_vector_get_nx(buffers, 0)])
    ELSE
      NULLIFY(WDFLG)
    END IF
    IF (SCALAR_INT_BUF(2) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 1)
      CALL c_f_pointer(curr_ptr, WDFLG_TMP, shape=[fstarpu_vector_get_nx(buffers, 1)])
    ELSE
      NULLIFY(WDFLG_TMP)
    END IF
    IF (SCALAR_INT_BUF(3) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 2)
      CALL c_f_pointer(curr_ptr, DOFW, shape=[fstarpu_vector_get_nx(buffers, 2)])
    ELSE
      NULLIFY(DOFW)
    END IF
    IF (SCALAR_INT_BUF(4) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 3)
      CALL c_f_pointer(curr_ptr, EL_UPDATED, shape=[fstarpu_vector_get_nx(buffers, 3)])
    ELSE
      NULLIFY(EL_UPDATED)
    END IF
    IF (SCALAR_INT_BUF(5) == 1) THEN
      curr_ptr = fstarpu_block_get_ptr(buffers, 4)
      CALL c_f_pointer(curr_ptr, LEDGE_NVEC, shape=[fstarpu_block_get_nx(buffers, 4), fstarpu_block_get_ny(buffers, 4), fstarpu_block_get_nz(buffers, 4)])
    ELSE
      NULLIFY(LEDGE_NVEC)
    END IF
    IF (SCALAR_INT_BUF(6) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 5)
      CALL c_f_pointer(curr_ptr, DOFS, shape=[fstarpu_vector_get_nx(buffers, 5)])
    ELSE
      NULLIFY(DOFS)
    END IF
    IF (SCALAR_INT_BUF(7) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 6)
      CALL c_f_pointer(curr_ptr, PCOUNT, shape=[fstarpu_vector_get_nx(buffers, 6)])
    ELSE
      NULLIFY(PCOUNT)
    END IF
    IF (SCALAR_INT_BUF(8) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 7)
      CALL c_f_pointer(curr_ptr, PDG, shape=[fstarpu_vector_get_nx(buffers, 7)])
    ELSE
      NULLIFY(PDG)
    END IF
    IF (SCALAR_INT_BUF(9) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 8)
      CALL c_f_pointer(curr_ptr, NCOUNT, shape=[fstarpu_vector_get_nx(buffers, 8)])
    ELSE
      NULLIFY(NCOUNT)
    END IF
    IF (SCALAR_INT_BUF(10) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 9)
      CALL c_f_pointer(curr_ptr, NEDEL, shape=[fstarpu_matrix_get_nx(buffers, 9), fstarpu_matrix_get_ny(buffers, 9)])
    ELSE
      NULLIFY(NEDEL)
    END IF
    IF (SCALAR_INT_BUF(11) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 10)
      CALL c_f_pointer(curr_ptr, NEDSD, shape=[fstarpu_matrix_get_nx(buffers, 10), fstarpu_matrix_get_ny(buffers, 10)])
    ELSE
      NULLIFY(NEDSD)
    END IF
    IF (SCALAR_INT_BUF(12) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 11)
      CALL c_f_pointer(curr_ptr, NEDNO, shape=[fstarpu_matrix_get_nx(buffers, 11), fstarpu_matrix_get_ny(buffers, 11)])
    ELSE
      NULLIFY(NEDNO)
    END IF
    IF (SCALAR_INT_BUF(13) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 12)
      CALL c_f_pointer(curr_ptr, NEDNO1, shape=[fstarpu_vector_get_nx(buffers, 12)])
    ELSE
      NULLIFY(NEDNO1)
    END IF
    IF (SCALAR_INT_BUF(14) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 13)
      CALL c_f_pointer(curr_ptr, NEDNO2, shape=[fstarpu_vector_get_nx(buffers, 13)])
    ELSE
      NULLIFY(NEDNO2)
    END IF
    IF (SCALAR_INT_BUF(15) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 14)
      CALL c_f_pointer(curr_ptr, NIEDN, shape=[fstarpu_vector_get_nx(buffers, 14)])
    ELSE
      NULLIFY(NIEDN)
    END IF
    IF (SCALAR_INT_BUF(16) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 15)
      CALL c_f_pointer(curr_ptr, NLEDN, shape=[fstarpu_vector_get_nx(buffers, 15)])
    ELSE
      NULLIFY(NLEDN)
    END IF
    IF (SCALAR_INT_BUF(17) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 16)
      CALL c_f_pointer(curr_ptr, NEEDN, shape=[fstarpu_vector_get_nx(buffers, 16)])
    ELSE
      NULLIFY(NEEDN)
    END IF
    IF (SCALAR_INT_BUF(18) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 17)
      CALL c_f_pointer(curr_ptr, NFEDN, shape=[fstarpu_vector_get_nx(buffers, 17)])
    ELSE
      NULLIFY(NFEDN)
    END IF
    IF (SCALAR_INT_BUF(19) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 18)
      CALL c_f_pointer(curr_ptr, NREDN, shape=[fstarpu_vector_get_nx(buffers, 18)])
    ELSE
      NULLIFY(NREDN)
    END IF
    IF (SCALAR_INT_BUF(20) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 19)
      CALL c_f_pointer(curr_ptr, NEBEDN, shape=[fstarpu_vector_get_nx(buffers, 19)])
    ELSE
      NULLIFY(NEBEDN)
    END IF
    IF (SCALAR_INT_BUF(21) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 20)
      CALL c_f_pointer(curr_ptr, NIBEDN, shape=[fstarpu_vector_get_nx(buffers, 20)])
    ELSE
      NULLIFY(NIBEDN)
    END IF
    IF (SCALAR_INT_BUF(22) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 21)
      CALL c_f_pointer(curr_ptr, NIBSEGN, shape=[fstarpu_matrix_get_nx(buffers, 21), fstarpu_matrix_get_ny(buffers, 21)])
    ELSE
      NULLIFY(NIBSEGN)
    END IF
    IF (SCALAR_INT_BUF(23) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 22)
      CALL c_f_pointer(curr_ptr, NEBSEGN, shape=[fstarpu_vector_get_nx(buffers, 22)])
    ELSE
      NULLIFY(NEBSEGN)
    END IF
    IF (SCALAR_INT_BUF(24) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 23)
      CALL c_f_pointer(curr_ptr, EL_NBORS, shape=[fstarpu_matrix_get_nx(buffers, 23), fstarpu_matrix_get_ny(buffers, 23)])
    ELSE
      NULLIFY(EL_NBORS)
    END IF
    IF (SCALAR_INT_BUF(25) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 24)
      CALL c_f_pointer(curr_ptr, BACKNODES, shape=[fstarpu_matrix_get_nx(buffers, 24), fstarpu_matrix_get_ny(buffers, 24)])
    ELSE
      NULLIFY(BACKNODES)
    END IF
    IF (SCALAR_INT_BUF(26) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 25)
      CALL c_f_pointer(curr_ptr, MARK, shape=[fstarpu_vector_get_nx(buffers, 25)])
    ELSE
      NULLIFY(MARK)
    END IF
    IF (SCALAR_INT_BUF(27) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 26)
      CALL c_f_pointer(curr_ptr, ATVD, shape=[fstarpu_matrix_get_nx(buffers, 26), fstarpu_matrix_get_ny(buffers, 26)])
    ELSE
      NULLIFY(ATVD)
    END IF
    IF (SCALAR_INT_BUF(28) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 27)
      CALL c_f_pointer(curr_ptr, BTVD, shape=[fstarpu_matrix_get_nx(buffers, 27), fstarpu_matrix_get_ny(buffers, 27)])
    ELSE
      NULLIFY(BTVD)
    END IF
    IF (SCALAR_INT_BUF(29) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 28)
      CALL c_f_pointer(curr_ptr, CTVD, shape=[fstarpu_matrix_get_nx(buffers, 28), fstarpu_matrix_get_ny(buffers, 28)])
    ELSE
      NULLIFY(CTVD)
    END IF
    IF (SCALAR_INT_BUF(30) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 29)
      CALL c_f_pointer(curr_ptr, DTVD, shape=[fstarpu_vector_get_nx(buffers, 29)])
    ELSE
      NULLIFY(DTVD)
    END IF
    IF (SCALAR_INT_BUF(31) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 30)
      CALL c_f_pointer(curr_ptr, MAX_BOA_DT, shape=[fstarpu_vector_get_nx(buffers, 30)])
    ELSE
      NULLIFY(MAX_BOA_DT)
    END IF
    IF (SCALAR_INT_BUF(32) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 31)
      CALL c_f_pointer(curr_ptr, e1, shape=[fstarpu_vector_get_nx(buffers, 31)])
    ELSE
      NULLIFY(e1)
    END IF
    IF (SCALAR_INT_BUF(33) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 32)
      CALL c_f_pointer(curr_ptr, balance, shape=[fstarpu_vector_get_nx(buffers, 32)])
    ELSE
      NULLIFY(balance)
    END IF
    IF (SCALAR_INT_BUF(34) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 33)
      CALL c_f_pointer(curr_ptr, RKC_T, shape=[fstarpu_vector_get_nx(buffers, 33)])
    ELSE
      NULLIFY(RKC_T)
    END IF
    IF (SCALAR_INT_BUF(35) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 34)
      CALL c_f_pointer(curr_ptr, RKC_U, shape=[fstarpu_vector_get_nx(buffers, 34)])
    ELSE
      NULLIFY(RKC_U)
    END IF
    IF (SCALAR_INT_BUF(36) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 35)
      CALL c_f_pointer(curr_ptr, RKC_Tprime, shape=[fstarpu_vector_get_nx(buffers, 35)])
    ELSE
      NULLIFY(RKC_Tprime)
    END IF
    IF (SCALAR_INT_BUF(37) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 36)
      CALL c_f_pointer(curr_ptr, RKC_Tdprime, shape=[fstarpu_vector_get_nx(buffers, 36)])
    ELSE
      NULLIFY(RKC_Tdprime)
    END IF
    IF (SCALAR_INT_BUF(38) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 37)
      CALL c_f_pointer(curr_ptr, RKC_a, shape=[fstarpu_vector_get_nx(buffers, 37)])
    ELSE
      NULLIFY(RKC_a)
    END IF
    IF (SCALAR_INT_BUF(39) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 38)
      CALL c_f_pointer(curr_ptr, RKC_b, shape=[fstarpu_vector_get_nx(buffers, 38)])
    ELSE
      NULLIFY(RKC_b)
    END IF
    IF (SCALAR_INT_BUF(40) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 39)
      CALL c_f_pointer(curr_ptr, RKC_c, shape=[fstarpu_vector_get_nx(buffers, 39)])
    ELSE
      NULLIFY(RKC_c)
    END IF
    IF (SCALAR_INT_BUF(41) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 40)
      CALL c_f_pointer(curr_ptr, RKC_mu, shape=[fstarpu_vector_get_nx(buffers, 40)])
    ELSE
      NULLIFY(RKC_mu)
    END IF
    IF (SCALAR_INT_BUF(42) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 41)
      CALL c_f_pointer(curr_ptr, RKC_tildemu, shape=[fstarpu_vector_get_nx(buffers, 41)])
    ELSE
      NULLIFY(RKC_tildemu)
    END IF
    IF (SCALAR_INT_BUF(43) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 42)
      CALL c_f_pointer(curr_ptr, RKC_nu, shape=[fstarpu_vector_get_nx(buffers, 42)])
    ELSE
      NULLIFY(RKC_nu)
    END IF
    IF (SCALAR_INT_BUF(44) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 43)
      CALL c_f_pointer(curr_ptr, RKC_gamma, shape=[fstarpu_vector_get_nx(buffers, 43)])
    ELSE
      NULLIFY(RKC_gamma)
    END IF
    IF (SCALAR_INT_BUF(45) == 1) THEN
      curr_ptr = fstarpu_block_get_ptr(buffers, 44)
      CALL c_f_pointer(curr_ptr, BATH, shape=[fstarpu_block_get_nx(buffers, 44), fstarpu_block_get_ny(buffers, 44), fstarpu_block_get_nz(buffers, 44)])
    ELSE
      NULLIFY(BATH)
    END IF
    IF (SCALAR_INT_BUF(46) == 1) THEN
      curr_ptr = fstarpu_block_get_ptr(buffers, 45)
      CALL c_f_pointer(curr_ptr, DBATHDX, shape=[fstarpu_block_get_nx(buffers, 45), fstarpu_block_get_ny(buffers, 45), fstarpu_block_get_nz(buffers, 45)])
    ELSE
      NULLIFY(DBATHDX)
    END IF
    IF (SCALAR_INT_BUF(47) == 1) THEN
      curr_ptr = fstarpu_block_get_ptr(buffers, 46)
      CALL c_f_pointer(curr_ptr, DBATHDY, shape=[fstarpu_block_get_nx(buffers, 46), fstarpu_block_get_ny(buffers, 46), fstarpu_block_get_nz(buffers, 46)])
    ELSE
      NULLIFY(DBATHDY)
    END IF
    IF (SCALAR_INT_BUF(48) == 1) THEN
      curr_ptr = fstarpu_block_get_ptr(buffers, 47)
      CALL c_f_pointer(curr_ptr, SFAC_ELEM, shape=[fstarpu_block_get_nx(buffers, 47), fstarpu_block_get_ny(buffers, 47), fstarpu_block_get_nz(buffers, 47)])
    ELSE
      NULLIFY(SFAC_ELEM)
    END IF
    IF (SCALAR_INT_BUF(49) == 1) THEN
      curr_ptr = fstarpu_tensor_get_ptr(buffers, 48)
      CALL c_f_pointer(curr_ptr, BATHED, shape=[fstarpu_tensor_get_nx(buffers, 48), fstarpu_tensor_get_ny(buffers, 48), fstarpu_tensor_get_nz(buffers, 48), fstarpu_tensor_get_nt(buffers, 48)])
    ELSE
      NULLIFY(BATHED)
    END IF
    IF (SCALAR_INT_BUF(50) == 1) THEN
      curr_ptr = fstarpu_tensor_get_ptr(buffers, 49)
      CALL c_f_pointer(curr_ptr, SFACED, shape=[fstarpu_tensor_get_nx(buffers, 49), fstarpu_tensor_get_ny(buffers, 49), fstarpu_tensor_get_nz(buffers, 49), fstarpu_tensor_get_nt(buffers, 49)])
    ELSE
      NULLIFY(SFACED)
    END IF
    IF (SCALAR_INT_BUF(51) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 50)
      CALL c_f_pointer(curr_ptr, COSNX, shape=[fstarpu_vector_get_nx(buffers, 50)])
    ELSE
      NULLIFY(COSNX)
    END IF
    IF (SCALAR_INT_BUF(52) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 51)
      CALL c_f_pointer(curr_ptr, SINNX, shape=[fstarpu_vector_get_nx(buffers, 51)])
    ELSE
      NULLIFY(SINNX)
    END IF
    IF (SCALAR_INT_BUF(53) == 1) THEN
      curr_ptr = fstarpu_block_get_ptr(buffers, 52)
      CALL c_f_pointer(curr_ptr, DP_NODE, shape=[fstarpu_block_get_nx(buffers, 52), fstarpu_block_get_ny(buffers, 52), fstarpu_block_get_nz(buffers, 52)])
    ELSE
      NULLIFY(DP_NODE)
    END IF
    IF (SCALAR_INT_BUF(54) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 53)
      CALL c_f_pointer(curr_ptr, DP_VOL, shape=[fstarpu_matrix_get_nx(buffers, 53), fstarpu_matrix_get_ny(buffers, 53)])
    ELSE
      NULLIFY(DP_VOL)
    END IF
    IF (SCALAR_INT_BUF(55) == 1) THEN
      curr_ptr = fstarpu_block_get_ptr(buffers, 54)
      CALL c_f_pointer(curr_ptr, DRPHI, shape=[fstarpu_block_get_nx(buffers, 54), fstarpu_block_get_ny(buffers, 54), fstarpu_block_get_nz(buffers, 54)])
    ELSE
      NULLIFY(DRPHI)
    END IF
    IF (SCALAR_INT_BUF(56) == 1) THEN
      curr_ptr = fstarpu_block_get_ptr(buffers, 55)
      CALL c_f_pointer(curr_ptr, DSPHI, shape=[fstarpu_block_get_nx(buffers, 55), fstarpu_block_get_ny(buffers, 55), fstarpu_block_get_nz(buffers, 55)])
    ELSE
      NULLIFY(DSPHI)
    END IF
    IF (SCALAR_INT_BUF(57) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 56)
      CALL c_f_pointer(curr_ptr, DRDX, shape=[fstarpu_vector_get_nx(buffers, 56)])
    ELSE
      NULLIFY(DRDX)
    END IF
    IF (SCALAR_INT_BUF(58) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 57)
      CALL c_f_pointer(curr_ptr, DSDX, shape=[fstarpu_vector_get_nx(buffers, 57)])
    ELSE
      NULLIFY(DSDX)
    END IF
    IF (SCALAR_INT_BUF(59) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 58)
      CALL c_f_pointer(curr_ptr, DRDY, shape=[fstarpu_vector_get_nx(buffers, 58)])
    ELSE
      NULLIFY(DRDY)
    END IF
    IF (SCALAR_INT_BUF(60) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 59)
      CALL c_f_pointer(curr_ptr, DSDY, shape=[fstarpu_vector_get_nx(buffers, 59)])
    ELSE
      NULLIFY(DSDY)
    END IF
    IF (SCALAR_INT_BUF(61) == 1) THEN
      curr_ptr = fstarpu_block_get_ptr(buffers, 60)
      CALL c_f_pointer(curr_ptr, DXPHI2, shape=[fstarpu_block_get_nx(buffers, 60), fstarpu_block_get_ny(buffers, 60), fstarpu_block_get_nz(buffers, 60)])
    ELSE
      NULLIFY(DXPHI2)
    END IF
    IF (SCALAR_INT_BUF(62) == 1) THEN
      curr_ptr = fstarpu_block_get_ptr(buffers, 61)
      CALL c_f_pointer(curr_ptr, DYPHI2, shape=[fstarpu_block_get_nx(buffers, 61), fstarpu_block_get_ny(buffers, 61), fstarpu_block_get_nz(buffers, 61)])
    ELSE
      NULLIFY(DYPHI2)
    END IF
    IF (SCALAR_INT_BUF(63) == 1) THEN
      curr_ptr = fstarpu_block_get_ptr(buffers, 62)
      CALL c_f_pointer(curr_ptr, PHI2, shape=[fstarpu_block_get_nx(buffers, 62), fstarpu_block_get_ny(buffers, 62), fstarpu_block_get_nz(buffers, 62)])
    ELSE
      NULLIFY(PHI2)
    END IF
    IF (SCALAR_INT_BUF(64) == 1) THEN
      curr_ptr = fstarpu_block_get_ptr(buffers, 63)
      CALL c_f_pointer(curr_ptr, EFA_DG, shape=[fstarpu_block_get_nx(buffers, 63), fstarpu_block_get_ny(buffers, 63), fstarpu_block_get_nz(buffers, 63)])
    ELSE
      NULLIFY(EFA_DG)
    END IF
    IF (SCALAR_INT_BUF(65) == 1) THEN
      curr_ptr = fstarpu_block_get_ptr(buffers, 64)
      CALL c_f_pointer(curr_ptr, EMO_DG, shape=[fstarpu_block_get_nx(buffers, 64), fstarpu_block_get_ny(buffers, 64), fstarpu_block_get_nz(buffers, 64)])
    ELSE
      NULLIFY(EMO_DG)
    END IF
    IF (SCALAR_INT_BUF(66) == 1) THEN
      curr_ptr = fstarpu_block_get_ptr(buffers, 65)
      CALL c_f_pointer(curr_ptr, UFA_DG, shape=[fstarpu_block_get_nx(buffers, 65), fstarpu_block_get_ny(buffers, 65), fstarpu_block_get_nz(buffers, 65)])
    ELSE
      NULLIFY(UFA_DG)
    END IF
    IF (SCALAR_INT_BUF(67) == 1) THEN
      curr_ptr = fstarpu_block_get_ptr(buffers, 66)
      CALL c_f_pointer(curr_ptr, UMO_DG, shape=[fstarpu_block_get_nx(buffers, 66), fstarpu_block_get_ny(buffers, 66), fstarpu_block_get_nz(buffers, 66)])
    ELSE
      NULLIFY(UMO_DG)
    END IF
    IF (SCALAR_INT_BUF(68) == 1) THEN
      curr_ptr = fstarpu_block_get_ptr(buffers, 67)
      CALL c_f_pointer(curr_ptr, VFA_DG, shape=[fstarpu_block_get_nx(buffers, 67), fstarpu_block_get_ny(buffers, 67), fstarpu_block_get_nz(buffers, 67)])
    ELSE
      NULLIFY(VFA_DG)
    END IF
    IF (SCALAR_INT_BUF(69) == 1) THEN
      curr_ptr = fstarpu_block_get_ptr(buffers, 68)
      CALL c_f_pointer(curr_ptr, VMO_DG, shape=[fstarpu_block_get_nx(buffers, 68), fstarpu_block_get_ny(buffers, 68), fstarpu_block_get_nz(buffers, 68)])
    ELSE
      NULLIFY(VMO_DG)
    END IF
    IF (SCALAR_INT_BUF(70) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 69)
      CALL c_f_pointer(curr_ptr, XLEN, shape=[fstarpu_vector_get_nx(buffers, 69)])
    ELSE
      NULLIFY(XLEN)
    END IF
    IF (SCALAR_INT_BUF(71) == 1) THEN
      curr_ptr = fstarpu_block_get_ptr(buffers, 70)
      CALL c_f_pointer(curr_ptr, HB, shape=[fstarpu_block_get_nx(buffers, 70), fstarpu_block_get_ny(buffers, 70), fstarpu_block_get_nz(buffers, 70)])
    ELSE
      NULLIFY(HB)
    END IF
    IF (SCALAR_INT_BUF(72) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 71)
      CALL c_f_pointer(curr_ptr, MANN, shape=[fstarpu_matrix_get_nx(buffers, 71), fstarpu_matrix_get_ny(buffers, 71)])
    ELSE
      NULLIFY(MANN)
    END IF
    IF (SCALAR_INT_BUF(73) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 72)
      CALL c_f_pointer(curr_ptr, IBHT, shape=[fstarpu_vector_get_nx(buffers, 72)])
    ELSE
      NULLIFY(IBHT)
    END IF
    IF (SCALAR_INT_BUF(74) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 73)
      CALL c_f_pointer(curr_ptr, EBHT, shape=[fstarpu_vector_get_nx(buffers, 73)])
    ELSE
      NULLIFY(EBHT)
    END IF
    IF (SCALAR_INT_BUF(75) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 74)
      CALL c_f_pointer(curr_ptr, EBCFSP, shape=[fstarpu_vector_get_nx(buffers, 74)])
    ELSE
      NULLIFY(EBCFSP)
    END IF
    IF (SCALAR_INT_BUF(76) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 75)
      CALL c_f_pointer(curr_ptr, IBCFSP, shape=[fstarpu_vector_get_nx(buffers, 75)])
    ELSE
      NULLIFY(IBCFSP)
    END IF
    IF (SCALAR_INT_BUF(77) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 76)
      CALL c_f_pointer(curr_ptr, IBCFSB, shape=[fstarpu_vector_get_nx(buffers, 76)])
    ELSE
      NULLIFY(IBCFSB)
    END IF
    IF (SCALAR_INT_BUF(78) == 1) THEN
      curr_ptr = fstarpu_tensor_get_ptr(buffers, 77)
      CALL c_f_pointer(curr_ptr, JACOBI, shape=[fstarpu_tensor_get_nx(buffers, 77), fstarpu_tensor_get_ny(buffers, 77), fstarpu_tensor_get_nz(buffers, 77), fstarpu_tensor_get_nt(buffers, 77)])
    ELSE
      NULLIFY(JACOBI)
    END IF
    IF (SCALAR_INT_BUF(79) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 78)
      CALL c_f_pointer(curr_ptr, M_INV, shape=[fstarpu_matrix_get_nx(buffers, 78), fstarpu_matrix_get_ny(buffers, 78)])
    ELSE
      NULLIFY(M_INV)
    END IF
    IF (SCALAR_INT_BUF(80) == 1) THEN
      curr_ptr = fstarpu_block_get_ptr(buffers, 79)
      CALL c_f_pointer(curr_ptr, phi_edge_fixed, shape=[fstarpu_block_get_nx(buffers, 79), fstarpu_block_get_ny(buffers, 79), fstarpu_block_get_nz(buffers, 79)])
    ELSE
      NULLIFY(phi_edge_fixed)
    END IF
    IF (SCALAR_INT_BUF(81) == 1) THEN
      curr_ptr = fstarpu_block_get_ptr(buffers, 80)
      CALL c_f_pointer(curr_ptr, PHI_AREA, shape=[fstarpu_block_get_nx(buffers, 80), fstarpu_block_get_ny(buffers, 80), fstarpu_block_get_nz(buffers, 80)])
    ELSE
      NULLIFY(PHI_AREA)
    END IF
    IF (SCALAR_INT_BUF(82) == 1) THEN
      curr_ptr = fstarpu_tensor_get_ptr(buffers, 81)
      CALL c_f_pointer(curr_ptr, PHI_EDGE, shape=[fstarpu_tensor_get_nx(buffers, 81), fstarpu_tensor_get_ny(buffers, 81), fstarpu_tensor_get_nz(buffers, 81), fstarpu_tensor_get_nt(buffers, 81)])
    ELSE
      NULLIFY(PHI_EDGE)
    END IF
    IF (SCALAR_INT_BUF(83) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 82)
      CALL c_f_pointer(curr_ptr, PHI_CENTER, shape=[fstarpu_matrix_get_nx(buffers, 82), fstarpu_matrix_get_ny(buffers, 82)])
    ELSE
      NULLIFY(PHI_CENTER)
    END IF
    IF (SCALAR_INT_BUF(84) == 1) THEN
      curr_ptr = fstarpu_block_get_ptr(buffers, 83)
      CALL c_f_pointer(curr_ptr, PHI_CORNER, shape=[fstarpu_block_get_nx(buffers, 83), fstarpu_block_get_ny(buffers, 83), fstarpu_block_get_nz(buffers, 83)])
    ELSE
      NULLIFY(PHI_CORNER)
    END IF
    IF (SCALAR_INT_BUF(85) == 1) THEN
      curr_ptr = fstarpu_block_get_ptr(buffers, 84)
      CALL c_f_pointer(curr_ptr, PHI_CHECK, shape=[fstarpu_block_get_nx(buffers, 84), fstarpu_block_get_ny(buffers, 84), fstarpu_block_get_nz(buffers, 84)])
    ELSE
      NULLIFY(PHI_CHECK)
    END IF
    IF (SCALAR_INT_BUF(86) == 1) THEN
      curr_ptr = fstarpu_tensor_get_ptr(buffers, 85)
      CALL c_f_pointer(curr_ptr, PHI_CORNER1, shape=[fstarpu_tensor_get_nx(buffers, 85), fstarpu_tensor_get_ny(buffers, 85), fstarpu_tensor_get_nz(buffers, 85), fstarpu_tensor_get_nt(buffers, 85)])
    ELSE
      NULLIFY(PHI_CORNER1)
    END IF
    IF (SCALAR_INT_BUF(87) == 1) THEN
      curr_ptr = fstarpu_block_get_ptr(buffers, 86)
      CALL c_f_pointer(curr_ptr, PHI_MID, shape=[fstarpu_block_get_nx(buffers, 86), fstarpu_block_get_ny(buffers, 86), fstarpu_block_get_nz(buffers, 86)])
    ELSE
      NULLIFY(PHI_MID)
    END IF
    IF (SCALAR_INT_BUF(88) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 87)
      CALL c_f_pointer(curr_ptr, PHI_INTEGRATED, shape=[fstarpu_matrix_get_nx(buffers, 87), fstarpu_matrix_get_ny(buffers, 87)])
    ELSE
      NULLIFY(PHI_INTEGRATED)
    END IF
    IF (SCALAR_INT_BUF(89) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 88)
      CALL c_f_pointer(curr_ptr, PSI_CHECK, shape=[fstarpu_matrix_get_nx(buffers, 88), fstarpu_matrix_get_ny(buffers, 88)])
    ELSE
      NULLIFY(PSI_CHECK)
    END IF
    IF (SCALAR_INT_BUF(90) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 89)
      CALL c_f_pointer(curr_ptr, PSI1, shape=[fstarpu_matrix_get_nx(buffers, 89), fstarpu_matrix_get_ny(buffers, 89)])
    ELSE
      NULLIFY(PSI1)
    END IF
    IF (SCALAR_INT_BUF(91) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 90)
      CALL c_f_pointer(curr_ptr, PSI2, shape=[fstarpu_matrix_get_nx(buffers, 90), fstarpu_matrix_get_ny(buffers, 90)])
    ELSE
      NULLIFY(PSI2)
    END IF
    IF (SCALAR_INT_BUF(92) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 91)
      CALL c_f_pointer(curr_ptr, PSI3, shape=[fstarpu_matrix_get_nx(buffers, 91), fstarpu_matrix_get_ny(buffers, 91)])
    ELSE
      NULLIFY(PSI3)
    END IF
    IF (SCALAR_INT_BUF(93) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 92)
      CALL c_f_pointer(curr_ptr, Q_HAT, shape=[fstarpu_vector_get_nx(buffers, 92)])
    ELSE
      NULLIFY(Q_HAT)
    END IF
    IF (SCALAR_INT_BUF(94) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 93)
      CALL c_f_pointer(curr_ptr, QIB, shape=[fstarpu_vector_get_nx(buffers, 93)])
    ELSE
      NULLIFY(QIB)
    END IF
    IF (SCALAR_INT_BUF(95) == 1) THEN
      curr_ptr = fstarpu_block_get_ptr(buffers, 94)
      CALL c_f_pointer(curr_ptr, QX, shape=[fstarpu_block_get_nx(buffers, 94), fstarpu_block_get_ny(buffers, 94), fstarpu_block_get_nz(buffers, 94)])
    ELSE
      NULLIFY(QX)
    END IF
    IF (SCALAR_INT_BUF(96) == 1) THEN
      curr_ptr = fstarpu_block_get_ptr(buffers, 95)
      CALL c_f_pointer(curr_ptr, QY, shape=[fstarpu_block_get_nx(buffers, 95), fstarpu_block_get_ny(buffers, 95), fstarpu_block_get_nz(buffers, 95)])
    ELSE
      NULLIFY(QY)
    END IF
    IF (SCALAR_INT_BUF(97) == 1) THEN
      curr_ptr = fstarpu_block_get_ptr(buffers, 96)
      CALL c_f_pointer(curr_ptr, ZE, shape=[fstarpu_block_get_nx(buffers, 96), fstarpu_block_get_ny(buffers, 96), fstarpu_block_get_nz(buffers, 96)])
    ELSE
      NULLIFY(ZE)
    END IF
    IF (SCALAR_INT_BUF(98) == 1) THEN
      curr_ptr = fstarpu_block_get_ptr(buffers, 97)
      CALL c_f_pointer(curr_ptr, ze_edge, shape=[fstarpu_block_get_nx(buffers, 97), fstarpu_block_get_ny(buffers, 97), fstarpu_block_get_nz(buffers, 97)])
    ELSE
      NULLIFY(ze_edge)
    END IF
    IF (SCALAR_INT_BUF(99) == 1) THEN
      curr_ptr = fstarpu_block_get_ptr(buffers, 98)
      CALL c_f_pointer(curr_ptr, qx_edge, shape=[fstarpu_block_get_nx(buffers, 98), fstarpu_block_get_ny(buffers, 98), fstarpu_block_get_nz(buffers, 98)])
    ELSE
      NULLIFY(qx_edge)
    END IF
    IF (SCALAR_INT_BUF(100) == 1) THEN
      curr_ptr = fstarpu_block_get_ptr(buffers, 99)
      CALL c_f_pointer(curr_ptr, qy_edge, shape=[fstarpu_block_get_nx(buffers, 99), fstarpu_block_get_ny(buffers, 99), fstarpu_block_get_nz(buffers, 99)])
    ELSE
      NULLIFY(qy_edge)
    END IF
    IF (SCALAR_INT_BUF(101) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 100)
      CALL c_f_pointer(curr_ptr, elem_edge, shape=[fstarpu_matrix_get_nx(buffers, 100), fstarpu_matrix_get_ny(buffers, 100)])
    ELSE
      NULLIFY(elem_edge)
    END IF
    IF (SCALAR_INT_BUF(102) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 101)
      CALL c_f_pointer(curr_ptr, nieds_count, shape=[fstarpu_vector_get_nx(buffers, 101)])
    ELSE
      NULLIFY(nieds_count)
    END IF
    IF (SCALAR_INT_BUF(103) == 1) THEN
      curr_ptr = fstarpu_tensor_get_ptr(buffers, 102)
      CALL c_f_pointer(curr_ptr, bed, shape=[fstarpu_tensor_get_nx(buffers, 102), fstarpu_tensor_get_ny(buffers, 102), fstarpu_tensor_get_nz(buffers, 102), fstarpu_tensor_get_nt(buffers, 102)])
    ELSE
      NULLIFY(bed)
    END IF
    IF (SCALAR_INT_BUF(104) == 1) THEN
      curr_ptr = fstarpu_block_get_ptr(buffers, 103)
      CALL c_f_pointer(curr_ptr, dynP, shape=[fstarpu_block_get_nx(buffers, 103), fstarpu_block_get_ny(buffers, 103), fstarpu_block_get_nz(buffers, 103)])
    ELSE
      NULLIFY(dynP)
    END IF
    IF (SCALAR_INT_BUF(105) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 104)
      CALL c_f_pointer(curr_ptr, dynP_MAX, shape=[fstarpu_vector_get_nx(buffers, 104)])
    ELSE
      NULLIFY(dynP_MAX)
    END IF
    IF (SCALAR_INT_BUF(106) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 105)
      CALL c_f_pointer(curr_ptr, dynP_MIN, shape=[fstarpu_vector_get_nx(buffers, 105)])
    ELSE
      NULLIFY(dynP_MIN)
    END IF
    IF (SCALAR_INT_BUF(107) == 1) THEN
      curr_ptr = fstarpu_block_get_ptr(buffers, 106)
      CALL c_f_pointer(curr_ptr, iota, shape=[fstarpu_block_get_nx(buffers, 106), fstarpu_block_get_ny(buffers, 106), fstarpu_block_get_nz(buffers, 106)])
    ELSE
      NULLIFY(iota)
    END IF
    IF (SCALAR_INT_BUF(108) == 1) THEN
      curr_ptr = fstarpu_block_get_ptr(buffers, 107)
      CALL c_f_pointer(curr_ptr, iotaa, shape=[fstarpu_block_get_nx(buffers, 107), fstarpu_block_get_ny(buffers, 107), fstarpu_block_get_nz(buffers, 107)])
    ELSE
      NULLIFY(iotaa)
    END IF
    IF (SCALAR_INT_BUF(109) == 1) THEN
      curr_ptr = fstarpu_block_get_ptr(buffers, 108)
      CALL c_f_pointer(curr_ptr, iota2, shape=[fstarpu_block_get_nx(buffers, 108), fstarpu_block_get_ny(buffers, 108), fstarpu_block_get_nz(buffers, 108)])
    ELSE
      NULLIFY(iota2)
    END IF
    IF (SCALAR_INT_BUF(110) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 109)
      CALL c_f_pointer(curr_ptr, iota_MAX, shape=[fstarpu_vector_get_nx(buffers, 109)])
    ELSE
      NULLIFY(iota_MAX)
    END IF
    IF (SCALAR_INT_BUF(111) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 110)
      CALL c_f_pointer(curr_ptr, iota_MIN, shape=[fstarpu_vector_get_nx(buffers, 110)])
    ELSE
      NULLIFY(iota_MIN)
    END IF
    IF (SCALAR_INT_BUF(112) == 1) THEN
      curr_ptr = fstarpu_block_get_ptr(buffers, 111)
      CALL c_f_pointer(curr_ptr, iotaa2, shape=[fstarpu_block_get_nx(buffers, 111), fstarpu_block_get_ny(buffers, 111), fstarpu_block_get_nz(buffers, 111)])
    ELSE
      NULLIFY(iotaa2)
    END IF
    IF (SCALAR_INT_BUF(113) == 1) THEN
      curr_ptr = fstarpu_block_get_ptr(buffers, 112)
      CALL c_f_pointer(curr_ptr, iotaa3, shape=[fstarpu_block_get_nx(buffers, 112), fstarpu_block_get_ny(buffers, 112), fstarpu_block_get_nz(buffers, 112)])
    ELSE
      NULLIFY(iotaa3)
    END IF
    IF (SCALAR_INT_BUF(114) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 113)
      CALL c_f_pointer(curr_ptr, iota2_MAX, shape=[fstarpu_vector_get_nx(buffers, 113)])
    ELSE
      NULLIFY(iota2_MAX)
    END IF
    IF (SCALAR_INT_BUF(115) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 114)
      CALL c_f_pointer(curr_ptr, iota2_MIN, shape=[fstarpu_vector_get_nx(buffers, 114)])
    ELSE
      NULLIFY(iota2_MIN)
    END IF
    IF (SCALAR_INT_BUF(116) == 1) THEN
      curr_ptr = fstarpu_block_get_ptr(buffers, 115)
      CALL c_f_pointer(curr_ptr, arrayfix, shape=[fstarpu_block_get_nx(buffers, 115), fstarpu_block_get_ny(buffers, 115), fstarpu_block_get_nz(buffers, 115)])
    ELSE
      NULLIFY(arrayfix)
    END IF
    IF (SCALAR_INT_BUF(117) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 116)
      CALL c_f_pointer(curr_ptr, CORI_EL, shape=[fstarpu_vector_get_nx(buffers, 116)])
    ELSE
      NULLIFY(CORI_EL)
    END IF
    IF (SCALAR_INT_BUF(118) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 117)
      CALL c_f_pointer(curr_ptr, FRIC_EL, shape=[fstarpu_vector_get_nx(buffers, 117)])
    ELSE
      NULLIFY(FRIC_EL)
    END IF
    IF (SCALAR_INT_BUF(119) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 118)
      CALL c_f_pointer(curr_ptr, ZE_MAX, shape=[fstarpu_vector_get_nx(buffers, 118)])
    ELSE
      NULLIFY(ZE_MAX)
    END IF
    IF (SCALAR_INT_BUF(120) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 119)
      CALL c_f_pointer(curr_ptr, ZE_MIN, shape=[fstarpu_vector_get_nx(buffers, 119)])
    ELSE
      NULLIFY(ZE_MIN)
    END IF
    IF (SCALAR_INT_BUF(121) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 120)
      CALL c_f_pointer(curr_ptr, DPE_MIN, shape=[fstarpu_vector_get_nx(buffers, 120)])
    ELSE
      NULLIFY(DPE_MIN)
    END IF
    IF (SCALAR_INT_BUF(122) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 121)
      CALL c_f_pointer(curr_ptr, WATER_DEPTH_OLD, shape=[fstarpu_matrix_get_nx(buffers, 121), fstarpu_matrix_get_ny(buffers, 121)])
    ELSE
      NULLIFY(WATER_DEPTH_OLD)
    END IF
    IF (SCALAR_INT_BUF(123) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 122)
      CALL c_f_pointer(curr_ptr, WATER_DEPTH, shape=[fstarpu_matrix_get_nx(buffers, 122), fstarpu_matrix_get_ny(buffers, 122)])
    ELSE
      NULLIFY(WATER_DEPTH)
    END IF
    IF (SCALAR_INT_BUF(124) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 123)
      CALL c_f_pointer(curr_ptr, ADVECTQX, shape=[fstarpu_vector_get_nx(buffers, 123)])
    ELSE
      NULLIFY(ADVECTQX)
    END IF
    IF (SCALAR_INT_BUF(125) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 124)
      CALL c_f_pointer(curr_ptr, ADVECTQY, shape=[fstarpu_vector_get_nx(buffers, 124)])
    ELSE
      NULLIFY(ADVECTQY)
    END IF
    IF (SCALAR_INT_BUF(126) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 125)
      CALL c_f_pointer(curr_ptr, SOURCEQX, shape=[fstarpu_vector_get_nx(buffers, 125)])
    ELSE
      NULLIFY(SOURCEQX)
    END IF
    IF (SCALAR_INT_BUF(127) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 126)
      CALL c_f_pointer(curr_ptr, SOURCEQY, shape=[fstarpu_vector_get_nx(buffers, 126)])
    ELSE
      NULLIFY(SOURCEQY)
    END IF
    IF (SCALAR_INT_BUF(128) == 1) THEN
      curr_ptr = fstarpu_tensor_get_ptr(buffers, 127)
      CALL c_f_pointer(curr_ptr, LZ, shape=[fstarpu_tensor_get_nx(buffers, 127), fstarpu_tensor_get_ny(buffers, 127), fstarpu_tensor_get_nz(buffers, 127), fstarpu_tensor_get_nt(buffers, 127)])
    ELSE
      NULLIFY(LZ)
    END IF
    IF (SCALAR_INT_BUF(129) == 1) THEN
      curr_ptr = fstarpu_tensor_get_ptr(buffers, 128)
      CALL c_f_pointer(curr_ptr, MZ, shape=[fstarpu_tensor_get_nx(buffers, 128), fstarpu_tensor_get_ny(buffers, 128), fstarpu_tensor_get_nz(buffers, 128), fstarpu_tensor_get_nt(buffers, 128)])
    ELSE
      NULLIFY(MZ)
    END IF
    IF (SCALAR_INT_BUF(130) == 1) THEN
      curr_ptr = fstarpu_tensor_get_ptr(buffers, 129)
      CALL c_f_pointer(curr_ptr, HZ, shape=[fstarpu_tensor_get_nx(buffers, 129), fstarpu_tensor_get_ny(buffers, 129), fstarpu_tensor_get_nz(buffers, 129), fstarpu_tensor_get_nt(buffers, 129)])
    ELSE
      NULLIFY(HZ)
    END IF
    IF (SCALAR_INT_BUF(131) == 1) THEN
      curr_ptr = fstarpu_tensor_get_ptr(buffers, 130)
      CALL c_f_pointer(curr_ptr, TZ, shape=[fstarpu_tensor_get_nx(buffers, 130), fstarpu_tensor_get_ny(buffers, 130), fstarpu_tensor_get_nz(buffers, 130), fstarpu_tensor_get_nt(buffers, 130)])
    ELSE
      NULLIFY(TZ)
    END IF
    IF (SCALAR_INT_BUF(132) == 1) THEN
      curr_ptr = fstarpu_block_get_ptr(buffers, 131)
      CALL c_f_pointer(curr_ptr, QNAM_DG, shape=[fstarpu_block_get_nx(buffers, 131), fstarpu_block_get_ny(buffers, 131), fstarpu_block_get_nz(buffers, 131)])
    ELSE
      NULLIFY(QNAM_DG)
    END IF
    IF (SCALAR_INT_BUF(133) == 1) THEN
      curr_ptr = fstarpu_block_get_ptr(buffers, 132)
      CALL c_f_pointer(curr_ptr, QNPH_DG, shape=[fstarpu_block_get_nx(buffers, 132), fstarpu_block_get_ny(buffers, 132), fstarpu_block_get_nz(buffers, 132)])
    ELSE
      NULLIFY(QNPH_DG)
    END IF
    IF (SCALAR_INT_BUF(134) == 1) THEN
      curr_ptr = fstarpu_block_get_ptr(buffers, 133)
      CALL c_f_pointer(curr_ptr, RHS_ZE, shape=[fstarpu_block_get_nx(buffers, 133), fstarpu_block_get_ny(buffers, 133), fstarpu_block_get_nz(buffers, 133)])
    ELSE
      NULLIFY(RHS_ZE)
    END IF
    IF (SCALAR_INT_BUF(135) == 1) THEN
      curr_ptr = fstarpu_tensor_get_ptr(buffers, 134)
      CALL c_f_pointer(curr_ptr, RHS_bed, shape=[fstarpu_tensor_get_nx(buffers, 134), fstarpu_tensor_get_ny(buffers, 134), fstarpu_tensor_get_nz(buffers, 134), fstarpu_tensor_get_nt(buffers, 134)])
    ELSE
      NULLIFY(RHS_bed)
    END IF
    IF (SCALAR_INT_BUF(136) == 1) THEN
      curr_ptr = fstarpu_block_get_ptr(buffers, 135)
      CALL c_f_pointer(curr_ptr, RHS_QX, shape=[fstarpu_block_get_nx(buffers, 135), fstarpu_block_get_ny(buffers, 135), fstarpu_block_get_nz(buffers, 135)])
    ELSE
      NULLIFY(RHS_QX)
    END IF
    IF (SCALAR_INT_BUF(137) == 1) THEN
      curr_ptr = fstarpu_block_get_ptr(buffers, 136)
      CALL c_f_pointer(curr_ptr, RHS_QY, shape=[fstarpu_block_get_nx(buffers, 136), fstarpu_block_get_ny(buffers, 136), fstarpu_block_get_nz(buffers, 136)])
    ELSE
      NULLIFY(RHS_QY)
    END IF
    IF (SCALAR_INT_BUF(138) == 1) THEN
      curr_ptr = fstarpu_block_get_ptr(buffers, 137)
      CALL c_f_pointer(curr_ptr, RHS_iota, shape=[fstarpu_block_get_nx(buffers, 137), fstarpu_block_get_ny(buffers, 137), fstarpu_block_get_nz(buffers, 137)])
    ELSE
      NULLIFY(RHS_iota)
    END IF
    IF (SCALAR_INT_BUF(139) == 1) THEN
      curr_ptr = fstarpu_block_get_ptr(buffers, 138)
      CALL c_f_pointer(curr_ptr, RHS_iota2, shape=[fstarpu_block_get_nx(buffers, 138), fstarpu_block_get_ny(buffers, 138), fstarpu_block_get_nz(buffers, 138)])
    ELSE
      NULLIFY(RHS_iota2)
    END IF
    IF (SCALAR_INT_BUF(140) == 1) THEN
      curr_ptr = fstarpu_block_get_ptr(buffers, 139)
      CALL c_f_pointer(curr_ptr, RHS_dynP, shape=[fstarpu_block_get_nx(buffers, 139), fstarpu_block_get_ny(buffers, 139), fstarpu_block_get_nz(buffers, 139)])
    ELSE
      NULLIFY(RHS_dynP)
    END IF
    IF (SCALAR_INT_BUF(141) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 140)
      CALL c_f_pointer(curr_ptr, RHS_bed_IN, shape=[fstarpu_matrix_get_nx(buffers, 140), fstarpu_matrix_get_ny(buffers, 140)])
    ELSE
      NULLIFY(RHS_bed_IN)
    END IF
    IF (SCALAR_INT_BUF(142) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 141)
      CALL c_f_pointer(curr_ptr, RHS_bed_EX, shape=[fstarpu_matrix_get_nx(buffers, 141), fstarpu_matrix_get_ny(buffers, 141)])
    ELSE
      NULLIFY(RHS_bed_EX)
    END IF
    IF (SCALAR_INT_BUF(143) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 142)
      CALL c_f_pointer(curr_ptr, bed_HAT_O, shape=[fstarpu_vector_get_nx(buffers, 142)])
    ELSE
      NULLIFY(bed_HAT_O)
    END IF
    IF (SCALAR_INT_BUF(144) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 143)
      CALL c_f_pointer(curr_ptr, XAGP, shape=[fstarpu_matrix_get_nx(buffers, 143), fstarpu_matrix_get_ny(buffers, 143)])
    ELSE
      NULLIFY(XAGP)
    END IF
    IF (SCALAR_INT_BUF(145) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 144)
      CALL c_f_pointer(curr_ptr, YAGP, shape=[fstarpu_matrix_get_nx(buffers, 144), fstarpu_matrix_get_ny(buffers, 144)])
    ELSE
      NULLIFY(YAGP)
    END IF
    IF (SCALAR_INT_BUF(146) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 145)
      CALL c_f_pointer(curr_ptr, WAGP, shape=[fstarpu_matrix_get_nx(buffers, 145), fstarpu_matrix_get_ny(buffers, 145)])
    ELSE
      NULLIFY(WAGP)
    END IF
    IF (SCALAR_INT_BUF(147) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 146)
      CALL c_f_pointer(curr_ptr, XEGP, shape=[fstarpu_matrix_get_nx(buffers, 146), fstarpu_matrix_get_ny(buffers, 146)])
    ELSE
      NULLIFY(XEGP)
    END IF
    IF (SCALAR_INT_BUF(148) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 147)
      CALL c_f_pointer(curr_ptr, YEGP, shape=[fstarpu_matrix_get_nx(buffers, 147), fstarpu_matrix_get_ny(buffers, 147)])
    ELSE
      NULLIFY(YEGP)
    END IF
    IF (SCALAR_INT_BUF(149) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 148)
      CALL c_f_pointer(curr_ptr, WEGP, shape=[fstarpu_matrix_get_nx(buffers, 148), fstarpu_matrix_get_ny(buffers, 148)])
    ELSE
      NULLIFY(WEGP)
    END IF
    IF (SCALAR_INT_BUF(150) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 149)
      CALL c_f_pointer(curr_ptr, SL3, shape=[fstarpu_matrix_get_nx(buffers, 149), fstarpu_matrix_get_ny(buffers, 149)])
    ELSE
      NULLIFY(SL3)
    END IF
    IF (SCALAR_INT_BUF(151) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 150)
      CALL c_f_pointer(curr_ptr, XBC, shape=[fstarpu_vector_get_nx(buffers, 150)])
    ELSE
      NULLIFY(XBC)
    END IF
    IF (SCALAR_INT_BUF(152) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 151)
      CALL c_f_pointer(curr_ptr, YBC, shape=[fstarpu_vector_get_nx(buffers, 151)])
    ELSE
      NULLIFY(YBC)
    END IF
    IF (SCALAR_INT_BUF(153) == 1) THEN
      curr_ptr = fstarpu_tensor_get_ptr(buffers, 152)
      CALL c_f_pointer(curr_ptr, XFAC, shape=[fstarpu_tensor_get_nx(buffers, 152), fstarpu_tensor_get_ny(buffers, 152), fstarpu_tensor_get_nz(buffers, 152), fstarpu_tensor_get_nt(buffers, 152)])
    ELSE
      NULLIFY(XFAC)
    END IF
    IF (SCALAR_INT_BUF(154) == 1) THEN
      curr_ptr = fstarpu_tensor_get_ptr(buffers, 153)
      CALL c_f_pointer(curr_ptr, YFAC, shape=[fstarpu_tensor_get_nx(buffers, 153), fstarpu_tensor_get_ny(buffers, 153), fstarpu_tensor_get_nz(buffers, 153), fstarpu_tensor_get_nt(buffers, 153)])
    ELSE
      NULLIFY(YFAC)
    END IF
    IF (SCALAR_INT_BUF(155) == 1) THEN
      curr_ptr = fstarpu_tensor_get_ptr(buffers, 154)
      CALL c_f_pointer(curr_ptr, SRFAC, shape=[fstarpu_tensor_get_nx(buffers, 154), fstarpu_tensor_get_ny(buffers, 154), fstarpu_tensor_get_nz(buffers, 154), fstarpu_tensor_get_nt(buffers, 154)])
    ELSE
      NULLIFY(SRFAC)
    END IF
    IF (SCALAR_INT_BUF(156) == 1) THEN
      curr_ptr = fstarpu_tensor_get_ptr(buffers, 155)
      CALL c_f_pointer(curr_ptr, EDGEQ, shape=[fstarpu_tensor_get_nx(buffers, 155), fstarpu_tensor_get_ny(buffers, 155), fstarpu_tensor_get_nz(buffers, 155), fstarpu_tensor_get_nt(buffers, 155)])
    ELSE
      NULLIFY(EDGEQ)
    END IF
    IF (SCALAR_INT_BUF(157) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 156)
      CALL c_f_pointer(curr_ptr, PHI, shape=[fstarpu_vector_get_nx(buffers, 156)])
    ELSE
      NULLIFY(PHI)
    END IF
    IF (SCALAR_INT_BUF(158) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 157)
      CALL c_f_pointer(curr_ptr, DPHIDZ1, shape=[fstarpu_vector_get_nx(buffers, 157)])
    ELSE
      NULLIFY(DPHIDZ1)
    END IF
    IF (SCALAR_INT_BUF(159) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 158)
      CALL c_f_pointer(curr_ptr, DPHIDZ2, shape=[fstarpu_vector_get_nx(buffers, 158)])
    ELSE
      NULLIFY(DPHIDZ2)
    END IF
    IF (SCALAR_INT_BUF(160) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 159)
      CALL c_f_pointer(curr_ptr, PHI_STAE, shape=[fstarpu_matrix_get_nx(buffers, 159), fstarpu_matrix_get_ny(buffers, 159)])
    ELSE
      NULLIFY(PHI_STAE)
    END IF
    IF (SCALAR_INT_BUF(161) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 160)
      CALL c_f_pointer(curr_ptr, PHI_STAV, shape=[fstarpu_matrix_get_nx(buffers, 160), fstarpu_matrix_get_ny(buffers, 160)])
    ELSE
      NULLIFY(PHI_STAV)
    END IF
    IF (SCALAR_INT_BUF(162) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 161)
      CALL c_f_pointer(curr_ptr, bed_IN, shape=[fstarpu_vector_get_nx(buffers, 161)])
    ELSE
      NULLIFY(bed_IN)
    END IF
    IF (SCALAR_INT_BUF(163) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 162)
      CALL c_f_pointer(curr_ptr, bed_EX, shape=[fstarpu_vector_get_nx(buffers, 162)])
    ELSE
      NULLIFY(bed_EX)
    END IF
    IF (SCALAR_INT_BUF(164) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 163)
      CALL c_f_pointer(curr_ptr, bed_HAT, shape=[fstarpu_vector_get_nx(buffers, 163)])
    ELSE
      NULLIFY(bed_HAT)
    END IF
    IF (SCALAR_INT_BUF(165) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 164)
      CALL c_f_pointer(curr_ptr, fact, shape=[fstarpu_vector_get_nx(buffers, 164)])
    ELSE
      NULLIFY(fact)
    END IF
    IF (SCALAR_INT_BUF(166) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 165)
      CALL c_f_pointer(curr_ptr, focal_neigh, shape=[fstarpu_matrix_get_nx(buffers, 165), fstarpu_matrix_get_ny(buffers, 165)])
    ELSE
      NULLIFY(focal_neigh)
    END IF
    IF (SCALAR_INT_BUF(167) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 166)
      CALL c_f_pointer(curr_ptr, focal_up, shape=[fstarpu_vector_get_nx(buffers, 166)])
    ELSE
      NULLIFY(focal_up)
    END IF
    IF (SCALAR_INT_BUF(168) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 167)
      CALL c_f_pointer(curr_ptr, bi, shape=[fstarpu_vector_get_nx(buffers, 167)])
    ELSE
      NULLIFY(bi)
    END IF
    IF (SCALAR_INT_BUF(169) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 168)
      CALL c_f_pointer(curr_ptr, bj, shape=[fstarpu_vector_get_nx(buffers, 168)])
    ELSE
      NULLIFY(bj)
    END IF
    IF (SCALAR_INT_BUF(170) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 169)
      CALL c_f_pointer(curr_ptr, XBCb, shape=[fstarpu_vector_get_nx(buffers, 169)])
    ELSE
      NULLIFY(XBCb)
    END IF
    IF (SCALAR_INT_BUF(171) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 170)
      CALL c_f_pointer(curr_ptr, YBCb, shape=[fstarpu_vector_get_nx(buffers, 170)])
    ELSE
      NULLIFY(YBCb)
    END IF
    IF (SCALAR_INT_BUF(172) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 171)
      CALL c_f_pointer(curr_ptr, xi1, shape=[fstarpu_matrix_get_nx(buffers, 171), fstarpu_matrix_get_ny(buffers, 171)])
    ELSE
      NULLIFY(xi1)
    END IF
    IF (SCALAR_INT_BUF(173) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 172)
      CALL c_f_pointer(curr_ptr, xi2, shape=[fstarpu_matrix_get_nx(buffers, 172), fstarpu_matrix_get_ny(buffers, 172)])
    ELSE
      NULLIFY(xi2)
    END IF
    IF (SCALAR_INT_BUF(174) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 173)
      CALL c_f_pointer(curr_ptr, xtransform, shape=[fstarpu_matrix_get_nx(buffers, 173), fstarpu_matrix_get_ny(buffers, 173)])
    ELSE
      NULLIFY(xtransform)
    END IF
    IF (SCALAR_INT_BUF(175) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 174)
      CALL c_f_pointer(curr_ptr, ytransform, shape=[fstarpu_matrix_get_nx(buffers, 174), fstarpu_matrix_get_ny(buffers, 174)])
    ELSE
      NULLIFY(ytransform)
    END IF
    IF (SCALAR_INT_BUF(176) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 175)
      CALL c_f_pointer(curr_ptr, xi1BCb, shape=[fstarpu_vector_get_nx(buffers, 175)])
    ELSE
      NULLIFY(xi1BCb)
    END IF
    IF (SCALAR_INT_BUF(177) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 176)
      CALL c_f_pointer(curr_ptr, xi2BCb, shape=[fstarpu_vector_get_nx(buffers, 176)])
    ELSE
      NULLIFY(xi2BCb)
    END IF
    IF (SCALAR_INT_BUF(178) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 177)
      CALL c_f_pointer(curr_ptr, xi1vert, shape=[fstarpu_matrix_get_nx(buffers, 177), fstarpu_matrix_get_ny(buffers, 177)])
    ELSE
      NULLIFY(xi1vert)
    END IF
    IF (SCALAR_INT_BUF(179) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 178)
      CALL c_f_pointer(curr_ptr, xi2vert, shape=[fstarpu_matrix_get_nx(buffers, 178), fstarpu_matrix_get_ny(buffers, 178)])
    ELSE
      NULLIFY(xi2vert)
    END IF
    IF (SCALAR_INT_BUF(180) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 179)
      CALL c_f_pointer(curr_ptr, xtransformv, shape=[fstarpu_matrix_get_nx(buffers, 179), fstarpu_matrix_get_ny(buffers, 179)])
    ELSE
      NULLIFY(xtransformv)
    END IF
    IF (SCALAR_INT_BUF(181) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 180)
      CALL c_f_pointer(curr_ptr, ytransformv, shape=[fstarpu_matrix_get_nx(buffers, 180), fstarpu_matrix_get_ny(buffers, 180)])
    ELSE
      NULLIFY(ytransformv)
    END IF
    IF (SCALAR_INT_BUF(182) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 181)
      CALL c_f_pointer(curr_ptr, XBCv, shape=[fstarpu_matrix_get_nx(buffers, 181), fstarpu_matrix_get_ny(buffers, 181)])
    ELSE
      NULLIFY(XBCv)
    END IF
    IF (SCALAR_INT_BUF(183) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 182)
      CALL c_f_pointer(curr_ptr, YBCv, shape=[fstarpu_matrix_get_nx(buffers, 182), fstarpu_matrix_get_ny(buffers, 182)])
    ELSE
      NULLIFY(YBCv)
    END IF
    IF (SCALAR_INT_BUF(184) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 183)
      CALL c_f_pointer(curr_ptr, xi1BCv, shape=[fstarpu_matrix_get_nx(buffers, 183), fstarpu_matrix_get_ny(buffers, 183)])
    ELSE
      NULLIFY(xi1BCv)
    END IF
    IF (SCALAR_INT_BUF(185) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 184)
      CALL c_f_pointer(curr_ptr, xi2BCv, shape=[fstarpu_matrix_get_nx(buffers, 184), fstarpu_matrix_get_ny(buffers, 184)])
    ELSE
      NULLIFY(xi2BCv)
    END IF
    IF (SCALAR_INT_BUF(186) == 1) THEN
      curr_ptr = fstarpu_block_get_ptr(buffers, 185)
      CALL c_f_pointer(curr_ptr, Area_integral, shape=[fstarpu_block_get_nx(buffers, 185), fstarpu_block_get_ny(buffers, 185), fstarpu_block_get_nz(buffers, 185)])
    ELSE
      NULLIFY(Area_integral)
    END IF
    IF (SCALAR_INT_BUF(187) == 1) THEN
      curr_ptr = fstarpu_tensor_get_ptr(buffers, 186)
      CALL c_f_pointer(curr_ptr, f, shape=[fstarpu_tensor_get_nx(buffers, 186), fstarpu_tensor_get_ny(buffers, 186), fstarpu_tensor_get_nz(buffers, 186), fstarpu_tensor_get_nt(buffers, 186)])
    ELSE
      NULLIFY(f)
    END IF
    IF (SCALAR_INT_BUF(188) == 1) THEN
      curr_ptr = fstarpu_tensor_get_ptr(buffers, 187)
      CALL c_f_pointer(curr_ptr, g0, shape=[fstarpu_tensor_get_nx(buffers, 187), fstarpu_tensor_get_ny(buffers, 187), fstarpu_tensor_get_nz(buffers, 187), fstarpu_tensor_get_nt(buffers, 187)])
    ELSE
      NULLIFY(g0)
    END IF
    IF (SCALAR_INT_BUF(189) == 1) THEN
      curr_ptr = fstarpu_tensor_get_ptr(buffers, 188)
      CALL c_f_pointer(curr_ptr, varsigma0, shape=[fstarpu_tensor_get_nx(buffers, 188), fstarpu_tensor_get_ny(buffers, 188), fstarpu_tensor_get_nz(buffers, 188), fstarpu_tensor_get_nt(buffers, 188)])
    ELSE
      NULLIFY(varsigma0)
    END IF
    IF (SCALAR_INT_BUF(190) == 1) THEN
      curr_ptr = fstarpu_tensor_get_ptr(buffers, 189)
      CALL c_f_pointer(curr_ptr, fv, shape=[fstarpu_tensor_get_nx(buffers, 189), fstarpu_tensor_get_ny(buffers, 189), fstarpu_tensor_get_nz(buffers, 189), fstarpu_tensor_get_nt(buffers, 189)])
    ELSE
      NULLIFY(fv)
    END IF
    IF (SCALAR_INT_BUF(191) == 1) THEN
      curr_ptr = fstarpu_tensor_get_ptr(buffers, 190)
      CALL c_f_pointer(curr_ptr, g0v, shape=[fstarpu_tensor_get_nx(buffers, 190), fstarpu_tensor_get_ny(buffers, 190), fstarpu_tensor_get_nz(buffers, 190), fstarpu_tensor_get_nt(buffers, 190)])
    ELSE
      NULLIFY(g0v)
    END IF
    IF (SCALAR_INT_BUF(192) == 1) THEN
      curr_ptr = fstarpu_tensor_get_ptr(buffers, 191)
      CALL c_f_pointer(curr_ptr, varsigma0v, shape=[fstarpu_tensor_get_nx(buffers, 191), fstarpu_tensor_get_ny(buffers, 191), fstarpu_tensor_get_nz(buffers, 191), fstarpu_tensor_get_nt(buffers, 191)])
    ELSE
      NULLIFY(varsigma0v)
    END IF
    IF (SCALAR_INT_BUF(193) == 1) THEN
      curr_ptr = fstarpu_block_get_ptr(buffers, 192)
      CALL c_f_pointer(curr_ptr, var2sigmag, shape=[fstarpu_block_get_nx(buffers, 192), fstarpu_block_get_ny(buffers, 192), fstarpu_block_get_nz(buffers, 192)])
    ELSE
      NULLIFY(var2sigmag)
    END IF
    IF (SCALAR_INT_BUF(194) == 1) THEN
      curr_ptr = fstarpu_block_get_ptr(buffers, 193)
      CALL c_f_pointer(curr_ptr, var2sigmav, shape=[fstarpu_block_get_nx(buffers, 193), fstarpu_block_get_ny(buffers, 193), fstarpu_block_get_nz(buffers, 193)])
    ELSE
      NULLIFY(var2sigmav)
    END IF
    IF (SCALAR_INT_BUF(195) == 1) THEN
      curr_ptr = fstarpu_tensor_get_ptr(buffers, 194)
      CALL c_f_pointer(curr_ptr, Nmatrix, shape=[fstarpu_tensor_get_nx(buffers, 194), fstarpu_tensor_get_ny(buffers, 194), fstarpu_tensor_get_nz(buffers, 194), fstarpu_tensor_get_nt(buffers, 194)])
    ELSE
      NULLIFY(Nmatrix)
    END IF
    IF (SCALAR_INT_BUF(196) == 1) THEN
      curr_ptr = fstarpu_tensor_get_ptr(buffers, 195)
      CALL c_f_pointer(curr_ptr, NmatrixInv, shape=[fstarpu_tensor_get_nx(buffers, 195), fstarpu_tensor_get_ny(buffers, 195), fstarpu_tensor_get_nz(buffers, 195), fstarpu_tensor_get_nt(buffers, 195)])
    ELSE
      NULLIFY(NmatrixInv)
    END IF
    IF (SCALAR_INT_BUF(197) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 196)
      CALL c_f_pointer(curr_ptr, deltx, shape=[fstarpu_vector_get_nx(buffers, 196)])
    ELSE
      NULLIFY(deltx)
    END IF
    IF (SCALAR_INT_BUF(198) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 197)
      CALL c_f_pointer(curr_ptr, delty, shape=[fstarpu_vector_get_nx(buffers, 197)])
    ELSE
      NULLIFY(delty)
    END IF
    IF (SCALAR_INT_BUF(199) == 1) THEN
      curr_ptr = fstarpu_block_get_ptr(buffers, 198)
      CALL c_f_pointer(curr_ptr, pmatrix, shape=[fstarpu_block_get_nx(buffers, 198), fstarpu_block_get_ny(buffers, 198), fstarpu_block_get_nz(buffers, 198)])
    ELSE
      NULLIFY(pmatrix)
    END IF
    IF (SCALAR_INT_BUF(200) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 199)
      CALL c_f_pointer(curr_ptr, ZEmin, shape=[fstarpu_matrix_get_nx(buffers, 199), fstarpu_matrix_get_ny(buffers, 199)])
    ELSE
      NULLIFY(ZEmin)
    END IF
    IF (SCALAR_INT_BUF(201) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 200)
      CALL c_f_pointer(curr_ptr, ZEmax, shape=[fstarpu_matrix_get_nx(buffers, 200), fstarpu_matrix_get_ny(buffers, 200)])
    ELSE
      NULLIFY(ZEmax)
    END IF
    IF (SCALAR_INT_BUF(202) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 201)
      CALL c_f_pointer(curr_ptr, QXmin, shape=[fstarpu_matrix_get_nx(buffers, 201), fstarpu_matrix_get_ny(buffers, 201)])
    ELSE
      NULLIFY(QXmin)
    END IF
    IF (SCALAR_INT_BUF(203) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 202)
      CALL c_f_pointer(curr_ptr, QXmax, shape=[fstarpu_matrix_get_nx(buffers, 202), fstarpu_matrix_get_ny(buffers, 202)])
    ELSE
      NULLIFY(QXmax)
    END IF
    IF (SCALAR_INT_BUF(204) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 203)
      CALL c_f_pointer(curr_ptr, QYmin, shape=[fstarpu_matrix_get_nx(buffers, 203), fstarpu_matrix_get_ny(buffers, 203)])
    ELSE
      NULLIFY(QYmin)
    END IF
    IF (SCALAR_INT_BUF(205) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 204)
      CALL c_f_pointer(curr_ptr, QYmax, shape=[fstarpu_matrix_get_nx(buffers, 204), fstarpu_matrix_get_ny(buffers, 204)])
    ELSE
      NULLIFY(QYmax)
    END IF
    IF (SCALAR_INT_BUF(206) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 205)
      CALL c_f_pointer(curr_ptr, iotamin, shape=[fstarpu_matrix_get_nx(buffers, 205), fstarpu_matrix_get_ny(buffers, 205)])
    ELSE
      NULLIFY(iotamin)
    END IF
    IF (SCALAR_INT_BUF(207) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 206)
      CALL c_f_pointer(curr_ptr, iotamax, shape=[fstarpu_matrix_get_nx(buffers, 206), fstarpu_matrix_get_ny(buffers, 206)])
    ELSE
      NULLIFY(iotamax)
    END IF
    IF (SCALAR_INT_BUF(208) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 207)
      CALL c_f_pointer(curr_ptr, iota2min, shape=[fstarpu_matrix_get_nx(buffers, 207), fstarpu_matrix_get_ny(buffers, 207)])
    ELSE
      NULLIFY(iota2min)
    END IF
    IF (SCALAR_INT_BUF(209) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 208)
      CALL c_f_pointer(curr_ptr, iota2max, shape=[fstarpu_matrix_get_nx(buffers, 208), fstarpu_matrix_get_ny(buffers, 208)])
    ELSE
      NULLIFY(iota2max)
    END IF
    IF (SCALAR_INT_BUF(210) == 1) THEN
      curr_ptr = fstarpu_block_get_ptr(buffers, 209)
      CALL c_f_pointer(curr_ptr, ZEtaylor, shape=[fstarpu_block_get_nx(buffers, 209), fstarpu_block_get_ny(buffers, 209), fstarpu_block_get_nz(buffers, 209)])
    ELSE
      NULLIFY(ZEtaylor)
    END IF
    IF (SCALAR_INT_BUF(211) == 1) THEN
      curr_ptr = fstarpu_block_get_ptr(buffers, 210)
      CALL c_f_pointer(curr_ptr, QXtaylor, shape=[fstarpu_block_get_nx(buffers, 210), fstarpu_block_get_ny(buffers, 210), fstarpu_block_get_nz(buffers, 210)])
    ELSE
      NULLIFY(QXtaylor)
    END IF
    IF (SCALAR_INT_BUF(212) == 1) THEN
      curr_ptr = fstarpu_block_get_ptr(buffers, 211)
      CALL c_f_pointer(curr_ptr, QYtaylor, shape=[fstarpu_block_get_nx(buffers, 211), fstarpu_block_get_ny(buffers, 211), fstarpu_block_get_nz(buffers, 211)])
    ELSE
      NULLIFY(QYtaylor)
    END IF
    IF (SCALAR_INT_BUF(213) == 1) THEN
      curr_ptr = fstarpu_block_get_ptr(buffers, 212)
      CALL c_f_pointer(curr_ptr, iotataylor, shape=[fstarpu_block_get_nx(buffers, 212), fstarpu_block_get_ny(buffers, 212), fstarpu_block_get_nz(buffers, 212)])
    ELSE
      NULLIFY(iotataylor)
    END IF
    IF (SCALAR_INT_BUF(214) == 1) THEN
      curr_ptr = fstarpu_block_get_ptr(buffers, 213)
      CALL c_f_pointer(curr_ptr, iota2taylor, shape=[fstarpu_block_get_nx(buffers, 213), fstarpu_block_get_ny(buffers, 213), fstarpu_block_get_nz(buffers, 213)])
    ELSE
      NULLIFY(iota2taylor)
    END IF
    IF (SCALAR_INT_BUF(215) == 1) THEN
      curr_ptr = fstarpu_block_get_ptr(buffers, 214)
      CALL c_f_pointer(curr_ptr, ZEtaylorvert, shape=[fstarpu_block_get_nx(buffers, 214), fstarpu_block_get_ny(buffers, 214), fstarpu_block_get_nz(buffers, 214)])
    ELSE
      NULLIFY(ZEtaylorvert)
    END IF
    IF (SCALAR_INT_BUF(216) == 1) THEN
      curr_ptr = fstarpu_block_get_ptr(buffers, 215)
      CALL c_f_pointer(curr_ptr, QXtaylorvert, shape=[fstarpu_block_get_nx(buffers, 215), fstarpu_block_get_ny(buffers, 215), fstarpu_block_get_nz(buffers, 215)])
    ELSE
      NULLIFY(QXtaylorvert)
    END IF
    IF (SCALAR_INT_BUF(217) == 1) THEN
      curr_ptr = fstarpu_block_get_ptr(buffers, 216)
      CALL c_f_pointer(curr_ptr, QYtaylorvert, shape=[fstarpu_block_get_nx(buffers, 216), fstarpu_block_get_ny(buffers, 216), fstarpu_block_get_nz(buffers, 216)])
    ELSE
      NULLIFY(QYtaylorvert)
    END IF
    IF (SCALAR_INT_BUF(218) == 1) THEN
      curr_ptr = fstarpu_block_get_ptr(buffers, 217)
      CALL c_f_pointer(curr_ptr, iotataylorvert, shape=[fstarpu_block_get_nx(buffers, 217), fstarpu_block_get_ny(buffers, 217), fstarpu_block_get_nz(buffers, 217)])
    ELSE
      NULLIFY(iotataylorvert)
    END IF
    IF (SCALAR_INT_BUF(219) == 1) THEN
      curr_ptr = fstarpu_block_get_ptr(buffers, 218)
      CALL c_f_pointer(curr_ptr, iota2taylorvert, shape=[fstarpu_block_get_nx(buffers, 218), fstarpu_block_get_ny(buffers, 218), fstarpu_block_get_nz(buffers, 218)])
    ELSE
      NULLIFY(iota2taylorvert)
    END IF
    IF (SCALAR_INT_BUF(220) == 1) THEN
      curr_ptr = fstarpu_block_get_ptr(buffers, 219)
      CALL c_f_pointer(curr_ptr, alphaZE0, shape=[fstarpu_block_get_nx(buffers, 219), fstarpu_block_get_ny(buffers, 219), fstarpu_block_get_nz(buffers, 219)])
    ELSE
      NULLIFY(alphaZE0)
    END IF
    IF (SCALAR_INT_BUF(221) == 1) THEN
      curr_ptr = fstarpu_block_get_ptr(buffers, 220)
      CALL c_f_pointer(curr_ptr, alphaQX0, shape=[fstarpu_block_get_nx(buffers, 220), fstarpu_block_get_ny(buffers, 220), fstarpu_block_get_nz(buffers, 220)])
    ELSE
      NULLIFY(alphaQX0)
    END IF
    IF (SCALAR_INT_BUF(222) == 1) THEN
      curr_ptr = fstarpu_block_get_ptr(buffers, 221)
      CALL c_f_pointer(curr_ptr, alphaQY0, shape=[fstarpu_block_get_nx(buffers, 221), fstarpu_block_get_ny(buffers, 221), fstarpu_block_get_nz(buffers, 221)])
    ELSE
      NULLIFY(alphaQY0)
    END IF
    IF (SCALAR_INT_BUF(223) == 1) THEN
      curr_ptr = fstarpu_block_get_ptr(buffers, 222)
      CALL c_f_pointer(curr_ptr, alphaiota0, shape=[fstarpu_block_get_nx(buffers, 222), fstarpu_block_get_ny(buffers, 222), fstarpu_block_get_nz(buffers, 222)])
    ELSE
      NULLIFY(alphaiota0)
    END IF
    IF (SCALAR_INT_BUF(224) == 1) THEN
      curr_ptr = fstarpu_block_get_ptr(buffers, 223)
      CALL c_f_pointer(curr_ptr, alphaiota20, shape=[fstarpu_block_get_nx(buffers, 223), fstarpu_block_get_ny(buffers, 223), fstarpu_block_get_nz(buffers, 223)])
    ELSE
      NULLIFY(alphaiota20)
    END IF
    IF (SCALAR_INT_BUF(225) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 224)
      CALL c_f_pointer(curr_ptr, alphaZE, shape=[fstarpu_matrix_get_nx(buffers, 224), fstarpu_matrix_get_ny(buffers, 224)])
    ELSE
      NULLIFY(alphaZE)
    END IF
    IF (SCALAR_INT_BUF(226) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 225)
      CALL c_f_pointer(curr_ptr, alphaQX, shape=[fstarpu_matrix_get_nx(buffers, 225), fstarpu_matrix_get_ny(buffers, 225)])
    ELSE
      NULLIFY(alphaQX)
    END IF
    IF (SCALAR_INT_BUF(227) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 226)
      CALL c_f_pointer(curr_ptr, alphaQY, shape=[fstarpu_matrix_get_nx(buffers, 226), fstarpu_matrix_get_ny(buffers, 226)])
    ELSE
      NULLIFY(alphaQY)
    END IF
    IF (SCALAR_INT_BUF(228) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 227)
      CALL c_f_pointer(curr_ptr, alphaiota, shape=[fstarpu_matrix_get_nx(buffers, 227), fstarpu_matrix_get_ny(buffers, 227)])
    ELSE
      NULLIFY(alphaiota)
    END IF
    IF (SCALAR_INT_BUF(229) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 228)
      CALL c_f_pointer(curr_ptr, alphaiota2, shape=[fstarpu_matrix_get_nx(buffers, 228), fstarpu_matrix_get_ny(buffers, 228)])
    ELSE
      NULLIFY(alphaiota2)
    END IF
    IF (SCALAR_INT_BUF(230) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 229)
      CALL c_f_pointer(curr_ptr, alphaZEm, shape=[fstarpu_matrix_get_nx(buffers, 229), fstarpu_matrix_get_ny(buffers, 229)])
    ELSE
      NULLIFY(alphaZEm)
    END IF
    IF (SCALAR_INT_BUF(231) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 230)
      CALL c_f_pointer(curr_ptr, alphaQXm, shape=[fstarpu_matrix_get_nx(buffers, 230), fstarpu_matrix_get_ny(buffers, 230)])
    ELSE
      NULLIFY(alphaQXm)
    END IF
    IF (SCALAR_INT_BUF(232) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 231)
      CALL c_f_pointer(curr_ptr, alphaQYm, shape=[fstarpu_matrix_get_nx(buffers, 231), fstarpu_matrix_get_ny(buffers, 231)])
    ELSE
      NULLIFY(alphaQYm)
    END IF
    IF (SCALAR_INT_BUF(233) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 232)
      CALL c_f_pointer(curr_ptr, alphaiotam, shape=[fstarpu_matrix_get_nx(buffers, 232), fstarpu_matrix_get_ny(buffers, 232)])
    ELSE
      NULLIFY(alphaiotam)
    END IF
    IF (SCALAR_INT_BUF(234) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 233)
      CALL c_f_pointer(curr_ptr, alphaiota2m, shape=[fstarpu_matrix_get_nx(buffers, 233), fstarpu_matrix_get_ny(buffers, 233)])
    ELSE
      NULLIFY(alphaiota2m)
    END IF
    IF (SCALAR_INT_BUF(235) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 234)
      CALL c_f_pointer(curr_ptr, alphaZE_max, shape=[fstarpu_matrix_get_nx(buffers, 234), fstarpu_matrix_get_ny(buffers, 234)])
    ELSE
      NULLIFY(alphaZE_max)
    END IF
    IF (SCALAR_INT_BUF(236) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 235)
      CALL c_f_pointer(curr_ptr, alphaQX_max, shape=[fstarpu_matrix_get_nx(buffers, 235), fstarpu_matrix_get_ny(buffers, 235)])
    ELSE
      NULLIFY(alphaQX_max)
    END IF
    IF (SCALAR_INT_BUF(237) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 236)
      CALL c_f_pointer(curr_ptr, alphaQY_max, shape=[fstarpu_matrix_get_nx(buffers, 236), fstarpu_matrix_get_ny(buffers, 236)])
    ELSE
      NULLIFY(alphaQY_max)
    END IF
    IF (SCALAR_INT_BUF(238) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 237)
      CALL c_f_pointer(curr_ptr, alphaiota_max, shape=[fstarpu_matrix_get_nx(buffers, 237), fstarpu_matrix_get_ny(buffers, 237)])
    ELSE
      NULLIFY(alphaiota_max)
    END IF
    IF (SCALAR_INT_BUF(239) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 238)
      CALL c_f_pointer(curr_ptr, alphaiota2_max, shape=[fstarpu_matrix_get_nx(buffers, 238), fstarpu_matrix_get_ny(buffers, 238)])
    ELSE
      NULLIFY(alphaiota2_max)
    END IF
    IF (SCALAR_INT_BUF(240) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 239)
      CALL c_f_pointer(curr_ptr, limitZE, shape=[fstarpu_matrix_get_nx(buffers, 239), fstarpu_matrix_get_ny(buffers, 239)])
    ELSE
      NULLIFY(limitZE)
    END IF
    IF (SCALAR_INT_BUF(241) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 240)
      CALL c_f_pointer(curr_ptr, limitQX, shape=[fstarpu_matrix_get_nx(buffers, 240), fstarpu_matrix_get_ny(buffers, 240)])
    ELSE
      NULLIFY(limitQX)
    END IF
    IF (SCALAR_INT_BUF(242) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 241)
      CALL c_f_pointer(curr_ptr, limitQY, shape=[fstarpu_matrix_get_nx(buffers, 241), fstarpu_matrix_get_ny(buffers, 241)])
    ELSE
      NULLIFY(limitQY)
    END IF
    IF (SCALAR_INT_BUF(243) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 242)
      CALL c_f_pointer(curr_ptr, limitiota, shape=[fstarpu_matrix_get_nx(buffers, 242), fstarpu_matrix_get_ny(buffers, 242)])
    ELSE
      NULLIFY(limitiota)
    END IF
    IF (SCALAR_INT_BUF(244) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 243)
      CALL c_f_pointer(curr_ptr, limitiota2, shape=[fstarpu_matrix_get_nx(buffers, 243), fstarpu_matrix_get_ny(buffers, 243)])
    ELSE
      NULLIFY(limitiota2)
    END IF
    IF (SCALAR_INT_BUF(245) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 244)
      CALL c_f_pointer(curr_ptr, ZEconst, shape=[fstarpu_matrix_get_nx(buffers, 244), fstarpu_matrix_get_ny(buffers, 244)])
    ELSE
      NULLIFY(ZEconst)
    END IF
    IF (SCALAR_INT_BUF(246) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 245)
      CALL c_f_pointer(curr_ptr, QXconst, shape=[fstarpu_matrix_get_nx(buffers, 245), fstarpu_matrix_get_ny(buffers, 245)])
    ELSE
      NULLIFY(QXconst)
    END IF
    IF (SCALAR_INT_BUF(247) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 246)
      CALL c_f_pointer(curr_ptr, QYconst, shape=[fstarpu_matrix_get_nx(buffers, 246), fstarpu_matrix_get_ny(buffers, 246)])
    ELSE
      NULLIFY(QYconst)
    END IF
    IF (SCALAR_INT_BUF(248) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 247)
      CALL c_f_pointer(curr_ptr, iotaconst, shape=[fstarpu_matrix_get_nx(buffers, 247), fstarpu_matrix_get_ny(buffers, 247)])
    ELSE
      NULLIFY(iotaconst)
    END IF
    IF (SCALAR_INT_BUF(249) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 248)
      CALL c_f_pointer(curr_ptr, iota2const, shape=[fstarpu_matrix_get_nx(buffers, 248), fstarpu_matrix_get_ny(buffers, 248)])
    ELSE
      NULLIFY(iota2const)
    END IF
    IF (SCALAR_INT_BUF(250) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 249)
      CALL c_f_pointer(curr_ptr, NODES_LG, shape=[fstarpu_vector_get_nx(buffers, 249)])
    ELSE
      NULLIFY(NODES_LG)
    END IF
    IF (SCALAR_INT_BUF(251) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 250)
      CALL c_f_pointer(curr_ptr, ANGTAB, shape=[fstarpu_matrix_get_nx(buffers, 250), fstarpu_matrix_get_ny(buffers, 250)])
    ELSE
      NULLIFY(ANGTAB)
    END IF
    IF (SCALAR_INT_BUF(252) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 251)
      CALL c_f_pointer(curr_ptr, CENTAB, shape=[fstarpu_matrix_get_nx(buffers, 251), fstarpu_matrix_get_ny(buffers, 251)])
    ELSE
      NULLIFY(CENTAB)
    END IF
    IF (SCALAR_INT_BUF(253) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 252)
      CALL c_f_pointer(curr_ptr, ELETAB, shape=[fstarpu_matrix_get_nx(buffers, 252), fstarpu_matrix_get_ny(buffers, 252)])
    ELSE
      NULLIFY(ELETAB)
    END IF
    IF (SCALAR_INT_BUF(254) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 253)
      CALL c_f_pointer(curr_ptr, DG_ANG, shape=[fstarpu_vector_get_nx(buffers, 253)])
    ELSE
      NULLIFY(DG_ANG)
    END IF
    IF (SCALAR_INT_BUF(255) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 254)
      CALL c_f_pointer(curr_ptr, DP_DG, shape=[fstarpu_vector_get_nx(buffers, 254)])
    ELSE
      NULLIFY(DP_DG)
    END IF
    IF (SCALAR_INT_BUF(256) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 255)
      CALL c_f_pointer(curr_ptr, EL_COUNT, shape=[fstarpu_vector_get_nx(buffers, 255)])
    ELSE
      NULLIFY(EL_COUNT)
    END IF
    IF (SCALAR_INT_BUF(257) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 256)
      CALL c_f_pointer(curr_ptr, NNOEL, shape=[fstarpu_matrix_get_nx(buffers, 256), fstarpu_matrix_get_ny(buffers, 256)])
    ELSE
      NULLIFY(NNOEL)
    END IF
    IF (SCALAR_INT_BUF(258) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 257)
      CALL c_f_pointer(curr_ptr, NNDEL, shape=[fstarpu_vector_get_nx(buffers, 257)])
    ELSE
      NULLIFY(NNDEL)
    END IF
    IF (SCALAR_INT_BUF(259) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 258)
      CALL c_f_pointer(curr_ptr, NDEL, shape=[fstarpu_matrix_get_nx(buffers, 258), fstarpu_matrix_get_ny(buffers, 258)])
    ELSE
      NULLIFY(NDEL)
    END IF
    IF (SCALAR_INT_BUF(260) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 259)
      CALL c_f_pointer(curr_ptr, FX_MID, shape=[fstarpu_matrix_get_nx(buffers, 259), fstarpu_matrix_get_ny(buffers, 259)])
    ELSE
      NULLIFY(FX_MID)
    END IF
    IF (SCALAR_INT_BUF(261) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 260)
      CALL c_f_pointer(curr_ptr, GX_MID, shape=[fstarpu_matrix_get_nx(buffers, 260), fstarpu_matrix_get_ny(buffers, 260)])
    ELSE
      NULLIFY(GX_MID)
    END IF
    IF (SCALAR_INT_BUF(262) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 261)
      CALL c_f_pointer(curr_ptr, HX_MID, shape=[fstarpu_matrix_get_nx(buffers, 261), fstarpu_matrix_get_ny(buffers, 261)])
    ELSE
      NULLIFY(HX_MID)
    END IF
    IF (SCALAR_INT_BUF(263) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 262)
      CALL c_f_pointer(curr_ptr, FY_MID, shape=[fstarpu_matrix_get_nx(buffers, 262), fstarpu_matrix_get_ny(buffers, 262)])
    ELSE
      NULLIFY(FY_MID)
    END IF
    IF (SCALAR_INT_BUF(264) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 263)
      CALL c_f_pointer(curr_ptr, GY_MID, shape=[fstarpu_matrix_get_nx(buffers, 263), fstarpu_matrix_get_ny(buffers, 263)])
    ELSE
      NULLIFY(GY_MID)
    END IF
    IF (SCALAR_INT_BUF(265) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 264)
      CALL c_f_pointer(curr_ptr, HY_MID, shape=[fstarpu_matrix_get_nx(buffers, 264), fstarpu_matrix_get_ny(buffers, 264)])
    ELSE
      NULLIFY(HY_MID)
    END IF
    IF (SCALAR_INT_BUF(266) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 265)
      CALL c_f_pointer(curr_ptr, ZE_C, shape=[fstarpu_vector_get_nx(buffers, 265)])
    ELSE
      NULLIFY(ZE_C)
    END IF
    IF (SCALAR_INT_BUF(267) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 266)
      CALL c_f_pointer(curr_ptr, QX_C, shape=[fstarpu_vector_get_nx(buffers, 266)])
    ELSE
      NULLIFY(QX_C)
    END IF
    IF (SCALAR_INT_BUF(268) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 267)
      CALL c_f_pointer(curr_ptr, QY_C, shape=[fstarpu_vector_get_nx(buffers, 267)])
    ELSE
      NULLIFY(QY_C)
    END IF
    IF (SCALAR_INT_BUF(269) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 268)
      CALL c_f_pointer(curr_ptr, dynP_DG, shape=[fstarpu_vector_get_nx(buffers, 268)])
    ELSE
      NULLIFY(dynP_DG)
    END IF
    IF (SCALAR_INT_BUF(270) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 269)
      CALL c_f_pointer(curr_ptr, iota2_DG, shape=[fstarpu_vector_get_nx(buffers, 269)])
    ELSE
      NULLIFY(iota2_DG)
    END IF
    IF (SCALAR_INT_BUF(271) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 270)
      CALL c_f_pointer(curr_ptr, iota_DG, shape=[fstarpu_vector_get_nx(buffers, 270)])
    ELSE
      NULLIFY(iota_DG)
    END IF
    IF (SCALAR_INT_BUF(272) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 271)
      CALL c_f_pointer(curr_ptr, iotaa_DG, shape=[fstarpu_vector_get_nx(buffers, 271)])
    ELSE
      NULLIFY(iotaa_DG)
    END IF
    IF (SCALAR_INT_BUF(273) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 272)
      CALL c_f_pointer(curr_ptr, bed_DG, shape=[fstarpu_matrix_get_nx(buffers, 272), fstarpu_matrix_get_ny(buffers, 272)])
    ELSE
      NULLIFY(bed_DG)
    END IF
    IF (SCALAR_INT_BUF(274) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 273)
      CALL c_f_pointer(curr_ptr, bed_N_int, shape=[fstarpu_vector_get_nx(buffers, 273)])
    ELSE
      NULLIFY(bed_N_int)
    END IF
    IF (SCALAR_INT_BUF(275) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 274)
      CALL c_f_pointer(curr_ptr, bed_N_ext, shape=[fstarpu_vector_get_nx(buffers, 274)])
    ELSE
      NULLIFY(bed_N_ext)
    END IF
    IF (SCALAR_INT_BUF(276) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 275)
      CALL c_f_pointer(curr_ptr, IDUMY, shape=[fstarpu_vector_get_nx(buffers, 275)])
    ELSE
      NULLIFY(IDUMY)
    END IF
    IF (SCALAR_INT_BUF(277) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 276)
      CALL c_f_pointer(curr_ptr, DUMY1, shape=[fstarpu_vector_get_nx(buffers, 276)])
    ELSE
      NULLIFY(DUMY1)
    END IF
    IF (SCALAR_INT_BUF(278) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 277)
      CALL c_f_pointer(curr_ptr, DUMY2, shape=[fstarpu_vector_get_nx(buffers, 277)])
    ELSE
      NULLIFY(DUMY2)
    END IF
    IF (SCALAR_INT_BUF(279) == 1) THEN
      curr_ptr = fstarpu_block_get_ptr(buffers, 278)
      CALL c_f_pointer(curr_ptr, DGDUMY1, shape=[fstarpu_block_get_nx(buffers, 278), fstarpu_block_get_ny(buffers, 278), fstarpu_block_get_nz(buffers, 278)])
    ELSE
      NULLIFY(DGDUMY1)
    END IF
    IF (SCALAR_INT_BUF(280) == 1) THEN
      curr_ptr = fstarpu_block_get_ptr(buffers, 279)
      CALL c_f_pointer(curr_ptr, DGDUMY2, shape=[fstarpu_block_get_nx(buffers, 279), fstarpu_block_get_ny(buffers, 279), fstarpu_block_get_nz(buffers, 279)])
    ELSE
      NULLIFY(DGDUMY2)
    END IF
    IF (SCALAR_INT_BUF(281) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 280)
      CALL c_f_pointer(curr_ptr, pdg_el, shape=[fstarpu_vector_get_nx(buffers, 280)])
    ELSE
      NULLIFY(pdg_el)
    END IF
    IF (SCALAR_INT_BUF(282) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 281)
      CALL c_f_pointer(curr_ptr, ETAS, shape=[fstarpu_vector_get_nx(buffers, 281)])
    ELSE
      NULLIFY(ETAS)
    END IF
    IF (SCALAR_INT_BUF(283) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 282)
      CALL c_f_pointer(curr_ptr, ETA1, shape=[fstarpu_vector_get_nx(buffers, 282)])
    ELSE
      NULLIFY(ETA1)
    END IF
    IF (SCALAR_INT_BUF(284) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 283)
      CALL c_f_pointer(curr_ptr, ETA2, shape=[fstarpu_vector_get_nx(buffers, 283)])
    ELSE
      NULLIFY(ETA2)
    END IF
    IF (SCALAR_INT_BUF(285) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 284)
      CALL c_f_pointer(curr_ptr, ETAMAX, shape=[fstarpu_vector_get_nx(buffers, 284)])
    ELSE
      NULLIFY(ETAMAX)
    END IF
    IF (SCALAR_INT_BUF(286) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 285)
      CALL c_f_pointer(curr_ptr, entrop, shape=[fstarpu_matrix_get_nx(buffers, 285), fstarpu_matrix_get_ny(buffers, 285)])
    ELSE
      NULLIFY(entrop)
    END IF
    IF (SCALAR_INT_BUF(287) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 286)
      CALL c_f_pointer(curr_ptr, tracer, shape=[fstarpu_vector_get_nx(buffers, 286)])
    ELSE
      NULLIFY(tracer)
    END IF
    IF (SCALAR_INT_BUF(288) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 287)
      CALL c_f_pointer(curr_ptr, tracer2, shape=[fstarpu_vector_get_nx(buffers, 287)])
    ELSE
      NULLIFY(tracer2)
    END IF
    IF (SCALAR_INT_BUF(289) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 288)
      CALL c_f_pointer(curr_ptr, MassMax, shape=[fstarpu_vector_get_nx(buffers, 288)])
    ELSE
      NULLIFY(MassMax)
    END IF
    IF (SCALAR_INT_BUF(290) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 289)
      CALL c_f_pointer(curr_ptr, bed_int, shape=[fstarpu_matrix_get_nx(buffers, 289), fstarpu_matrix_get_ny(buffers, 289)])
    ELSE
      NULLIFY(bed_int)
    END IF
    IF (SCALAR_INT_BUF(291) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 290)
      CALL c_f_pointer(curr_ptr, UU1, shape=[fstarpu_vector_get_nx(buffers, 290)])
    ELSE
      NULLIFY(UU1)
    END IF
    IF (SCALAR_INT_BUF(292) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 291)
      CALL c_f_pointer(curr_ptr, UU2, shape=[fstarpu_vector_get_nx(buffers, 291)])
    ELSE
      NULLIFY(UU2)
    END IF
    IF (SCALAR_INT_BUF(293) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 292)
      CALL c_f_pointer(curr_ptr, VV1, shape=[fstarpu_vector_get_nx(buffers, 292)])
    ELSE
      NULLIFY(VV1)
    END IF
    IF (SCALAR_INT_BUF(294) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 293)
      CALL c_f_pointer(curr_ptr, VV2, shape=[fstarpu_vector_get_nx(buffers, 293)])
    ELSE
      NULLIFY(VV2)
    END IF
    IF (SCALAR_INT_BUF(295) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 294)
      CALL c_f_pointer(curr_ptr, dyn_P, shape=[fstarpu_vector_get_nx(buffers, 294)])
    ELSE
      NULLIFY(dyn_P)
    END IF
    IF (SCALAR_INT_BUF(296) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 295)
      CALL c_f_pointer(curr_ptr, DP, shape=[fstarpu_vector_get_nx(buffers, 295)])
    ELSE
      NULLIFY(DP)
    END IF
    IF (SCALAR_INT_BUF(297) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 296)
      CALL c_f_pointer(curr_ptr, DP0, shape=[fstarpu_vector_get_nx(buffers, 296)])
    ELSE
      NULLIFY(DP0)
    END IF
    IF (SCALAR_INT_BUF(298) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 297)
      CALL c_f_pointer(curr_ptr, DPe, shape=[fstarpu_vector_get_nx(buffers, 297)])
    ELSE
      NULLIFY(DPe)
    END IF
    IF (SCALAR_INT_BUF(299) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 298)
      CALL c_f_pointer(curr_ptr, SFAC, shape=[fstarpu_vector_get_nx(buffers, 298)])
    ELSE
      NULLIFY(SFAC)
    END IF
    IF (SCALAR_INT_BUF(300) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 299)
      CALL c_f_pointer(curr_ptr, QU, shape=[fstarpu_vector_get_nx(buffers, 299)])
    ELSE
      NULLIFY(QU)
    END IF
    IF (SCALAR_INT_BUF(301) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 300)
      CALL c_f_pointer(curr_ptr, QV, shape=[fstarpu_vector_get_nx(buffers, 300)])
    ELSE
      NULLIFY(QV)
    END IF
    IF (SCALAR_INT_BUF(302) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 301)
      CALL c_f_pointer(curr_ptr, QW, shape=[fstarpu_vector_get_nx(buffers, 301)])
    ELSE
      NULLIFY(QW)
    END IF
    IF (SCALAR_INT_BUF(303) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 302)
      CALL c_f_pointer(curr_ptr, CORIF, shape=[fstarpu_vector_get_nx(buffers, 302)])
    ELSE
      NULLIFY(CORIF)
    END IF
    IF (SCALAR_INT_BUF(304) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 303)
      CALL c_f_pointer(curr_ptr, TPK, shape=[fstarpu_vector_get_nx(buffers, 303)])
    ELSE
      NULLIFY(TPK)
    END IF
    IF (SCALAR_INT_BUF(305) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 304)
      CALL c_f_pointer(curr_ptr, FFT, shape=[fstarpu_vector_get_nx(buffers, 304)])
    ELSE
      NULLIFY(FFT)
    END IF
    IF (SCALAR_INT_BUF(306) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 305)
      CALL c_f_pointer(curr_ptr, FACET, shape=[fstarpu_vector_get_nx(buffers, 305)])
    ELSE
      NULLIFY(FACET)
    END IF
    IF (SCALAR_INT_BUF(307) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 306)
      CALL c_f_pointer(curr_ptr, ETRF, shape=[fstarpu_vector_get_nx(buffers, 306)])
    ELSE
      NULLIFY(ETRF)
    END IF
    IF (SCALAR_INT_BUF(308) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 307)
      CALL c_f_pointer(curr_ptr, ESBIN1, shape=[fstarpu_vector_get_nx(buffers, 307)])
    ELSE
      NULLIFY(ESBIN1)
    END IF
    IF (SCALAR_INT_BUF(309) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 308)
      CALL c_f_pointer(curr_ptr, ESBIN2, shape=[fstarpu_vector_get_nx(buffers, 308)])
    ELSE
      NULLIFY(ESBIN2)
    END IF
    IF (SCALAR_INT_BUF(310) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 309)
      CALL c_f_pointer(curr_ptr, QTEMA, shape=[fstarpu_matrix_get_nx(buffers, 309), fstarpu_matrix_get_ny(buffers, 309)])
    ELSE
      NULLIFY(QTEMA)
    END IF
    IF (SCALAR_INT_BUF(311) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 310)
      CALL c_f_pointer(curr_ptr, QTEMB, shape=[fstarpu_matrix_get_nx(buffers, 310), fstarpu_matrix_get_ny(buffers, 310)])
    ELSE
      NULLIFY(QTEMB)
    END IF
    IF (SCALAR_INT_BUF(312) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 311)
      CALL c_f_pointer(curr_ptr, QN2, shape=[fstarpu_vector_get_nx(buffers, 311)])
    ELSE
      NULLIFY(QN2)
    END IF
    IF (SCALAR_INT_BUF(313) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 312)
      CALL c_f_pointer(curr_ptr, BNDLEN2O3, shape=[fstarpu_vector_get_nx(buffers, 312)])
    ELSE
      NULLIFY(BNDLEN2O3)
    END IF
    IF (SCALAR_INT_BUF(314) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 313)
      CALL c_f_pointer(curr_ptr, CSII, shape=[fstarpu_vector_get_nx(buffers, 313)])
    ELSE
      NULLIFY(CSII)
    END IF
    IF (SCALAR_INT_BUF(315) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 314)
      CALL c_f_pointer(curr_ptr, SIII, shape=[fstarpu_vector_get_nx(buffers, 314)])
    ELSE
      NULLIFY(SIII)
    END IF
    IF (SCALAR_INT_BUF(316) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 315)
      CALL c_f_pointer(curr_ptr, QNAM, shape=[fstarpu_matrix_get_nx(buffers, 315), fstarpu_matrix_get_ny(buffers, 315)])
    ELSE
      NULLIFY(QNAM)
    END IF
    IF (SCALAR_INT_BUF(317) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 316)
      CALL c_f_pointer(curr_ptr, QNPH, shape=[fstarpu_matrix_get_nx(buffers, 316), fstarpu_matrix_get_ny(buffers, 316)])
    ELSE
      NULLIFY(QNPH)
    END IF
    IF (SCALAR_INT_BUF(318) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 317)
      CALL c_f_pointer(curr_ptr, QNIN1, shape=[fstarpu_vector_get_nx(buffers, 317)])
    ELSE
      NULLIFY(QNIN1)
    END IF
    IF (SCALAR_INT_BUF(319) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 318)
      CALL c_f_pointer(curr_ptr, QNIN2, shape=[fstarpu_vector_get_nx(buffers, 318)])
    ELSE
      NULLIFY(QNIN2)
    END IF
    IF (SCALAR_INT_BUF(320) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 319)
      CALL c_f_pointer(curr_ptr, CSI, shape=[fstarpu_vector_get_nx(buffers, 319)])
    ELSE
      NULLIFY(CSI)
    END IF
    IF (SCALAR_INT_BUF(321) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 320)
      CALL c_f_pointer(curr_ptr, SII, shape=[fstarpu_vector_get_nx(buffers, 320)])
    ELSE
      NULLIFY(SII)
    END IF
    IF (SCALAR_INT_BUF(322) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 321)
      CALL c_f_pointer(curr_ptr, ET00, shape=[fstarpu_vector_get_nx(buffers, 321)])
    ELSE
      NULLIFY(ET00)
    END IF
    IF (SCALAR_INT_BUF(323) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 322)
      CALL c_f_pointer(curr_ptr, BT00, shape=[fstarpu_vector_get_nx(buffers, 322)])
    ELSE
      NULLIFY(BT00)
    END IF
    IF (SCALAR_INT_BUF(324) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 323)
      CALL c_f_pointer(curr_ptr, STAIE1, shape=[fstarpu_vector_get_nx(buffers, 323)])
    ELSE
      NULLIFY(STAIE1)
    END IF
    IF (SCALAR_INT_BUF(325) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 324)
      CALL c_f_pointer(curr_ptr, STAIE2, shape=[fstarpu_vector_get_nx(buffers, 324)])
    ELSE
      NULLIFY(STAIE2)
    END IF
    IF (SCALAR_INT_BUF(326) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 325)
      CALL c_f_pointer(curr_ptr, STAIE3, shape=[fstarpu_vector_get_nx(buffers, 325)])
    ELSE
      NULLIFY(STAIE3)
    END IF
    IF (SCALAR_INT_BUF(327) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 326)
      CALL c_f_pointer(curr_ptr, XEV, shape=[fstarpu_vector_get_nx(buffers, 326)])
    ELSE
      NULLIFY(XEV)
    END IF
    IF (SCALAR_INT_BUF(328) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 327)
      CALL c_f_pointer(curr_ptr, YEV, shape=[fstarpu_vector_get_nx(buffers, 327)])
    ELSE
      NULLIFY(YEV)
    END IF
    IF (SCALAR_INT_BUF(329) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 328)
      CALL c_f_pointer(curr_ptr, SLEV, shape=[fstarpu_vector_get_nx(buffers, 328)])
    ELSE
      NULLIFY(SLEV)
    END IF
    IF (SCALAR_INT_BUF(330) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 329)
      CALL c_f_pointer(curr_ptr, SFEV, shape=[fstarpu_vector_get_nx(buffers, 329)])
    ELSE
      NULLIFY(SFEV)
    END IF
    IF (SCALAR_INT_BUF(331) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 330)
      CALL c_f_pointer(curr_ptr, UU00, shape=[fstarpu_vector_get_nx(buffers, 330)])
    ELSE
      NULLIFY(UU00)
    END IF
    IF (SCALAR_INT_BUF(332) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 331)
      CALL c_f_pointer(curr_ptr, VV00, shape=[fstarpu_vector_get_nx(buffers, 331)])
    ELSE
      NULLIFY(VV00)
    END IF
    IF (SCALAR_INT_BUF(333) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 332)
      CALL c_f_pointer(curr_ptr, STAIV1, shape=[fstarpu_vector_get_nx(buffers, 332)])
    ELSE
      NULLIFY(STAIV1)
    END IF
    IF (SCALAR_INT_BUF(334) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 333)
      CALL c_f_pointer(curr_ptr, STAIV2, shape=[fstarpu_vector_get_nx(buffers, 333)])
    ELSE
      NULLIFY(STAIV2)
    END IF
    IF (SCALAR_INT_BUF(335) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 334)
      CALL c_f_pointer(curr_ptr, STAIV3, shape=[fstarpu_vector_get_nx(buffers, 334)])
    ELSE
      NULLIFY(STAIV3)
    END IF
    IF (SCALAR_INT_BUF(336) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 335)
      CALL c_f_pointer(curr_ptr, XEC, shape=[fstarpu_vector_get_nx(buffers, 335)])
    ELSE
      NULLIFY(XEC)
    END IF
    IF (SCALAR_INT_BUF(337) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 336)
      CALL c_f_pointer(curr_ptr, YEC, shape=[fstarpu_vector_get_nx(buffers, 336)])
    ELSE
      NULLIFY(YEC)
    END IF
    IF (SCALAR_INT_BUF(338) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 337)
      CALL c_f_pointer(curr_ptr, SLEC, shape=[fstarpu_vector_get_nx(buffers, 337)])
    ELSE
      NULLIFY(SLEC)
    END IF
    IF (SCALAR_INT_BUF(339) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 338)
      CALL c_f_pointer(curr_ptr, SFEC, shape=[fstarpu_vector_get_nx(buffers, 338)])
    ELSE
      NULLIFY(SFEC)
    END IF
    IF (SCALAR_INT_BUF(340) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 339)
      CALL c_f_pointer(curr_ptr, CC00, shape=[fstarpu_vector_get_nx(buffers, 339)])
    ELSE
      NULLIFY(CC00)
    END IF
    IF (SCALAR_INT_BUF(341) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 340)
      CALL c_f_pointer(curr_ptr, STAIC1, shape=[fstarpu_vector_get_nx(buffers, 340)])
    ELSE
      NULLIFY(STAIC1)
    END IF
    IF (SCALAR_INT_BUF(342) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 341)
      CALL c_f_pointer(curr_ptr, STAIC2, shape=[fstarpu_vector_get_nx(buffers, 341)])
    ELSE
      NULLIFY(STAIC2)
    END IF
    IF (SCALAR_INT_BUF(343) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 342)
      CALL c_f_pointer(curr_ptr, STAIC3, shape=[fstarpu_vector_get_nx(buffers, 342)])
    ELSE
      NULLIFY(STAIC3)
    END IF
    IF (SCALAR_INT_BUF(344) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 343)
      CALL c_f_pointer(curr_ptr, XEM, shape=[fstarpu_vector_get_nx(buffers, 343)])
    ELSE
      NULLIFY(XEM)
    END IF
    IF (SCALAR_INT_BUF(345) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 344)
      CALL c_f_pointer(curr_ptr, YEM, shape=[fstarpu_vector_get_nx(buffers, 344)])
    ELSE
      NULLIFY(YEM)
    END IF
    IF (SCALAR_INT_BUF(346) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 345)
      CALL c_f_pointer(curr_ptr, SLEM, shape=[fstarpu_vector_get_nx(buffers, 345)])
    ELSE
      NULLIFY(SLEM)
    END IF
    IF (SCALAR_INT_BUF(347) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 346)
      CALL c_f_pointer(curr_ptr, SFEM, shape=[fstarpu_vector_get_nx(buffers, 346)])
    ELSE
      NULLIFY(SFEM)
    END IF
    IF (SCALAR_INT_BUF(348) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 347)
      CALL c_f_pointer(curr_ptr, RMU00, shape=[fstarpu_vector_get_nx(buffers, 347)])
    ELSE
      NULLIFY(RMU00)
    END IF
    IF (SCALAR_INT_BUF(349) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 348)
      CALL c_f_pointer(curr_ptr, RMV00, shape=[fstarpu_vector_get_nx(buffers, 348)])
    ELSE
      NULLIFY(RMV00)
    END IF
    IF (SCALAR_INT_BUF(350) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 349)
      CALL c_f_pointer(curr_ptr, RMP00, shape=[fstarpu_vector_get_nx(buffers, 349)])
    ELSE
      NULLIFY(RMP00)
    END IF
    IF (SCALAR_INT_BUF(351) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 350)
      CALL c_f_pointer(curr_ptr, STAIM1, shape=[fstarpu_vector_get_nx(buffers, 350)])
    ELSE
      NULLIFY(STAIM1)
    END IF
    IF (SCALAR_INT_BUF(352) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 351)
      CALL c_f_pointer(curr_ptr, STAIM2, shape=[fstarpu_vector_get_nx(buffers, 351)])
    ELSE
      NULLIFY(STAIM2)
    END IF
    IF (SCALAR_INT_BUF(353) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 352)
      CALL c_f_pointer(curr_ptr, STAIM3, shape=[fstarpu_vector_get_nx(buffers, 352)])
    ELSE
      NULLIFY(STAIM3)
    END IF
    IF (SCALAR_INT_BUF(354) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 353)
      CALL c_f_pointer(curr_ptr, CH1, shape=[fstarpu_vector_get_nx(buffers, 353)])
    ELSE
      NULLIFY(CH1)
    END IF
    IF (SCALAR_INT_BUF(355) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 354)
      CALL c_f_pointer(curr_ptr, QB, shape=[fstarpu_vector_get_nx(buffers, 354)])
    ELSE
      NULLIFY(QB)
    END IF
    IF (SCALAR_INT_BUF(356) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 355)
      CALL c_f_pointer(curr_ptr, QA, shape=[fstarpu_vector_get_nx(buffers, 355)])
    ELSE
      NULLIFY(QA)
    END IF
    IF (SCALAR_INT_BUF(357) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 356)
      CALL c_f_pointer(curr_ptr, SOURSIN, shape=[fstarpu_vector_get_nx(buffers, 356)])
    ELSE
      NULLIFY(SOURSIN)
    END IF
    IF (SCALAR_INT_BUF(358) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 357)
      CALL c_f_pointer(curr_ptr, WSX1, shape=[fstarpu_vector_get_nx(buffers, 357)])
    ELSE
      NULLIFY(WSX1)
    END IF
    IF (SCALAR_INT_BUF(359) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 358)
      CALL c_f_pointer(curr_ptr, WSY1, shape=[fstarpu_vector_get_nx(buffers, 358)])
    ELSE
      NULLIFY(WSY1)
    END IF
    IF (SCALAR_INT_BUF(360) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 359)
      CALL c_f_pointer(curr_ptr, PR1, shape=[fstarpu_vector_get_nx(buffers, 359)])
    ELSE
      NULLIFY(PR1)
    END IF
    IF (SCALAR_INT_BUF(361) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 360)
      CALL c_f_pointer(curr_ptr, WSX2, shape=[fstarpu_vector_get_nx(buffers, 360)])
    ELSE
      NULLIFY(WSX2)
    END IF
    IF (SCALAR_INT_BUF(362) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 361)
      CALL c_f_pointer(curr_ptr, WSY2, shape=[fstarpu_vector_get_nx(buffers, 361)])
    ELSE
      NULLIFY(WSY2)
    END IF
    IF (SCALAR_INT_BUF(363) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 362)
      CALL c_f_pointer(curr_ptr, PR2, shape=[fstarpu_vector_get_nx(buffers, 362)])
    ELSE
      NULLIFY(PR2)
    END IF
    IF (SCALAR_INT_BUF(364) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 363)
      CALL c_f_pointer(curr_ptr, WVNX1, shape=[fstarpu_vector_get_nx(buffers, 363)])
    ELSE
      NULLIFY(WVNX1)
    END IF
    IF (SCALAR_INT_BUF(365) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 364)
      CALL c_f_pointer(curr_ptr, WVNY1, shape=[fstarpu_vector_get_nx(buffers, 364)])
    ELSE
      NULLIFY(WVNY1)
    END IF
    IF (SCALAR_INT_BUF(366) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 365)
      CALL c_f_pointer(curr_ptr, PRN1, shape=[fstarpu_vector_get_nx(buffers, 365)])
    ELSE
      NULLIFY(PRN1)
    END IF
    IF (SCALAR_INT_BUF(367) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 366)
      CALL c_f_pointer(curr_ptr, WVNX2, shape=[fstarpu_vector_get_nx(buffers, 366)])
    ELSE
      NULLIFY(WVNX2)
    END IF
    IF (SCALAR_INT_BUF(368) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 367)
      CALL c_f_pointer(curr_ptr, WVNY2, shape=[fstarpu_vector_get_nx(buffers, 367)])
    ELSE
      NULLIFY(WVNY2)
    END IF
    IF (SCALAR_INT_BUF(369) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 368)
      CALL c_f_pointer(curr_ptr, PRN2, shape=[fstarpu_vector_get_nx(buffers, 368)])
    ELSE
      NULLIFY(PRN2)
    END IF
    IF (SCALAR_INT_BUF(370) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 369)
      CALL c_f_pointer(curr_ptr, RSNX1, shape=[fstarpu_vector_get_nx(buffers, 369)])
    ELSE
      NULLIFY(RSNX1)
    END IF
    IF (SCALAR_INT_BUF(371) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 370)
      CALL c_f_pointer(curr_ptr, RSNY1, shape=[fstarpu_vector_get_nx(buffers, 370)])
    ELSE
      NULLIFY(RSNY1)
    END IF
    IF (SCALAR_INT_BUF(372) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 371)
      CALL c_f_pointer(curr_ptr, RSNX2, shape=[fstarpu_vector_get_nx(buffers, 371)])
    ELSE
      NULLIFY(RSNX2)
    END IF
    IF (SCALAR_INT_BUF(373) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 372)
      CALL c_f_pointer(curr_ptr, RSNY2, shape=[fstarpu_vector_get_nx(buffers, 372)])
    ELSE
      NULLIFY(RSNY2)
    END IF
    IF (SCALAR_INT_BUF(374) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 373)
      CALL c_f_pointer(curr_ptr, RSNXOUT, shape=[fstarpu_vector_get_nx(buffers, 373)])
    ELSE
      NULLIFY(RSNXOUT)
    END IF
    IF (SCALAR_INT_BUF(375) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 374)
      CALL c_f_pointer(curr_ptr, RSNYOUT, shape=[fstarpu_vector_get_nx(buffers, 374)])
    ELSE
      NULLIFY(RSNYOUT)
    END IF
    IF (SCALAR_INT_BUF(376) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 375)
      CALL c_f_pointer(curr_ptr, WAVE_T1, shape=[fstarpu_vector_get_nx(buffers, 375)])
    ELSE
      NULLIFY(WAVE_T1)
    END IF
    IF (SCALAR_INT_BUF(377) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 376)
      CALL c_f_pointer(curr_ptr, WAVE_H1, shape=[fstarpu_vector_get_nx(buffers, 376)])
    ELSE
      NULLIFY(WAVE_H1)
    END IF
    IF (SCALAR_INT_BUF(378) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 377)
      CALL c_f_pointer(curr_ptr, WAVE_A1, shape=[fstarpu_vector_get_nx(buffers, 377)])
    ELSE
      NULLIFY(WAVE_A1)
    END IF
    IF (SCALAR_INT_BUF(379) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 378)
      CALL c_f_pointer(curr_ptr, WAVE_D1, shape=[fstarpu_vector_get_nx(buffers, 378)])
    ELSE
      NULLIFY(WAVE_D1)
    END IF
    IF (SCALAR_INT_BUF(380) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 379)
      CALL c_f_pointer(curr_ptr, WAVE_T2, shape=[fstarpu_vector_get_nx(buffers, 379)])
    ELSE
      NULLIFY(WAVE_T2)
    END IF
    IF (SCALAR_INT_BUF(381) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 380)
      CALL c_f_pointer(curr_ptr, WAVE_H2, shape=[fstarpu_vector_get_nx(buffers, 380)])
    ELSE
      NULLIFY(WAVE_H2)
    END IF
    IF (SCALAR_INT_BUF(382) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 381)
      CALL c_f_pointer(curr_ptr, WAVE_A2, shape=[fstarpu_vector_get_nx(buffers, 381)])
    ELSE
      NULLIFY(WAVE_A2)
    END IF
    IF (SCALAR_INT_BUF(383) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 382)
      CALL c_f_pointer(curr_ptr, WAVE_D2, shape=[fstarpu_vector_get_nx(buffers, 382)])
    ELSE
      NULLIFY(WAVE_D2)
    END IF
    IF (SCALAR_INT_BUF(384) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 383)
      CALL c_f_pointer(curr_ptr, WAVE_T, shape=[fstarpu_vector_get_nx(buffers, 383)])
    ELSE
      NULLIFY(WAVE_T)
    END IF
    IF (SCALAR_INT_BUF(385) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 384)
      CALL c_f_pointer(curr_ptr, WAVE_H, shape=[fstarpu_vector_get_nx(buffers, 384)])
    ELSE
      NULLIFY(WAVE_H)
    END IF
    IF (SCALAR_INT_BUF(386) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 385)
      CALL c_f_pointer(curr_ptr, WAVE_A, shape=[fstarpu_vector_get_nx(buffers, 385)])
    ELSE
      NULLIFY(WAVE_A)
    END IF
    IF (SCALAR_INT_BUF(387) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 386)
      CALL c_f_pointer(curr_ptr, WAVE_D, shape=[fstarpu_vector_get_nx(buffers, 386)])
    ELSE
      NULLIFY(WAVE_D)
    END IF
    IF (SCALAR_INT_BUF(388) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 387)
      CALL c_f_pointer(curr_ptr, WB, shape=[fstarpu_vector_get_nx(buffers, 387)])
    ELSE
      NULLIFY(WB)
    END IF
    IF (SCALAR_INT_BUF(389) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 388)
      CALL c_f_pointer(curr_ptr, WVNXOUT, shape=[fstarpu_vector_get_nx(buffers, 388)])
    ELSE
      NULLIFY(WVNXOUT)
    END IF
    IF (SCALAR_INT_BUF(390) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 389)
      CALL c_f_pointer(curr_ptr, WVNYOUT, shape=[fstarpu_vector_get_nx(buffers, 389)])
    ELSE
      NULLIFY(WVNYOUT)
    END IF
    IF (SCALAR_INT_BUF(391) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 390)
      CALL c_f_pointer(curr_ptr, TKXX, shape=[fstarpu_vector_get_nx(buffers, 390)])
    ELSE
      NULLIFY(TKXX)
    END IF
    IF (SCALAR_INT_BUF(392) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 391)
      CALL c_f_pointer(curr_ptr, TKYY, shape=[fstarpu_vector_get_nx(buffers, 391)])
    ELSE
      NULLIFY(TKYY)
    END IF
    IF (SCALAR_INT_BUF(393) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 392)
      CALL c_f_pointer(curr_ptr, TKXY, shape=[fstarpu_vector_get_nx(buffers, 392)])
    ELSE
      NULLIFY(TKXY)
    END IF
    IF (SCALAR_INT_BUF(394) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 393)
      CALL c_f_pointer(curr_ptr, EMO, shape=[fstarpu_matrix_get_nx(buffers, 393), fstarpu_matrix_get_ny(buffers, 393)])
    ELSE
      NULLIFY(EMO)
    END IF
    IF (SCALAR_INT_BUF(395) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 394)
      CALL c_f_pointer(curr_ptr, EFA, shape=[fstarpu_matrix_get_nx(buffers, 394), fstarpu_matrix_get_ny(buffers, 394)])
    ELSE
      NULLIFY(EFA)
    END IF
    IF (SCALAR_INT_BUF(396) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 395)
      CALL c_f_pointer(curr_ptr, UMO, shape=[fstarpu_matrix_get_nx(buffers, 395), fstarpu_matrix_get_ny(buffers, 395)])
    ELSE
      NULLIFY(UMO)
    END IF
    IF (SCALAR_INT_BUF(397) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 396)
      CALL c_f_pointer(curr_ptr, UFA, shape=[fstarpu_matrix_get_nx(buffers, 396), fstarpu_matrix_get_ny(buffers, 396)])
    ELSE
      NULLIFY(UFA)
    END IF
    IF (SCALAR_INT_BUF(398) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 397)
      CALL c_f_pointer(curr_ptr, VMO, shape=[fstarpu_matrix_get_nx(buffers, 397), fstarpu_matrix_get_ny(buffers, 397)])
    ELSE
      NULLIFY(VMO)
    END IF
    IF (SCALAR_INT_BUF(399) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 398)
      CALL c_f_pointer(curr_ptr, VFA, shape=[fstarpu_matrix_get_nx(buffers, 398), fstarpu_matrix_get_ny(buffers, 398)])
    ELSE
      NULLIFY(VFA)
    END IF
    IF (SCALAR_INT_BUF(400) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 399)
      CALL c_f_pointer(curr_ptr, XEL, shape=[fstarpu_vector_get_nx(buffers, 399)])
    ELSE
      NULLIFY(XEL)
    END IF
    IF (SCALAR_INT_BUF(401) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 400)
      CALL c_f_pointer(curr_ptr, YEL, shape=[fstarpu_vector_get_nx(buffers, 400)])
    ELSE
      NULLIFY(YEL)
    END IF
    IF (SCALAR_INT_BUF(402) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 401)
      CALL c_f_pointer(curr_ptr, SLEL, shape=[fstarpu_vector_get_nx(buffers, 401)])
    ELSE
      NULLIFY(SLEL)
    END IF
    IF (SCALAR_INT_BUF(403) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 402)
      CALL c_f_pointer(curr_ptr, SFEL, shape=[fstarpu_vector_get_nx(buffers, 402)])
    ELSE
      NULLIFY(SFEL)
    END IF
    IF (SCALAR_INT_BUF(404) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 403)
      CALL c_f_pointer(curr_ptr, AREAS, shape=[fstarpu_vector_get_nx(buffers, 403)])
    ELSE
      NULLIFY(AREAS)
    END IF
    IF (SCALAR_INT_BUF(405) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 404)
      CALL c_f_pointer(curr_ptr, SFACDUB, shape=[fstarpu_matrix_get_nx(buffers, 404), fstarpu_matrix_get_ny(buffers, 404)])
    ELSE
      NULLIFY(SFACDUB)
    END IF
    IF (SCALAR_INT_BUF(406) == 1) THEN
      curr_ptr = fstarpu_block_get_ptr(buffers, 405)
      CALL c_f_pointer(curr_ptr, YDUB, shape=[fstarpu_block_get_nx(buffers, 405), fstarpu_block_get_ny(buffers, 405), fstarpu_block_get_nz(buffers, 405)])
    ELSE
      NULLIFY(YDUB)
    END IF
    IF (SCALAR_INT_BUF(407) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 406)
      CALL c_f_pointer(curr_ptr, RTEMP2, shape=[fstarpu_vector_get_nx(buffers, 406)])
    ELSE
      NULLIFY(RTEMP2)
    END IF
    IF (SCALAR_INT_BUF(408) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 407)
      CALL c_f_pointer(curr_ptr, AUV11, shape=[fstarpu_vector_get_nx(buffers, 407)])
    ELSE
      NULLIFY(AUV11)
    END IF
    IF (SCALAR_INT_BUF(409) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 408)
      CALL c_f_pointer(curr_ptr, AUV12, shape=[fstarpu_vector_get_nx(buffers, 408)])
    ELSE
      NULLIFY(AUV12)
    END IF
    IF (SCALAR_INT_BUF(410) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 409)
      CALL c_f_pointer(curr_ptr, AUV13, shape=[fstarpu_vector_get_nx(buffers, 409)])
    ELSE
      NULLIFY(AUV13)
    END IF
    IF (SCALAR_INT_BUF(411) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 410)
      CALL c_f_pointer(curr_ptr, AUV14, shape=[fstarpu_vector_get_nx(buffers, 410)])
    ELSE
      NULLIFY(AUV14)
    END IF
    IF (SCALAR_INT_BUF(412) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 411)
      CALL c_f_pointer(curr_ptr, AUVXX, shape=[fstarpu_vector_get_nx(buffers, 411)])
    ELSE
      NULLIFY(AUVXX)
    END IF
    IF (SCALAR_INT_BUF(413) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 412)
      CALL c_f_pointer(curr_ptr, AUVYY, shape=[fstarpu_vector_get_nx(buffers, 412)])
    ELSE
      NULLIFY(AUVYY)
    END IF
    IF (SCALAR_INT_BUF(414) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 413)
      CALL c_f_pointer(curr_ptr, AUVXY, shape=[fstarpu_vector_get_nx(buffers, 413)])
    ELSE
      NULLIFY(AUVXY)
    END IF
    IF (SCALAR_INT_BUF(415) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 414)
      CALL c_f_pointer(curr_ptr, AUVYX, shape=[fstarpu_vector_get_nx(buffers, 414)])
    ELSE
      NULLIFY(AUVYX)
    END IF
    IF (SCALAR_INT_BUF(416) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 415)
      CALL c_f_pointer(curr_ptr, DUU1, shape=[fstarpu_vector_get_nx(buffers, 415)])
    ELSE
      NULLIFY(DUU1)
    END IF
    IF (SCALAR_INT_BUF(417) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 416)
      CALL c_f_pointer(curr_ptr, DUV1, shape=[fstarpu_vector_get_nx(buffers, 416)])
    ELSE
      NULLIFY(DUV1)
    END IF
    IF (SCALAR_INT_BUF(418) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 417)
      CALL c_f_pointer(curr_ptr, DVV1, shape=[fstarpu_vector_get_nx(buffers, 417)])
    ELSE
      NULLIFY(DVV1)
    END IF
    IF (SCALAR_INT_BUF(419) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 418)
      CALL c_f_pointer(curr_ptr, BSX1, shape=[fstarpu_vector_get_nx(buffers, 418)])
    ELSE
      NULLIFY(BSX1)
    END IF
    IF (SCALAR_INT_BUF(420) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 419)
      CALL c_f_pointer(curr_ptr, BSY1, shape=[fstarpu_vector_get_nx(buffers, 419)])
    ELSE
      NULLIFY(BSY1)
    END IF
    IF (SCALAR_INT_BUF(421) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 420)
      CALL c_f_pointer(curr_ptr, TIP1, shape=[fstarpu_vector_get_nx(buffers, 420)])
    ELSE
      NULLIFY(TIP1)
    END IF
    IF (SCALAR_INT_BUF(422) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 421)
      CALL c_f_pointer(curr_ptr, TIP2, shape=[fstarpu_vector_get_nx(buffers, 421)])
    ELSE
      NULLIFY(TIP2)
    END IF
    IF (SCALAR_INT_BUF(423) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 422)
      CALL c_f_pointer(curr_ptr, SALTAMP, shape=[fstarpu_matrix_get_nx(buffers, 422), fstarpu_matrix_get_ny(buffers, 422)])
    ELSE
      NULLIFY(SALTAMP)
    END IF
    IF (SCALAR_INT_BUF(424) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 423)
      CALL c_f_pointer(curr_ptr, SALTPHA, shape=[fstarpu_matrix_get_nx(buffers, 423), fstarpu_matrix_get_ny(buffers, 423)])
    ELSE
      NULLIFY(SALTPHA)
    END IF
    IF (SCALAR_INT_BUF(425) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 424)
      CALL c_f_pointer(curr_ptr, OBCCOEF, shape=[fstarpu_matrix_get_nx(buffers, 424), fstarpu_matrix_get_ny(buffers, 424)])
    ELSE
      NULLIFY(OBCCOEF)
    END IF
    IF (SCALAR_INT_BUF(426) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 425)
      CALL c_f_pointer(curr_ptr, COEF, shape=[fstarpu_matrix_get_nx(buffers, 425), fstarpu_matrix_get_ny(buffers, 425)])
    ELSE
      NULLIFY(COEF)
    END IF
    IF (SCALAR_INT_BUF(427) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 426)
      CALL c_f_pointer(curr_ptr, WKSP, shape=[fstarpu_vector_get_nx(buffers, 426)])
    ELSE
      NULLIFY(WKSP)
    END IF
    IF (SCALAR_INT_BUF(428) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 427)
      CALL c_f_pointer(curr_ptr, RPARM, shape=[fstarpu_vector_get_nx(buffers, 427)])
    ELSE
      NULLIFY(RPARM)
    END IF
    IF (SCALAR_INT_BUF(429) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 428)
      CALL c_f_pointer(curr_ptr, ABD, shape=[fstarpu_matrix_get_nx(buffers, 428), fstarpu_matrix_get_ny(buffers, 428)])
    ELSE
      NULLIFY(ABD)
    END IF
    IF (SCALAR_INT_BUF(430) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 429)
      CALL c_f_pointer(curr_ptr, ZX, shape=[fstarpu_vector_get_nx(buffers, 429)])
    ELSE
      NULLIFY(ZX)
    END IF
    IF (SCALAR_INT_BUF(431) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 430)
      CALL c_f_pointer(curr_ptr, GRAVX, shape=[fstarpu_vector_get_nx(buffers, 430)])
    ELSE
      NULLIFY(GRAVX)
    END IF
    IF (SCALAR_INT_BUF(432) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 431)
      CALL c_f_pointer(curr_ptr, GRAVY, shape=[fstarpu_vector_get_nx(buffers, 431)])
    ELSE
      NULLIFY(GRAVY)
    END IF
    IF (SCALAR_INT_BUF(433) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 432)
      CALL c_f_pointer(curr_ptr, ME2GW, shape=[fstarpu_vector_get_nx(buffers, 432)])
    ELSE
      NULLIFY(ME2GW)
    END IF
    IF (SCALAR_INT_BUF(434) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 433)
      CALL c_f_pointer(curr_ptr, NBV, shape=[fstarpu_vector_get_nx(buffers, 433)])
    ELSE
      NULLIFY(NBV)
    END IF
    IF (SCALAR_INT_BUF(435) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 434)
      CALL c_f_pointer(curr_ptr, LBCODEI, shape=[fstarpu_vector_get_nx(buffers, 434)])
    ELSE
      NULLIFY(LBCODEI)
    END IF
    IF (SCALAR_INT_BUF(436) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 435)
      CALL c_f_pointer(curr_ptr, NNODECODE, shape=[fstarpu_vector_get_nx(buffers, 435)])
    ELSE
      NULLIFY(NNODECODE)
    END IF
    IF (SCALAR_INT_BUF(437) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 436)
      CALL c_f_pointer(curr_ptr, NODECODE, shape=[fstarpu_vector_get_nx(buffers, 436)])
    ELSE
      NULLIFY(NODECODE)
    END IF
    IF (SCALAR_INT_BUF(438) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 437)
      CALL c_f_pointer(curr_ptr, NODEREP, shape=[fstarpu_vector_get_nx(buffers, 437)])
    ELSE
      NULLIFY(NODEREP)
    END IF
    IF (SCALAR_INT_BUF(439) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 438)
      CALL c_f_pointer(curr_ptr, NIBCNT, shape=[fstarpu_vector_get_nx(buffers, 438)])
    ELSE
      NULLIFY(NIBCNT)
    END IF
    IF (SCALAR_INT_BUF(440) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 439)
      CALL c_f_pointer(curr_ptr, NM, shape=[fstarpu_matrix_get_nx(buffers, 439), fstarpu_matrix_get_ny(buffers, 439)])
    ELSE
      NULLIFY(NM)
    END IF
    IF (SCALAR_INT_BUF(441) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 440)
      CALL c_f_pointer(curr_ptr, NNEIGH, shape=[fstarpu_vector_get_nx(buffers, 440)])
    ELSE
      NULLIFY(NNEIGH)
    END IF
    IF (SCALAR_INT_BUF(442) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 441)
      CALL c_f_pointer(curr_ptr, MJU, shape=[fstarpu_vector_get_nx(buffers, 441)])
    ELSE
      NULLIFY(MJU)
    END IF
    IF (SCALAR_INT_BUF(443) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 442)
      CALL c_f_pointer(curr_ptr, NODELE, shape=[fstarpu_vector_get_nx(buffers, 442)])
    ELSE
      NULLIFY(NODELE)
    END IF
    IF (SCALAR_INT_BUF(444) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 443)
      CALL c_f_pointer(curr_ptr, NEITAB, shape=[fstarpu_matrix_get_nx(buffers, 443), fstarpu_matrix_get_ny(buffers, 443)])
    ELSE
      NULLIFY(NEITAB)
    END IF
    IF (SCALAR_INT_BUF(445) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 444)
      CALL c_f_pointer(curr_ptr, NNEIGH_ELEM, shape=[fstarpu_vector_get_nx(buffers, 444)])
    ELSE
      NULLIFY(NNEIGH_ELEM)
    END IF
    IF (SCALAR_INT_BUF(446) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 445)
      CALL c_f_pointer(curr_ptr, NIBNODECODE, shape=[fstarpu_vector_get_nx(buffers, 445)])
    ELSE
      NULLIFY(NIBNODECODE)
    END IF
    IF (SCALAR_INT_BUF(447) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 446)
      CALL c_f_pointer(curr_ptr, NEIGH_ELEM, shape=[fstarpu_matrix_get_nx(buffers, 446), fstarpu_matrix_get_ny(buffers, 446)])
    ELSE
      NULLIFY(NEIGH_ELEM)
    END IF
    IF (SCALAR_INT_BUF(448) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 447)
      CALL c_f_pointer(curr_ptr, LBCODE, shape=[fstarpu_vector_get_nx(buffers, 447)])
    ELSE
      NULLIFY(LBCODE)
    END IF
    IF (SCALAR_INT_BUF(449) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 448)
      CALL c_f_pointer(curr_ptr, NNC, shape=[fstarpu_vector_get_nx(buffers, 448)])
    ELSE
      NULLIFY(NNC)
    END IF
    IF (SCALAR_INT_BUF(450) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 449)
      CALL c_f_pointer(curr_ptr, NNE, shape=[fstarpu_vector_get_nx(buffers, 449)])
    ELSE
      NULLIFY(NNE)
    END IF
    IF (SCALAR_INT_BUF(451) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 450)
      CALL c_f_pointer(curr_ptr, NNV, shape=[fstarpu_vector_get_nx(buffers, 450)])
    ELSE
      NULLIFY(NNV)
    END IF
    IF (SCALAR_INT_BUF(452) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 451)
      CALL c_f_pointer(curr_ptr, NNM, shape=[fstarpu_vector_get_nx(buffers, 451)])
    ELSE
      NULLIFY(NNM)
    END IF
    IF (SCALAR_INT_BUF(453) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 452)
      CALL c_f_pointer(curr_ptr, IWKSP, shape=[fstarpu_vector_get_nx(buffers, 452)])
    ELSE
      NULLIFY(IWKSP)
    END IF
    IF (SCALAR_INT_BUF(454) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 453)
      CALL c_f_pointer(curr_ptr, IPARM, shape=[fstarpu_vector_get_nx(buffers, 453)])
    ELSE
      NULLIFY(IPARM)
    END IF
    IF (SCALAR_INT_BUF(455) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 454)
      CALL c_f_pointer(curr_ptr, IPV, shape=[fstarpu_vector_get_nx(buffers, 454)])
    ELSE
      NULLIFY(IPV)
    END IF
    IF (SCALAR_INT_BUF(456) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 455)
      CALL c_f_pointer(curr_ptr, NVDLL, shape=[fstarpu_vector_get_nx(buffers, 455)])
    ELSE
      NULLIFY(NVDLL)
    END IF
    IF (SCALAR_INT_BUF(457) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 456)
      CALL c_f_pointer(curr_ptr, NBD, shape=[fstarpu_vector_get_nx(buffers, 456)])
    ELSE
      NULLIFY(NBD)
    END IF
    IF (SCALAR_INT_BUF(458) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 457)
      CALL c_f_pointer(curr_ptr, NBDV, shape=[fstarpu_matrix_get_nx(buffers, 457), fstarpu_matrix_get_ny(buffers, 457)])
    ELSE
      NULLIFY(NBDV)
    END IF
    IF (SCALAR_INT_BUF(459) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 458)
      CALL c_f_pointer(curr_ptr, NVELL, shape=[fstarpu_vector_get_nx(buffers, 458)])
    ELSE
      NULLIFY(NVELL)
    END IF
    IF (SCALAR_INT_BUF(460) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 459)
      CALL c_f_pointer(curr_ptr, NBVV, shape=[fstarpu_matrix_get_nx(buffers, 459), fstarpu_matrix_get_ny(buffers, 459)])
    ELSE
      NULLIFY(NBVV)
    END IF
    IF (SCALAR_INT_BUF(461) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 460)
      CALL c_f_pointer(curr_ptr, NELED, shape=[fstarpu_matrix_get_nx(buffers, 460), fstarpu_matrix_get_ny(buffers, 460)])
    ELSE
      NULLIFY(NELED)
    END IF
    IF (SCALAR_INT_BUF(462) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 461)
      CALL c_f_pointer(curr_ptr, SEGTYPE, shape=[fstarpu_vector_get_nx(buffers, 461)])
    ELSE
      NULLIFY(SEGTYPE)
    END IF
    IF (SCALAR_INT_BUF(463) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 462)
      CALL c_f_pointer(curr_ptr, NOT_AN_EDGE, shape=[fstarpu_vector_get_nx(buffers, 462)])
    ELSE
      NULLIFY(NOT_AN_EDGE)
    END IF
    IF (SCALAR_INT_BUF(464) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 463)
      CALL c_f_pointer(curr_ptr, WEIR_BUDDY_NODE, shape=[fstarpu_matrix_get_nx(buffers, 463), fstarpu_matrix_get_ny(buffers, 463)])
    ELSE
      NULLIFY(WEIR_BUDDY_NODE)
    END IF
    IF (SCALAR_INT_BUF(465) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 464)
      CALL c_f_pointer(curr_ptr, ONE_OR_TWO, shape=[fstarpu_vector_get_nx(buffers, 464)])
    ELSE
      NULLIFY(ONE_OR_TWO)
    END IF
    IF (SCALAR_INT_BUF(466) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 465)
      CALL c_f_pointer(curr_ptr, EDFLG, shape=[fstarpu_matrix_get_nx(buffers, 465), fstarpu_matrix_get_ny(buffers, 465)])
    ELSE
      NULLIFY(EDFLG)
    END IF
    IF (SCALAR_INT_BUF(467) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 466)
      CALL c_f_pointer(curr_ptr, BARLANHTR, shape=[fstarpu_vector_get_nx(buffers, 466)])
    ELSE
      NULLIFY(BARLANHTR)
    END IF
    IF (SCALAR_INT_BUF(468) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 467)
      CALL c_f_pointer(curr_ptr, BARLANCFSPR, shape=[fstarpu_vector_get_nx(buffers, 467)])
    ELSE
      NULLIFY(BARLANCFSPR)
    END IF
    IF (SCALAR_INT_BUF(469) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 468)
      CALL c_f_pointer(curr_ptr, BARINHTR, shape=[fstarpu_vector_get_nx(buffers, 468)])
    ELSE
      NULLIFY(BARINHTR)
    END IF
    IF (SCALAR_INT_BUF(470) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 469)
      CALL c_f_pointer(curr_ptr, BARINCFSBR, shape=[fstarpu_vector_get_nx(buffers, 469)])
    ELSE
      NULLIFY(BARINCFSBR)
    END IF
    IF (SCALAR_INT_BUF(471) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 470)
      CALL c_f_pointer(curr_ptr, BARINCFSPR, shape=[fstarpu_vector_get_nx(buffers, 470)])
    ELSE
      NULLIFY(BARINCFSPR)
    END IF
    IF (SCALAR_INT_BUF(472) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 471)
      CALL c_f_pointer(curr_ptr, PIPEHTR, shape=[fstarpu_vector_get_nx(buffers, 471)])
    ELSE
      NULLIFY(PIPEHTR)
    END IF
    IF (SCALAR_INT_BUF(473) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 472)
      CALL c_f_pointer(curr_ptr, PIPECOEFR, shape=[fstarpu_vector_get_nx(buffers, 472)])
    ELSE
      NULLIFY(PIPECOEFR)
    END IF
    IF (SCALAR_INT_BUF(474) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 473)
      CALL c_f_pointer(curr_ptr, PIPEDIAMR, shape=[fstarpu_vector_get_nx(buffers, 473)])
    ELSE
      NULLIFY(PIPEDIAMR)
    END IF
    IF (SCALAR_INT_BUF(475) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 474)
      CALL c_f_pointer(curr_ptr, BARLANHT, shape=[fstarpu_vector_get_nx(buffers, 474)])
    ELSE
      NULLIFY(BARLANHT)
    END IF
    IF (SCALAR_INT_BUF(476) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 475)
      CALL c_f_pointer(curr_ptr, BARLANCFSP, shape=[fstarpu_vector_get_nx(buffers, 475)])
    ELSE
      NULLIFY(BARLANCFSP)
    END IF
    IF (SCALAR_INT_BUF(477) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 476)
      CALL c_f_pointer(curr_ptr, FFF, shape=[fstarpu_vector_get_nx(buffers, 476)])
    ELSE
      NULLIFY(FFF)
    END IF
    IF (SCALAR_INT_BUF(478) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 477)
      CALL c_f_pointer(curr_ptr, FFACE, shape=[fstarpu_vector_get_nx(buffers, 477)])
    ELSE
      NULLIFY(FFACE)
    END IF
    IF (SCALAR_INT_BUF(479) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 478)
      CALL c_f_pointer(curr_ptr, BTRAN3, shape=[fstarpu_vector_get_nx(buffers, 478)])
    ELSE
      NULLIFY(BTRAN3)
    END IF
    IF (SCALAR_INT_BUF(480) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 479)
      CALL c_f_pointer(curr_ptr, BTRAN4, shape=[fstarpu_vector_get_nx(buffers, 479)])
    ELSE
      NULLIFY(BTRAN4)
    END IF
    IF (SCALAR_INT_BUF(481) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 480)
      CALL c_f_pointer(curr_ptr, BTRAN5, shape=[fstarpu_vector_get_nx(buffers, 480)])
    ELSE
      NULLIFY(BTRAN5)
    END IF
    IF (SCALAR_INT_BUF(482) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 481)
      CALL c_f_pointer(curr_ptr, BTRAN6, shape=[fstarpu_vector_get_nx(buffers, 481)])
    ELSE
      NULLIFY(BTRAN6)
    END IF
    IF (SCALAR_INT_BUF(483) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 482)
      CALL c_f_pointer(curr_ptr, BTRAN7, shape=[fstarpu_vector_get_nx(buffers, 482)])
    ELSE
      NULLIFY(BTRAN7)
    END IF
    IF (SCALAR_INT_BUF(484) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 483)
      CALL c_f_pointer(curr_ptr, BTRAN8, shape=[fstarpu_vector_get_nx(buffers, 483)])
    ELSE
      NULLIFY(BTRAN8)
    END IF
    IF (SCALAR_INT_BUF(485) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 484)
      CALL c_f_pointer(curr_ptr, BARINHT, shape=[fstarpu_vector_get_nx(buffers, 484)])
    ELSE
      NULLIFY(BARINHT)
    END IF
    IF (SCALAR_INT_BUF(486) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 485)
      CALL c_f_pointer(curr_ptr, BARINCFSB, shape=[fstarpu_vector_get_nx(buffers, 485)])
    ELSE
      NULLIFY(BARINCFSB)
    END IF
    IF (SCALAR_INT_BUF(487) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 486)
      CALL c_f_pointer(curr_ptr, BARINCFSP, shape=[fstarpu_vector_get_nx(buffers, 486)])
    ELSE
      NULLIFY(BARINCFSP)
    END IF
    IF (SCALAR_INT_BUF(488) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 487)
      CALL c_f_pointer(curr_ptr, PIPEHT, shape=[fstarpu_vector_get_nx(buffers, 487)])
    ELSE
      NULLIFY(PIPEHT)
    END IF
    IF (SCALAR_INT_BUF(489) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 488)
      CALL c_f_pointer(curr_ptr, PIPECOEF, shape=[fstarpu_vector_get_nx(buffers, 488)])
    ELSE
      NULLIFY(PIPECOEF)
    END IF
    IF (SCALAR_INT_BUF(490) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 489)
      CALL c_f_pointer(curr_ptr, PIPEDIAM, shape=[fstarpu_vector_get_nx(buffers, 489)])
    ELSE
      NULLIFY(PIPEDIAM)
    END IF
    IF (SCALAR_INT_BUF(491) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 490)
      CALL c_f_pointer(curr_ptr, RBARWL1AVG, shape=[fstarpu_vector_get_nx(buffers, 490)])
    ELSE
      NULLIFY(RBARWL1AVG)
    END IF
    IF (SCALAR_INT_BUF(492) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 491)
      CALL c_f_pointer(curr_ptr, RBARWL2AVG, shape=[fstarpu_vector_get_nx(buffers, 491)])
    ELSE
      NULLIFY(RBARWL2AVG)
    END IF
    IF (SCALAR_INT_BUF(493) == 1) THEN
      curr_ptr = fstarpu_matrix_get_ptr(buffers, 492)
      CALL c_f_pointer(curr_ptr, ELEXLEN, shape=[fstarpu_matrix_get_nx(buffers, 492), fstarpu_matrix_get_ny(buffers, 492)])
    ELSE
      NULLIFY(ELEXLEN)
    END IF
    IF (SCALAR_INT_BUF(494) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 493)
      CALL c_f_pointer(curr_ptr, IBCONN, shape=[fstarpu_vector_get_nx(buffers, 493)])
    ELSE
      NULLIFY(IBCONN)
    END IF
    IF (SCALAR_INT_BUF(495) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 494)
      CALL c_f_pointer(curr_ptr, IBCONNR, shape=[fstarpu_vector_get_nx(buffers, 494)])
    ELSE
      NULLIFY(IBCONNR)
    END IF
    IF (SCALAR_INT_BUF(496) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 495)
      CALL c_f_pointer(curr_ptr, NTRAN1, shape=[fstarpu_vector_get_nx(buffers, 495)])
    ELSE
      NULLIFY(NTRAN1)
    END IF
    IF (SCALAR_INT_BUF(497) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 496)
      CALL c_f_pointer(curr_ptr, NTRAN2, shape=[fstarpu_vector_get_nx(buffers, 496)])
    ELSE
      NULLIFY(NTRAN2)
    END IF
    IF (SCALAR_INT_BUF(498) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 497)
      CALL c_f_pointer(curr_ptr, BK, shape=[fstarpu_vector_get_nx(buffers, 497)])
    ELSE
      NULLIFY(BK)
    END IF
    IF (SCALAR_INT_BUF(499) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 498)
      CALL c_f_pointer(curr_ptr, BALPHA, shape=[fstarpu_vector_get_nx(buffers, 498)])
    ELSE
      NULLIFY(BALPHA)
    END IF
    IF (SCALAR_INT_BUF(500) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 499)
      CALL c_f_pointer(curr_ptr, BDELX, shape=[fstarpu_vector_get_nx(buffers, 499)])
    ELSE
      NULLIFY(BDELX)
    END IF
    IF (SCALAR_INT_BUF(501) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 500)
      CALL c_f_pointer(curr_ptr, NBNNUM, shape=[fstarpu_vector_get_nx(buffers, 500)])
    ELSE
      NULLIFY(NBNNUM)
    END IF
    IF (SCALAR_INT_BUF(502) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 501)
      CALL c_f_pointer(curr_ptr, AMIG, shape=[fstarpu_vector_get_nx(buffers, 501)])
    ELSE
      NULLIFY(AMIG)
    END IF
    IF (SCALAR_INT_BUF(503) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 502)
      CALL c_f_pointer(curr_ptr, AMIGT, shape=[fstarpu_vector_get_nx(buffers, 502)])
    ELSE
      NULLIFY(AMIGT)
    END IF
    IF (SCALAR_INT_BUF(504) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 503)
      CALL c_f_pointer(curr_ptr, FAMIG, shape=[fstarpu_vector_get_nx(buffers, 503)])
    ELSE
      NULLIFY(FAMIG)
    END IF
    IF (SCALAR_INT_BUF(505) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 504)
      CALL c_f_pointer(curr_ptr, PER, shape=[fstarpu_vector_get_nx(buffers, 504)])
    ELSE
      NULLIFY(PER)
    END IF
    IF (SCALAR_INT_BUF(506) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 505)
      CALL c_f_pointer(curr_ptr, PERT, shape=[fstarpu_vector_get_nx(buffers, 505)])
    ELSE
      NULLIFY(PERT)
    END IF
    IF (SCALAR_INT_BUF(507) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 506)
      CALL c_f_pointer(curr_ptr, FPER, shape=[fstarpu_vector_get_nx(buffers, 506)])
    ELSE
      NULLIFY(FPER)
    END IF
    IF (SCALAR_INT_BUF(508) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 507)
      CALL c_f_pointer(curr_ptr, FREQ, shape=[fstarpu_vector_get_nx(buffers, 507)])
    ELSE
      NULLIFY(FREQ)
    END IF
    IF (SCALAR_INT_BUF(509) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 508)
      CALL c_f_pointer(curr_ptr, FF, shape=[fstarpu_vector_get_nx(buffers, 508)])
    ELSE
      NULLIFY(FF)
    END IF
    IF (SCALAR_INT_BUF(510) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 509)
      CALL c_f_pointer(curr_ptr, FACE, shape=[fstarpu_vector_get_nx(buffers, 509)])
    ELSE
      NULLIFY(FACE)
    END IF
    IF (SCALAR_INT_BUF(511) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 510)
      CALL c_f_pointer(curr_ptr, SLAM, shape=[fstarpu_vector_get_nx(buffers, 510)])
    ELSE
      NULLIFY(SLAM)
    END IF
    IF (SCALAR_INT_BUF(512) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 511)
      CALL c_f_pointer(curr_ptr, SFEA, shape=[fstarpu_vector_get_nx(buffers, 511)])
    ELSE
      NULLIFY(SFEA)
    END IF
    IF (SCALAR_INT_BUF(513) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 512)
      CALL c_f_pointer(curr_ptr, X, shape=[fstarpu_vector_get_nx(buffers, 512)])
    ELSE
      NULLIFY(X)
    END IF
    IF (SCALAR_INT_BUF(514) == 1) THEN
      curr_ptr = fstarpu_vector_get_ptr(buffers, 513)
      CALL c_f_pointer(curr_ptr, Y, shape=[fstarpu_vector_get_nx(buffers, 513)])
    ELSE
      NULLIFY(Y)
    END IF
  END SUBROUTINE DGSWEM_STATE_ACTIVATE

  SUBROUTINE DGSWEM_DEALLOC_ALLOCATABLES()
#ifdef  REAL4
#endif
#ifdef REAL16
#endif
#ifndef REAL4
#ifndef REAL16
#endif
#endif
#if defined(CMPI) || defined(TASKPAR)
#else
#endif
#ifdef SWAN
    IF (ALLOCATED(iarray)) DEALLOCATE(iarray)
    IF (ALLOCATED(iarray_g)) DEALLOCATE(iarray_g)
    IF (ALLOCATED(imap)) DEALLOCATE(imap)
    IF (ALLOCATED(array)) DEALLOCATE(array)
    IF (ALLOCATED(array2)) DEALLOCATE(array2)
    IF (ALLOCATED(array3)) DEALLOCATE(array3)
    IF (ALLOCATED(array_g)) DEALLOCATE(array_g)
    IF (ALLOCATED(array2_g)) DEALLOCATE(array2_g)
    IF (ALLOCATED(array3_g)) DEALLOCATE(array3_g)
    IF (ALLOCATED(hotstart)) DEALLOCATE(hotstart)
    IF (ALLOCATED(hotstart_g)) DEALLOCATE(hotstart_g)
#endif
    IF (ALLOCATED(TIPOTAG)) DEALLOCATE(TIPOTAG)
    IF (ALLOCATED(BOUNTAG)) DEALLOCATE(BOUNTAG)
    IF (ALLOCATED(FBOUNTAG)) DEALLOCATE(FBOUNTAG)
    IF (ALLOCATED(SwanWaveRefrac)) DEALLOCATE(SwanWaveRefrac)
    IF (ALLOCATED(STARTDRY)) DEALLOCATE(STARTDRY)
    IF (ALLOCATED(FRIC)) DEALLOCATE(FRIC)
    IF (ALLOCATED(TAU0VAR)) DEALLOCATE(TAU0VAR)
    IF (ALLOCATED(TAU0BASE)) DEALLOCATE(TAU0BASE)
    IF (ALLOCATED(z0land)) DEALLOCATE(z0land)
    IF (ALLOCATED(vcanopy)) DEALLOCATE(vcanopy)
    IF (ALLOCATED(BridgePilings)) DEALLOCATE(BridgePilings)
    IF (ALLOCATED(Chezy)) DEALLOCATE(Chezy)
    IF (ALLOCATED(ManningsN)) DEALLOCATE(ManningsN)
    IF (ALLOCATED(GeoidOffset)) DEALLOCATE(GeoidOffset)
    IF (ALLOCATED(EVM)) DEALLOCATE(EVM)
    IF (ALLOCATED(EVC)) DEALLOCATE(EVC)
  END SUBROUTINE DGSWEM_DEALLOC_ALLOCATABLES

END MODULE DAGSWEM_STATE
