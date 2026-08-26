MODULE DAGSWEM_STATE
  USE SIZES
  USE GLOBAL
  USE DG
  USE FSTARPU_MOD
  USE ISO_C_BINDING
  IMPLICIT NONE

  INTEGER, PARAMETER :: NUM_STATE_HANDLES = 250

  CONTAINS

  SUBROUTINE DGSWEM_STATE_REGISTER(handles)
    TYPE(C_PTR), INTENT(OUT) :: handles(NUM_STATE_HANDLES)
    IF (ASSOCIATED(WDFLG)) THEN
      CALL fstarpu_vector_data_register(handles(1), 0, C_LOC(WDFLG(LBOUND(WDFLG,1))), SIZE(WDFLG,1), C_SIZEOF(WDFLG(LBOUND(WDFLG,1))))
    ELSE
      handles(1) = C_NULL_PTR
    END IF
    NULLIFY(WDFLG)
    IF (ASSOCIATED(DOFS)) THEN
      CALL fstarpu_vector_data_register(handles(2), 0, C_LOC(DOFS(LBOUND(DOFS,1))), SIZE(DOFS,1), C_SIZEOF(DOFS(LBOUND(DOFS,1))))
    ELSE
      handles(2) = C_NULL_PTR
    END IF
    NULLIFY(DOFS)
    IF (ASSOCIATED(PCOUNT)) THEN
      CALL fstarpu_vector_data_register(handles(3), 0, C_LOC(PCOUNT(LBOUND(PCOUNT,1))), SIZE(PCOUNT,1), C_SIZEOF(PCOUNT(LBOUND(PCOUNT,1))))
    ELSE
      handles(3) = C_NULL_PTR
    END IF
    NULLIFY(PCOUNT)
    IF (ASSOCIATED(PDG)) THEN
      CALL fstarpu_vector_data_register(handles(4), 0, C_LOC(PDG(LBOUND(PDG,1))), SIZE(PDG,1), C_SIZEOF(PDG(LBOUND(PDG,1))))
    ELSE
      handles(4) = C_NULL_PTR
    END IF
    NULLIFY(PDG)
    IF (ASSOCIATED(NCOUNT)) THEN
      CALL fstarpu_vector_data_register(handles(5), 0, C_LOC(NCOUNT(LBOUND(NCOUNT,1))), SIZE(NCOUNT,1), C_SIZEOF(NCOUNT(LBOUND(NCOUNT,1))))
    ELSE
      handles(5) = C_NULL_PTR
    END IF
    NULLIFY(NCOUNT)
    IF (ASSOCIATED(NEDEL)) THEN
      CALL fstarpu_matrix_data_register(handles(6), 0, C_LOC(NEDEL(LBOUND(NEDEL,1),LBOUND(NEDEL,2))), SIZE(NEDEL,1), SIZE(NEDEL,1), SIZE(NEDEL,2), C_SIZEOF(NEDEL(LBOUND(NEDEL,1),LBOUND(NEDEL,2))))
    ELSE
      handles(6) = C_NULL_PTR
    END IF
    NULLIFY(NEDEL)
    IF (ASSOCIATED(NEDSD)) THEN
      CALL fstarpu_matrix_data_register(handles(7), 0, C_LOC(NEDSD(LBOUND(NEDSD,1),LBOUND(NEDSD,2))), SIZE(NEDSD,1), SIZE(NEDSD,1), SIZE(NEDSD,2), C_SIZEOF(NEDSD(LBOUND(NEDSD,1),LBOUND(NEDSD,2))))
    ELSE
      handles(7) = C_NULL_PTR
    END IF
    NULLIFY(NEDSD)
    IF (ASSOCIATED(NEDNO)) THEN
      CALL fstarpu_matrix_data_register(handles(8), 0, C_LOC(NEDNO(LBOUND(NEDNO,1),LBOUND(NEDNO,2))), SIZE(NEDNO,1), SIZE(NEDNO,1), SIZE(NEDNO,2), C_SIZEOF(NEDNO(LBOUND(NEDNO,1),LBOUND(NEDNO,2))))
    ELSE
      handles(8) = C_NULL_PTR
    END IF
    NULLIFY(NEDNO)
    IF (ASSOCIATED(NEDNO1)) THEN
      CALL fstarpu_vector_data_register(handles(9), 0, C_LOC(NEDNO1(LBOUND(NEDNO1,1))), SIZE(NEDNO1,1), C_SIZEOF(NEDNO1(LBOUND(NEDNO1,1))))
    ELSE
      handles(9) = C_NULL_PTR
    END IF
    NULLIFY(NEDNO1)
    IF (ASSOCIATED(NEDNO2)) THEN
      CALL fstarpu_vector_data_register(handles(10), 0, C_LOC(NEDNO2(LBOUND(NEDNO2,1))), SIZE(NEDNO2,1), C_SIZEOF(NEDNO2(LBOUND(NEDNO2,1))))
    ELSE
      handles(10) = C_NULL_PTR
    END IF
    NULLIFY(NEDNO2)
    IF (ASSOCIATED(NIEDN)) THEN
      CALL fstarpu_vector_data_register(handles(11), 0, C_LOC(NIEDN(LBOUND(NIEDN,1))), SIZE(NIEDN,1), C_SIZEOF(NIEDN(LBOUND(NIEDN,1))))
    ELSE
      handles(11) = C_NULL_PTR
    END IF
    NULLIFY(NIEDN)
    IF (ASSOCIATED(NLEDN)) THEN
      CALL fstarpu_vector_data_register(handles(12), 0, C_LOC(NLEDN(LBOUND(NLEDN,1))), SIZE(NLEDN,1), C_SIZEOF(NLEDN(LBOUND(NLEDN,1))))
    ELSE
      handles(12) = C_NULL_PTR
    END IF
    NULLIFY(NLEDN)
    IF (ASSOCIATED(NEEDN)) THEN
      CALL fstarpu_vector_data_register(handles(13), 0, C_LOC(NEEDN(LBOUND(NEEDN,1))), SIZE(NEEDN,1), C_SIZEOF(NEEDN(LBOUND(NEEDN,1))))
    ELSE
      handles(13) = C_NULL_PTR
    END IF
    NULLIFY(NEEDN)
    IF (ASSOCIATED(NFEDN)) THEN
      CALL fstarpu_vector_data_register(handles(14), 0, C_LOC(NFEDN(LBOUND(NFEDN,1))), SIZE(NFEDN,1), C_SIZEOF(NFEDN(LBOUND(NFEDN,1))))
    ELSE
      handles(14) = C_NULL_PTR
    END IF
    NULLIFY(NFEDN)
    IF (ASSOCIATED(NREDN)) THEN
      CALL fstarpu_vector_data_register(handles(15), 0, C_LOC(NREDN(LBOUND(NREDN,1))), SIZE(NREDN,1), C_SIZEOF(NREDN(LBOUND(NREDN,1))))
    ELSE
      handles(15) = C_NULL_PTR
    END IF
    NULLIFY(NREDN)
    IF (ASSOCIATED(NEBEDN)) THEN
      CALL fstarpu_vector_data_register(handles(16), 0, C_LOC(NEBEDN(LBOUND(NEBEDN,1))), SIZE(NEBEDN,1), C_SIZEOF(NEBEDN(LBOUND(NEBEDN,1))))
    ELSE
      handles(16) = C_NULL_PTR
    END IF
    NULLIFY(NEBEDN)
    IF (ASSOCIATED(NIBEDN)) THEN
      CALL fstarpu_vector_data_register(handles(17), 0, C_LOC(NIBEDN(LBOUND(NIBEDN,1))), SIZE(NIBEDN,1), C_SIZEOF(NIBEDN(LBOUND(NIBEDN,1))))
    ELSE
      handles(17) = C_NULL_PTR
    END IF
    NULLIFY(NIBEDN)
    IF (ASSOCIATED(NIBSEGN)) THEN
      CALL fstarpu_matrix_data_register(handles(18), 0, C_LOC(NIBSEGN(LBOUND(NIBSEGN,1),LBOUND(NIBSEGN,2))), SIZE(NIBSEGN,1), SIZE(NIBSEGN,1), SIZE(NIBSEGN,2), C_SIZEOF(NIBSEGN(LBOUND(NIBSEGN,1),LBOUND(NIBSEGN,2))))
    ELSE
      handles(18) = C_NULL_PTR
    END IF
    NULLIFY(NIBSEGN)
    IF (ASSOCIATED(NEBSEGN)) THEN
      CALL fstarpu_vector_data_register(handles(19), 0, C_LOC(NEBSEGN(LBOUND(NEBSEGN,1))), SIZE(NEBSEGN,1), C_SIZEOF(NEBSEGN(LBOUND(NEBSEGN,1))))
    ELSE
      handles(19) = C_NULL_PTR
    END IF
    NULLIFY(NEBSEGN)
    IF (ASSOCIATED(EL_NBORS)) THEN
      CALL fstarpu_matrix_data_register(handles(20), 0, C_LOC(EL_NBORS(LBOUND(EL_NBORS,1),LBOUND(EL_NBORS,2))), SIZE(EL_NBORS,1), SIZE(EL_NBORS,1), SIZE(EL_NBORS,2), C_SIZEOF(EL_NBORS(LBOUND(EL_NBORS,1),LBOUND(EL_NBORS,2))))
    ELSE
      handles(20) = C_NULL_PTR
    END IF
    NULLIFY(EL_NBORS)
    IF (ASSOCIATED(BACKNODES)) THEN
      CALL fstarpu_matrix_data_register(handles(21), 0, C_LOC(BACKNODES(LBOUND(BACKNODES,1),LBOUND(BACKNODES,2))), SIZE(BACKNODES,1), SIZE(BACKNODES,1), SIZE(BACKNODES,2), C_SIZEOF(BACKNODES(LBOUND(BACKNODES,1),LBOUND(BACKNODES,2))))
    ELSE
      handles(21) = C_NULL_PTR
    END IF
    NULLIFY(BACKNODES)
    IF (ASSOCIATED(ATVD)) THEN
      CALL fstarpu_matrix_data_register(handles(22), 0, C_LOC(ATVD(LBOUND(ATVD,1),LBOUND(ATVD,2))), SIZE(ATVD,1), SIZE(ATVD,1), SIZE(ATVD,2), C_SIZEOF(ATVD(LBOUND(ATVD,1),LBOUND(ATVD,2))))
    ELSE
      handles(22) = C_NULL_PTR
    END IF
    NULLIFY(ATVD)
    IF (ASSOCIATED(BTVD)) THEN
      CALL fstarpu_matrix_data_register(handles(23), 0, C_LOC(BTVD(LBOUND(BTVD,1),LBOUND(BTVD,2))), SIZE(BTVD,1), SIZE(BTVD,1), SIZE(BTVD,2), C_SIZEOF(BTVD(LBOUND(BTVD,1),LBOUND(BTVD,2))))
    ELSE
      handles(23) = C_NULL_PTR
    END IF
    NULLIFY(BTVD)
    IF (ASSOCIATED(CTVD)) THEN
      CALL fstarpu_matrix_data_register(handles(24), 0, C_LOC(CTVD(LBOUND(CTVD,1),LBOUND(CTVD,2))), SIZE(CTVD,1), SIZE(CTVD,1), SIZE(CTVD,2), C_SIZEOF(CTVD(LBOUND(CTVD,1),LBOUND(CTVD,2))))
    ELSE
      handles(24) = C_NULL_PTR
    END IF
    NULLIFY(CTVD)
    IF (ASSOCIATED(DTVD)) THEN
      CALL fstarpu_vector_data_register(handles(25), 0, C_LOC(DTVD(LBOUND(DTVD,1))), SIZE(DTVD,1), C_SIZEOF(DTVD(LBOUND(DTVD,1))))
    ELSE
      handles(25) = C_NULL_PTR
    END IF
    NULLIFY(DTVD)
    IF (ASSOCIATED(MAX_BOA_DT)) THEN
      CALL fstarpu_vector_data_register(handles(26), 0, C_LOC(MAX_BOA_DT(LBOUND(MAX_BOA_DT,1))), SIZE(MAX_BOA_DT,1), C_SIZEOF(MAX_BOA_DT(LBOUND(MAX_BOA_DT,1))))
    ELSE
      handles(26) = C_NULL_PTR
    END IF
    NULLIFY(MAX_BOA_DT)
    IF (ASSOCIATED(e1)) THEN
      CALL fstarpu_vector_data_register(handles(27), 0, C_LOC(e1(LBOUND(e1,1))), SIZE(e1,1), C_SIZEOF(e1(LBOUND(e1,1))))
    ELSE
      handles(27) = C_NULL_PTR
    END IF
    NULLIFY(e1)
    IF (ASSOCIATED(balance)) THEN
      CALL fstarpu_vector_data_register(handles(28), 0, C_LOC(balance(LBOUND(balance,1))), SIZE(balance,1), C_SIZEOF(balance(LBOUND(balance,1))))
    ELSE
      handles(28) = C_NULL_PTR
    END IF
    NULLIFY(balance)
    IF (ASSOCIATED(RKC_Tdprime)) THEN
      CALL fstarpu_vector_data_register(handles(29), 0, C_LOC(RKC_Tdprime(LBOUND(RKC_Tdprime,1))), SIZE(RKC_Tdprime,1), C_SIZEOF(RKC_Tdprime(LBOUND(RKC_Tdprime,1))))
    ELSE
      handles(29) = C_NULL_PTR
    END IF
    NULLIFY(RKC_Tdprime)
    IF (ASSOCIATED(RKC_a)) THEN
      CALL fstarpu_vector_data_register(handles(30), 0, C_LOC(RKC_a(LBOUND(RKC_a,1))), SIZE(RKC_a,1), C_SIZEOF(RKC_a(LBOUND(RKC_a,1))))
    ELSE
      handles(30) = C_NULL_PTR
    END IF
    NULLIFY(RKC_a)
    IF (ASSOCIATED(RKC_b)) THEN
      CALL fstarpu_vector_data_register(handles(31), 0, C_LOC(RKC_b(LBOUND(RKC_b,1))), SIZE(RKC_b,1), C_SIZEOF(RKC_b(LBOUND(RKC_b,1))))
    ELSE
      handles(31) = C_NULL_PTR
    END IF
    NULLIFY(RKC_b)
    IF (ASSOCIATED(RKC_c)) THEN
      CALL fstarpu_vector_data_register(handles(32), 0, C_LOC(RKC_c(LBOUND(RKC_c,1))), SIZE(RKC_c,1), C_SIZEOF(RKC_c(LBOUND(RKC_c,1))))
    ELSE
      handles(32) = C_NULL_PTR
    END IF
    NULLIFY(RKC_c)
    IF (ASSOCIATED(RKC_mu)) THEN
      CALL fstarpu_vector_data_register(handles(33), 0, C_LOC(RKC_mu(LBOUND(RKC_mu,1))), SIZE(RKC_mu,1), C_SIZEOF(RKC_mu(LBOUND(RKC_mu,1))))
    ELSE
      handles(33) = C_NULL_PTR
    END IF
    NULLIFY(RKC_mu)
    IF (ASSOCIATED(RKC_tildemu)) THEN
      CALL fstarpu_vector_data_register(handles(34), 0, C_LOC(RKC_tildemu(LBOUND(RKC_tildemu,1))), SIZE(RKC_tildemu,1), C_SIZEOF(RKC_tildemu(LBOUND(RKC_tildemu,1))))
    ELSE
      handles(34) = C_NULL_PTR
    END IF
    NULLIFY(RKC_tildemu)
    IF (ASSOCIATED(RKC_nu)) THEN
      CALL fstarpu_vector_data_register(handles(35), 0, C_LOC(RKC_nu(LBOUND(RKC_nu,1))), SIZE(RKC_nu,1), C_SIZEOF(RKC_nu(LBOUND(RKC_nu,1))))
    ELSE
      handles(35) = C_NULL_PTR
    END IF
    NULLIFY(RKC_nu)
    IF (ASSOCIATED(BATH)) THEN
      CALL fstarpu_block_data_register(handles(36), 0, C_LOC(BATH(LBOUND(BATH,1),LBOUND(BATH,2),LBOUND(BATH,3))), SIZE(BATH,1), SIZE(BATH,1)*SIZE(BATH,2), SIZE(BATH,1), SIZE(BATH,2), SIZE(BATH,3), C_SIZEOF(BATH(LBOUND(BATH,1),LBOUND(BATH,2),LBOUND(BATH,3))))
    ELSE
      handles(36) = C_NULL_PTR
    END IF
    NULLIFY(BATH)
    IF (ASSOCIATED(DBATHDX)) THEN
      CALL fstarpu_block_data_register(handles(37), 0, C_LOC(DBATHDX(LBOUND(DBATHDX,1),LBOUND(DBATHDX,2),LBOUND(DBATHDX,3))), SIZE(DBATHDX,1), SIZE(DBATHDX,1)*SIZE(DBATHDX,2), SIZE(DBATHDX,1), SIZE(DBATHDX,2), SIZE(DBATHDX,3), C_SIZEOF(DBATHDX(LBOUND(DBATHDX,1),LBOUND(DBATHDX,2),LBOUND(DBATHDX,3))))
    ELSE
      handles(37) = C_NULL_PTR
    END IF
    NULLIFY(DBATHDX)
    IF (ASSOCIATED(DBATHDY)) THEN
      CALL fstarpu_block_data_register(handles(38), 0, C_LOC(DBATHDY(LBOUND(DBATHDY,1),LBOUND(DBATHDY,2),LBOUND(DBATHDY,3))), SIZE(DBATHDY,1), SIZE(DBATHDY,1)*SIZE(DBATHDY,2), SIZE(DBATHDY,1), SIZE(DBATHDY,2), SIZE(DBATHDY,3), C_SIZEOF(DBATHDY(LBOUND(DBATHDY,1),LBOUND(DBATHDY,2),LBOUND(DBATHDY,3))))
    ELSE
      handles(38) = C_NULL_PTR
    END IF
    NULLIFY(DBATHDY)
    IF (ASSOCIATED(SFAC_ELEM)) THEN
      CALL fstarpu_block_data_register(handles(39), 0, C_LOC(SFAC_ELEM(LBOUND(SFAC_ELEM,1),LBOUND(SFAC_ELEM,2),LBOUND(SFAC_ELEM,3))), SIZE(SFAC_ELEM,1), SIZE(SFAC_ELEM,1)*SIZE(SFAC_ELEM,2), SIZE(SFAC_ELEM,1), SIZE(SFAC_ELEM,2), SIZE(SFAC_ELEM,3), C_SIZEOF(SFAC_ELEM(LBOUND(SFAC_ELEM,1),LBOUND(SFAC_ELEM,2),LBOUND(SFAC_ELEM,3))))
    ELSE
      handles(39) = C_NULL_PTR
    END IF
    NULLIFY(SFAC_ELEM)
    IF (ASSOCIATED(BATHED)) THEN
      CALL fstarpu_tensor_data_register(handles(40), 0, C_LOC(BATHED(LBOUND(BATHED,1),LBOUND(BATHED,2),LBOUND(BATHED,3),LBOUND(BATHED,4))), SIZE(BATHED,1), SIZE(BATHED,1)*SIZE(BATHED,2), SIZE(BATHED,1)*SIZE(BATHED,2)*SIZE(BATHED,3), SIZE(BATHED,1), SIZE(BATHED,2), SIZE(BATHED,3), SIZE(BATHED,4), C_SIZEOF(BATHED(LBOUND(BATHED,1),LBOUND(BATHED,2),LBOUND(BATHED,3),LBOUND(BATHED,4))))
    ELSE
      handles(40) = C_NULL_PTR
    END IF
    NULLIFY(BATHED)
    IF (ASSOCIATED(SFACED)) THEN
      CALL fstarpu_tensor_data_register(handles(41), 0, C_LOC(SFACED(LBOUND(SFACED,1),LBOUND(SFACED,2),LBOUND(SFACED,3),LBOUND(SFACED,4))), SIZE(SFACED,1), SIZE(SFACED,1)*SIZE(SFACED,2), SIZE(SFACED,1)*SIZE(SFACED,2)*SIZE(SFACED,3), SIZE(SFACED,1), SIZE(SFACED,2), SIZE(SFACED,3), SIZE(SFACED,4), C_SIZEOF(SFACED(LBOUND(SFACED,1),LBOUND(SFACED,2),LBOUND(SFACED,3),LBOUND(SFACED,4))))
    ELSE
      handles(41) = C_NULL_PTR
    END IF
    NULLIFY(SFACED)
    IF (ASSOCIATED(COSNX)) THEN
      CALL fstarpu_vector_data_register(handles(42), 0, C_LOC(COSNX(LBOUND(COSNX,1))), SIZE(COSNX,1), C_SIZEOF(COSNX(LBOUND(COSNX,1))))
    ELSE
      handles(42) = C_NULL_PTR
    END IF
    NULLIFY(COSNX)
    IF (ASSOCIATED(SINNX)) THEN
      CALL fstarpu_vector_data_register(handles(43), 0, C_LOC(SINNX(LBOUND(SINNX,1))), SIZE(SINNX,1), C_SIZEOF(SINNX(LBOUND(SINNX,1))))
    ELSE
      handles(43) = C_NULL_PTR
    END IF
    NULLIFY(SINNX)
    IF (ASSOCIATED(DP_NODE)) THEN
      CALL fstarpu_block_data_register(handles(44), 0, C_LOC(DP_NODE(LBOUND(DP_NODE,1),LBOUND(DP_NODE,2),LBOUND(DP_NODE,3))), SIZE(DP_NODE,1), SIZE(DP_NODE,1)*SIZE(DP_NODE,2), SIZE(DP_NODE,1), SIZE(DP_NODE,2), SIZE(DP_NODE,3), C_SIZEOF(DP_NODE(LBOUND(DP_NODE,1),LBOUND(DP_NODE,2),LBOUND(DP_NODE,3))))
    ELSE
      handles(44) = C_NULL_PTR
    END IF
    NULLIFY(DP_NODE)
    IF (ASSOCIATED(DP_VOL)) THEN
      CALL fstarpu_matrix_data_register(handles(45), 0, C_LOC(DP_VOL(LBOUND(DP_VOL,1),LBOUND(DP_VOL,2))), SIZE(DP_VOL,1), SIZE(DP_VOL,1), SIZE(DP_VOL,2), C_SIZEOF(DP_VOL(LBOUND(DP_VOL,1),LBOUND(DP_VOL,2))))
    ELSE
      handles(45) = C_NULL_PTR
    END IF
    NULLIFY(DP_VOL)
    IF (ASSOCIATED(DRPHI)) THEN
      CALL fstarpu_block_data_register(handles(46), 0, C_LOC(DRPHI(LBOUND(DRPHI,1),LBOUND(DRPHI,2),LBOUND(DRPHI,3))), SIZE(DRPHI,1), SIZE(DRPHI,1)*SIZE(DRPHI,2), SIZE(DRPHI,1), SIZE(DRPHI,2), SIZE(DRPHI,3), C_SIZEOF(DRPHI(LBOUND(DRPHI,1),LBOUND(DRPHI,2),LBOUND(DRPHI,3))))
    ELSE
      handles(46) = C_NULL_PTR
    END IF
    NULLIFY(DRPHI)
    IF (ASSOCIATED(DSPHI)) THEN
      CALL fstarpu_block_data_register(handles(47), 0, C_LOC(DSPHI(LBOUND(DSPHI,1),LBOUND(DSPHI,2),LBOUND(DSPHI,3))), SIZE(DSPHI,1), SIZE(DSPHI,1)*SIZE(DSPHI,2), SIZE(DSPHI,1), SIZE(DSPHI,2), SIZE(DSPHI,3), C_SIZEOF(DSPHI(LBOUND(DSPHI,1),LBOUND(DSPHI,2),LBOUND(DSPHI,3))))
    ELSE
      handles(47) = C_NULL_PTR
    END IF
    NULLIFY(DSPHI)
    IF (ASSOCIATED(DRDX)) THEN
      CALL fstarpu_vector_data_register(handles(48), 0, C_LOC(DRDX(LBOUND(DRDX,1))), SIZE(DRDX,1), C_SIZEOF(DRDX(LBOUND(DRDX,1))))
    ELSE
      handles(48) = C_NULL_PTR
    END IF
    NULLIFY(DRDX)
    IF (ASSOCIATED(DSDX)) THEN
      CALL fstarpu_vector_data_register(handles(49), 0, C_LOC(DSDX(LBOUND(DSDX,1))), SIZE(DSDX,1), C_SIZEOF(DSDX(LBOUND(DSDX,1))))
    ELSE
      handles(49) = C_NULL_PTR
    END IF
    NULLIFY(DSDX)
    IF (ASSOCIATED(DRDY)) THEN
      CALL fstarpu_vector_data_register(handles(50), 0, C_LOC(DRDY(LBOUND(DRDY,1))), SIZE(DRDY,1), C_SIZEOF(DRDY(LBOUND(DRDY,1))))
    ELSE
      handles(50) = C_NULL_PTR
    END IF
    NULLIFY(DRDY)
    IF (ASSOCIATED(DSDY)) THEN
      CALL fstarpu_vector_data_register(handles(51), 0, C_LOC(DSDY(LBOUND(DSDY,1))), SIZE(DSDY,1), C_SIZEOF(DSDY(LBOUND(DSDY,1))))
    ELSE
      handles(51) = C_NULL_PTR
    END IF
    NULLIFY(DSDY)
    IF (ASSOCIATED(EFA_DG)) THEN
      CALL fstarpu_block_data_register(handles(52), 0, C_LOC(EFA_DG(LBOUND(EFA_DG,1),LBOUND(EFA_DG,2),LBOUND(EFA_DG,3))), SIZE(EFA_DG,1), SIZE(EFA_DG,1)*SIZE(EFA_DG,2), SIZE(EFA_DG,1), SIZE(EFA_DG,2), SIZE(EFA_DG,3), C_SIZEOF(EFA_DG(LBOUND(EFA_DG,1),LBOUND(EFA_DG,2),LBOUND(EFA_DG,3))))
    ELSE
      handles(52) = C_NULL_PTR
    END IF
    NULLIFY(EFA_DG)
    IF (ASSOCIATED(EMO_DG)) THEN
      CALL fstarpu_block_data_register(handles(53), 0, C_LOC(EMO_DG(LBOUND(EMO_DG,1),LBOUND(EMO_DG,2),LBOUND(EMO_DG,3))), SIZE(EMO_DG,1), SIZE(EMO_DG,1)*SIZE(EMO_DG,2), SIZE(EMO_DG,1), SIZE(EMO_DG,2), SIZE(EMO_DG,3), C_SIZEOF(EMO_DG(LBOUND(EMO_DG,1),LBOUND(EMO_DG,2),LBOUND(EMO_DG,3))))
    ELSE
      handles(53) = C_NULL_PTR
    END IF
    NULLIFY(EMO_DG)
    IF (ASSOCIATED(XLEN)) THEN
      CALL fstarpu_vector_data_register(handles(54), 0, C_LOC(XLEN(LBOUND(XLEN,1))), SIZE(XLEN,1), C_SIZEOF(XLEN(LBOUND(XLEN,1))))
    ELSE
      handles(54) = C_NULL_PTR
    END IF
    NULLIFY(XLEN)
    IF (ASSOCIATED(HB)) THEN
      CALL fstarpu_block_data_register(handles(55), 0, C_LOC(HB(LBOUND(HB,1),LBOUND(HB,2),LBOUND(HB,3))), SIZE(HB,1), SIZE(HB,1)*SIZE(HB,2), SIZE(HB,1), SIZE(HB,2), SIZE(HB,3), C_SIZEOF(HB(LBOUND(HB,1),LBOUND(HB,2),LBOUND(HB,3))))
    ELSE
      handles(55) = C_NULL_PTR
    END IF
    NULLIFY(HB)
    IF (ASSOCIATED(IBHT)) THEN
      CALL fstarpu_vector_data_register(handles(56), 0, C_LOC(IBHT(LBOUND(IBHT,1))), SIZE(IBHT,1), C_SIZEOF(IBHT(LBOUND(IBHT,1))))
    ELSE
      handles(56) = C_NULL_PTR
    END IF
    NULLIFY(IBHT)
    IF (ASSOCIATED(EBHT)) THEN
      CALL fstarpu_vector_data_register(handles(57), 0, C_LOC(EBHT(LBOUND(EBHT,1))), SIZE(EBHT,1), C_SIZEOF(EBHT(LBOUND(EBHT,1))))
    ELSE
      handles(57) = C_NULL_PTR
    END IF
    NULLIFY(EBHT)
    IF (ASSOCIATED(EBCFSP)) THEN
      CALL fstarpu_vector_data_register(handles(58), 0, C_LOC(EBCFSP(LBOUND(EBCFSP,1))), SIZE(EBCFSP,1), C_SIZEOF(EBCFSP(LBOUND(EBCFSP,1))))
    ELSE
      handles(58) = C_NULL_PTR
    END IF
    NULLIFY(EBCFSP)
    IF (ASSOCIATED(IBCFSP)) THEN
      CALL fstarpu_vector_data_register(handles(59), 0, C_LOC(IBCFSP(LBOUND(IBCFSP,1))), SIZE(IBCFSP,1), C_SIZEOF(IBCFSP(LBOUND(IBCFSP,1))))
    ELSE
      handles(59) = C_NULL_PTR
    END IF
    NULLIFY(IBCFSP)
    IF (ASSOCIATED(IBCFSB)) THEN
      CALL fstarpu_vector_data_register(handles(60), 0, C_LOC(IBCFSB(LBOUND(IBCFSB,1))), SIZE(IBCFSB,1), C_SIZEOF(IBCFSB(LBOUND(IBCFSB,1))))
    ELSE
      handles(60) = C_NULL_PTR
    END IF
    NULLIFY(IBCFSB)
    IF (ASSOCIATED(M_INV)) THEN
      CALL fstarpu_matrix_data_register(handles(61), 0, C_LOC(M_INV(LBOUND(M_INV,1),LBOUND(M_INV,2))), SIZE(M_INV,1), SIZE(M_INV,1), SIZE(M_INV,2), C_SIZEOF(M_INV(LBOUND(M_INV,1),LBOUND(M_INV,2))))
    ELSE
      handles(61) = C_NULL_PTR
    END IF
    NULLIFY(M_INV)
    IF (ASSOCIATED(phi_edge_fixed)) THEN
      CALL fstarpu_block_data_register(handles(62), 0, C_LOC(phi_edge_fixed(LBOUND(phi_edge_fixed,1),LBOUND(phi_edge_fixed,2),LBOUND(phi_edge_fixed,3))), SIZE(phi_edge_fixed,1), SIZE(phi_edge_fixed,1)*SIZE(phi_edge_fixed,2), SIZE(phi_edge_fixed,1), SIZE(phi_edge_fixed,2), SIZE(phi_edge_fixed,3), C_SIZEOF(phi_edge_fixed(LBOUND(phi_edge_fixed,1),LBOUND(phi_edge_fixed,2),LBOUND(phi_edge_fixed,3))))
    ELSE
      handles(62) = C_NULL_PTR
    END IF
    NULLIFY(phi_edge_fixed)
    IF (ASSOCIATED(PHI_AREA)) THEN
      CALL fstarpu_block_data_register(handles(63), 0, C_LOC(PHI_AREA(LBOUND(PHI_AREA,1),LBOUND(PHI_AREA,2),LBOUND(PHI_AREA,3))), SIZE(PHI_AREA,1), SIZE(PHI_AREA,1)*SIZE(PHI_AREA,2), SIZE(PHI_AREA,1), SIZE(PHI_AREA,2), SIZE(PHI_AREA,3), C_SIZEOF(PHI_AREA(LBOUND(PHI_AREA,1),LBOUND(PHI_AREA,2),LBOUND(PHI_AREA,3))))
    ELSE
      handles(63) = C_NULL_PTR
    END IF
    NULLIFY(PHI_AREA)
    IF (ASSOCIATED(PHI_EDGE)) THEN
      CALL fstarpu_tensor_data_register(handles(64), 0, C_LOC(PHI_EDGE(LBOUND(PHI_EDGE,1),LBOUND(PHI_EDGE,2),LBOUND(PHI_EDGE,3),LBOUND(PHI_EDGE,4))), SIZE(PHI_EDGE,1), SIZE(PHI_EDGE,1)*SIZE(PHI_EDGE,2), SIZE(PHI_EDGE,1)*SIZE(PHI_EDGE,2)*SIZE(PHI_EDGE,3), SIZE(PHI_EDGE,1), SIZE(PHI_EDGE,2), SIZE(PHI_EDGE,3), SIZE(PHI_EDGE,4), C_SIZEOF(PHI_EDGE(LBOUND(PHI_EDGE,1),LBOUND(PHI_EDGE,2),LBOUND(PHI_EDGE,3),LBOUND(PHI_EDGE,4))))
    ELSE
      handles(64) = C_NULL_PTR
    END IF
    NULLIFY(PHI_EDGE)
    IF (ASSOCIATED(PHI_CENTER)) THEN
      CALL fstarpu_matrix_data_register(handles(65), 0, C_LOC(PHI_CENTER(LBOUND(PHI_CENTER,1),LBOUND(PHI_CENTER,2))), SIZE(PHI_CENTER,1), SIZE(PHI_CENTER,1), SIZE(PHI_CENTER,2), C_SIZEOF(PHI_CENTER(LBOUND(PHI_CENTER,1),LBOUND(PHI_CENTER,2))))
    ELSE
      handles(65) = C_NULL_PTR
    END IF
    NULLIFY(PHI_CENTER)
    IF (ASSOCIATED(PHI_CORNER)) THEN
      CALL fstarpu_block_data_register(handles(66), 0, C_LOC(PHI_CORNER(LBOUND(PHI_CORNER,1),LBOUND(PHI_CORNER,2),LBOUND(PHI_CORNER,3))), SIZE(PHI_CORNER,1), SIZE(PHI_CORNER,1)*SIZE(PHI_CORNER,2), SIZE(PHI_CORNER,1), SIZE(PHI_CORNER,2), SIZE(PHI_CORNER,3), C_SIZEOF(PHI_CORNER(LBOUND(PHI_CORNER,1),LBOUND(PHI_CORNER,2),LBOUND(PHI_CORNER,3))))
    ELSE
      handles(66) = C_NULL_PTR
    END IF
    NULLIFY(PHI_CORNER)
    IF (ASSOCIATED(PHI_CHECK)) THEN
      CALL fstarpu_block_data_register(handles(67), 0, C_LOC(PHI_CHECK(LBOUND(PHI_CHECK,1),LBOUND(PHI_CHECK,2),LBOUND(PHI_CHECK,3))), SIZE(PHI_CHECK,1), SIZE(PHI_CHECK,1)*SIZE(PHI_CHECK,2), SIZE(PHI_CHECK,1), SIZE(PHI_CHECK,2), SIZE(PHI_CHECK,3), C_SIZEOF(PHI_CHECK(LBOUND(PHI_CHECK,1),LBOUND(PHI_CHECK,2),LBOUND(PHI_CHECK,3))))
    ELSE
      handles(67) = C_NULL_PTR
    END IF
    NULLIFY(PHI_CHECK)
    IF (ASSOCIATED(PHI_INTEGRATED)) THEN
      CALL fstarpu_matrix_data_register(handles(68), 0, C_LOC(PHI_INTEGRATED(LBOUND(PHI_INTEGRATED,1),LBOUND(PHI_INTEGRATED,2))), SIZE(PHI_INTEGRATED,1), SIZE(PHI_INTEGRATED,1), SIZE(PHI_INTEGRATED,2), C_SIZEOF(PHI_INTEGRATED(LBOUND(PHI_INTEGRATED,1),LBOUND(PHI_INTEGRATED,2))))
    ELSE
      handles(68) = C_NULL_PTR
    END IF
    NULLIFY(PHI_INTEGRATED)
    IF (ASSOCIATED(PSI1)) THEN
      CALL fstarpu_matrix_data_register(handles(69), 0, C_LOC(PSI1(LBOUND(PSI1,1),LBOUND(PSI1,2))), SIZE(PSI1,1), SIZE(PSI1,1), SIZE(PSI1,2), C_SIZEOF(PSI1(LBOUND(PSI1,1),LBOUND(PSI1,2))))
    ELSE
      handles(69) = C_NULL_PTR
    END IF
    NULLIFY(PSI1)
    IF (ASSOCIATED(PSI2)) THEN
      CALL fstarpu_matrix_data_register(handles(70), 0, C_LOC(PSI2(LBOUND(PSI2,1),LBOUND(PSI2,2))), SIZE(PSI2,1), SIZE(PSI2,1), SIZE(PSI2,2), C_SIZEOF(PSI2(LBOUND(PSI2,1),LBOUND(PSI2,2))))
    ELSE
      handles(70) = C_NULL_PTR
    END IF
    NULLIFY(PSI2)
    IF (ASSOCIATED(PSI3)) THEN
      CALL fstarpu_matrix_data_register(handles(71), 0, C_LOC(PSI3(LBOUND(PSI3,1),LBOUND(PSI3,2))), SIZE(PSI3,1), SIZE(PSI3,1), SIZE(PSI3,2), C_SIZEOF(PSI3(LBOUND(PSI3,1),LBOUND(PSI3,2))))
    ELSE
      handles(71) = C_NULL_PTR
    END IF
    NULLIFY(PSI3)
    IF (ASSOCIATED(QIB)) THEN
      CALL fstarpu_vector_data_register(handles(72), 0, C_LOC(QIB(LBOUND(QIB,1))), SIZE(QIB,1), C_SIZEOF(QIB(LBOUND(QIB,1))))
    ELSE
      handles(72) = C_NULL_PTR
    END IF
    NULLIFY(QIB)
    IF (ASSOCIATED(QX)) THEN
      CALL fstarpu_block_data_register(handles(73), 0, C_LOC(QX(LBOUND(QX,1),LBOUND(QX,2),LBOUND(QX,3))), SIZE(QX,1), SIZE(QX,1)*SIZE(QX,2), SIZE(QX,1), SIZE(QX,2), SIZE(QX,3), C_SIZEOF(QX(LBOUND(QX,1),LBOUND(QX,2),LBOUND(QX,3))))
    ELSE
      handles(73) = C_NULL_PTR
    END IF
    NULLIFY(QX)
    IF (ASSOCIATED(QY)) THEN
      CALL fstarpu_block_data_register(handles(74), 0, C_LOC(QY(LBOUND(QY,1),LBOUND(QY,2),LBOUND(QY,3))), SIZE(QY,1), SIZE(QY,1)*SIZE(QY,2), SIZE(QY,1), SIZE(QY,2), SIZE(QY,3), C_SIZEOF(QY(LBOUND(QY,1),LBOUND(QY,2),LBOUND(QY,3))))
    ELSE
      handles(74) = C_NULL_PTR
    END IF
    NULLIFY(QY)
    IF (ASSOCIATED(ZE)) THEN
      CALL fstarpu_block_data_register(handles(75), 0, C_LOC(ZE(LBOUND(ZE,1),LBOUND(ZE,2),LBOUND(ZE,3))), SIZE(ZE,1), SIZE(ZE,1)*SIZE(ZE,2), SIZE(ZE,1), SIZE(ZE,2), SIZE(ZE,3), C_SIZEOF(ZE(LBOUND(ZE,1),LBOUND(ZE,2),LBOUND(ZE,3))))
    ELSE
      handles(75) = C_NULL_PTR
    END IF
    NULLIFY(ZE)
    IF (ASSOCIATED(ze_edge)) THEN
      CALL fstarpu_block_data_register(handles(76), 0, C_LOC(ze_edge(LBOUND(ze_edge,1),LBOUND(ze_edge,2),LBOUND(ze_edge,3))), SIZE(ze_edge,1), SIZE(ze_edge,1)*SIZE(ze_edge,2), SIZE(ze_edge,1), SIZE(ze_edge,2), SIZE(ze_edge,3), C_SIZEOF(ze_edge(LBOUND(ze_edge,1),LBOUND(ze_edge,2),LBOUND(ze_edge,3))))
    ELSE
      handles(76) = C_NULL_PTR
    END IF
    NULLIFY(ze_edge)
    IF (ASSOCIATED(qx_edge)) THEN
      CALL fstarpu_block_data_register(handles(77), 0, C_LOC(qx_edge(LBOUND(qx_edge,1),LBOUND(qx_edge,2),LBOUND(qx_edge,3))), SIZE(qx_edge,1), SIZE(qx_edge,1)*SIZE(qx_edge,2), SIZE(qx_edge,1), SIZE(qx_edge,2), SIZE(qx_edge,3), C_SIZEOF(qx_edge(LBOUND(qx_edge,1),LBOUND(qx_edge,2),LBOUND(qx_edge,3))))
    ELSE
      handles(77) = C_NULL_PTR
    END IF
    NULLIFY(qx_edge)
    IF (ASSOCIATED(qy_edge)) THEN
      CALL fstarpu_block_data_register(handles(78), 0, C_LOC(qy_edge(LBOUND(qy_edge,1),LBOUND(qy_edge,2),LBOUND(qy_edge,3))), SIZE(qy_edge,1), SIZE(qy_edge,1)*SIZE(qy_edge,2), SIZE(qy_edge,1), SIZE(qy_edge,2), SIZE(qy_edge,3), C_SIZEOF(qy_edge(LBOUND(qy_edge,1),LBOUND(qy_edge,2),LBOUND(qy_edge,3))))
    ELSE
      handles(78) = C_NULL_PTR
    END IF
    NULLIFY(qy_edge)
    IF (ASSOCIATED(elem_edge)) THEN
      CALL fstarpu_matrix_data_register(handles(79), 0, C_LOC(elem_edge(LBOUND(elem_edge,1),LBOUND(elem_edge,2))), SIZE(elem_edge,1), SIZE(elem_edge,1), SIZE(elem_edge,2), C_SIZEOF(elem_edge(LBOUND(elem_edge,1),LBOUND(elem_edge,2))))
    ELSE
      handles(79) = C_NULL_PTR
    END IF
    NULLIFY(elem_edge)
    IF (ASSOCIATED(nieds_count)) THEN
      CALL fstarpu_vector_data_register(handles(80), 0, C_LOC(nieds_count(LBOUND(nieds_count,1))), SIZE(nieds_count,1), C_SIZEOF(nieds_count(LBOUND(nieds_count,1))))
    ELSE
      handles(80) = C_NULL_PTR
    END IF
    NULLIFY(nieds_count)
    IF (ASSOCIATED(bed)) THEN
      CALL fstarpu_tensor_data_register(handles(81), 0, C_LOC(bed(LBOUND(bed,1),LBOUND(bed,2),LBOUND(bed,3),LBOUND(bed,4))), SIZE(bed,1), SIZE(bed,1)*SIZE(bed,2), SIZE(bed,1)*SIZE(bed,2)*SIZE(bed,3), SIZE(bed,1), SIZE(bed,2), SIZE(bed,3), SIZE(bed,4), C_SIZEOF(bed(LBOUND(bed,1),LBOUND(bed,2),LBOUND(bed,3),LBOUND(bed,4))))
    ELSE
      handles(81) = C_NULL_PTR
    END IF
    NULLIFY(bed)
    IF (ASSOCIATED(dynP)) THEN
      CALL fstarpu_block_data_register(handles(82), 0, C_LOC(dynP(LBOUND(dynP,1),LBOUND(dynP,2),LBOUND(dynP,3))), SIZE(dynP,1), SIZE(dynP,1)*SIZE(dynP,2), SIZE(dynP,1), SIZE(dynP,2), SIZE(dynP,3), C_SIZEOF(dynP(LBOUND(dynP,1),LBOUND(dynP,2),LBOUND(dynP,3))))
    ELSE
      handles(82) = C_NULL_PTR
    END IF
    NULLIFY(dynP)
    IF (ASSOCIATED(dynP_MAX)) THEN
      CALL fstarpu_vector_data_register(handles(83), 0, C_LOC(dynP_MAX(LBOUND(dynP_MAX,1))), SIZE(dynP_MAX,1), C_SIZEOF(dynP_MAX(LBOUND(dynP_MAX,1))))
    ELSE
      handles(83) = C_NULL_PTR
    END IF
    NULLIFY(dynP_MAX)
    IF (ASSOCIATED(dynP_MIN)) THEN
      CALL fstarpu_vector_data_register(handles(84), 0, C_LOC(dynP_MIN(LBOUND(dynP_MIN,1))), SIZE(dynP_MIN,1), C_SIZEOF(dynP_MIN(LBOUND(dynP_MIN,1))))
    ELSE
      handles(84) = C_NULL_PTR
    END IF
    NULLIFY(dynP_MIN)
    IF (ASSOCIATED(iota)) THEN
      CALL fstarpu_block_data_register(handles(85), 0, C_LOC(iota(LBOUND(iota,1),LBOUND(iota,2),LBOUND(iota,3))), SIZE(iota,1), SIZE(iota,1)*SIZE(iota,2), SIZE(iota,1), SIZE(iota,2), SIZE(iota,3), C_SIZEOF(iota(LBOUND(iota,1),LBOUND(iota,2),LBOUND(iota,3))))
    ELSE
      handles(85) = C_NULL_PTR
    END IF
    NULLIFY(iota)
    IF (ASSOCIATED(iotaa)) THEN
      CALL fstarpu_block_data_register(handles(86), 0, C_LOC(iotaa(LBOUND(iotaa,1),LBOUND(iotaa,2),LBOUND(iotaa,3))), SIZE(iotaa,1), SIZE(iotaa,1)*SIZE(iotaa,2), SIZE(iotaa,1), SIZE(iotaa,2), SIZE(iotaa,3), C_SIZEOF(iotaa(LBOUND(iotaa,1),LBOUND(iotaa,2),LBOUND(iotaa,3))))
    ELSE
      handles(86) = C_NULL_PTR
    END IF
    NULLIFY(iotaa)
    IF (ASSOCIATED(iota2)) THEN
      CALL fstarpu_block_data_register(handles(87), 0, C_LOC(iota2(LBOUND(iota2,1),LBOUND(iota2,2),LBOUND(iota2,3))), SIZE(iota2,1), SIZE(iota2,1)*SIZE(iota2,2), SIZE(iota2,1), SIZE(iota2,2), SIZE(iota2,3), C_SIZEOF(iota2(LBOUND(iota2,1),LBOUND(iota2,2),LBOUND(iota2,3))))
    ELSE
      handles(87) = C_NULL_PTR
    END IF
    NULLIFY(iota2)
    IF (ASSOCIATED(arrayfix)) THEN
      CALL fstarpu_block_data_register(handles(88), 0, C_LOC(arrayfix(LBOUND(arrayfix,1),LBOUND(arrayfix,2),LBOUND(arrayfix,3))), SIZE(arrayfix,1), SIZE(arrayfix,1)*SIZE(arrayfix,2), SIZE(arrayfix,1), SIZE(arrayfix,2), SIZE(arrayfix,3), C_SIZEOF(arrayfix(LBOUND(arrayfix,1),LBOUND(arrayfix,2),LBOUND(arrayfix,3))))
    ELSE
      handles(88) = C_NULL_PTR
    END IF
    NULLIFY(arrayfix)
    IF (ASSOCIATED(CORI_EL)) THEN
      CALL fstarpu_vector_data_register(handles(89), 0, C_LOC(CORI_EL(LBOUND(CORI_EL,1))), SIZE(CORI_EL,1), C_SIZEOF(CORI_EL(LBOUND(CORI_EL,1))))
    ELSE
      handles(89) = C_NULL_PTR
    END IF
    NULLIFY(CORI_EL)
    IF (ASSOCIATED(FRIC_EL)) THEN
      CALL fstarpu_vector_data_register(handles(90), 0, C_LOC(FRIC_EL(LBOUND(FRIC_EL,1))), SIZE(FRIC_EL,1), C_SIZEOF(FRIC_EL(LBOUND(FRIC_EL,1))))
    ELSE
      handles(90) = C_NULL_PTR
    END IF
    NULLIFY(FRIC_EL)
    IF (ASSOCIATED(ZE_MAX)) THEN
      CALL fstarpu_vector_data_register(handles(91), 0, C_LOC(ZE_MAX(LBOUND(ZE_MAX,1))), SIZE(ZE_MAX,1), C_SIZEOF(ZE_MAX(LBOUND(ZE_MAX,1))))
    ELSE
      handles(91) = C_NULL_PTR
    END IF
    NULLIFY(ZE_MAX)
    IF (ASSOCIATED(ZE_MIN)) THEN
      CALL fstarpu_vector_data_register(handles(92), 0, C_LOC(ZE_MIN(LBOUND(ZE_MIN,1))), SIZE(ZE_MIN,1), C_SIZEOF(ZE_MIN(LBOUND(ZE_MIN,1))))
    ELSE
      handles(92) = C_NULL_PTR
    END IF
    NULLIFY(ZE_MIN)
    IF (ASSOCIATED(DPE_MIN)) THEN
      CALL fstarpu_vector_data_register(handles(93), 0, C_LOC(DPE_MIN(LBOUND(DPE_MIN,1))), SIZE(DPE_MIN,1), C_SIZEOF(DPE_MIN(LBOUND(DPE_MIN,1))))
    ELSE
      handles(93) = C_NULL_PTR
    END IF
    NULLIFY(DPE_MIN)
    IF (ASSOCIATED(ADVECTQX)) THEN
      CALL fstarpu_vector_data_register(handles(94), 0, C_LOC(ADVECTQX(LBOUND(ADVECTQX,1))), SIZE(ADVECTQX,1), C_SIZEOF(ADVECTQX(LBOUND(ADVECTQX,1))))
    ELSE
      handles(94) = C_NULL_PTR
    END IF
    NULLIFY(ADVECTQX)
    IF (ASSOCIATED(ADVECTQY)) THEN
      CALL fstarpu_vector_data_register(handles(95), 0, C_LOC(ADVECTQY(LBOUND(ADVECTQY,1))), SIZE(ADVECTQY,1), C_SIZEOF(ADVECTQY(LBOUND(ADVECTQY,1))))
    ELSE
      handles(95) = C_NULL_PTR
    END IF
    NULLIFY(ADVECTQY)
    IF (ASSOCIATED(SOURCEQX)) THEN
      CALL fstarpu_vector_data_register(handles(96), 0, C_LOC(SOURCEQX(LBOUND(SOURCEQX,1))), SIZE(SOURCEQX,1), C_SIZEOF(SOURCEQX(LBOUND(SOURCEQX,1))))
    ELSE
      handles(96) = C_NULL_PTR
    END IF
    NULLIFY(SOURCEQX)
    IF (ASSOCIATED(SOURCEQY)) THEN
      CALL fstarpu_vector_data_register(handles(97), 0, C_LOC(SOURCEQY(LBOUND(SOURCEQY,1))), SIZE(SOURCEQY,1), C_SIZEOF(SOURCEQY(LBOUND(SOURCEQY,1))))
    ELSE
      handles(97) = C_NULL_PTR
    END IF
    NULLIFY(SOURCEQY)
    IF (ASSOCIATED(LZ)) THEN
      CALL fstarpu_tensor_data_register(handles(98), 0, C_LOC(LZ(LBOUND(LZ,1),LBOUND(LZ,2),LBOUND(LZ,3),LBOUND(LZ,4))), SIZE(LZ,1), SIZE(LZ,1)*SIZE(LZ,2), SIZE(LZ,1)*SIZE(LZ,2)*SIZE(LZ,3), SIZE(LZ,1), SIZE(LZ,2), SIZE(LZ,3), SIZE(LZ,4), C_SIZEOF(LZ(LBOUND(LZ,1),LBOUND(LZ,2),LBOUND(LZ,3),LBOUND(LZ,4))))
    ELSE
      handles(98) = C_NULL_PTR
    END IF
    NULLIFY(LZ)
    IF (ASSOCIATED(MZ)) THEN
      CALL fstarpu_tensor_data_register(handles(99), 0, C_LOC(MZ(LBOUND(MZ,1),LBOUND(MZ,2),LBOUND(MZ,3),LBOUND(MZ,4))), SIZE(MZ,1), SIZE(MZ,1)*SIZE(MZ,2), SIZE(MZ,1)*SIZE(MZ,2)*SIZE(MZ,3), SIZE(MZ,1), SIZE(MZ,2), SIZE(MZ,3), SIZE(MZ,4), C_SIZEOF(MZ(LBOUND(MZ,1),LBOUND(MZ,2),LBOUND(MZ,3),LBOUND(MZ,4))))
    ELSE
      handles(99) = C_NULL_PTR
    END IF
    NULLIFY(MZ)
    IF (ASSOCIATED(HZ)) THEN
      CALL fstarpu_tensor_data_register(handles(100), 0, C_LOC(HZ(LBOUND(HZ,1),LBOUND(HZ,2),LBOUND(HZ,3),LBOUND(HZ,4))), SIZE(HZ,1), SIZE(HZ,1)*SIZE(HZ,2), SIZE(HZ,1)*SIZE(HZ,2)*SIZE(HZ,3), SIZE(HZ,1), SIZE(HZ,2), SIZE(HZ,3), SIZE(HZ,4), C_SIZEOF(HZ(LBOUND(HZ,1),LBOUND(HZ,2),LBOUND(HZ,3),LBOUND(HZ,4))))
    ELSE
      handles(100) = C_NULL_PTR
    END IF
    NULLIFY(HZ)
    IF (ASSOCIATED(TZ)) THEN
      CALL fstarpu_tensor_data_register(handles(101), 0, C_LOC(TZ(LBOUND(TZ,1),LBOUND(TZ,2),LBOUND(TZ,3),LBOUND(TZ,4))), SIZE(TZ,1), SIZE(TZ,1)*SIZE(TZ,2), SIZE(TZ,1)*SIZE(TZ,2)*SIZE(TZ,3), SIZE(TZ,1), SIZE(TZ,2), SIZE(TZ,3), SIZE(TZ,4), C_SIZEOF(TZ(LBOUND(TZ,1),LBOUND(TZ,2),LBOUND(TZ,3),LBOUND(TZ,4))))
    ELSE
      handles(101) = C_NULL_PTR
    END IF
    NULLIFY(TZ)
    IF (ASSOCIATED(QNAM_DG)) THEN
      CALL fstarpu_block_data_register(handles(102), 0, C_LOC(QNAM_DG(LBOUND(QNAM_DG,1),LBOUND(QNAM_DG,2),LBOUND(QNAM_DG,3))), SIZE(QNAM_DG,1), SIZE(QNAM_DG,1)*SIZE(QNAM_DG,2), SIZE(QNAM_DG,1), SIZE(QNAM_DG,2), SIZE(QNAM_DG,3), C_SIZEOF(QNAM_DG(LBOUND(QNAM_DG,1),LBOUND(QNAM_DG,2),LBOUND(QNAM_DG,3))))
    ELSE
      handles(102) = C_NULL_PTR
    END IF
    NULLIFY(QNAM_DG)
    IF (ASSOCIATED(QNPH_DG)) THEN
      CALL fstarpu_block_data_register(handles(103), 0, C_LOC(QNPH_DG(LBOUND(QNPH_DG,1),LBOUND(QNPH_DG,2),LBOUND(QNPH_DG,3))), SIZE(QNPH_DG,1), SIZE(QNPH_DG,1)*SIZE(QNPH_DG,2), SIZE(QNPH_DG,1), SIZE(QNPH_DG,2), SIZE(QNPH_DG,3), C_SIZEOF(QNPH_DG(LBOUND(QNPH_DG,1),LBOUND(QNPH_DG,2),LBOUND(QNPH_DG,3))))
    ELSE
      handles(103) = C_NULL_PTR
    END IF
    NULLIFY(QNPH_DG)
    IF (ASSOCIATED(RHS_ZE)) THEN
      CALL fstarpu_block_data_register(handles(104), 0, C_LOC(RHS_ZE(LBOUND(RHS_ZE,1),LBOUND(RHS_ZE,2),LBOUND(RHS_ZE,3))), SIZE(RHS_ZE,1), SIZE(RHS_ZE,1)*SIZE(RHS_ZE,2), SIZE(RHS_ZE,1), SIZE(RHS_ZE,2), SIZE(RHS_ZE,3), C_SIZEOF(RHS_ZE(LBOUND(RHS_ZE,1),LBOUND(RHS_ZE,2),LBOUND(RHS_ZE,3))))
    ELSE
      handles(104) = C_NULL_PTR
    END IF
    NULLIFY(RHS_ZE)
    IF (ASSOCIATED(RHS_bed)) THEN
      CALL fstarpu_tensor_data_register(handles(105), 0, C_LOC(RHS_bed(LBOUND(RHS_bed,1),LBOUND(RHS_bed,2),LBOUND(RHS_bed,3),LBOUND(RHS_bed,4))), SIZE(RHS_bed,1), SIZE(RHS_bed,1)*SIZE(RHS_bed,2), SIZE(RHS_bed,1)*SIZE(RHS_bed,2)*SIZE(RHS_bed,3), SIZE(RHS_bed,1), SIZE(RHS_bed,2), SIZE(RHS_bed,3), SIZE(RHS_bed,4), C_SIZEOF(RHS_bed(LBOUND(RHS_bed,1),LBOUND(RHS_bed,2),LBOUND(RHS_bed,3),LBOUND(RHS_bed,4))))
    ELSE
      handles(105) = C_NULL_PTR
    END IF
    NULLIFY(RHS_bed)
    IF (ASSOCIATED(RHS_QX)) THEN
      CALL fstarpu_block_data_register(handles(106), 0, C_LOC(RHS_QX(LBOUND(RHS_QX,1),LBOUND(RHS_QX,2),LBOUND(RHS_QX,3))), SIZE(RHS_QX,1), SIZE(RHS_QX,1)*SIZE(RHS_QX,2), SIZE(RHS_QX,1), SIZE(RHS_QX,2), SIZE(RHS_QX,3), C_SIZEOF(RHS_QX(LBOUND(RHS_QX,1),LBOUND(RHS_QX,2),LBOUND(RHS_QX,3))))
    ELSE
      handles(106) = C_NULL_PTR
    END IF
    NULLIFY(RHS_QX)
    IF (ASSOCIATED(RHS_QY)) THEN
      CALL fstarpu_block_data_register(handles(107), 0, C_LOC(RHS_QY(LBOUND(RHS_QY,1),LBOUND(RHS_QY,2),LBOUND(RHS_QY,3))), SIZE(RHS_QY,1), SIZE(RHS_QY,1)*SIZE(RHS_QY,2), SIZE(RHS_QY,1), SIZE(RHS_QY,2), SIZE(RHS_QY,3), C_SIZEOF(RHS_QY(LBOUND(RHS_QY,1),LBOUND(RHS_QY,2),LBOUND(RHS_QY,3))))
    ELSE
      handles(107) = C_NULL_PTR
    END IF
    NULLIFY(RHS_QY)
    IF (ASSOCIATED(RHS_iota)) THEN
      CALL fstarpu_block_data_register(handles(108), 0, C_LOC(RHS_iota(LBOUND(RHS_iota,1),LBOUND(RHS_iota,2),LBOUND(RHS_iota,3))), SIZE(RHS_iota,1), SIZE(RHS_iota,1)*SIZE(RHS_iota,2), SIZE(RHS_iota,1), SIZE(RHS_iota,2), SIZE(RHS_iota,3), C_SIZEOF(RHS_iota(LBOUND(RHS_iota,1),LBOUND(RHS_iota,2),LBOUND(RHS_iota,3))))
    ELSE
      handles(108) = C_NULL_PTR
    END IF
    NULLIFY(RHS_iota)
    IF (ASSOCIATED(RHS_iota2)) THEN
      CALL fstarpu_block_data_register(handles(109), 0, C_LOC(RHS_iota2(LBOUND(RHS_iota2,1),LBOUND(RHS_iota2,2),LBOUND(RHS_iota2,3))), SIZE(RHS_iota2,1), SIZE(RHS_iota2,1)*SIZE(RHS_iota2,2), SIZE(RHS_iota2,1), SIZE(RHS_iota2,2), SIZE(RHS_iota2,3), C_SIZEOF(RHS_iota2(LBOUND(RHS_iota2,1),LBOUND(RHS_iota2,2),LBOUND(RHS_iota2,3))))
    ELSE
      handles(109) = C_NULL_PTR
    END IF
    NULLIFY(RHS_iota2)
    IF (ASSOCIATED(XAGP)) THEN
      CALL fstarpu_matrix_data_register(handles(110), 0, C_LOC(XAGP(LBOUND(XAGP,1),LBOUND(XAGP,2))), SIZE(XAGP,1), SIZE(XAGP,1), SIZE(XAGP,2), C_SIZEOF(XAGP(LBOUND(XAGP,1),LBOUND(XAGP,2))))
    ELSE
      handles(110) = C_NULL_PTR
    END IF
    NULLIFY(XAGP)
    IF (ASSOCIATED(YAGP)) THEN
      CALL fstarpu_matrix_data_register(handles(111), 0, C_LOC(YAGP(LBOUND(YAGP,1),LBOUND(YAGP,2))), SIZE(YAGP,1), SIZE(YAGP,1), SIZE(YAGP,2), C_SIZEOF(YAGP(LBOUND(YAGP,1),LBOUND(YAGP,2))))
    ELSE
      handles(111) = C_NULL_PTR
    END IF
    NULLIFY(YAGP)
    IF (ASSOCIATED(WAGP)) THEN
      CALL fstarpu_matrix_data_register(handles(112), 0, C_LOC(WAGP(LBOUND(WAGP,1),LBOUND(WAGP,2))), SIZE(WAGP,1), SIZE(WAGP,1), SIZE(WAGP,2), C_SIZEOF(WAGP(LBOUND(WAGP,1),LBOUND(WAGP,2))))
    ELSE
      handles(112) = C_NULL_PTR
    END IF
    NULLIFY(WAGP)
    IF (ASSOCIATED(XEGP)) THEN
      CALL fstarpu_matrix_data_register(handles(113), 0, C_LOC(XEGP(LBOUND(XEGP,1),LBOUND(XEGP,2))), SIZE(XEGP,1), SIZE(XEGP,1), SIZE(XEGP,2), C_SIZEOF(XEGP(LBOUND(XEGP,1),LBOUND(XEGP,2))))
    ELSE
      handles(113) = C_NULL_PTR
    END IF
    NULLIFY(XEGP)
    IF (ASSOCIATED(YEGP)) THEN
      CALL fstarpu_matrix_data_register(handles(114), 0, C_LOC(YEGP(LBOUND(YEGP,1),LBOUND(YEGP,2))), SIZE(YEGP,1), SIZE(YEGP,1), SIZE(YEGP,2), C_SIZEOF(YEGP(LBOUND(YEGP,1),LBOUND(YEGP,2))))
    ELSE
      handles(114) = C_NULL_PTR
    END IF
    NULLIFY(YEGP)
    IF (ASSOCIATED(WEGP)) THEN
      CALL fstarpu_matrix_data_register(handles(115), 0, C_LOC(WEGP(LBOUND(WEGP,1),LBOUND(WEGP,2))), SIZE(WEGP,1), SIZE(WEGP,1), SIZE(WEGP,2), C_SIZEOF(WEGP(LBOUND(WEGP,1),LBOUND(WEGP,2))))
    ELSE
      handles(115) = C_NULL_PTR
    END IF
    NULLIFY(WEGP)
    IF (ASSOCIATED(SL3)) THEN
      CALL fstarpu_matrix_data_register(handles(116), 0, C_LOC(SL3(LBOUND(SL3,1),LBOUND(SL3,2))), SIZE(SL3,1), SIZE(SL3,1), SIZE(SL3,2), C_SIZEOF(SL3(LBOUND(SL3,1),LBOUND(SL3,2))))
    ELSE
      handles(116) = C_NULL_PTR
    END IF
    NULLIFY(SL3)
    IF (ASSOCIATED(XBC)) THEN
      CALL fstarpu_vector_data_register(handles(117), 0, C_LOC(XBC(LBOUND(XBC,1))), SIZE(XBC,1), C_SIZEOF(XBC(LBOUND(XBC,1))))
    ELSE
      handles(117) = C_NULL_PTR
    END IF
    NULLIFY(XBC)
    IF (ASSOCIATED(YBC)) THEN
      CALL fstarpu_vector_data_register(handles(118), 0, C_LOC(YBC(LBOUND(YBC,1))), SIZE(YBC,1), C_SIZEOF(YBC(LBOUND(YBC,1))))
    ELSE
      handles(118) = C_NULL_PTR
    END IF
    NULLIFY(YBC)
    IF (ASSOCIATED(XFAC)) THEN
      CALL fstarpu_tensor_data_register(handles(119), 0, C_LOC(XFAC(LBOUND(XFAC,1),LBOUND(XFAC,2),LBOUND(XFAC,3),LBOUND(XFAC,4))), SIZE(XFAC,1), SIZE(XFAC,1)*SIZE(XFAC,2), SIZE(XFAC,1)*SIZE(XFAC,2)*SIZE(XFAC,3), SIZE(XFAC,1), SIZE(XFAC,2), SIZE(XFAC,3), SIZE(XFAC,4), C_SIZEOF(XFAC(LBOUND(XFAC,1),LBOUND(XFAC,2),LBOUND(XFAC,3),LBOUND(XFAC,4))))
    ELSE
      handles(119) = C_NULL_PTR
    END IF
    NULLIFY(XFAC)
    IF (ASSOCIATED(YFAC)) THEN
      CALL fstarpu_tensor_data_register(handles(120), 0, C_LOC(YFAC(LBOUND(YFAC,1),LBOUND(YFAC,2),LBOUND(YFAC,3),LBOUND(YFAC,4))), SIZE(YFAC,1), SIZE(YFAC,1)*SIZE(YFAC,2), SIZE(YFAC,1)*SIZE(YFAC,2)*SIZE(YFAC,3), SIZE(YFAC,1), SIZE(YFAC,2), SIZE(YFAC,3), SIZE(YFAC,4), C_SIZEOF(YFAC(LBOUND(YFAC,1),LBOUND(YFAC,2),LBOUND(YFAC,3),LBOUND(YFAC,4))))
    ELSE
      handles(120) = C_NULL_PTR
    END IF
    NULLIFY(YFAC)
    IF (ASSOCIATED(EDGEQ)) THEN
      CALL fstarpu_tensor_data_register(handles(121), 0, C_LOC(EDGEQ(LBOUND(EDGEQ,1),LBOUND(EDGEQ,2),LBOUND(EDGEQ,3),LBOUND(EDGEQ,4))), SIZE(EDGEQ,1), SIZE(EDGEQ,1)*SIZE(EDGEQ,2), SIZE(EDGEQ,1)*SIZE(EDGEQ,2)*SIZE(EDGEQ,3), SIZE(EDGEQ,1), SIZE(EDGEQ,2), SIZE(EDGEQ,3), SIZE(EDGEQ,4), C_SIZEOF(EDGEQ(LBOUND(EDGEQ,1),LBOUND(EDGEQ,2),LBOUND(EDGEQ,3),LBOUND(EDGEQ,4))))
    ELSE
      handles(121) = C_NULL_PTR
    END IF
    NULLIFY(EDGEQ)
    IF (ASSOCIATED(bed_IN)) THEN
      CALL fstarpu_vector_data_register(handles(122), 0, C_LOC(bed_IN(LBOUND(bed_IN,1))), SIZE(bed_IN,1), C_SIZEOF(bed_IN(LBOUND(bed_IN,1))))
    ELSE
      handles(122) = C_NULL_PTR
    END IF
    NULLIFY(bed_IN)
    IF (ASSOCIATED(bed_EX)) THEN
      CALL fstarpu_vector_data_register(handles(123), 0, C_LOC(bed_EX(LBOUND(bed_EX,1))), SIZE(bed_EX,1), C_SIZEOF(bed_EX(LBOUND(bed_EX,1))))
    ELSE
      handles(123) = C_NULL_PTR
    END IF
    NULLIFY(bed_EX)
    IF (ASSOCIATED(bed_HAT)) THEN
      CALL fstarpu_vector_data_register(handles(124), 0, C_LOC(bed_HAT(LBOUND(bed_HAT,1))), SIZE(bed_HAT,1), C_SIZEOF(bed_HAT(LBOUND(bed_HAT,1))))
    ELSE
      handles(124) = C_NULL_PTR
    END IF
    NULLIFY(bed_HAT)
    IF (ASSOCIATED(fact)) THEN
      CALL fstarpu_vector_data_register(handles(125), 0, C_LOC(fact(LBOUND(fact,1))), SIZE(fact,1), C_SIZEOF(fact(LBOUND(fact,1))))
    ELSE
      handles(125) = C_NULL_PTR
    END IF
    NULLIFY(fact)
    IF (ASSOCIATED(focal_neigh)) THEN
      CALL fstarpu_matrix_data_register(handles(126), 0, C_LOC(focal_neigh(LBOUND(focal_neigh,1),LBOUND(focal_neigh,2))), SIZE(focal_neigh,1), SIZE(focal_neigh,1), SIZE(focal_neigh,2), C_SIZEOF(focal_neigh(LBOUND(focal_neigh,1),LBOUND(focal_neigh,2))))
    ELSE
      handles(126) = C_NULL_PTR
    END IF
    NULLIFY(focal_neigh)
    IF (ASSOCIATED(focal_up)) THEN
      CALL fstarpu_vector_data_register(handles(127), 0, C_LOC(focal_up(LBOUND(focal_up,1))), SIZE(focal_up,1), C_SIZEOF(focal_up(LBOUND(focal_up,1))))
    ELSE
      handles(127) = C_NULL_PTR
    END IF
    NULLIFY(focal_up)
    IF (ASSOCIATED(bi)) THEN
      CALL fstarpu_vector_data_register(handles(128), 0, C_LOC(bi(LBOUND(bi,1))), SIZE(bi,1), C_SIZEOF(bi(LBOUND(bi,1))))
    ELSE
      handles(128) = C_NULL_PTR
    END IF
    NULLIFY(bi)
    IF (ASSOCIATED(XBCb)) THEN
      CALL fstarpu_vector_data_register(handles(129), 0, C_LOC(XBCb(LBOUND(XBCb,1))), SIZE(XBCb,1), C_SIZEOF(XBCb(LBOUND(XBCb,1))))
    ELSE
      handles(129) = C_NULL_PTR
    END IF
    NULLIFY(XBCb)
    IF (ASSOCIATED(YBCb)) THEN
      CALL fstarpu_vector_data_register(handles(130), 0, C_LOC(YBCb(LBOUND(YBCb,1))), SIZE(YBCb,1), C_SIZEOF(YBCb(LBOUND(YBCb,1))))
    ELSE
      handles(130) = C_NULL_PTR
    END IF
    NULLIFY(YBCb)
    IF (ASSOCIATED(xi1)) THEN
      CALL fstarpu_matrix_data_register(handles(131), 0, C_LOC(xi1(LBOUND(xi1,1),LBOUND(xi1,2))), SIZE(xi1,1), SIZE(xi1,1), SIZE(xi1,2), C_SIZEOF(xi1(LBOUND(xi1,1),LBOUND(xi1,2))))
    ELSE
      handles(131) = C_NULL_PTR
    END IF
    NULLIFY(xi1)
    IF (ASSOCIATED(xi2)) THEN
      CALL fstarpu_matrix_data_register(handles(132), 0, C_LOC(xi2(LBOUND(xi2,1),LBOUND(xi2,2))), SIZE(xi2,1), SIZE(xi2,1), SIZE(xi2,2), C_SIZEOF(xi2(LBOUND(xi2,1),LBOUND(xi2,2))))
    ELSE
      handles(132) = C_NULL_PTR
    END IF
    NULLIFY(xi2)
    IF (ASSOCIATED(xtransform)) THEN
      CALL fstarpu_matrix_data_register(handles(133), 0, C_LOC(xtransform(LBOUND(xtransform,1),LBOUND(xtransform,2))), SIZE(xtransform,1), SIZE(xtransform,1), SIZE(xtransform,2), C_SIZEOF(xtransform(LBOUND(xtransform,1),LBOUND(xtransform,2))))
    ELSE
      handles(133) = C_NULL_PTR
    END IF
    NULLIFY(xtransform)
    IF (ASSOCIATED(ytransform)) THEN
      CALL fstarpu_matrix_data_register(handles(134), 0, C_LOC(ytransform(LBOUND(ytransform,1),LBOUND(ytransform,2))), SIZE(ytransform,1), SIZE(ytransform,1), SIZE(ytransform,2), C_SIZEOF(ytransform(LBOUND(ytransform,1),LBOUND(ytransform,2))))
    ELSE
      handles(134) = C_NULL_PTR
    END IF
    NULLIFY(ytransform)
    IF (ASSOCIATED(xi1BCb)) THEN
      CALL fstarpu_vector_data_register(handles(135), 0, C_LOC(xi1BCb(LBOUND(xi1BCb,1))), SIZE(xi1BCb,1), C_SIZEOF(xi1BCb(LBOUND(xi1BCb,1))))
    ELSE
      handles(135) = C_NULL_PTR
    END IF
    NULLIFY(xi1BCb)
    IF (ASSOCIATED(xi2BCb)) THEN
      CALL fstarpu_vector_data_register(handles(136), 0, C_LOC(xi2BCb(LBOUND(xi2BCb,1))), SIZE(xi2BCb,1), C_SIZEOF(xi2BCb(LBOUND(xi2BCb,1))))
    ELSE
      handles(136) = C_NULL_PTR
    END IF
    NULLIFY(xi2BCb)
    IF (ASSOCIATED(xi1vert)) THEN
      CALL fstarpu_matrix_data_register(handles(137), 0, C_LOC(xi1vert(LBOUND(xi1vert,1),LBOUND(xi1vert,2))), SIZE(xi1vert,1), SIZE(xi1vert,1), SIZE(xi1vert,2), C_SIZEOF(xi1vert(LBOUND(xi1vert,1),LBOUND(xi1vert,2))))
    ELSE
      handles(137) = C_NULL_PTR
    END IF
    NULLIFY(xi1vert)
    IF (ASSOCIATED(xi2vert)) THEN
      CALL fstarpu_matrix_data_register(handles(138), 0, C_LOC(xi2vert(LBOUND(xi2vert,1),LBOUND(xi2vert,2))), SIZE(xi2vert,1), SIZE(xi2vert,1), SIZE(xi2vert,2), C_SIZEOF(xi2vert(LBOUND(xi2vert,1),LBOUND(xi2vert,2))))
    ELSE
      handles(138) = C_NULL_PTR
    END IF
    NULLIFY(xi2vert)
    IF (ASSOCIATED(xtransformv)) THEN
      CALL fstarpu_matrix_data_register(handles(139), 0, C_LOC(xtransformv(LBOUND(xtransformv,1),LBOUND(xtransformv,2))), SIZE(xtransformv,1), SIZE(xtransformv,1), SIZE(xtransformv,2), C_SIZEOF(xtransformv(LBOUND(xtransformv,1),LBOUND(xtransformv,2))))
    ELSE
      handles(139) = C_NULL_PTR
    END IF
    NULLIFY(xtransformv)
    IF (ASSOCIATED(XBCv)) THEN
      CALL fstarpu_matrix_data_register(handles(140), 0, C_LOC(XBCv(LBOUND(XBCv,1),LBOUND(XBCv,2))), SIZE(XBCv,1), SIZE(XBCv,1), SIZE(XBCv,2), C_SIZEOF(XBCv(LBOUND(XBCv,1),LBOUND(XBCv,2))))
    ELSE
      handles(140) = C_NULL_PTR
    END IF
    NULLIFY(XBCv)
    IF (ASSOCIATED(YBCv)) THEN
      CALL fstarpu_matrix_data_register(handles(141), 0, C_LOC(YBCv(LBOUND(YBCv,1),LBOUND(YBCv,2))), SIZE(YBCv,1), SIZE(YBCv,1), SIZE(YBCv,2), C_SIZEOF(YBCv(LBOUND(YBCv,1),LBOUND(YBCv,2))))
    ELSE
      handles(141) = C_NULL_PTR
    END IF
    NULLIFY(YBCv)
    IF (ASSOCIATED(xi1BCv)) THEN
      CALL fstarpu_matrix_data_register(handles(142), 0, C_LOC(xi1BCv(LBOUND(xi1BCv,1),LBOUND(xi1BCv,2))), SIZE(xi1BCv,1), SIZE(xi1BCv,1), SIZE(xi1BCv,2), C_SIZEOF(xi1BCv(LBOUND(xi1BCv,1),LBOUND(xi1BCv,2))))
    ELSE
      handles(142) = C_NULL_PTR
    END IF
    NULLIFY(xi1BCv)
    IF (ASSOCIATED(xi2BCv)) THEN
      CALL fstarpu_matrix_data_register(handles(143), 0, C_LOC(xi2BCv(LBOUND(xi2BCv,1),LBOUND(xi2BCv,2))), SIZE(xi2BCv,1), SIZE(xi2BCv,1), SIZE(xi2BCv,2), C_SIZEOF(xi2BCv(LBOUND(xi2BCv,1),LBOUND(xi2BCv,2))))
    ELSE
      handles(143) = C_NULL_PTR
    END IF
    NULLIFY(xi2BCv)
    IF (ASSOCIATED(Area_integral)) THEN
      CALL fstarpu_block_data_register(handles(144), 0, C_LOC(Area_integral(LBOUND(Area_integral,1),LBOUND(Area_integral,2),LBOUND(Area_integral,3))), SIZE(Area_integral,1), SIZE(Area_integral,1)*SIZE(Area_integral,2), SIZE(Area_integral,1), SIZE(Area_integral,2), SIZE(Area_integral,3), C_SIZEOF(Area_integral(LBOUND(Area_integral,1),LBOUND(Area_integral,2),LBOUND(Area_integral,3))))
    ELSE
      handles(144) = C_NULL_PTR
    END IF
    NULLIFY(Area_integral)
    IF (ASSOCIATED(f)) THEN
      CALL fstarpu_tensor_data_register(handles(145), 0, C_LOC(f(LBOUND(f,1),LBOUND(f,2),LBOUND(f,3),LBOUND(f,4))), SIZE(f,1), SIZE(f,1)*SIZE(f,2), SIZE(f,1)*SIZE(f,2)*SIZE(f,3), SIZE(f,1), SIZE(f,2), SIZE(f,3), SIZE(f,4), C_SIZEOF(f(LBOUND(f,1),LBOUND(f,2),LBOUND(f,3),LBOUND(f,4))))
    ELSE
      handles(145) = C_NULL_PTR
    END IF
    NULLIFY(f)
    IF (ASSOCIATED(g0)) THEN
      CALL fstarpu_tensor_data_register(handles(146), 0, C_LOC(g0(LBOUND(g0,1),LBOUND(g0,2),LBOUND(g0,3),LBOUND(g0,4))), SIZE(g0,1), SIZE(g0,1)*SIZE(g0,2), SIZE(g0,1)*SIZE(g0,2)*SIZE(g0,3), SIZE(g0,1), SIZE(g0,2), SIZE(g0,3), SIZE(g0,4), C_SIZEOF(g0(LBOUND(g0,1),LBOUND(g0,2),LBOUND(g0,3),LBOUND(g0,4))))
    ELSE
      handles(146) = C_NULL_PTR
    END IF
    NULLIFY(g0)
    IF (ASSOCIATED(varsigma0)) THEN
      CALL fstarpu_tensor_data_register(handles(147), 0, C_LOC(varsigma0(LBOUND(varsigma0,1),LBOUND(varsigma0,2),LBOUND(varsigma0,3),LBOUND(varsigma0,4))), SIZE(varsigma0,1), SIZE(varsigma0,1)*SIZE(varsigma0,2), SIZE(varsigma0,1)*SIZE(varsigma0,2)*SIZE(varsigma0,3), SIZE(varsigma0,1), SIZE(varsigma0,2), SIZE(varsigma0,3), SIZE(varsigma0,4), C_SIZEOF(varsigma0(LBOUND(varsigma0,1),LBOUND(varsigma0,2),LBOUND(varsigma0,3),LBOUND(varsigma0,4))))
    ELSE
      handles(147) = C_NULL_PTR
    END IF
    NULLIFY(varsigma0)
    IF (ASSOCIATED(fv)) THEN
      CALL fstarpu_tensor_data_register(handles(148), 0, C_LOC(fv(LBOUND(fv,1),LBOUND(fv,2),LBOUND(fv,3),LBOUND(fv,4))), SIZE(fv,1), SIZE(fv,1)*SIZE(fv,2), SIZE(fv,1)*SIZE(fv,2)*SIZE(fv,3), SIZE(fv,1), SIZE(fv,2), SIZE(fv,3), SIZE(fv,4), C_SIZEOF(fv(LBOUND(fv,1),LBOUND(fv,2),LBOUND(fv,3),LBOUND(fv,4))))
    ELSE
      handles(148) = C_NULL_PTR
    END IF
    NULLIFY(fv)
    IF (ASSOCIATED(g0v)) THEN
      CALL fstarpu_tensor_data_register(handles(149), 0, C_LOC(g0v(LBOUND(g0v,1),LBOUND(g0v,2),LBOUND(g0v,3),LBOUND(g0v,4))), SIZE(g0v,1), SIZE(g0v,1)*SIZE(g0v,2), SIZE(g0v,1)*SIZE(g0v,2)*SIZE(g0v,3), SIZE(g0v,1), SIZE(g0v,2), SIZE(g0v,3), SIZE(g0v,4), C_SIZEOF(g0v(LBOUND(g0v,1),LBOUND(g0v,2),LBOUND(g0v,3),LBOUND(g0v,4))))
    ELSE
      handles(149) = C_NULL_PTR
    END IF
    NULLIFY(g0v)
    IF (ASSOCIATED(var2sigmag)) THEN
      CALL fstarpu_block_data_register(handles(150), 0, C_LOC(var2sigmag(LBOUND(var2sigmag,1),LBOUND(var2sigmag,2),LBOUND(var2sigmag,3))), SIZE(var2sigmag,1), SIZE(var2sigmag,1)*SIZE(var2sigmag,2), SIZE(var2sigmag,1), SIZE(var2sigmag,2), SIZE(var2sigmag,3), C_SIZEOF(var2sigmag(LBOUND(var2sigmag,1),LBOUND(var2sigmag,2),LBOUND(var2sigmag,3))))
    ELSE
      handles(150) = C_NULL_PTR
    END IF
    NULLIFY(var2sigmag)
    IF (ASSOCIATED(var2sigmav)) THEN
      CALL fstarpu_block_data_register(handles(151), 0, C_LOC(var2sigmav(LBOUND(var2sigmav,1),LBOUND(var2sigmav,2),LBOUND(var2sigmav,3))), SIZE(var2sigmav,1), SIZE(var2sigmav,1)*SIZE(var2sigmav,2), SIZE(var2sigmav,1), SIZE(var2sigmav,2), SIZE(var2sigmav,3), C_SIZEOF(var2sigmav(LBOUND(var2sigmav,1),LBOUND(var2sigmav,2),LBOUND(var2sigmav,3))))
    ELSE
      handles(151) = C_NULL_PTR
    END IF
    NULLIFY(var2sigmav)
    IF (ASSOCIATED(Nmatrix)) THEN
      CALL fstarpu_tensor_data_register(handles(152), 0, C_LOC(Nmatrix(LBOUND(Nmatrix,1),LBOUND(Nmatrix,2),LBOUND(Nmatrix,3),LBOUND(Nmatrix,4))), SIZE(Nmatrix,1), SIZE(Nmatrix,1)*SIZE(Nmatrix,2), SIZE(Nmatrix,1)*SIZE(Nmatrix,2)*SIZE(Nmatrix,3), SIZE(Nmatrix,1), SIZE(Nmatrix,2), SIZE(Nmatrix,3), SIZE(Nmatrix,4), C_SIZEOF(Nmatrix(LBOUND(Nmatrix,1),LBOUND(Nmatrix,2),LBOUND(Nmatrix,3),LBOUND(Nmatrix,4))))
    ELSE
      handles(152) = C_NULL_PTR
    END IF
    NULLIFY(Nmatrix)
    IF (ASSOCIATED(NmatrixInv)) THEN
      CALL fstarpu_tensor_data_register(handles(153), 0, C_LOC(NmatrixInv(LBOUND(NmatrixInv,1),LBOUND(NmatrixInv,2),LBOUND(NmatrixInv,3),LBOUND(NmatrixInv,4))), SIZE(NmatrixInv,1), SIZE(NmatrixInv,1)*SIZE(NmatrixInv,2), SIZE(NmatrixInv,1)*SIZE(NmatrixInv,2)*SIZE(NmatrixInv,3), SIZE(NmatrixInv,1), SIZE(NmatrixInv,2), SIZE(NmatrixInv,3), SIZE(NmatrixInv,4), C_SIZEOF(NmatrixInv(LBOUND(NmatrixInv,1),LBOUND(NmatrixInv,2),LBOUND(NmatrixInv,3),LBOUND(NmatrixInv,4))))
    ELSE
      handles(153) = C_NULL_PTR
    END IF
    NULLIFY(NmatrixInv)
    IF (ASSOCIATED(deltx)) THEN
      CALL fstarpu_vector_data_register(handles(154), 0, C_LOC(deltx(LBOUND(deltx,1))), SIZE(deltx,1), C_SIZEOF(deltx(LBOUND(deltx,1))))
    ELSE
      handles(154) = C_NULL_PTR
    END IF
    NULLIFY(deltx)
    IF (ASSOCIATED(delty)) THEN
      CALL fstarpu_vector_data_register(handles(155), 0, C_LOC(delty(LBOUND(delty,1))), SIZE(delty,1), C_SIZEOF(delty(LBOUND(delty,1))))
    ELSE
      handles(155) = C_NULL_PTR
    END IF
    NULLIFY(delty)
    IF (ASSOCIATED(pmatrix)) THEN
      CALL fstarpu_block_data_register(handles(156), 0, C_LOC(pmatrix(LBOUND(pmatrix,1),LBOUND(pmatrix,2),LBOUND(pmatrix,3))), SIZE(pmatrix,1), SIZE(pmatrix,1)*SIZE(pmatrix,2), SIZE(pmatrix,1), SIZE(pmatrix,2), SIZE(pmatrix,3), C_SIZEOF(pmatrix(LBOUND(pmatrix,1),LBOUND(pmatrix,2),LBOUND(pmatrix,3))))
    ELSE
      handles(156) = C_NULL_PTR
    END IF
    NULLIFY(pmatrix)
    IF (ASSOCIATED(ZEmin)) THEN
      CALL fstarpu_matrix_data_register(handles(157), 0, C_LOC(ZEmin(LBOUND(ZEmin,1),LBOUND(ZEmin,2))), SIZE(ZEmin,1), SIZE(ZEmin,1), SIZE(ZEmin,2), C_SIZEOF(ZEmin(LBOUND(ZEmin,1),LBOUND(ZEmin,2))))
    ELSE
      handles(157) = C_NULL_PTR
    END IF
    NULLIFY(ZEmin)
    IF (ASSOCIATED(ZEmax)) THEN
      CALL fstarpu_matrix_data_register(handles(158), 0, C_LOC(ZEmax(LBOUND(ZEmax,1),LBOUND(ZEmax,2))), SIZE(ZEmax,1), SIZE(ZEmax,1), SIZE(ZEmax,2), C_SIZEOF(ZEmax(LBOUND(ZEmax,1),LBOUND(ZEmax,2))))
    ELSE
      handles(158) = C_NULL_PTR
    END IF
    NULLIFY(ZEmax)
    IF (ASSOCIATED(QXmin)) THEN
      CALL fstarpu_matrix_data_register(handles(159), 0, C_LOC(QXmin(LBOUND(QXmin,1),LBOUND(QXmin,2))), SIZE(QXmin,1), SIZE(QXmin,1), SIZE(QXmin,2), C_SIZEOF(QXmin(LBOUND(QXmin,1),LBOUND(QXmin,2))))
    ELSE
      handles(159) = C_NULL_PTR
    END IF
    NULLIFY(QXmin)
    IF (ASSOCIATED(QXmax)) THEN
      CALL fstarpu_matrix_data_register(handles(160), 0, C_LOC(QXmax(LBOUND(QXmax,1),LBOUND(QXmax,2))), SIZE(QXmax,1), SIZE(QXmax,1), SIZE(QXmax,2), C_SIZEOF(QXmax(LBOUND(QXmax,1),LBOUND(QXmax,2))))
    ELSE
      handles(160) = C_NULL_PTR
    END IF
    NULLIFY(QXmax)
    IF (ASSOCIATED(QYmin)) THEN
      CALL fstarpu_matrix_data_register(handles(161), 0, C_LOC(QYmin(LBOUND(QYmin,1),LBOUND(QYmin,2))), SIZE(QYmin,1), SIZE(QYmin,1), SIZE(QYmin,2), C_SIZEOF(QYmin(LBOUND(QYmin,1),LBOUND(QYmin,2))))
    ELSE
      handles(161) = C_NULL_PTR
    END IF
    NULLIFY(QYmin)
    IF (ASSOCIATED(QYmax)) THEN
      CALL fstarpu_matrix_data_register(handles(162), 0, C_LOC(QYmax(LBOUND(QYmax,1),LBOUND(QYmax,2))), SIZE(QYmax,1), SIZE(QYmax,1), SIZE(QYmax,2), C_SIZEOF(QYmax(LBOUND(QYmax,1),LBOUND(QYmax,2))))
    ELSE
      handles(162) = C_NULL_PTR
    END IF
    NULLIFY(QYmax)
    IF (ASSOCIATED(iotamin)) THEN
      CALL fstarpu_matrix_data_register(handles(163), 0, C_LOC(iotamin(LBOUND(iotamin,1),LBOUND(iotamin,2))), SIZE(iotamin,1), SIZE(iotamin,1), SIZE(iotamin,2), C_SIZEOF(iotamin(LBOUND(iotamin,1),LBOUND(iotamin,2))))
    ELSE
      handles(163) = C_NULL_PTR
    END IF
    NULLIFY(iotamin)
    IF (ASSOCIATED(iotamax)) THEN
      CALL fstarpu_matrix_data_register(handles(164), 0, C_LOC(iotamax(LBOUND(iotamax,1),LBOUND(iotamax,2))), SIZE(iotamax,1), SIZE(iotamax,1), SIZE(iotamax,2), C_SIZEOF(iotamax(LBOUND(iotamax,1),LBOUND(iotamax,2))))
    ELSE
      handles(164) = C_NULL_PTR
    END IF
    NULLIFY(iotamax)
    IF (ASSOCIATED(iota2min)) THEN
      CALL fstarpu_matrix_data_register(handles(165), 0, C_LOC(iota2min(LBOUND(iota2min,1),LBOUND(iota2min,2))), SIZE(iota2min,1), SIZE(iota2min,1), SIZE(iota2min,2), C_SIZEOF(iota2min(LBOUND(iota2min,1),LBOUND(iota2min,2))))
    ELSE
      handles(165) = C_NULL_PTR
    END IF
    NULLIFY(iota2min)
    IF (ASSOCIATED(iota2max)) THEN
      CALL fstarpu_matrix_data_register(handles(166), 0, C_LOC(iota2max(LBOUND(iota2max,1),LBOUND(iota2max,2))), SIZE(iota2max,1), SIZE(iota2max,1), SIZE(iota2max,2), C_SIZEOF(iota2max(LBOUND(iota2max,1),LBOUND(iota2max,2))))
    ELSE
      handles(166) = C_NULL_PTR
    END IF
    NULLIFY(iota2max)
    IF (ASSOCIATED(ANGTAB)) THEN
      CALL fstarpu_matrix_data_register(handles(167), 0, C_LOC(ANGTAB(LBOUND(ANGTAB,1),LBOUND(ANGTAB,2))), SIZE(ANGTAB,1), SIZE(ANGTAB,1), SIZE(ANGTAB,2), C_SIZEOF(ANGTAB(LBOUND(ANGTAB,1),LBOUND(ANGTAB,2))))
    ELSE
      handles(167) = C_NULL_PTR
    END IF
    NULLIFY(ANGTAB)
    IF (ASSOCIATED(CENTAB)) THEN
      CALL fstarpu_matrix_data_register(handles(168), 0, C_LOC(CENTAB(LBOUND(CENTAB,1),LBOUND(CENTAB,2))), SIZE(CENTAB,1), SIZE(CENTAB,1), SIZE(CENTAB,2), C_SIZEOF(CENTAB(LBOUND(CENTAB,1),LBOUND(CENTAB,2))))
    ELSE
      handles(168) = C_NULL_PTR
    END IF
    NULLIFY(CENTAB)
    IF (ASSOCIATED(ELETAB)) THEN
      CALL fstarpu_matrix_data_register(handles(169), 0, C_LOC(ELETAB(LBOUND(ELETAB,1),LBOUND(ELETAB,2))), SIZE(ELETAB,1), SIZE(ELETAB,1), SIZE(ELETAB,2), C_SIZEOF(ELETAB(LBOUND(ELETAB,1),LBOUND(ELETAB,2))))
    ELSE
      handles(169) = C_NULL_PTR
    END IF
    NULLIFY(ELETAB)
    IF (ASSOCIATED(EL_COUNT)) THEN
      CALL fstarpu_vector_data_register(handles(170), 0, C_LOC(EL_COUNT(LBOUND(EL_COUNT,1))), SIZE(EL_COUNT,1), C_SIZEOF(EL_COUNT(LBOUND(EL_COUNT,1))))
    ELSE
      handles(170) = C_NULL_PTR
    END IF
    NULLIFY(EL_COUNT)
    IF (ASSOCIATED(NNDEL)) THEN
      CALL fstarpu_vector_data_register(handles(171), 0, C_LOC(NNDEL(LBOUND(NNDEL,1))), SIZE(NNDEL,1), C_SIZEOF(NNDEL(LBOUND(NNDEL,1))))
    ELSE
      handles(171) = C_NULL_PTR
    END IF
    NULLIFY(NNDEL)
    IF (ASSOCIATED(NDEL)) THEN
      CALL fstarpu_matrix_data_register(handles(172), 0, C_LOC(NDEL(LBOUND(NDEL,1),LBOUND(NDEL,2))), SIZE(NDEL,1), SIZE(NDEL,1), SIZE(NDEL,2), C_SIZEOF(NDEL(LBOUND(NDEL,1),LBOUND(NDEL,2))))
    ELSE
      handles(172) = C_NULL_PTR
    END IF
    NULLIFY(NDEL)
    IF (ASSOCIATED(iota2_DG)) THEN
      CALL fstarpu_vector_data_register(handles(173), 0, C_LOC(iota2_DG(LBOUND(iota2_DG,1))), SIZE(iota2_DG,1), C_SIZEOF(iota2_DG(LBOUND(iota2_DG,1))))
    ELSE
      handles(173) = C_NULL_PTR
    END IF
    NULLIFY(iota2_DG)
    IF (ASSOCIATED(iota_DG)) THEN
      CALL fstarpu_vector_data_register(handles(174), 0, C_LOC(iota_DG(LBOUND(iota_DG,1))), SIZE(iota_DG,1), C_SIZEOF(iota_DG(LBOUND(iota_DG,1))))
    ELSE
      handles(174) = C_NULL_PTR
    END IF
    NULLIFY(iota_DG)
    IF (ASSOCIATED(iotaa_DG)) THEN
      CALL fstarpu_vector_data_register(handles(175), 0, C_LOC(iotaa_DG(LBOUND(iotaa_DG,1))), SIZE(iotaa_DG,1), C_SIZEOF(iotaa_DG(LBOUND(iotaa_DG,1))))
    ELSE
      handles(175) = C_NULL_PTR
    END IF
    NULLIFY(iotaa_DG)
    IF (ASSOCIATED(bed_DG)) THEN
      CALL fstarpu_matrix_data_register(handles(176), 0, C_LOC(bed_DG(LBOUND(bed_DG,1),LBOUND(bed_DG,2))), SIZE(bed_DG,1), SIZE(bed_DG,1), SIZE(bed_DG,2), C_SIZEOF(bed_DG(LBOUND(bed_DG,1),LBOUND(bed_DG,2))))
    ELSE
      handles(176) = C_NULL_PTR
    END IF
    NULLIFY(bed_DG)
    IF (ASSOCIATED(bed_N_int)) THEN
      CALL fstarpu_vector_data_register(handles(177), 0, C_LOC(bed_N_int(LBOUND(bed_N_int,1))), SIZE(bed_N_int,1), C_SIZEOF(bed_N_int(LBOUND(bed_N_int,1))))
    ELSE
      handles(177) = C_NULL_PTR
    END IF
    NULLIFY(bed_N_int)
    IF (ASSOCIATED(bed_N_ext)) THEN
      CALL fstarpu_vector_data_register(handles(178), 0, C_LOC(bed_N_ext(LBOUND(bed_N_ext,1))), SIZE(bed_N_ext,1), C_SIZEOF(bed_N_ext(LBOUND(bed_N_ext,1))))
    ELSE
      handles(178) = C_NULL_PTR
    END IF
    NULLIFY(bed_N_ext)
    IF (ASSOCIATED(pdg_el)) THEN
      CALL fstarpu_vector_data_register(handles(179), 0, C_LOC(pdg_el(LBOUND(pdg_el,1))), SIZE(pdg_el,1), C_SIZEOF(pdg_el(LBOUND(pdg_el,1))))
    ELSE
      handles(179) = C_NULL_PTR
    END IF
    NULLIFY(pdg_el)
    IF (ASSOCIATED(ETAS)) THEN
      CALL fstarpu_vector_data_register(handles(180), 0, C_LOC(ETAS(LBOUND(ETAS,1))), SIZE(ETAS,1), C_SIZEOF(ETAS(LBOUND(ETAS,1))))
    ELSE
      handles(180) = C_NULL_PTR
    END IF
    NULLIFY(ETAS)
    IF (ASSOCIATED(ETA1)) THEN
      CALL fstarpu_vector_data_register(handles(181), 0, C_LOC(ETA1(LBOUND(ETA1,1))), SIZE(ETA1,1), C_SIZEOF(ETA1(LBOUND(ETA1,1))))
    ELSE
      handles(181) = C_NULL_PTR
    END IF
    NULLIFY(ETA1)
    IF (ASSOCIATED(ETA2)) THEN
      CALL fstarpu_vector_data_register(handles(182), 0, C_LOC(ETA2(LBOUND(ETA2,1))), SIZE(ETA2,1), C_SIZEOF(ETA2(LBOUND(ETA2,1))))
    ELSE
      handles(182) = C_NULL_PTR
    END IF
    NULLIFY(ETA2)
    IF (ASSOCIATED(ETAMAX)) THEN
      CALL fstarpu_vector_data_register(handles(183), 0, C_LOC(ETAMAX(LBOUND(ETAMAX,1))), SIZE(ETAMAX,1), C_SIZEOF(ETAMAX(LBOUND(ETAMAX,1))))
    ELSE
      handles(183) = C_NULL_PTR
    END IF
    NULLIFY(ETAMAX)
    IF (ASSOCIATED(entrop)) THEN
      CALL fstarpu_matrix_data_register(handles(184), 0, C_LOC(entrop(LBOUND(entrop,1),LBOUND(entrop,2))), SIZE(entrop,1), SIZE(entrop,1), SIZE(entrop,2), C_SIZEOF(entrop(LBOUND(entrop,1),LBOUND(entrop,2))))
    ELSE
      handles(184) = C_NULL_PTR
    END IF
    NULLIFY(entrop)
    IF (ASSOCIATED(tracer)) THEN
      CALL fstarpu_vector_data_register(handles(185), 0, C_LOC(tracer(LBOUND(tracer,1))), SIZE(tracer,1), C_SIZEOF(tracer(LBOUND(tracer,1))))
    ELSE
      handles(185) = C_NULL_PTR
    END IF
    NULLIFY(tracer)
    IF (ASSOCIATED(tracer2)) THEN
      CALL fstarpu_vector_data_register(handles(186), 0, C_LOC(tracer2(LBOUND(tracer2,1))), SIZE(tracer2,1), C_SIZEOF(tracer2(LBOUND(tracer2,1))))
    ELSE
      handles(186) = C_NULL_PTR
    END IF
    NULLIFY(tracer2)
    IF (ASSOCIATED(MassMax)) THEN
      CALL fstarpu_vector_data_register(handles(187), 0, C_LOC(MassMax(LBOUND(MassMax,1))), SIZE(MassMax,1), C_SIZEOF(MassMax(LBOUND(MassMax,1))))
    ELSE
      handles(187) = C_NULL_PTR
    END IF
    NULLIFY(MassMax)
    IF (ASSOCIATED(bed_int)) THEN
      CALL fstarpu_matrix_data_register(handles(188), 0, C_LOC(bed_int(LBOUND(bed_int,1),LBOUND(bed_int,2))), SIZE(bed_int,1), SIZE(bed_int,1), SIZE(bed_int,2), C_SIZEOF(bed_int(LBOUND(bed_int,1),LBOUND(bed_int,2))))
    ELSE
      handles(188) = C_NULL_PTR
    END IF
    NULLIFY(bed_int)
    IF (ASSOCIATED(DP)) THEN
      CALL fstarpu_vector_data_register(handles(189), 0, C_LOC(DP(LBOUND(DP,1))), SIZE(DP,1), C_SIZEOF(DP(LBOUND(DP,1))))
    ELSE
      handles(189) = C_NULL_PTR
    END IF
    NULLIFY(DP)
    IF (ASSOCIATED(DP0)) THEN
      CALL fstarpu_vector_data_register(handles(190), 0, C_LOC(DP0(LBOUND(DP0,1))), SIZE(DP0,1), C_SIZEOF(DP0(LBOUND(DP0,1))))
    ELSE
      handles(190) = C_NULL_PTR
    END IF
    NULLIFY(DP0)
    IF (ASSOCIATED(DPe)) THEN
      CALL fstarpu_vector_data_register(handles(191), 0, C_LOC(DPe(LBOUND(DPe,1))), SIZE(DPe,1), C_SIZEOF(DPe(LBOUND(DPe,1))))
    ELSE
      handles(191) = C_NULL_PTR
    END IF
    NULLIFY(DPe)
    IF (ASSOCIATED(SFAC)) THEN
      CALL fstarpu_vector_data_register(handles(192), 0, C_LOC(SFAC(LBOUND(SFAC,1))), SIZE(SFAC,1), C_SIZEOF(SFAC(LBOUND(SFAC,1))))
    ELSE
      handles(192) = C_NULL_PTR
    END IF
    NULLIFY(SFAC)
    IF (ASSOCIATED(CORIF)) THEN
      CALL fstarpu_vector_data_register(handles(193), 0, C_LOC(CORIF(LBOUND(CORIF,1))), SIZE(CORIF,1), C_SIZEOF(CORIF(LBOUND(CORIF,1))))
    ELSE
      handles(193) = C_NULL_PTR
    END IF
    NULLIFY(CORIF)
    IF (ASSOCIATED(ESBIN1)) THEN
      CALL fstarpu_vector_data_register(handles(194), 0, C_LOC(ESBIN1(LBOUND(ESBIN1,1))), SIZE(ESBIN1,1), C_SIZEOF(ESBIN1(LBOUND(ESBIN1,1))))
    ELSE
      handles(194) = C_NULL_PTR
    END IF
    NULLIFY(ESBIN1)
    IF (ASSOCIATED(ESBIN2)) THEN
      CALL fstarpu_vector_data_register(handles(195), 0, C_LOC(ESBIN2(LBOUND(ESBIN2,1))), SIZE(ESBIN2,1), C_SIZEOF(ESBIN2(LBOUND(ESBIN2,1))))
    ELSE
      handles(195) = C_NULL_PTR
    END IF
    NULLIFY(ESBIN2)
    IF (ASSOCIATED(QNIN1)) THEN
      CALL fstarpu_vector_data_register(handles(196), 0, C_LOC(QNIN1(LBOUND(QNIN1,1))), SIZE(QNIN1,1), C_SIZEOF(QNIN1(LBOUND(QNIN1,1))))
    ELSE
      handles(196) = C_NULL_PTR
    END IF
    NULLIFY(QNIN1)
    IF (ASSOCIATED(QNIN2)) THEN
      CALL fstarpu_vector_data_register(handles(197), 0, C_LOC(QNIN2(LBOUND(QNIN2,1))), SIZE(QNIN2,1), C_SIZEOF(QNIN2(LBOUND(QNIN2,1))))
    ELSE
      handles(197) = C_NULL_PTR
    END IF
    NULLIFY(QNIN2)
    IF (ASSOCIATED(WSX2)) THEN
      CALL fstarpu_vector_data_register(handles(198), 0, C_LOC(WSX2(LBOUND(WSX2,1))), SIZE(WSX2,1), C_SIZEOF(WSX2(LBOUND(WSX2,1))))
    ELSE
      handles(198) = C_NULL_PTR
    END IF
    NULLIFY(WSX2)
    IF (ASSOCIATED(WSY2)) THEN
      CALL fstarpu_vector_data_register(handles(199), 0, C_LOC(WSY2(LBOUND(WSY2,1))), SIZE(WSY2,1), C_SIZEOF(WSY2(LBOUND(WSY2,1))))
    ELSE
      handles(199) = C_NULL_PTR
    END IF
    NULLIFY(WSY2)
    IF (ASSOCIATED(PR2)) THEN
      CALL fstarpu_vector_data_register(handles(200), 0, C_LOC(PR2(LBOUND(PR2,1))), SIZE(PR2,1), C_SIZEOF(PR2(LBOUND(PR2,1))))
    ELSE
      handles(200) = C_NULL_PTR
    END IF
    NULLIFY(PR2)
    IF (ASSOCIATED(AREAS)) THEN
      CALL fstarpu_vector_data_register(handles(201), 0, C_LOC(AREAS(LBOUND(AREAS,1))), SIZE(AREAS,1), C_SIZEOF(AREAS(LBOUND(AREAS,1))))
    ELSE
      handles(201) = C_NULL_PTR
    END IF
    NULLIFY(AREAS)
    IF (ASSOCIATED(SFACDUB)) THEN
      CALL fstarpu_matrix_data_register(handles(202), 0, C_LOC(SFACDUB(LBOUND(SFACDUB,1),LBOUND(SFACDUB,2))), SIZE(SFACDUB,1), SIZE(SFACDUB,1), SIZE(SFACDUB,2), C_SIZEOF(SFACDUB(LBOUND(SFACDUB,1),LBOUND(SFACDUB,2))))
    ELSE
      handles(202) = C_NULL_PTR
    END IF
    NULLIFY(SFACDUB)
    IF (ASSOCIATED(YDUB)) THEN
      CALL fstarpu_block_data_register(handles(203), 0, C_LOC(YDUB(LBOUND(YDUB,1),LBOUND(YDUB,2),LBOUND(YDUB,3))), SIZE(YDUB,1), SIZE(YDUB,1)*SIZE(YDUB,2), SIZE(YDUB,1), SIZE(YDUB,2), SIZE(YDUB,3), C_SIZEOF(YDUB(LBOUND(YDUB,1),LBOUND(YDUB,2),LBOUND(YDUB,3))))
    ELSE
      handles(203) = C_NULL_PTR
    END IF
    NULLIFY(YDUB)
    IF (ASSOCIATED(TIP1)) THEN
      CALL fstarpu_vector_data_register(handles(204), 0, C_LOC(TIP1(LBOUND(TIP1,1))), SIZE(TIP1,1), C_SIZEOF(TIP1(LBOUND(TIP1,1))))
    ELSE
      handles(204) = C_NULL_PTR
    END IF
    NULLIFY(TIP1)
    IF (ASSOCIATED(TIP2)) THEN
      CALL fstarpu_vector_data_register(handles(205), 0, C_LOC(TIP2(LBOUND(TIP2,1))), SIZE(TIP2,1), C_SIZEOF(TIP2(LBOUND(TIP2,1))))
    ELSE
      handles(205) = C_NULL_PTR
    END IF
    NULLIFY(TIP2)
    IF (ASSOCIATED(NBV)) THEN
      CALL fstarpu_vector_data_register(handles(206), 0, C_LOC(NBV(LBOUND(NBV,1))), SIZE(NBV,1), C_SIZEOF(NBV(LBOUND(NBV,1))))
    ELSE
      handles(206) = C_NULL_PTR
    END IF
    NULLIFY(NBV)
    IF (ASSOCIATED(LBCODEI)) THEN
      CALL fstarpu_vector_data_register(handles(207), 0, C_LOC(LBCODEI(LBOUND(LBCODEI,1))), SIZE(LBCODEI,1), C_SIZEOF(LBCODEI(LBOUND(LBCODEI,1))))
    ELSE
      handles(207) = C_NULL_PTR
    END IF
    NULLIFY(LBCODEI)
    IF (ASSOCIATED(NNODECODE)) THEN
      CALL fstarpu_vector_data_register(handles(208), 0, C_LOC(NNODECODE(LBOUND(NNODECODE,1))), SIZE(NNODECODE,1), C_SIZEOF(NNODECODE(LBOUND(NNODECODE,1))))
    ELSE
      handles(208) = C_NULL_PTR
    END IF
    NULLIFY(NNODECODE)
    IF (ASSOCIATED(NODECODE)) THEN
      CALL fstarpu_vector_data_register(handles(209), 0, C_LOC(NODECODE(LBOUND(NODECODE,1))), SIZE(NODECODE,1), C_SIZEOF(NODECODE(LBOUND(NODECODE,1))))
    ELSE
      handles(209) = C_NULL_PTR
    END IF
    NULLIFY(NODECODE)
    IF (ASSOCIATED(NODEREP)) THEN
      CALL fstarpu_vector_data_register(handles(210), 0, C_LOC(NODEREP(LBOUND(NODEREP,1))), SIZE(NODEREP,1), C_SIZEOF(NODEREP(LBOUND(NODEREP,1))))
    ELSE
      handles(210) = C_NULL_PTR
    END IF
    NULLIFY(NODEREP)
    IF (ASSOCIATED(NM)) THEN
      CALL fstarpu_matrix_data_register(handles(211), 0, C_LOC(NM(LBOUND(NM,1),LBOUND(NM,2))), SIZE(NM,1), SIZE(NM,1), SIZE(NM,2), C_SIZEOF(NM(LBOUND(NM,1),LBOUND(NM,2))))
    ELSE
      handles(211) = C_NULL_PTR
    END IF
    NULLIFY(NM)
    IF (ASSOCIATED(NEITAB)) THEN
      CALL fstarpu_matrix_data_register(handles(212), 0, C_LOC(NEITAB(LBOUND(NEITAB,1),LBOUND(NEITAB,2))), SIZE(NEITAB,1), SIZE(NEITAB,1), SIZE(NEITAB,2), C_SIZEOF(NEITAB(LBOUND(NEITAB,1),LBOUND(NEITAB,2))))
    ELSE
      handles(212) = C_NULL_PTR
    END IF
    NULLIFY(NEITAB)
    IF (ASSOCIATED(NNEIGH_ELEM)) THEN
      CALL fstarpu_vector_data_register(handles(213), 0, C_LOC(NNEIGH_ELEM(LBOUND(NNEIGH_ELEM,1))), SIZE(NNEIGH_ELEM,1), C_SIZEOF(NNEIGH_ELEM(LBOUND(NNEIGH_ELEM,1))))
    ELSE
      handles(213) = C_NULL_PTR
    END IF
    NULLIFY(NNEIGH_ELEM)
    IF (ASSOCIATED(NIBNODECODE)) THEN
      CALL fstarpu_vector_data_register(handles(214), 0, C_LOC(NIBNODECODE(LBOUND(NIBNODECODE,1))), SIZE(NIBNODECODE,1), C_SIZEOF(NIBNODECODE(LBOUND(NIBNODECODE,1))))
    ELSE
      handles(214) = C_NULL_PTR
    END IF
    NULLIFY(NIBNODECODE)
    IF (ASSOCIATED(NEIGH_ELEM)) THEN
      CALL fstarpu_matrix_data_register(handles(215), 0, C_LOC(NEIGH_ELEM(LBOUND(NEIGH_ELEM,1),LBOUND(NEIGH_ELEM,2))), SIZE(NEIGH_ELEM,1), SIZE(NEIGH_ELEM,1), SIZE(NEIGH_ELEM,2), C_SIZEOF(NEIGH_ELEM(LBOUND(NEIGH_ELEM,1),LBOUND(NEIGH_ELEM,2))))
    ELSE
      handles(215) = C_NULL_PTR
    END IF
    NULLIFY(NEIGH_ELEM)
    IF (ASSOCIATED(NVDLL)) THEN
      CALL fstarpu_vector_data_register(handles(216), 0, C_LOC(NVDLL(LBOUND(NVDLL,1))), SIZE(NVDLL,1), C_SIZEOF(NVDLL(LBOUND(NVDLL,1))))
    ELSE
      handles(216) = C_NULL_PTR
    END IF
    NULLIFY(NVDLL)
    IF (ASSOCIATED(NBD)) THEN
      CALL fstarpu_vector_data_register(handles(217), 0, C_LOC(NBD(LBOUND(NBD,1))), SIZE(NBD,1), C_SIZEOF(NBD(LBOUND(NBD,1))))
    ELSE
      handles(217) = C_NULL_PTR
    END IF
    NULLIFY(NBD)
    IF (ASSOCIATED(NBDV)) THEN
      CALL fstarpu_matrix_data_register(handles(218), 0, C_LOC(NBDV(LBOUND(NBDV,1),LBOUND(NBDV,2))), SIZE(NBDV,1), SIZE(NBDV,1), SIZE(NBDV,2), C_SIZEOF(NBDV(LBOUND(NBDV,1),LBOUND(NBDV,2))))
    ELSE
      handles(218) = C_NULL_PTR
    END IF
    NULLIFY(NBDV)
    IF (ASSOCIATED(NVELL)) THEN
      CALL fstarpu_vector_data_register(handles(219), 0, C_LOC(NVELL(LBOUND(NVELL,1))), SIZE(NVELL,1), C_SIZEOF(NVELL(LBOUND(NVELL,1))))
    ELSE
      handles(219) = C_NULL_PTR
    END IF
    NULLIFY(NVELL)
    IF (ASSOCIATED(NBVV)) THEN
      CALL fstarpu_matrix_data_register(handles(220), 0, C_LOC(NBVV(LBOUND(NBVV,1),LBOUND(NBVV,2))), SIZE(NBVV,1), SIZE(NBVV,1), SIZE(NBVV,2), C_SIZEOF(NBVV(LBOUND(NBVV,1),LBOUND(NBVV,2))))
    ELSE
      handles(220) = C_NULL_PTR
    END IF
    NULLIFY(NBVV)
    IF (ASSOCIATED(NELED)) THEN
      CALL fstarpu_matrix_data_register(handles(221), 0, C_LOC(NELED(LBOUND(NELED,1),LBOUND(NELED,2))), SIZE(NELED,1), SIZE(NELED,1), SIZE(NELED,2), C_SIZEOF(NELED(LBOUND(NELED,1),LBOUND(NELED,2))))
    ELSE
      handles(221) = C_NULL_PTR
    END IF
    NULLIFY(NELED)
    IF (ASSOCIATED(SEGTYPE)) THEN
      CALL fstarpu_vector_data_register(handles(222), 0, C_LOC(SEGTYPE(LBOUND(SEGTYPE,1))), SIZE(SEGTYPE,1), C_SIZEOF(SEGTYPE(LBOUND(SEGTYPE,1))))
    ELSE
      handles(222) = C_NULL_PTR
    END IF
    NULLIFY(SEGTYPE)
    IF (ASSOCIATED(NOT_AN_EDGE)) THEN
      CALL fstarpu_vector_data_register(handles(223), 0, C_LOC(NOT_AN_EDGE(LBOUND(NOT_AN_EDGE,1))), SIZE(NOT_AN_EDGE,1), C_SIZEOF(NOT_AN_EDGE(LBOUND(NOT_AN_EDGE,1))))
    ELSE
      handles(223) = C_NULL_PTR
    END IF
    NULLIFY(NOT_AN_EDGE)
    IF (ASSOCIATED(WEIR_BUDDY_NODE)) THEN
      CALL fstarpu_matrix_data_register(handles(224), 0, C_LOC(WEIR_BUDDY_NODE(LBOUND(WEIR_BUDDY_NODE,1),LBOUND(WEIR_BUDDY_NODE,2))), SIZE(WEIR_BUDDY_NODE,1), SIZE(WEIR_BUDDY_NODE,1), SIZE(WEIR_BUDDY_NODE,2), C_SIZEOF(WEIR_BUDDY_NODE(LBOUND(WEIR_BUDDY_NODE,1),LBOUND(WEIR_BUDDY_NODE,2))))
    ELSE
      handles(224) = C_NULL_PTR
    END IF
    NULLIFY(WEIR_BUDDY_NODE)
    IF (ASSOCIATED(ONE_OR_TWO)) THEN
      CALL fstarpu_vector_data_register(handles(225), 0, C_LOC(ONE_OR_TWO(LBOUND(ONE_OR_TWO,1))), SIZE(ONE_OR_TWO,1), C_SIZEOF(ONE_OR_TWO(LBOUND(ONE_OR_TWO,1))))
    ELSE
      handles(225) = C_NULL_PTR
    END IF
    NULLIFY(ONE_OR_TWO)
    IF (ASSOCIATED(EDFLG)) THEN
      CALL fstarpu_matrix_data_register(handles(226), 0, C_LOC(EDFLG(LBOUND(EDFLG,1),LBOUND(EDFLG,2))), SIZE(EDFLG,1), SIZE(EDFLG,1), SIZE(EDFLG,2), C_SIZEOF(EDFLG(LBOUND(EDFLG,1),LBOUND(EDFLG,2))))
    ELSE
      handles(226) = C_NULL_PTR
    END IF
    NULLIFY(EDFLG)
    IF (ASSOCIATED(FFF)) THEN
      CALL fstarpu_vector_data_register(handles(227), 0, C_LOC(FFF(LBOUND(FFF,1))), SIZE(FFF,1), C_SIZEOF(FFF(LBOUND(FFF,1))))
    ELSE
      handles(227) = C_NULL_PTR
    END IF
    NULLIFY(FFF)
    IF (ASSOCIATED(FFACE)) THEN
      CALL fstarpu_vector_data_register(handles(228), 0, C_LOC(FFACE(LBOUND(FFACE,1))), SIZE(FFACE,1), C_SIZEOF(FFACE(LBOUND(FFACE,1))))
    ELSE
      handles(228) = C_NULL_PTR
    END IF
    NULLIFY(FFACE)
    IF (ASSOCIATED(BARINHT)) THEN
      CALL fstarpu_vector_data_register(handles(229), 0, C_LOC(BARINHT(LBOUND(BARINHT,1))), SIZE(BARINHT,1), C_SIZEOF(BARINHT(LBOUND(BARINHT,1))))
    ELSE
      handles(229) = C_NULL_PTR
    END IF
    NULLIFY(BARINHT)
    IF (ASSOCIATED(BARINCFSB)) THEN
      CALL fstarpu_vector_data_register(handles(230), 0, C_LOC(BARINCFSB(LBOUND(BARINCFSB,1))), SIZE(BARINCFSB,1), C_SIZEOF(BARINCFSB(LBOUND(BARINCFSB,1))))
    ELSE
      handles(230) = C_NULL_PTR
    END IF
    NULLIFY(BARINCFSB)
    IF (ASSOCIATED(BARINCFSP)) THEN
      CALL fstarpu_vector_data_register(handles(231), 0, C_LOC(BARINCFSP(LBOUND(BARINCFSP,1))), SIZE(BARINCFSP,1), C_SIZEOF(BARINCFSP(LBOUND(BARINCFSP,1))))
    ELSE
      handles(231) = C_NULL_PTR
    END IF
    NULLIFY(BARINCFSP)
    IF (ASSOCIATED(RBARWL1AVG)) THEN
      CALL fstarpu_vector_data_register(handles(232), 0, C_LOC(RBARWL1AVG(LBOUND(RBARWL1AVG,1))), SIZE(RBARWL1AVG,1), C_SIZEOF(RBARWL1AVG(LBOUND(RBARWL1AVG,1))))
    ELSE
      handles(232) = C_NULL_PTR
    END IF
    NULLIFY(RBARWL1AVG)
    IF (ASSOCIATED(RBARWL2AVG)) THEN
      CALL fstarpu_vector_data_register(handles(233), 0, C_LOC(RBARWL2AVG(LBOUND(RBARWL2AVG,1))), SIZE(RBARWL2AVG,1), C_SIZEOF(RBARWL2AVG(LBOUND(RBARWL2AVG,1))))
    ELSE
      handles(233) = C_NULL_PTR
    END IF
    NULLIFY(RBARWL2AVG)
    IF (ASSOCIATED(IBCONN)) THEN
      CALL fstarpu_vector_data_register(handles(234), 0, C_LOC(IBCONN(LBOUND(IBCONN,1))), SIZE(IBCONN,1), C_SIZEOF(IBCONN(LBOUND(IBCONN,1))))
    ELSE
      handles(234) = C_NULL_PTR
    END IF
    NULLIFY(IBCONN)
    IF (ASSOCIATED(IBCONNR)) THEN
      CALL fstarpu_vector_data_register(handles(235), 0, C_LOC(IBCONNR(LBOUND(IBCONNR,1))), SIZE(IBCONNR,1), C_SIZEOF(IBCONNR(LBOUND(IBCONNR,1))))
    ELSE
      handles(235) = C_NULL_PTR
    END IF
    NULLIFY(IBCONNR)
    IF (ASSOCIATED(NTRAN1)) THEN
      CALL fstarpu_vector_data_register(handles(236), 0, C_LOC(NTRAN1(LBOUND(NTRAN1,1))), SIZE(NTRAN1,1), C_SIZEOF(NTRAN1(LBOUND(NTRAN1,1))))
    ELSE
      handles(236) = C_NULL_PTR
    END IF
    NULLIFY(NTRAN1)
    IF (ASSOCIATED(NTRAN2)) THEN
      CALL fstarpu_vector_data_register(handles(237), 0, C_LOC(NTRAN2(LBOUND(NTRAN2,1))), SIZE(NTRAN2,1), C_SIZEOF(NTRAN2(LBOUND(NTRAN2,1))))
    ELSE
      handles(237) = C_NULL_PTR
    END IF
    NULLIFY(NTRAN2)
    IF (ASSOCIATED(AMIG)) THEN
      CALL fstarpu_vector_data_register(handles(238), 0, C_LOC(AMIG(LBOUND(AMIG,1))), SIZE(AMIG,1), C_SIZEOF(AMIG(LBOUND(AMIG,1))))
    ELSE
      handles(238) = C_NULL_PTR
    END IF
    NULLIFY(AMIG)
    IF (ASSOCIATED(AMIGT)) THEN
      CALL fstarpu_vector_data_register(handles(239), 0, C_LOC(AMIGT(LBOUND(AMIGT,1))), SIZE(AMIGT,1), C_SIZEOF(AMIGT(LBOUND(AMIGT,1))))
    ELSE
      handles(239) = C_NULL_PTR
    END IF
    NULLIFY(AMIGT)
    IF (ASSOCIATED(FAMIG)) THEN
      CALL fstarpu_vector_data_register(handles(240), 0, C_LOC(FAMIG(LBOUND(FAMIG,1))), SIZE(FAMIG,1), C_SIZEOF(FAMIG(LBOUND(FAMIG,1))))
    ELSE
      handles(240) = C_NULL_PTR
    END IF
    NULLIFY(FAMIG)
    IF (ASSOCIATED(PER)) THEN
      CALL fstarpu_vector_data_register(handles(241), 0, C_LOC(PER(LBOUND(PER,1))), SIZE(PER,1), C_SIZEOF(PER(LBOUND(PER,1))))
    ELSE
      handles(241) = C_NULL_PTR
    END IF
    NULLIFY(PER)
    IF (ASSOCIATED(PERT)) THEN
      CALL fstarpu_vector_data_register(handles(242), 0, C_LOC(PERT(LBOUND(PERT,1))), SIZE(PERT,1), C_SIZEOF(PERT(LBOUND(PERT,1))))
    ELSE
      handles(242) = C_NULL_PTR
    END IF
    NULLIFY(PERT)
    IF (ASSOCIATED(FPER)) THEN
      CALL fstarpu_vector_data_register(handles(243), 0, C_LOC(FPER(LBOUND(FPER,1))), SIZE(FPER,1), C_SIZEOF(FPER(LBOUND(FPER,1))))
    ELSE
      handles(243) = C_NULL_PTR
    END IF
    NULLIFY(FPER)
    IF (ASSOCIATED(FREQ)) THEN
      CALL fstarpu_vector_data_register(handles(244), 0, C_LOC(FREQ(LBOUND(FREQ,1))), SIZE(FREQ,1), C_SIZEOF(FREQ(LBOUND(FREQ,1))))
    ELSE
      handles(244) = C_NULL_PTR
    END IF
    NULLIFY(FREQ)
    IF (ASSOCIATED(FF)) THEN
      CALL fstarpu_vector_data_register(handles(245), 0, C_LOC(FF(LBOUND(FF,1))), SIZE(FF,1), C_SIZEOF(FF(LBOUND(FF,1))))
    ELSE
      handles(245) = C_NULL_PTR
    END IF
    NULLIFY(FF)
    IF (ASSOCIATED(FACE)) THEN
      CALL fstarpu_vector_data_register(handles(246), 0, C_LOC(FACE(LBOUND(FACE,1))), SIZE(FACE,1), C_SIZEOF(FACE(LBOUND(FACE,1))))
    ELSE
      handles(246) = C_NULL_PTR
    END IF
    NULLIFY(FACE)
    IF (ASSOCIATED(SLAM)) THEN
      CALL fstarpu_vector_data_register(handles(247), 0, C_LOC(SLAM(LBOUND(SLAM,1))), SIZE(SLAM,1), C_SIZEOF(SLAM(LBOUND(SLAM,1))))
    ELSE
      handles(247) = C_NULL_PTR
    END IF
    NULLIFY(SLAM)
    IF (ASSOCIATED(SFEA)) THEN
      CALL fstarpu_vector_data_register(handles(248), 0, C_LOC(SFEA(LBOUND(SFEA,1))), SIZE(SFEA,1), C_SIZEOF(SFEA(LBOUND(SFEA,1))))
    ELSE
      handles(248) = C_NULL_PTR
    END IF
    NULLIFY(SFEA)
    IF (ASSOCIATED(X)) THEN
      CALL fstarpu_vector_data_register(handles(249), 0, C_LOC(X(LBOUND(X,1))), SIZE(X,1), C_SIZEOF(X(LBOUND(X,1))))
    ELSE
      handles(249) = C_NULL_PTR
    END IF
    NULLIFY(X)
    IF (ASSOCIATED(Y)) THEN
      CALL fstarpu_vector_data_register(handles(250), 0, C_LOC(Y(LBOUND(Y,1))), SIZE(Y,1), C_SIZEOF(Y(LBOUND(Y,1))))
    ELSE
      handles(250) = C_NULL_PTR
    END IF
    NULLIFY(Y)
  END SUBROUTINE DGSWEM_STATE_REGISTER

  SUBROUTINE DGSWEM_STATE_ACTIVATE(buffers)
    TYPE(C_PTR), VALUE, INTENT(IN) :: buffers
    TYPE(C_PTR) :: curr_ptr
    curr_ptr = fstarpu_vector_get_ptr(buffers, 0)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, WDFLG, shape=[fstarpu_vector_get_nx(buffers, 0)])
    ELSE
      NULLIFY(WDFLG)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 1)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, DOFS, shape=[fstarpu_vector_get_nx(buffers, 1)])
    ELSE
      NULLIFY(DOFS)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 2)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, PCOUNT, shape=[fstarpu_vector_get_nx(buffers, 2)])
    ELSE
      NULLIFY(PCOUNT)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 3)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, PDG, shape=[fstarpu_vector_get_nx(buffers, 3)])
    ELSE
      NULLIFY(PDG)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 4)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, NCOUNT, shape=[fstarpu_vector_get_nx(buffers, 4)])
    ELSE
      NULLIFY(NCOUNT)
    END IF
    curr_ptr = fstarpu_matrix_get_ptr(buffers, 5)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, NEDEL, shape=[fstarpu_matrix_get_nx(buffers, 5), fstarpu_matrix_get_ny(buffers, 5)])
    ELSE
      NULLIFY(NEDEL)
    END IF
    curr_ptr = fstarpu_matrix_get_ptr(buffers, 6)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, NEDSD, shape=[fstarpu_matrix_get_nx(buffers, 6), fstarpu_matrix_get_ny(buffers, 6)])
    ELSE
      NULLIFY(NEDSD)
    END IF
    curr_ptr = fstarpu_matrix_get_ptr(buffers, 7)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, NEDNO, shape=[fstarpu_matrix_get_nx(buffers, 7), fstarpu_matrix_get_ny(buffers, 7)])
    ELSE
      NULLIFY(NEDNO)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 8)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, NEDNO1, shape=[fstarpu_vector_get_nx(buffers, 8)])
    ELSE
      NULLIFY(NEDNO1)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 9)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, NEDNO2, shape=[fstarpu_vector_get_nx(buffers, 9)])
    ELSE
      NULLIFY(NEDNO2)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 10)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, NIEDN, shape=[fstarpu_vector_get_nx(buffers, 10)])
    ELSE
      NULLIFY(NIEDN)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 11)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, NLEDN, shape=[fstarpu_vector_get_nx(buffers, 11)])
    ELSE
      NULLIFY(NLEDN)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 12)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, NEEDN, shape=[fstarpu_vector_get_nx(buffers, 12)])
    ELSE
      NULLIFY(NEEDN)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 13)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, NFEDN, shape=[fstarpu_vector_get_nx(buffers, 13)])
    ELSE
      NULLIFY(NFEDN)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 14)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, NREDN, shape=[fstarpu_vector_get_nx(buffers, 14)])
    ELSE
      NULLIFY(NREDN)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 15)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, NEBEDN, shape=[fstarpu_vector_get_nx(buffers, 15)])
    ELSE
      NULLIFY(NEBEDN)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 16)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, NIBEDN, shape=[fstarpu_vector_get_nx(buffers, 16)])
    ELSE
      NULLIFY(NIBEDN)
    END IF
    curr_ptr = fstarpu_matrix_get_ptr(buffers, 17)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, NIBSEGN, shape=[fstarpu_matrix_get_nx(buffers, 17), fstarpu_matrix_get_ny(buffers, 17)])
    ELSE
      NULLIFY(NIBSEGN)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 18)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, NEBSEGN, shape=[fstarpu_vector_get_nx(buffers, 18)])
    ELSE
      NULLIFY(NEBSEGN)
    END IF
    curr_ptr = fstarpu_matrix_get_ptr(buffers, 19)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, EL_NBORS, shape=[fstarpu_matrix_get_nx(buffers, 19), fstarpu_matrix_get_ny(buffers, 19)])
    ELSE
      NULLIFY(EL_NBORS)
    END IF
    curr_ptr = fstarpu_matrix_get_ptr(buffers, 20)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, BACKNODES, shape=[fstarpu_matrix_get_nx(buffers, 20), fstarpu_matrix_get_ny(buffers, 20)])
    ELSE
      NULLIFY(BACKNODES)
    END IF
    curr_ptr = fstarpu_matrix_get_ptr(buffers, 21)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, ATVD, shape=[fstarpu_matrix_get_nx(buffers, 21), fstarpu_matrix_get_ny(buffers, 21)])
    ELSE
      NULLIFY(ATVD)
    END IF
    curr_ptr = fstarpu_matrix_get_ptr(buffers, 22)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, BTVD, shape=[fstarpu_matrix_get_nx(buffers, 22), fstarpu_matrix_get_ny(buffers, 22)])
    ELSE
      NULLIFY(BTVD)
    END IF
    curr_ptr = fstarpu_matrix_get_ptr(buffers, 23)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, CTVD, shape=[fstarpu_matrix_get_nx(buffers, 23), fstarpu_matrix_get_ny(buffers, 23)])
    ELSE
      NULLIFY(CTVD)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 24)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, DTVD, shape=[fstarpu_vector_get_nx(buffers, 24)])
    ELSE
      NULLIFY(DTVD)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 25)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, MAX_BOA_DT, shape=[fstarpu_vector_get_nx(buffers, 25)])
    ELSE
      NULLIFY(MAX_BOA_DT)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 26)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, e1, shape=[fstarpu_vector_get_nx(buffers, 26)])
    ELSE
      NULLIFY(e1)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 27)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, balance, shape=[fstarpu_vector_get_nx(buffers, 27)])
    ELSE
      NULLIFY(balance)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 28)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, RKC_Tdprime, shape=[fstarpu_vector_get_nx(buffers, 28)])
    ELSE
      NULLIFY(RKC_Tdprime)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 29)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, RKC_a, shape=[fstarpu_vector_get_nx(buffers, 29)])
    ELSE
      NULLIFY(RKC_a)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 30)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, RKC_b, shape=[fstarpu_vector_get_nx(buffers, 30)])
    ELSE
      NULLIFY(RKC_b)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 31)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, RKC_c, shape=[fstarpu_vector_get_nx(buffers, 31)])
    ELSE
      NULLIFY(RKC_c)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 32)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, RKC_mu, shape=[fstarpu_vector_get_nx(buffers, 32)])
    ELSE
      NULLIFY(RKC_mu)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 33)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, RKC_tildemu, shape=[fstarpu_vector_get_nx(buffers, 33)])
    ELSE
      NULLIFY(RKC_tildemu)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 34)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, RKC_nu, shape=[fstarpu_vector_get_nx(buffers, 34)])
    ELSE
      NULLIFY(RKC_nu)
    END IF
    curr_ptr = fstarpu_block_get_ptr(buffers, 35)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, BATH, shape=[fstarpu_block_get_nx(buffers, 35), fstarpu_block_get_ny(buffers, 35), fstarpu_block_get_nz(buffers, 35)])
    ELSE
      NULLIFY(BATH)
    END IF
    curr_ptr = fstarpu_block_get_ptr(buffers, 36)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, DBATHDX, shape=[fstarpu_block_get_nx(buffers, 36), fstarpu_block_get_ny(buffers, 36), fstarpu_block_get_nz(buffers, 36)])
    ELSE
      NULLIFY(DBATHDX)
    END IF
    curr_ptr = fstarpu_block_get_ptr(buffers, 37)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, DBATHDY, shape=[fstarpu_block_get_nx(buffers, 37), fstarpu_block_get_ny(buffers, 37), fstarpu_block_get_nz(buffers, 37)])
    ELSE
      NULLIFY(DBATHDY)
    END IF
    curr_ptr = fstarpu_block_get_ptr(buffers, 38)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, SFAC_ELEM, shape=[fstarpu_block_get_nx(buffers, 38), fstarpu_block_get_ny(buffers, 38), fstarpu_block_get_nz(buffers, 38)])
    ELSE
      NULLIFY(SFAC_ELEM)
    END IF
    curr_ptr = fstarpu_tensor_get_ptr(buffers, 39)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, BATHED, shape=[fstarpu_tensor_get_nx(buffers, 39), fstarpu_tensor_get_ny(buffers, 39), fstarpu_tensor_get_nz(buffers, 39), fstarpu_tensor_get_nt(buffers, 39)])
    ELSE
      NULLIFY(BATHED)
    END IF
    curr_ptr = fstarpu_tensor_get_ptr(buffers, 40)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, SFACED, shape=[fstarpu_tensor_get_nx(buffers, 40), fstarpu_tensor_get_ny(buffers, 40), fstarpu_tensor_get_nz(buffers, 40), fstarpu_tensor_get_nt(buffers, 40)])
    ELSE
      NULLIFY(SFACED)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 41)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, COSNX, shape=[fstarpu_vector_get_nx(buffers, 41)])
    ELSE
      NULLIFY(COSNX)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 42)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, SINNX, shape=[fstarpu_vector_get_nx(buffers, 42)])
    ELSE
      NULLIFY(SINNX)
    END IF
    curr_ptr = fstarpu_block_get_ptr(buffers, 43)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, DP_NODE, shape=[fstarpu_block_get_nx(buffers, 43), fstarpu_block_get_ny(buffers, 43), fstarpu_block_get_nz(buffers, 43)])
    ELSE
      NULLIFY(DP_NODE)
    END IF
    curr_ptr = fstarpu_matrix_get_ptr(buffers, 44)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, DP_VOL, shape=[fstarpu_matrix_get_nx(buffers, 44), fstarpu_matrix_get_ny(buffers, 44)])
    ELSE
      NULLIFY(DP_VOL)
    END IF
    curr_ptr = fstarpu_block_get_ptr(buffers, 45)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, DRPHI, shape=[fstarpu_block_get_nx(buffers, 45), fstarpu_block_get_ny(buffers, 45), fstarpu_block_get_nz(buffers, 45)])
    ELSE
      NULLIFY(DRPHI)
    END IF
    curr_ptr = fstarpu_block_get_ptr(buffers, 46)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, DSPHI, shape=[fstarpu_block_get_nx(buffers, 46), fstarpu_block_get_ny(buffers, 46), fstarpu_block_get_nz(buffers, 46)])
    ELSE
      NULLIFY(DSPHI)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 47)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, DRDX, shape=[fstarpu_vector_get_nx(buffers, 47)])
    ELSE
      NULLIFY(DRDX)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 48)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, DSDX, shape=[fstarpu_vector_get_nx(buffers, 48)])
    ELSE
      NULLIFY(DSDX)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 49)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, DRDY, shape=[fstarpu_vector_get_nx(buffers, 49)])
    ELSE
      NULLIFY(DRDY)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 50)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, DSDY, shape=[fstarpu_vector_get_nx(buffers, 50)])
    ELSE
      NULLIFY(DSDY)
    END IF
    curr_ptr = fstarpu_block_get_ptr(buffers, 51)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, EFA_DG, shape=[fstarpu_block_get_nx(buffers, 51), fstarpu_block_get_ny(buffers, 51), fstarpu_block_get_nz(buffers, 51)])
    ELSE
      NULLIFY(EFA_DG)
    END IF
    curr_ptr = fstarpu_block_get_ptr(buffers, 52)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, EMO_DG, shape=[fstarpu_block_get_nx(buffers, 52), fstarpu_block_get_ny(buffers, 52), fstarpu_block_get_nz(buffers, 52)])
    ELSE
      NULLIFY(EMO_DG)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 53)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, XLEN, shape=[fstarpu_vector_get_nx(buffers, 53)])
    ELSE
      NULLIFY(XLEN)
    END IF
    curr_ptr = fstarpu_block_get_ptr(buffers, 54)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, HB, shape=[fstarpu_block_get_nx(buffers, 54), fstarpu_block_get_ny(buffers, 54), fstarpu_block_get_nz(buffers, 54)])
    ELSE
      NULLIFY(HB)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 55)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, IBHT, shape=[fstarpu_vector_get_nx(buffers, 55)])
    ELSE
      NULLIFY(IBHT)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 56)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, EBHT, shape=[fstarpu_vector_get_nx(buffers, 56)])
    ELSE
      NULLIFY(EBHT)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 57)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, EBCFSP, shape=[fstarpu_vector_get_nx(buffers, 57)])
    ELSE
      NULLIFY(EBCFSP)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 58)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, IBCFSP, shape=[fstarpu_vector_get_nx(buffers, 58)])
    ELSE
      NULLIFY(IBCFSP)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 59)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, IBCFSB, shape=[fstarpu_vector_get_nx(buffers, 59)])
    ELSE
      NULLIFY(IBCFSB)
    END IF
    curr_ptr = fstarpu_matrix_get_ptr(buffers, 60)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, M_INV, shape=[fstarpu_matrix_get_nx(buffers, 60), fstarpu_matrix_get_ny(buffers, 60)])
    ELSE
      NULLIFY(M_INV)
    END IF
    curr_ptr = fstarpu_block_get_ptr(buffers, 61)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, phi_edge_fixed, shape=[fstarpu_block_get_nx(buffers, 61), fstarpu_block_get_ny(buffers, 61), fstarpu_block_get_nz(buffers, 61)])
    ELSE
      NULLIFY(phi_edge_fixed)
    END IF
    curr_ptr = fstarpu_block_get_ptr(buffers, 62)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, PHI_AREA, shape=[fstarpu_block_get_nx(buffers, 62), fstarpu_block_get_ny(buffers, 62), fstarpu_block_get_nz(buffers, 62)])
    ELSE
      NULLIFY(PHI_AREA)
    END IF
    curr_ptr = fstarpu_tensor_get_ptr(buffers, 63)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, PHI_EDGE, shape=[fstarpu_tensor_get_nx(buffers, 63), fstarpu_tensor_get_ny(buffers, 63), fstarpu_tensor_get_nz(buffers, 63), fstarpu_tensor_get_nt(buffers, 63)])
    ELSE
      NULLIFY(PHI_EDGE)
    END IF
    curr_ptr = fstarpu_matrix_get_ptr(buffers, 64)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, PHI_CENTER, shape=[fstarpu_matrix_get_nx(buffers, 64), fstarpu_matrix_get_ny(buffers, 64)])
    ELSE
      NULLIFY(PHI_CENTER)
    END IF
    curr_ptr = fstarpu_block_get_ptr(buffers, 65)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, PHI_CORNER, shape=[fstarpu_block_get_nx(buffers, 65), fstarpu_block_get_ny(buffers, 65), fstarpu_block_get_nz(buffers, 65)])
    ELSE
      NULLIFY(PHI_CORNER)
    END IF
    curr_ptr = fstarpu_block_get_ptr(buffers, 66)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, PHI_CHECK, shape=[fstarpu_block_get_nx(buffers, 66), fstarpu_block_get_ny(buffers, 66), fstarpu_block_get_nz(buffers, 66)])
    ELSE
      NULLIFY(PHI_CHECK)
    END IF
    curr_ptr = fstarpu_matrix_get_ptr(buffers, 67)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, PHI_INTEGRATED, shape=[fstarpu_matrix_get_nx(buffers, 67), fstarpu_matrix_get_ny(buffers, 67)])
    ELSE
      NULLIFY(PHI_INTEGRATED)
    END IF
    curr_ptr = fstarpu_matrix_get_ptr(buffers, 68)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, PSI1, shape=[fstarpu_matrix_get_nx(buffers, 68), fstarpu_matrix_get_ny(buffers, 68)])
    ELSE
      NULLIFY(PSI1)
    END IF
    curr_ptr = fstarpu_matrix_get_ptr(buffers, 69)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, PSI2, shape=[fstarpu_matrix_get_nx(buffers, 69), fstarpu_matrix_get_ny(buffers, 69)])
    ELSE
      NULLIFY(PSI2)
    END IF
    curr_ptr = fstarpu_matrix_get_ptr(buffers, 70)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, PSI3, shape=[fstarpu_matrix_get_nx(buffers, 70), fstarpu_matrix_get_ny(buffers, 70)])
    ELSE
      NULLIFY(PSI3)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 71)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, QIB, shape=[fstarpu_vector_get_nx(buffers, 71)])
    ELSE
      NULLIFY(QIB)
    END IF
    curr_ptr = fstarpu_block_get_ptr(buffers, 72)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, QX, shape=[fstarpu_block_get_nx(buffers, 72), fstarpu_block_get_ny(buffers, 72), fstarpu_block_get_nz(buffers, 72)])
    ELSE
      NULLIFY(QX)
    END IF
    curr_ptr = fstarpu_block_get_ptr(buffers, 73)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, QY, shape=[fstarpu_block_get_nx(buffers, 73), fstarpu_block_get_ny(buffers, 73), fstarpu_block_get_nz(buffers, 73)])
    ELSE
      NULLIFY(QY)
    END IF
    curr_ptr = fstarpu_block_get_ptr(buffers, 74)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, ZE, shape=[fstarpu_block_get_nx(buffers, 74), fstarpu_block_get_ny(buffers, 74), fstarpu_block_get_nz(buffers, 74)])
    ELSE
      NULLIFY(ZE)
    END IF
    curr_ptr = fstarpu_block_get_ptr(buffers, 75)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, ze_edge, shape=[fstarpu_block_get_nx(buffers, 75), fstarpu_block_get_ny(buffers, 75), fstarpu_block_get_nz(buffers, 75)])
    ELSE
      NULLIFY(ze_edge)
    END IF
    curr_ptr = fstarpu_block_get_ptr(buffers, 76)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, qx_edge, shape=[fstarpu_block_get_nx(buffers, 76), fstarpu_block_get_ny(buffers, 76), fstarpu_block_get_nz(buffers, 76)])
    ELSE
      NULLIFY(qx_edge)
    END IF
    curr_ptr = fstarpu_block_get_ptr(buffers, 77)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, qy_edge, shape=[fstarpu_block_get_nx(buffers, 77), fstarpu_block_get_ny(buffers, 77), fstarpu_block_get_nz(buffers, 77)])
    ELSE
      NULLIFY(qy_edge)
    END IF
    curr_ptr = fstarpu_matrix_get_ptr(buffers, 78)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, elem_edge, shape=[fstarpu_matrix_get_nx(buffers, 78), fstarpu_matrix_get_ny(buffers, 78)])
    ELSE
      NULLIFY(elem_edge)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 79)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, nieds_count, shape=[fstarpu_vector_get_nx(buffers, 79)])
    ELSE
      NULLIFY(nieds_count)
    END IF
    curr_ptr = fstarpu_tensor_get_ptr(buffers, 80)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, bed, shape=[fstarpu_tensor_get_nx(buffers, 80), fstarpu_tensor_get_ny(buffers, 80), fstarpu_tensor_get_nz(buffers, 80), fstarpu_tensor_get_nt(buffers, 80)])
    ELSE
      NULLIFY(bed)
    END IF
    curr_ptr = fstarpu_block_get_ptr(buffers, 81)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, dynP, shape=[fstarpu_block_get_nx(buffers, 81), fstarpu_block_get_ny(buffers, 81), fstarpu_block_get_nz(buffers, 81)])
    ELSE
      NULLIFY(dynP)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 82)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, dynP_MAX, shape=[fstarpu_vector_get_nx(buffers, 82)])
    ELSE
      NULLIFY(dynP_MAX)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 83)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, dynP_MIN, shape=[fstarpu_vector_get_nx(buffers, 83)])
    ELSE
      NULLIFY(dynP_MIN)
    END IF
    curr_ptr = fstarpu_block_get_ptr(buffers, 84)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, iota, shape=[fstarpu_block_get_nx(buffers, 84), fstarpu_block_get_ny(buffers, 84), fstarpu_block_get_nz(buffers, 84)])
    ELSE
      NULLIFY(iota)
    END IF
    curr_ptr = fstarpu_block_get_ptr(buffers, 85)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, iotaa, shape=[fstarpu_block_get_nx(buffers, 85), fstarpu_block_get_ny(buffers, 85), fstarpu_block_get_nz(buffers, 85)])
    ELSE
      NULLIFY(iotaa)
    END IF
    curr_ptr = fstarpu_block_get_ptr(buffers, 86)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, iota2, shape=[fstarpu_block_get_nx(buffers, 86), fstarpu_block_get_ny(buffers, 86), fstarpu_block_get_nz(buffers, 86)])
    ELSE
      NULLIFY(iota2)
    END IF
    curr_ptr = fstarpu_block_get_ptr(buffers, 87)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, arrayfix, shape=[fstarpu_block_get_nx(buffers, 87), fstarpu_block_get_ny(buffers, 87), fstarpu_block_get_nz(buffers, 87)])
    ELSE
      NULLIFY(arrayfix)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 88)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, CORI_EL, shape=[fstarpu_vector_get_nx(buffers, 88)])
    ELSE
      NULLIFY(CORI_EL)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 89)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, FRIC_EL, shape=[fstarpu_vector_get_nx(buffers, 89)])
    ELSE
      NULLIFY(FRIC_EL)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 90)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, ZE_MAX, shape=[fstarpu_vector_get_nx(buffers, 90)])
    ELSE
      NULLIFY(ZE_MAX)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 91)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, ZE_MIN, shape=[fstarpu_vector_get_nx(buffers, 91)])
    ELSE
      NULLIFY(ZE_MIN)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 92)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, DPE_MIN, shape=[fstarpu_vector_get_nx(buffers, 92)])
    ELSE
      NULLIFY(DPE_MIN)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 93)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, ADVECTQX, shape=[fstarpu_vector_get_nx(buffers, 93)])
    ELSE
      NULLIFY(ADVECTQX)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 94)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, ADVECTQY, shape=[fstarpu_vector_get_nx(buffers, 94)])
    ELSE
      NULLIFY(ADVECTQY)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 95)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, SOURCEQX, shape=[fstarpu_vector_get_nx(buffers, 95)])
    ELSE
      NULLIFY(SOURCEQX)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 96)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, SOURCEQY, shape=[fstarpu_vector_get_nx(buffers, 96)])
    ELSE
      NULLIFY(SOURCEQY)
    END IF
    curr_ptr = fstarpu_tensor_get_ptr(buffers, 97)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, LZ, shape=[fstarpu_tensor_get_nx(buffers, 97), fstarpu_tensor_get_ny(buffers, 97), fstarpu_tensor_get_nz(buffers, 97), fstarpu_tensor_get_nt(buffers, 97)])
    ELSE
      NULLIFY(LZ)
    END IF
    curr_ptr = fstarpu_tensor_get_ptr(buffers, 98)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, MZ, shape=[fstarpu_tensor_get_nx(buffers, 98), fstarpu_tensor_get_ny(buffers, 98), fstarpu_tensor_get_nz(buffers, 98), fstarpu_tensor_get_nt(buffers, 98)])
    ELSE
      NULLIFY(MZ)
    END IF
    curr_ptr = fstarpu_tensor_get_ptr(buffers, 99)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, HZ, shape=[fstarpu_tensor_get_nx(buffers, 99), fstarpu_tensor_get_ny(buffers, 99), fstarpu_tensor_get_nz(buffers, 99), fstarpu_tensor_get_nt(buffers, 99)])
    ELSE
      NULLIFY(HZ)
    END IF
    curr_ptr = fstarpu_tensor_get_ptr(buffers, 100)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, TZ, shape=[fstarpu_tensor_get_nx(buffers, 100), fstarpu_tensor_get_ny(buffers, 100), fstarpu_tensor_get_nz(buffers, 100), fstarpu_tensor_get_nt(buffers, 100)])
    ELSE
      NULLIFY(TZ)
    END IF
    curr_ptr = fstarpu_block_get_ptr(buffers, 101)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, QNAM_DG, shape=[fstarpu_block_get_nx(buffers, 101), fstarpu_block_get_ny(buffers, 101), fstarpu_block_get_nz(buffers, 101)])
    ELSE
      NULLIFY(QNAM_DG)
    END IF
    curr_ptr = fstarpu_block_get_ptr(buffers, 102)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, QNPH_DG, shape=[fstarpu_block_get_nx(buffers, 102), fstarpu_block_get_ny(buffers, 102), fstarpu_block_get_nz(buffers, 102)])
    ELSE
      NULLIFY(QNPH_DG)
    END IF
    curr_ptr = fstarpu_block_get_ptr(buffers, 103)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, RHS_ZE, shape=[fstarpu_block_get_nx(buffers, 103), fstarpu_block_get_ny(buffers, 103), fstarpu_block_get_nz(buffers, 103)])
    ELSE
      NULLIFY(RHS_ZE)
    END IF
    curr_ptr = fstarpu_tensor_get_ptr(buffers, 104)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, RHS_bed, shape=[fstarpu_tensor_get_nx(buffers, 104), fstarpu_tensor_get_ny(buffers, 104), fstarpu_tensor_get_nz(buffers, 104), fstarpu_tensor_get_nt(buffers, 104)])
    ELSE
      NULLIFY(RHS_bed)
    END IF
    curr_ptr = fstarpu_block_get_ptr(buffers, 105)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, RHS_QX, shape=[fstarpu_block_get_nx(buffers, 105), fstarpu_block_get_ny(buffers, 105), fstarpu_block_get_nz(buffers, 105)])
    ELSE
      NULLIFY(RHS_QX)
    END IF
    curr_ptr = fstarpu_block_get_ptr(buffers, 106)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, RHS_QY, shape=[fstarpu_block_get_nx(buffers, 106), fstarpu_block_get_ny(buffers, 106), fstarpu_block_get_nz(buffers, 106)])
    ELSE
      NULLIFY(RHS_QY)
    END IF
    curr_ptr = fstarpu_block_get_ptr(buffers, 107)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, RHS_iota, shape=[fstarpu_block_get_nx(buffers, 107), fstarpu_block_get_ny(buffers, 107), fstarpu_block_get_nz(buffers, 107)])
    ELSE
      NULLIFY(RHS_iota)
    END IF
    curr_ptr = fstarpu_block_get_ptr(buffers, 108)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, RHS_iota2, shape=[fstarpu_block_get_nx(buffers, 108), fstarpu_block_get_ny(buffers, 108), fstarpu_block_get_nz(buffers, 108)])
    ELSE
      NULLIFY(RHS_iota2)
    END IF
    curr_ptr = fstarpu_matrix_get_ptr(buffers, 109)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, XAGP, shape=[fstarpu_matrix_get_nx(buffers, 109), fstarpu_matrix_get_ny(buffers, 109)])
    ELSE
      NULLIFY(XAGP)
    END IF
    curr_ptr = fstarpu_matrix_get_ptr(buffers, 110)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, YAGP, shape=[fstarpu_matrix_get_nx(buffers, 110), fstarpu_matrix_get_ny(buffers, 110)])
    ELSE
      NULLIFY(YAGP)
    END IF
    curr_ptr = fstarpu_matrix_get_ptr(buffers, 111)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, WAGP, shape=[fstarpu_matrix_get_nx(buffers, 111), fstarpu_matrix_get_ny(buffers, 111)])
    ELSE
      NULLIFY(WAGP)
    END IF
    curr_ptr = fstarpu_matrix_get_ptr(buffers, 112)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, XEGP, shape=[fstarpu_matrix_get_nx(buffers, 112), fstarpu_matrix_get_ny(buffers, 112)])
    ELSE
      NULLIFY(XEGP)
    END IF
    curr_ptr = fstarpu_matrix_get_ptr(buffers, 113)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, YEGP, shape=[fstarpu_matrix_get_nx(buffers, 113), fstarpu_matrix_get_ny(buffers, 113)])
    ELSE
      NULLIFY(YEGP)
    END IF
    curr_ptr = fstarpu_matrix_get_ptr(buffers, 114)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, WEGP, shape=[fstarpu_matrix_get_nx(buffers, 114), fstarpu_matrix_get_ny(buffers, 114)])
    ELSE
      NULLIFY(WEGP)
    END IF
    curr_ptr = fstarpu_matrix_get_ptr(buffers, 115)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, SL3, shape=[fstarpu_matrix_get_nx(buffers, 115), fstarpu_matrix_get_ny(buffers, 115)])
    ELSE
      NULLIFY(SL3)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 116)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, XBC, shape=[fstarpu_vector_get_nx(buffers, 116)])
    ELSE
      NULLIFY(XBC)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 117)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, YBC, shape=[fstarpu_vector_get_nx(buffers, 117)])
    ELSE
      NULLIFY(YBC)
    END IF
    curr_ptr = fstarpu_tensor_get_ptr(buffers, 118)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, XFAC, shape=[fstarpu_tensor_get_nx(buffers, 118), fstarpu_tensor_get_ny(buffers, 118), fstarpu_tensor_get_nz(buffers, 118), fstarpu_tensor_get_nt(buffers, 118)])
    ELSE
      NULLIFY(XFAC)
    END IF
    curr_ptr = fstarpu_tensor_get_ptr(buffers, 119)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, YFAC, shape=[fstarpu_tensor_get_nx(buffers, 119), fstarpu_tensor_get_ny(buffers, 119), fstarpu_tensor_get_nz(buffers, 119), fstarpu_tensor_get_nt(buffers, 119)])
    ELSE
      NULLIFY(YFAC)
    END IF
    curr_ptr = fstarpu_tensor_get_ptr(buffers, 120)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, EDGEQ, shape=[fstarpu_tensor_get_nx(buffers, 120), fstarpu_tensor_get_ny(buffers, 120), fstarpu_tensor_get_nz(buffers, 120), fstarpu_tensor_get_nt(buffers, 120)])
    ELSE
      NULLIFY(EDGEQ)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 121)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, bed_IN, shape=[fstarpu_vector_get_nx(buffers, 121)])
    ELSE
      NULLIFY(bed_IN)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 122)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, bed_EX, shape=[fstarpu_vector_get_nx(buffers, 122)])
    ELSE
      NULLIFY(bed_EX)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 123)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, bed_HAT, shape=[fstarpu_vector_get_nx(buffers, 123)])
    ELSE
      NULLIFY(bed_HAT)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 124)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, fact, shape=[fstarpu_vector_get_nx(buffers, 124)])
    ELSE
      NULLIFY(fact)
    END IF
    curr_ptr = fstarpu_matrix_get_ptr(buffers, 125)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, focal_neigh, shape=[fstarpu_matrix_get_nx(buffers, 125), fstarpu_matrix_get_ny(buffers, 125)])
    ELSE
      NULLIFY(focal_neigh)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 126)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, focal_up, shape=[fstarpu_vector_get_nx(buffers, 126)])
    ELSE
      NULLIFY(focal_up)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 127)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, bi, shape=[fstarpu_vector_get_nx(buffers, 127)])
    ELSE
      NULLIFY(bi)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 128)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, XBCb, shape=[fstarpu_vector_get_nx(buffers, 128)])
    ELSE
      NULLIFY(XBCb)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 129)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, YBCb, shape=[fstarpu_vector_get_nx(buffers, 129)])
    ELSE
      NULLIFY(YBCb)
    END IF
    curr_ptr = fstarpu_matrix_get_ptr(buffers, 130)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, xi1, shape=[fstarpu_matrix_get_nx(buffers, 130), fstarpu_matrix_get_ny(buffers, 130)])
    ELSE
      NULLIFY(xi1)
    END IF
    curr_ptr = fstarpu_matrix_get_ptr(buffers, 131)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, xi2, shape=[fstarpu_matrix_get_nx(buffers, 131), fstarpu_matrix_get_ny(buffers, 131)])
    ELSE
      NULLIFY(xi2)
    END IF
    curr_ptr = fstarpu_matrix_get_ptr(buffers, 132)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, xtransform, shape=[fstarpu_matrix_get_nx(buffers, 132), fstarpu_matrix_get_ny(buffers, 132)])
    ELSE
      NULLIFY(xtransform)
    END IF
    curr_ptr = fstarpu_matrix_get_ptr(buffers, 133)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, ytransform, shape=[fstarpu_matrix_get_nx(buffers, 133), fstarpu_matrix_get_ny(buffers, 133)])
    ELSE
      NULLIFY(ytransform)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 134)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, xi1BCb, shape=[fstarpu_vector_get_nx(buffers, 134)])
    ELSE
      NULLIFY(xi1BCb)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 135)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, xi2BCb, shape=[fstarpu_vector_get_nx(buffers, 135)])
    ELSE
      NULLIFY(xi2BCb)
    END IF
    curr_ptr = fstarpu_matrix_get_ptr(buffers, 136)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, xi1vert, shape=[fstarpu_matrix_get_nx(buffers, 136), fstarpu_matrix_get_ny(buffers, 136)])
    ELSE
      NULLIFY(xi1vert)
    END IF
    curr_ptr = fstarpu_matrix_get_ptr(buffers, 137)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, xi2vert, shape=[fstarpu_matrix_get_nx(buffers, 137), fstarpu_matrix_get_ny(buffers, 137)])
    ELSE
      NULLIFY(xi2vert)
    END IF
    curr_ptr = fstarpu_matrix_get_ptr(buffers, 138)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, xtransformv, shape=[fstarpu_matrix_get_nx(buffers, 138), fstarpu_matrix_get_ny(buffers, 138)])
    ELSE
      NULLIFY(xtransformv)
    END IF
    curr_ptr = fstarpu_matrix_get_ptr(buffers, 139)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, XBCv, shape=[fstarpu_matrix_get_nx(buffers, 139), fstarpu_matrix_get_ny(buffers, 139)])
    ELSE
      NULLIFY(XBCv)
    END IF
    curr_ptr = fstarpu_matrix_get_ptr(buffers, 140)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, YBCv, shape=[fstarpu_matrix_get_nx(buffers, 140), fstarpu_matrix_get_ny(buffers, 140)])
    ELSE
      NULLIFY(YBCv)
    END IF
    curr_ptr = fstarpu_matrix_get_ptr(buffers, 141)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, xi1BCv, shape=[fstarpu_matrix_get_nx(buffers, 141), fstarpu_matrix_get_ny(buffers, 141)])
    ELSE
      NULLIFY(xi1BCv)
    END IF
    curr_ptr = fstarpu_matrix_get_ptr(buffers, 142)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, xi2BCv, shape=[fstarpu_matrix_get_nx(buffers, 142), fstarpu_matrix_get_ny(buffers, 142)])
    ELSE
      NULLIFY(xi2BCv)
    END IF
    curr_ptr = fstarpu_block_get_ptr(buffers, 143)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, Area_integral, shape=[fstarpu_block_get_nx(buffers, 143), fstarpu_block_get_ny(buffers, 143), fstarpu_block_get_nz(buffers, 143)])
    ELSE
      NULLIFY(Area_integral)
    END IF
    curr_ptr = fstarpu_tensor_get_ptr(buffers, 144)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, f, shape=[fstarpu_tensor_get_nx(buffers, 144), fstarpu_tensor_get_ny(buffers, 144), fstarpu_tensor_get_nz(buffers, 144), fstarpu_tensor_get_nt(buffers, 144)])
    ELSE
      NULLIFY(f)
    END IF
    curr_ptr = fstarpu_tensor_get_ptr(buffers, 145)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, g0, shape=[fstarpu_tensor_get_nx(buffers, 145), fstarpu_tensor_get_ny(buffers, 145), fstarpu_tensor_get_nz(buffers, 145), fstarpu_tensor_get_nt(buffers, 145)])
    ELSE
      NULLIFY(g0)
    END IF
    curr_ptr = fstarpu_tensor_get_ptr(buffers, 146)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, varsigma0, shape=[fstarpu_tensor_get_nx(buffers, 146), fstarpu_tensor_get_ny(buffers, 146), fstarpu_tensor_get_nz(buffers, 146), fstarpu_tensor_get_nt(buffers, 146)])
    ELSE
      NULLIFY(varsigma0)
    END IF
    curr_ptr = fstarpu_tensor_get_ptr(buffers, 147)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, fv, shape=[fstarpu_tensor_get_nx(buffers, 147), fstarpu_tensor_get_ny(buffers, 147), fstarpu_tensor_get_nz(buffers, 147), fstarpu_tensor_get_nt(buffers, 147)])
    ELSE
      NULLIFY(fv)
    END IF
    curr_ptr = fstarpu_tensor_get_ptr(buffers, 148)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, g0v, shape=[fstarpu_tensor_get_nx(buffers, 148), fstarpu_tensor_get_ny(buffers, 148), fstarpu_tensor_get_nz(buffers, 148), fstarpu_tensor_get_nt(buffers, 148)])
    ELSE
      NULLIFY(g0v)
    END IF
    curr_ptr = fstarpu_block_get_ptr(buffers, 149)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, var2sigmag, shape=[fstarpu_block_get_nx(buffers, 149), fstarpu_block_get_ny(buffers, 149), fstarpu_block_get_nz(buffers, 149)])
    ELSE
      NULLIFY(var2sigmag)
    END IF
    curr_ptr = fstarpu_block_get_ptr(buffers, 150)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, var2sigmav, shape=[fstarpu_block_get_nx(buffers, 150), fstarpu_block_get_ny(buffers, 150), fstarpu_block_get_nz(buffers, 150)])
    ELSE
      NULLIFY(var2sigmav)
    END IF
    curr_ptr = fstarpu_tensor_get_ptr(buffers, 151)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, Nmatrix, shape=[fstarpu_tensor_get_nx(buffers, 151), fstarpu_tensor_get_ny(buffers, 151), fstarpu_tensor_get_nz(buffers, 151), fstarpu_tensor_get_nt(buffers, 151)])
    ELSE
      NULLIFY(Nmatrix)
    END IF
    curr_ptr = fstarpu_tensor_get_ptr(buffers, 152)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, NmatrixInv, shape=[fstarpu_tensor_get_nx(buffers, 152), fstarpu_tensor_get_ny(buffers, 152), fstarpu_tensor_get_nz(buffers, 152), fstarpu_tensor_get_nt(buffers, 152)])
    ELSE
      NULLIFY(NmatrixInv)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 153)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, deltx, shape=[fstarpu_vector_get_nx(buffers, 153)])
    ELSE
      NULLIFY(deltx)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 154)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, delty, shape=[fstarpu_vector_get_nx(buffers, 154)])
    ELSE
      NULLIFY(delty)
    END IF
    curr_ptr = fstarpu_block_get_ptr(buffers, 155)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, pmatrix, shape=[fstarpu_block_get_nx(buffers, 155), fstarpu_block_get_ny(buffers, 155), fstarpu_block_get_nz(buffers, 155)])
    ELSE
      NULLIFY(pmatrix)
    END IF
    curr_ptr = fstarpu_matrix_get_ptr(buffers, 156)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, ZEmin, shape=[fstarpu_matrix_get_nx(buffers, 156), fstarpu_matrix_get_ny(buffers, 156)])
    ELSE
      NULLIFY(ZEmin)
    END IF
    curr_ptr = fstarpu_matrix_get_ptr(buffers, 157)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, ZEmax, shape=[fstarpu_matrix_get_nx(buffers, 157), fstarpu_matrix_get_ny(buffers, 157)])
    ELSE
      NULLIFY(ZEmax)
    END IF
    curr_ptr = fstarpu_matrix_get_ptr(buffers, 158)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, QXmin, shape=[fstarpu_matrix_get_nx(buffers, 158), fstarpu_matrix_get_ny(buffers, 158)])
    ELSE
      NULLIFY(QXmin)
    END IF
    curr_ptr = fstarpu_matrix_get_ptr(buffers, 159)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, QXmax, shape=[fstarpu_matrix_get_nx(buffers, 159), fstarpu_matrix_get_ny(buffers, 159)])
    ELSE
      NULLIFY(QXmax)
    END IF
    curr_ptr = fstarpu_matrix_get_ptr(buffers, 160)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, QYmin, shape=[fstarpu_matrix_get_nx(buffers, 160), fstarpu_matrix_get_ny(buffers, 160)])
    ELSE
      NULLIFY(QYmin)
    END IF
    curr_ptr = fstarpu_matrix_get_ptr(buffers, 161)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, QYmax, shape=[fstarpu_matrix_get_nx(buffers, 161), fstarpu_matrix_get_ny(buffers, 161)])
    ELSE
      NULLIFY(QYmax)
    END IF
    curr_ptr = fstarpu_matrix_get_ptr(buffers, 162)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, iotamin, shape=[fstarpu_matrix_get_nx(buffers, 162), fstarpu_matrix_get_ny(buffers, 162)])
    ELSE
      NULLIFY(iotamin)
    END IF
    curr_ptr = fstarpu_matrix_get_ptr(buffers, 163)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, iotamax, shape=[fstarpu_matrix_get_nx(buffers, 163), fstarpu_matrix_get_ny(buffers, 163)])
    ELSE
      NULLIFY(iotamax)
    END IF
    curr_ptr = fstarpu_matrix_get_ptr(buffers, 164)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, iota2min, shape=[fstarpu_matrix_get_nx(buffers, 164), fstarpu_matrix_get_ny(buffers, 164)])
    ELSE
      NULLIFY(iota2min)
    END IF
    curr_ptr = fstarpu_matrix_get_ptr(buffers, 165)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, iota2max, shape=[fstarpu_matrix_get_nx(buffers, 165), fstarpu_matrix_get_ny(buffers, 165)])
    ELSE
      NULLIFY(iota2max)
    END IF
    curr_ptr = fstarpu_matrix_get_ptr(buffers, 166)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, ANGTAB, shape=[fstarpu_matrix_get_nx(buffers, 166), fstarpu_matrix_get_ny(buffers, 166)])
    ELSE
      NULLIFY(ANGTAB)
    END IF
    curr_ptr = fstarpu_matrix_get_ptr(buffers, 167)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, CENTAB, shape=[fstarpu_matrix_get_nx(buffers, 167), fstarpu_matrix_get_ny(buffers, 167)])
    ELSE
      NULLIFY(CENTAB)
    END IF
    curr_ptr = fstarpu_matrix_get_ptr(buffers, 168)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, ELETAB, shape=[fstarpu_matrix_get_nx(buffers, 168), fstarpu_matrix_get_ny(buffers, 168)])
    ELSE
      NULLIFY(ELETAB)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 169)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, EL_COUNT, shape=[fstarpu_vector_get_nx(buffers, 169)])
    ELSE
      NULLIFY(EL_COUNT)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 170)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, NNDEL, shape=[fstarpu_vector_get_nx(buffers, 170)])
    ELSE
      NULLIFY(NNDEL)
    END IF
    curr_ptr = fstarpu_matrix_get_ptr(buffers, 171)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, NDEL, shape=[fstarpu_matrix_get_nx(buffers, 171), fstarpu_matrix_get_ny(buffers, 171)])
    ELSE
      NULLIFY(NDEL)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 172)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, iota2_DG, shape=[fstarpu_vector_get_nx(buffers, 172)])
    ELSE
      NULLIFY(iota2_DG)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 173)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, iota_DG, shape=[fstarpu_vector_get_nx(buffers, 173)])
    ELSE
      NULLIFY(iota_DG)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 174)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, iotaa_DG, shape=[fstarpu_vector_get_nx(buffers, 174)])
    ELSE
      NULLIFY(iotaa_DG)
    END IF
    curr_ptr = fstarpu_matrix_get_ptr(buffers, 175)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, bed_DG, shape=[fstarpu_matrix_get_nx(buffers, 175), fstarpu_matrix_get_ny(buffers, 175)])
    ELSE
      NULLIFY(bed_DG)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 176)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, bed_N_int, shape=[fstarpu_vector_get_nx(buffers, 176)])
    ELSE
      NULLIFY(bed_N_int)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 177)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, bed_N_ext, shape=[fstarpu_vector_get_nx(buffers, 177)])
    ELSE
      NULLIFY(bed_N_ext)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 178)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, pdg_el, shape=[fstarpu_vector_get_nx(buffers, 178)])
    ELSE
      NULLIFY(pdg_el)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 179)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, ETAS, shape=[fstarpu_vector_get_nx(buffers, 179)])
    ELSE
      NULLIFY(ETAS)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 180)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, ETA1, shape=[fstarpu_vector_get_nx(buffers, 180)])
    ELSE
      NULLIFY(ETA1)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 181)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, ETA2, shape=[fstarpu_vector_get_nx(buffers, 181)])
    ELSE
      NULLIFY(ETA2)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 182)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, ETAMAX, shape=[fstarpu_vector_get_nx(buffers, 182)])
    ELSE
      NULLIFY(ETAMAX)
    END IF
    curr_ptr = fstarpu_matrix_get_ptr(buffers, 183)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, entrop, shape=[fstarpu_matrix_get_nx(buffers, 183), fstarpu_matrix_get_ny(buffers, 183)])
    ELSE
      NULLIFY(entrop)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 184)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, tracer, shape=[fstarpu_vector_get_nx(buffers, 184)])
    ELSE
      NULLIFY(tracer)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 185)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, tracer2, shape=[fstarpu_vector_get_nx(buffers, 185)])
    ELSE
      NULLIFY(tracer2)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 186)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, MassMax, shape=[fstarpu_vector_get_nx(buffers, 186)])
    ELSE
      NULLIFY(MassMax)
    END IF
    curr_ptr = fstarpu_matrix_get_ptr(buffers, 187)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, bed_int, shape=[fstarpu_matrix_get_nx(buffers, 187), fstarpu_matrix_get_ny(buffers, 187)])
    ELSE
      NULLIFY(bed_int)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 188)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, DP, shape=[fstarpu_vector_get_nx(buffers, 188)])
    ELSE
      NULLIFY(DP)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 189)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, DP0, shape=[fstarpu_vector_get_nx(buffers, 189)])
    ELSE
      NULLIFY(DP0)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 190)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, DPe, shape=[fstarpu_vector_get_nx(buffers, 190)])
    ELSE
      NULLIFY(DPe)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 191)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, SFAC, shape=[fstarpu_vector_get_nx(buffers, 191)])
    ELSE
      NULLIFY(SFAC)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 192)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, CORIF, shape=[fstarpu_vector_get_nx(buffers, 192)])
    ELSE
      NULLIFY(CORIF)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 193)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, ESBIN1, shape=[fstarpu_vector_get_nx(buffers, 193)])
    ELSE
      NULLIFY(ESBIN1)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 194)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, ESBIN2, shape=[fstarpu_vector_get_nx(buffers, 194)])
    ELSE
      NULLIFY(ESBIN2)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 195)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, QNIN1, shape=[fstarpu_vector_get_nx(buffers, 195)])
    ELSE
      NULLIFY(QNIN1)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 196)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, QNIN2, shape=[fstarpu_vector_get_nx(buffers, 196)])
    ELSE
      NULLIFY(QNIN2)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 197)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, WSX2, shape=[fstarpu_vector_get_nx(buffers, 197)])
    ELSE
      NULLIFY(WSX2)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 198)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, WSY2, shape=[fstarpu_vector_get_nx(buffers, 198)])
    ELSE
      NULLIFY(WSY2)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 199)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, PR2, shape=[fstarpu_vector_get_nx(buffers, 199)])
    ELSE
      NULLIFY(PR2)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 200)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, AREAS, shape=[fstarpu_vector_get_nx(buffers, 200)])
    ELSE
      NULLIFY(AREAS)
    END IF
    curr_ptr = fstarpu_matrix_get_ptr(buffers, 201)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, SFACDUB, shape=[fstarpu_matrix_get_nx(buffers, 201), fstarpu_matrix_get_ny(buffers, 201)])
    ELSE
      NULLIFY(SFACDUB)
    END IF
    curr_ptr = fstarpu_block_get_ptr(buffers, 202)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, YDUB, shape=[fstarpu_block_get_nx(buffers, 202), fstarpu_block_get_ny(buffers, 202), fstarpu_block_get_nz(buffers, 202)])
    ELSE
      NULLIFY(YDUB)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 203)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, TIP1, shape=[fstarpu_vector_get_nx(buffers, 203)])
    ELSE
      NULLIFY(TIP1)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 204)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, TIP2, shape=[fstarpu_vector_get_nx(buffers, 204)])
    ELSE
      NULLIFY(TIP2)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 205)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, NBV, shape=[fstarpu_vector_get_nx(buffers, 205)])
    ELSE
      NULLIFY(NBV)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 206)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, LBCODEI, shape=[fstarpu_vector_get_nx(buffers, 206)])
    ELSE
      NULLIFY(LBCODEI)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 207)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, NNODECODE, shape=[fstarpu_vector_get_nx(buffers, 207)])
    ELSE
      NULLIFY(NNODECODE)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 208)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, NODECODE, shape=[fstarpu_vector_get_nx(buffers, 208)])
    ELSE
      NULLIFY(NODECODE)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 209)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, NODEREP, shape=[fstarpu_vector_get_nx(buffers, 209)])
    ELSE
      NULLIFY(NODEREP)
    END IF
    curr_ptr = fstarpu_matrix_get_ptr(buffers, 210)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, NM, shape=[fstarpu_matrix_get_nx(buffers, 210), fstarpu_matrix_get_ny(buffers, 210)])
    ELSE
      NULLIFY(NM)
    END IF
    curr_ptr = fstarpu_matrix_get_ptr(buffers, 211)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, NEITAB, shape=[fstarpu_matrix_get_nx(buffers, 211), fstarpu_matrix_get_ny(buffers, 211)])
    ELSE
      NULLIFY(NEITAB)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 212)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, NNEIGH_ELEM, shape=[fstarpu_vector_get_nx(buffers, 212)])
    ELSE
      NULLIFY(NNEIGH_ELEM)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 213)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, NIBNODECODE, shape=[fstarpu_vector_get_nx(buffers, 213)])
    ELSE
      NULLIFY(NIBNODECODE)
    END IF
    curr_ptr = fstarpu_matrix_get_ptr(buffers, 214)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, NEIGH_ELEM, shape=[fstarpu_matrix_get_nx(buffers, 214), fstarpu_matrix_get_ny(buffers, 214)])
    ELSE
      NULLIFY(NEIGH_ELEM)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 215)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, NVDLL, shape=[fstarpu_vector_get_nx(buffers, 215)])
    ELSE
      NULLIFY(NVDLL)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 216)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, NBD, shape=[fstarpu_vector_get_nx(buffers, 216)])
    ELSE
      NULLIFY(NBD)
    END IF
    curr_ptr = fstarpu_matrix_get_ptr(buffers, 217)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, NBDV, shape=[fstarpu_matrix_get_nx(buffers, 217), fstarpu_matrix_get_ny(buffers, 217)])
    ELSE
      NULLIFY(NBDV)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 218)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, NVELL, shape=[fstarpu_vector_get_nx(buffers, 218)])
    ELSE
      NULLIFY(NVELL)
    END IF
    curr_ptr = fstarpu_matrix_get_ptr(buffers, 219)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, NBVV, shape=[fstarpu_matrix_get_nx(buffers, 219), fstarpu_matrix_get_ny(buffers, 219)])
    ELSE
      NULLIFY(NBVV)
    END IF
    curr_ptr = fstarpu_matrix_get_ptr(buffers, 220)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, NELED, shape=[fstarpu_matrix_get_nx(buffers, 220), fstarpu_matrix_get_ny(buffers, 220)])
    ELSE
      NULLIFY(NELED)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 221)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, SEGTYPE, shape=[fstarpu_vector_get_nx(buffers, 221)])
    ELSE
      NULLIFY(SEGTYPE)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 222)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, NOT_AN_EDGE, shape=[fstarpu_vector_get_nx(buffers, 222)])
    ELSE
      NULLIFY(NOT_AN_EDGE)
    END IF
    curr_ptr = fstarpu_matrix_get_ptr(buffers, 223)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, WEIR_BUDDY_NODE, shape=[fstarpu_matrix_get_nx(buffers, 223), fstarpu_matrix_get_ny(buffers, 223)])
    ELSE
      NULLIFY(WEIR_BUDDY_NODE)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 224)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, ONE_OR_TWO, shape=[fstarpu_vector_get_nx(buffers, 224)])
    ELSE
      NULLIFY(ONE_OR_TWO)
    END IF
    curr_ptr = fstarpu_matrix_get_ptr(buffers, 225)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, EDFLG, shape=[fstarpu_matrix_get_nx(buffers, 225), fstarpu_matrix_get_ny(buffers, 225)])
    ELSE
      NULLIFY(EDFLG)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 226)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, FFF, shape=[fstarpu_vector_get_nx(buffers, 226)])
    ELSE
      NULLIFY(FFF)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 227)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, FFACE, shape=[fstarpu_vector_get_nx(buffers, 227)])
    ELSE
      NULLIFY(FFACE)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 228)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, BARINHT, shape=[fstarpu_vector_get_nx(buffers, 228)])
    ELSE
      NULLIFY(BARINHT)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 229)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, BARINCFSB, shape=[fstarpu_vector_get_nx(buffers, 229)])
    ELSE
      NULLIFY(BARINCFSB)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 230)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, BARINCFSP, shape=[fstarpu_vector_get_nx(buffers, 230)])
    ELSE
      NULLIFY(BARINCFSP)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 231)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, RBARWL1AVG, shape=[fstarpu_vector_get_nx(buffers, 231)])
    ELSE
      NULLIFY(RBARWL1AVG)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 232)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, RBARWL2AVG, shape=[fstarpu_vector_get_nx(buffers, 232)])
    ELSE
      NULLIFY(RBARWL2AVG)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 233)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, IBCONN, shape=[fstarpu_vector_get_nx(buffers, 233)])
    ELSE
      NULLIFY(IBCONN)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 234)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, IBCONNR, shape=[fstarpu_vector_get_nx(buffers, 234)])
    ELSE
      NULLIFY(IBCONNR)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 235)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, NTRAN1, shape=[fstarpu_vector_get_nx(buffers, 235)])
    ELSE
      NULLIFY(NTRAN1)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 236)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, NTRAN2, shape=[fstarpu_vector_get_nx(buffers, 236)])
    ELSE
      NULLIFY(NTRAN2)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 237)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, AMIG, shape=[fstarpu_vector_get_nx(buffers, 237)])
    ELSE
      NULLIFY(AMIG)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 238)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, AMIGT, shape=[fstarpu_vector_get_nx(buffers, 238)])
    ELSE
      NULLIFY(AMIGT)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 239)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, FAMIG, shape=[fstarpu_vector_get_nx(buffers, 239)])
    ELSE
      NULLIFY(FAMIG)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 240)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, PER, shape=[fstarpu_vector_get_nx(buffers, 240)])
    ELSE
      NULLIFY(PER)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 241)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, PERT, shape=[fstarpu_vector_get_nx(buffers, 241)])
    ELSE
      NULLIFY(PERT)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 242)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, FPER, shape=[fstarpu_vector_get_nx(buffers, 242)])
    ELSE
      NULLIFY(FPER)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 243)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, FREQ, shape=[fstarpu_vector_get_nx(buffers, 243)])
    ELSE
      NULLIFY(FREQ)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 244)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, FF, shape=[fstarpu_vector_get_nx(buffers, 244)])
    ELSE
      NULLIFY(FF)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 245)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, FACE, shape=[fstarpu_vector_get_nx(buffers, 245)])
    ELSE
      NULLIFY(FACE)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 246)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, SLAM, shape=[fstarpu_vector_get_nx(buffers, 246)])
    ELSE
      NULLIFY(SLAM)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 247)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, SFEA, shape=[fstarpu_vector_get_nx(buffers, 247)])
    ELSE
      NULLIFY(SFEA)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 248)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, X, shape=[fstarpu_vector_get_nx(buffers, 248)])
    ELSE
      NULLIFY(X)
    END IF
    curr_ptr = fstarpu_vector_get_ptr(buffers, 249)
    IF (C_ASSOCIATED(curr_ptr)) THEN
      CALL c_f_pointer(curr_ptr, Y, shape=[fstarpu_vector_get_nx(buffers, 249)])
    ELSE
      NULLIFY(Y)
    END IF
  END SUBROUTINE DGSWEM_STATE_ACTIVATE
END MODULE DAGSWEM_STATE
