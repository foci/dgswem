      SUBROUTINE first_loop(IT,lvl)
      USE DG

      IMPLICIT NONE
      INTEGER(4) IT, lvl, kk, kkk, lvl_next

      IF (lvl==1) THEN

	! WORKING AT SMALLEST TS LEVEL
	DO kk=1,Msteps(1)
		! Take M-steps on level=1
		tnkp1_a(1) = tn_a(1) + DtDom(1)
		tnkp1_h(1) = tn_h(1) + DtDom(1)
		CALL DG_HYDRO_TIMESTEP(IT,lvl)
		tn_a(1) = tnkp1_a(1)
		tn_h(1) = tnkp1_h(1)
	ENDDO

      ELSE

	! TIME M-steps on the level==lvl
	DO kkk=1,Msteps(lvl)
		tnkp1_a(lvl) = tn_a(lvl) + DtDom(lvl)
		tnkp1_h(lvl) = tn_h(lvl) + DtDom(lvl)
		lvl_next = lvl - 1

		! Level next is the smallest ts
		IF (lvl_next==1) THEN
			DO kk=1,Msteps(1)
				tnkp1_a(1) = tn_a(1) + DtDom(1)
				tnkp1_h(1) = tn_h(1) + DtDom(1)
				CALL DG_HYDRO_TIMESTEP(IT,lvl_next)
				tn_a(1) = tnkp1_a(1)
				tn_h(1) = tnkp1_h(1)
			ENDDO
		ELSE
			CALL SECOND_LOOP(IT,lvl_next)
		ENDIF

		! Compute ts at level==lvl
		CALL DG_HYDRO_TIMESTEP(IT,lvl)
		tn_a(lvl) = tnkp1_a(lvl)
                tn_h(lvl) = tnkp1_h(lvl)
	ENDDO
      ENDIF

      RETURN
      END SUBROUTINE first_loop


