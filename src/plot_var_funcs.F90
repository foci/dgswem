module plot_var_funcs
  use sizes, only : sz
  implicit none
  
contains
  pure function ze_o(ze,b,qx,qy)
    real(sz) ze_o
    real(sz), intent(in) :: ze,b,qx,qy
    
    ze_o = ze

    return
  end function ze_o

  pure function b_o(ze,b,qx,qy)
    real(sz) b_o
    real(sz), intent(in) :: ze,b,qx,qy
    
    b_o = b

    return
  end function b_o

  pure function qx_o(ze,b,qx,qy)
    real(sz) qx_o
    real(sz), intent(in) :: ze,b,qx,qy
      
    qx_o = qx

    return
  end function qx_o

  pure function qy_o(ze,b,qx,qy)
    real(sz) qy_o
    real(sz), intent(in) :: ze,b,qx,qy
      
    qy_o = qy

    return
  end function qy_o

  pure function u_o(ze,b,qx,qy)
    real(sz) u_o
    real(sz), intent(in) :: ze,b,qx,qy
      
    IF ( ABS(ZE + B) < 10*MAX(ABS(ZE),ABS(B))*epsilon() ) THEN
       u_o = 0.D0
    ELSE
       u_o = qx/(ze + b)
    ENDIF

    return
  end function u_o

  pure function v_o(ze,b,qx,qy)
    real(sz) v_o
    real(sz), intent(in) :: ze,b,qx,qy
      
    IF ( ABS(ZE + B) < 10*MAX(ABS(ZE),ABS(B))*epsilon() ) THEN
       v_o = 0.D0
    ELSE
       v_o = qy/(ze + b)
    ENDIF

    return
  end function v_o

end module plot_var_funcs
