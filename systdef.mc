/*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< S A C   P R O J E C T >>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/
/* Assigns a name to the system equations, which could be written, as was explained
   in Section 1.3, either in general, affine or differential form. When a system is
   defined the software attempts to compute the other representations and stored them
   temporarily during the Maxima session.
   Input  : sys  = name to be used for the dynamic system.
            eq   = matrix ( system dx/dt = f(x,u) )
                   array of matrices [f,h] ( system dx/dt = f(x,u), y = h(x) )
            vars = list with the state, control and output variables:
                                vars=[state,control,output].
   Output : sys_F,sys_G,sys_H [matrix]    = variables created to store data of the
                                            affine form of representation.
            sys_Fg,sys_H [matrix]        = variables created to store data of the
                                            general form of representation.
            sys_dF,sys_dG,sys_dH [matrix] = variables created to store data of the
                                            differential form of representation.
            sys_var_sta,sys_var_con,sys_var_out [vector] = state, control and output
                                                           variables of the system.  */
/*===================================================================================*/

systdef(_name,_eq,_var) := block([_fg,_h,_n,_p,_df,_dh,_sv,_cv],
  if (listp(_eq)) then ( /*  */
     if length(_eq) # 2 then error("two arguments are expected"),
     _fg:pop(_eq),
    _h:pop(_eq),
    if not matrixp(_h) then _h:matrix([_h])
    )
   else (
    _fg:_eq,
    _h:zeromatrix(1,1)
    ),
  _n:length(_fg),
  _p:length(_h),
  concat(_name,_Fg)::_fg,
  concat(_name,_H)::_h,

  _m:length(vars(_fg,0,_var[2])),
  _sv:transpose(matrix(makelist(_var[1][i],i,1,_n))),
  _cv:transpose(matrix(makelist(_var[2][i],i,1,_m))),
  concat(_name,'_sv):: _sv,
  concat(_name,'_cv):: _cv,
  _df:gen2dif(_fg,_sv,_cv),
  _dh:gen2dif(_h,_sv,_cv),
  concat(_name,_dF)::submat(1,_n,_df,1,_n),
  concat(_name,_dG)::submat(1,_n,_df,_n+2,_n+_m+1),
  concat(_name,_dH)::submat(1,_p,_dh,1,_n),
  _af:gen2aff(_fg,_var[1],_var[2]),
  /* TODO FIX bug non affine */
  if (_af#done) then(
  concat(_name,_F)::submat(1,_n,_af,1,1),
  concat(_name,_P)::submat(1,_n,_af,2,2),
  concat(_name,_G)::submat(1,_n,_af,3,length(transpose(_af))),
  concat(_name,_H)::_gh
  ),
    disp("System Defined")
);

/*===================================================================================*/
/* Last modification date:  19/05/11*/
