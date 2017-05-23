/**
 * @file system-utils.mc
 * @author A. Garate-Garcia, R. Cuesta-Garcia, and L.A. Marquez-Martinez
 * @date May 2017
 * @brief Routines for definition of dynamical systems.
 *
 */
 /*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< S A C   P R O J E C T >>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/
 /* Computes the differential form of representation, sys_dF and sys_dG, starting from
    the general form, sys_Fg.
    Input  : Fg [matrix]     = matrix sys F g of dimension n x 1
             varsta [vector] = vector of state variables, its of dimension n x 1
             varcon [vector] = vector of control variables, its of dimension m x 1
    Output : dFG [matrix] = variable of dimension n x (n+m). It has the form [sys_dF|sys_dG],
                            where sys_dF is n x n and sys_dG is n x m                 */
 /*===================================================================================*/

 gen2dif(_fg,_varsta,_varcon) := block([_l,_r,_m,_n,_df,_dg,_j,_k,_sta,_con],
     _r:maxd(_fg),   /* maximum delay of Fg */
     _l:length(_fg),
     _m:length(_varsta),
     _n:length(_varcon),
     _df:zeromatrix(_m,_m),
     _dg:zeromatrix(_m,_n),
     _dd:zeromatrix(_m,1),
     for _j:1 thru _m do(
         for _k:1 thru _l do(
             _df[_k,_j]:sum(diff(_fg[_k,1],_varsta[_j,1](t-i),1)*_D^i,i,0,_r)
         )
     ),
     for _j:1 thru _n do(
         for _k:1 thru _l do(
             _dg[_k,_j]:sum(diff(_fg[_k,1],_varcon[_j,1](t-i),1)*_D^i,i,0,_r)
         )
     ),
     for _k:1 thru _l do(
         _dd[_k,1]:sum(diff(_fg[_k,1],q(t-i),1)*_D^i,i,0,_r)
     ),
     addcol(_df,_dd,_dg)
 )$
 
/**
 * @brief System definition
 * @author L.A. Marquez-Martinez
 *
 * Given a system of the form
 \f{eqnarray*}{
    \dot x &=& f(x_\tau,\ u_\tau)\\
    y &=& h(x_\tau)\\
    \mbox{where}\\
    x_\tau &=& x(t),x(t-1),\ldots,x(t-s)\\
    u_\tau &=& u(t),u(t-1),\ldots,u(t-s)\\
 \f}
 * assigns variables to describe it.
 *
 * <b>Usage</b>
 * @code
 * (%i1) load("sac.mc")$
 *
 * @endcode
 *
 * @param name Base name for the variables that will describe the system.
 * @param eq It can be
 * - Matrix f(x(*),u(*))
 * - List of matrices [f(x(*),u(*)), h(x(*))]
 * @param vars List of symbols (typically [x, u, y] ) that represent the state, control,
 * and output variables.
 * @return system
 * @warning Affine and differential form still pending!
 */
systdef(
  /*v symbol */ name,
  /*v var */ eq,
  /*v list */ vars
  ) := block([_fg,_h,_n,_p,_df,_dh,_sv,_cv],
  if (listp(eq)) then (
     if length(eq) # 2 then error("two arguments are expected"),
     _fg:pop(eq),
    _h:pop(eq),
    if not matrixp(_h) then _h:matrix([_h])
    )
   else (
    _fg:eq,
    _h:zeromatrix(1,1)
    ),
  _n:length(_fg),
  _p:length(_h),
  concat(name,_Fg)::_fg,
  concat(name,_H)::_h,
  _m:length(vars(_fg,0,vars[2])),
  _sv:transpose(matrix(makelist(vars[1][i],i,1,_n))),
  _cv:transpose(matrix(makelist(vars[2][i],i,1,_m))),
  concat(name,'_sv):: _sv,
  concat(name,'_cv):: _cv,
  _df:gen2dif(_fg,_sv,_cv),
  _dh:gen2dif(_h,_sv,_cv),
  concat(name,_dF)::submat(1,_n,_df,1,_n),
  concat(name,_dG)::submat(1,_n,_df,_n+2,_n+_m+1),
  concat(name,_dH)::submat(1,_p,_dh,1,_n),
  _af:gen2aff(_fg,vars[1],vars[2]),
  /* TODO FIX bug non affine */
  if (_af#done) then(
  concat(name,_F)::submat(1,_n,_af,1,1),
  concat(name,_P)::submat(1,_n,_af,2,2),
  concat(name,_G)::submat(1,_n,_af,3,length(transpose(_af))),
  concat(name,_H)::_gh
  ),
    disp("System Defined")
);
