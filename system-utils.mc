/**
 * @file system-utils.mc
 * @author L.A. Marquez-Martinez
 * @date May 2017
 * @brief Routines for definition of dynamical systems.
 *
 */

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
 * creates a structure to store it.
 *
 * <b>Usage</b>
 * @code
 * (%i1) load("sac.mc")$
 *
 * @endcode
 *
 * @param eq It can be
 * - Matrix f(x(*),u(*))
 * - List of matrices [f(x(*),u(*)), h(x(*))]
 * @param vars List of symbols (typically [x, u, y] ) that represent the state, control,
 * and output variables.
 * @return system
 */
/*v system */ systdef(
  /*v var    */     eq,
  /*v list   */     vars
  ) := block([name,tmp,varlist,s,ss],
  name:new(sys),  /* create new system structure */
  if (listp(eq)) then (
     if length(eq) # 2 then error("two arguments are expected"), /* [f,h] */
     name@fg:pop(eq),
     name@h:pop(eq),
     name@p:length(name@h),
    if not matrixp(name@h) then name@h:matrix([name@h])
    )
   else (
    name@fg:eq
    ),
  name@n:length(name@fg),
  tmp:pop(vars),
  name@statevar:makelist(tmp[i](t),i,1,name@n),
  tmp:pop(vars),
  if vars # [] then name@outputvar:pop(vars),
  /* find size of control input */
  varlist:showratvars(name@fg),
  name@m:0,
  for s in varlist do
    ( ss:inpart(s,0),
      if ((ss=tmp) and (name@m=0))
        then (
          name@m:1,
          name@controlvar:[tmp(t)]
         )
        else (
          if not symbolp(ss)
            then if inpart(ss,0)=tmp then name@m:max(name@m,inpart(ss,1))
            )
    ), /* for s*/
  if (name@m >0 ) then
        if member(tmp(t),varlist) /* check if we used u[1](t) with m=1 */
            then name@controlvar:[tmp(t)]
            else name@controlvar:makelist(tmp[i](t),i,1,name@m),
  name@taumax:maxd(name@fg),
  name@dF:sum(grad(name@fg,tshift(name@statevar,i))*_D^i,i,0,name@taumax),
  if name@m>0 then
     name@g:sum(grad(name@fg,tshift(name@controlvar,i))*_D^i,i,0,name@taumax),
  return(name)
  );
