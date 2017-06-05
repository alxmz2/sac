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
 (%i1) load("sac.mc")$

 (%i2) f:matrix([x[2](t-2)*u[1](t)],[u[2](t-3)]);
                               [ x (t - 2) u (t) ]
                               [  2         1    ]
 (%o2)                         [                 ]
                               [    u (t - 3)    ]
                               [     2           ]
 (%i3) h:x[1](t-1)$

 (%i4) S:systdef([f,h],[x,u,y]);
                           [            2 ]      [ x (t - 2)   0  ]
                           [ 0  u (t) _D  ]      [  2             ]
 (%o4) sys(affine, f, dF = [     1        ], g = [                ],
                           [              ]      [              3 ]
                           [ 0      0     ]      [     0      _D  ]
      [ x (t - 2) u (t) ]
      [  2         1    ]
 fg = [                 ], h = [ x (t - 1) ], n = 2, m = 2, p = 1,
      [    u (t - 3)    ]      [  1        ]
      [     2           ]
 statevar = [x (t), x (t)], controlvar = [u (t), u (t)], outputvar = [y (t)],
              1      2                     1      2                    1
 taumax = 3)
 * @endcode
 *
 * @param eq It can be
 * - Matrix \f$f(x_\tau,\ u_\tau)\f$
 * - List of matrices \f$[f(x_\tau,\ u_\tau),\ h(x_\tau)]\f$.
 * @param vars (optional) List of symbols that represent the state, control,
 * and output variables, in that order. If not given, it will be set to [x,u,y].
 * @return system
 * @todo
 * - Convert to affine form
 * - Set the affine flag
 * - compute g(_D)
 * - accept use of u(t) instead of u[1](t) when m=1
 * @bug g(_D) is not correct when system is not affine
 * @note Even if there is only one control input, it has to be noted as u[1](t).
 */
/*v sys systdef (expr eq, list vars)  */
/*v // */ systdef( [pars]
  ) := block([vars,name,tmp,varlist,s,ss],
  name:new(sys),           /* create new system structure        */
  if length(pars)=1
     then (                /* if only one parameter is given,    */
          eq:pop(pars),    /* it is f or [f,h]                   */
          vars:[x,u,y]     /* set default variables              */
     )
     else (
          eq:pop(pars),    /* f, or [f,h]                        */
          vars:pop(pars)   /* [state, control, output]           */
     ),
  if (listp(eq)) then (
     if length(eq) # 2 then error("two arguments are expected"), /* [f,h] */
     name@fg:pop(eq),
     name@h:pop(eq),
     name@p:length(name@h),
    if not matrixp(name@h) then name@h:matrix([name@h])
    )
   else (
    name@fg:eq,
    name@h:0,
    name@p:0
    ),
  name@n:length(name@fg),
  tmp:pop(vars),
  name@statevar:makelist(tmp[i](t),i,1,name@n),
  tmp:pop(vars),
  if vars # []
    then
      (
       name@outputvar:pop(vars),
       if name@h #0 then name@outputvar:makelist(name@outputvar[i](t),i,1,length(name@h))
      ),
  /* find maximal delay */
  name@taumax:maxd(name@fg),
  name@fg:subst(makelist(tmp(t-i)=tmp[1](t-i),i,0,name@taumax),name@fg),
  /* find size of control input */
  varlist:showratvars(name@fg),
  name@m:0,
  for s in varlist do  /* looks for u(t-j) or u[i](t-j) */
    ( if not atom(s) then ss:inpart(s,0) else ss:s,
      if ((ss=tmp) and (name@m=0))
       	   then (
       	 	   name@m:1,
          	 name@controlvar:[tmp[1](t)]
             )
        else (
          if not atom(ss)
            then if inpart(ss,0)=tmp then name@m:max(name@m,inpart(ss,1))
            )
    ), /* for s*/

/* compute dF */
  name@dF:sum(grad(name@fg,tshift(name@statevar,i))*_D^i,i,0,name@taumax),

/* compute g(_D) */
  if (name@m >0 ) then (
     name@controlvar:makelist(tmp[i](t),i,1,name@m),
     name@g:sum(grad(name@fg,tshift(name@controlvar,i))*_D^i,i,0,name@taumax)
  ),
  return(name)
)$
/**
* @brief
* @author A. Garate-Garcia, R. Cuesta-Garcia, and L.A. Marquez-Martinez
*
* Computes the submodules \f$H_k\f$
*
*
* <b>Usage</b>
* @code
(%i1) load("sac.mc")$
(%i2) f:matrix([s*(x[2](t)-x[1](t))+u(t)],[x[1](t)*(b-x[3](t))-x[2](t)],[x[1](t)*x[2](t)-a*x[3](t)]);
                        [ u(t) + s (x (t) - x (t))  ]
                        [            2       1      ]
                        [                           ]
(%o2)                    [ x (t) (b - x (t)) - x (t) ]
                        [  1          3        2    ]
                        [                           ]
                        [   x (t) x (t) - a x (t)   ]
                        [    1     2         3      ]
(%i3) lorenz:systdef(f)$
(%i4) hk(lorenz);
(%o4) [[2, [del(x (t)), del(x (t))]], [3,
                2           3
              [x (t) del(x (t)) - b del(x (t)) + x (t) del(x (t))]], [inf, 0]]
                3         3              3        2         2
(%i5) lorenz@hk:hk(lorenz);
(%o5) [[2, [del(x (t)), del(x (t))]], [3,
                2           3
              [x (t) del(x (t)) - b del(x (t)) + x (t) del(x (t))]], [inf, 0]]
                3         3              3        2         2

* @endcode
*
* @param s system structure
* @return list of pairs [k,{hk}], where {hk} is a basis for submodule Hk.
*/

/*v list     */ hk(
/*v system s */
             ):=block([g,Hi,dx,hk],
  g:copy(s@g),
  gi:copy(g),
  dx:transpose(matrix(map(del,s@statevar))),
  hk:[],
  dimhk:s@m,
  for i:1 thru s@n do (
     Hi:left_kernel(g)*^dx,
     if Hi = 0 then (hk:append(hk,[[inf,0]]),return()),
     if matrixp(Hi) then hk:append(hk,[[i+1,flatten(args(matrixmap(num,Hi)))]])
                    else hk:append(hk,[[i+1,[num(Hi)]]]),
     gi:s@dF*^gi-d_dt(gi,s),
     g:addcol(g,gi)
  ),
return(hk)
)$
