/**
 * @file system-utils.mc
 * @author L.A. Marquez-Martinez
 * @date May 2017
 * @brief Routines for definition of dynamical systems.
 *
 */
 /**
 * @brief
 * @author L.A. Marquez-Martinez
 *
 * given an expression f, and a symbol s, returns max k such that \f$\partial f/\partial s[k]\neq0\f$
 *
 * <b>Usage</b>
 * @code
 * (%i1) load("sac.mc")$
(%i2) find_max_idx(matrix([sin(u[3](t-3)+u(t))+1],[x[3](t)]),u);
(%o2)                                  3
 * @endcode
 *
 * @param f expression
 * @param s name of a variable
 * @return max k such that \f$exists\,i\in I\!\!N,\ \partial f / \partial s[k](t-i)\neq 0\f$.
 * @see systdef
 */

/*v int */ find_max_idx(
/*v expr */ f,
/*v symbol */ s
      ):=block([l,e],
if numberp(f) then return(f),
if atom(f) then return(-1),
e:inpart(f,0),
if not atom(e)
   then (if inpart(e,0)=s
       then (/* if it is u(t) then return 1 else it is u[k](t) and return k */
             if inpart(e,1)=t then return(1) else return(inpart(e,1))
            )
        )
    else if e=s then if inpart(f,1)=t then return(1) else return(inpart(f,1)),
e:sublist(args(f),lambda([u],not atom(u))),
return(apply(max,map(lambda([u],find_max_idx(u,s)),e)))
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
 * - accept use of u(t) instead of u[1](t) when m=1
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
  tmp:pop(vars), /* control var name */
  if vars # []
    then
      (
       name@outputvar:pop(vars),
       if name@h #0 then name@outputvar:makelist(name@outputvar[i](t),i,1,length(name@h))
      ),
  /* find maximal delay */
  name@taumax:maxd(name@fg),
  /* substitute u(t-i) by u[1](t-i) */
  name@fg:subst(makelist(tmp(t-i)=tmp[1](t-i),i,0,name@taumax),name@fg),
  /* find size of control input */
  name@m:find_max_idx(name@fg,tmp),

/* compute dF */
  name@dF:sum(grad(name@fg,tshift(name@statevar,i))*_D^i,i,0,name@taumax),

/* compute g(_D) */
  if (name@m >0 ) then (
     name@controlvar:makelist(tmp[i](t),i,1,name@m),
     name@g:sum(grad(name@fg,tshift(name@controlvar,i))*_D^i,i,0,name@taumax)
  ) else name@g:zeromatrix(name@n,1),
  /* system is not affine if g depends on u */
  name@affine: is(find_max_idx(name@g,tmp)<0),
  if name@affine then
      name@f:subst(map(lambda([u],u=0),flatten(makelist(tshift(name@controlvar,i),i,0,name@taumax))),name@fg),
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
(%i5) lorenz@hk;
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
/*v system */ s
             ):=block([g,Hi,dx,hk,j],
  g:copy(s@g),
  gi:copy(g),
  dx:transpose(matrix(map(del,s@statevar))),
  hk:[],
  dimhk:s@m,
  for i:1 thru s@n do (
     if i=s@n then j:inf else j:i+1,
     Hi:left_kernel(g)*^dx,
     if Hi = 0 then (hk:append(hk,[[inf,0]]),return()),
     if matrixp(Hi) then hk:append(hk,[[j,flatten(args(matrixmap(num,Hi)))]])
                    else hk:append(hk,[[j,[num(Hi)]]]),
     gi:s@dF*^gi-d_dt(gi,s),
     g:addcol(g,gi)
  ),
  s@hk:hk
)$
