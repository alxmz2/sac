/**
 * @file system-utils.mc
 * @author L.A. Marquez-Martinez
 * @date May 2017
 * @brief Routines for definition of dynamical systems.
 *
 */

/* List all variables that depend on t */

makealltvarslist(e):=if not atom(e)
     then ((if matchvar(e) then push(e,mylist)
                           else map(makealltvarslist,args(e))))$

showalltvars(e) := block ([mylist : [] ], makealltvarslist(e), unique(mylist))$

/* This function is similar to showalltvars, but keeping del() or diff() operators. */

 maketvarslist(e):=if not atom(e)
     then (if op(e) = nounify(diff) or op(e) = nounify(del)
               then push(e,mylist)	
               else (if matchvar(e) 
                         then push(e,mylist)
                         else map(maketvarslist,args(e))))$
matchvar(e):=not atom(e)
        and ((diff(args(e),t)=[1]) or subvarp(op(e)))
        and not freeof('t,args(e))$

showtvars(e) := block ([mylist : [] ], maketvarslist(e), unique(mylist))$

/* find max index of a subscripted variable s in expression f */
find_max_idx(f,s):=block([l,maxidx:minf], 
   l:showalltvars(f),
   for e in l do
      ( if subvarp(op(e)) 
           then (if op(op(e))=s then maxidx:max(maxidx,part(op(e),1)))
           else (if op(e)=s then maxidx:max(maxidx,1))
       ),
   return(maxidx)
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
 *
 * There are three main forms of defining a system:
 *
 * A) S:systdef(f);              ( no output x=state, u=control)
 *
 * @code
(%i1) load("sac.mc")$

(%i2) f:matrix([s*(x[2](t)-x[1](t))+u(t)],[x[1](t)*(b-x[3](t))-x[2](t)],[x[1](t)*x[2](t)-a*x[3](t)])$

(%i3) lorenz:systdef(f);
                             [     s (x (t) - x (t))     ]
                             [         2       1         ]
                             [                           ]
(%o3) sys(affine = true, f = [ x (t) (b - x (t)) - x (t) ],
                             [  1          3        2    ]
                             [                           ]
                             [   x (t) x (t) - a x (t)   ]
                             [    1     2         3      ]
     [    - s       s       0    ]
     [                           ]      [ 1 ]
     [ b - x (t)   - 1   - x (t) ]      [   ]
dF = [      3               1    ], g = [ 0 ],
     [                           ]      [   ]
     [   x (t)    x (t)    - a   ]      [ 0 ]
     [    2        1             ]
     [ s (x (t) - x (t)) + u (t) ]
     [     2       1        1    ]
     [                           ]
fg = [ x (t) (b - x (t)) - x (t) ], h = 0, n = 3, m = 1, p = 0,
     [  1          3        2    ]
     [                           ]
     [   x (t) x (t) - a x (t)   ]
     [    1     2         3      ]
statevar = [x (t), x (t), x (t)], controlvar = [u (t)], outputvar = y,
             1      2      3                     1
taumax = 0, hk)
 * @endcode
 *
 *
 * B) S:systdef([f,h]);          ( output \f$y=h(x)\f$ )
 *
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
 * C) S:systdef([f,h],[n,v,z]);  ( n=state, v=control, z=output variables)
 *
 * @code
(%i1) load("sac.mc")$

(%i2) f:matrix([n[2](t)],[v[1](t)])$

(%i3) h:matrix([n[1](t)],[n[2](t-1)]);
                                 [   n (t)   ]
                                 [    1      ]
(%o3)                            [           ]
                                 [ n (t - 1) ]
                                 [  2        ]
(%i4) S:systdef([f,h],[n,v,z]);
                             [ n (t) ]       [ 0  1 ]      [ 0 ]
(%o4) sys(affine = true, f = [  2    ], dF = [      ], g = [   ],
                             [       ]       [ 0  0 ]      [ 1 ]
                             [   0   ]
     [ n (t) ]      [   n (t)   ]
     [  2    ]      [    1      ]
fg = [       ], h = [           ], n = 2, m = 1, p = 2,
     [ v (t) ]      [ n (t - 1) ]
     [  1    ]      [  2        ]
statevar = [n (t), n (t)], controlvar = [v (t)], outputvar = [z (t), z (t)],
             1      2                     1                    1      2
taumax = 0, hk)
 * @endcode
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
  name@dF:sum(gradfnc(name@fg,tshift(name@statevar,i))*_D^i,i,0,name@taumax),

/* compute g(_D) */
  if (name@m >0 ) then (
     name@controlvar:makelist(tmp[i](t),i,1,name@m),
     name@g:sum(gradfnc(name@fg,tshift(name@controlvar,i))*_D^i,i,0,name@taumax)
  ) else name@g:zeromatrix(name@n,1),
  /* system is not affine if g depends on u */
  name@affine: is(find_max_idx(name@g,tmp)<0),
  if name@affine then
      name@f:subst(map(lambda([u],u=0),flatten(makelist(tshift(name@controlvar,i),i,0,name@taumax))),name@fg),
 /* compute hk's  */
  hk(name),
  return(name)
)$

/**
* @brief
* @author A. Garate-Garcia, R. Cuesta-Garcia, and L.A. Marquez-Martinez
*
* Computes the submodules \f$H_k\f$. After calling this routine, a specific \f$H_k\f$ can be recovered using assoc(k,S@@hk).
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
(%i6) assoc(2,lorenz@hk);
(%o6)         [del(x (t)), del(x (t))]
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
