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
   and (istvar(e) or subvarp(op(e)))  /* istvar replaces (diff(args(e),t)=[1])  */
   and not freeof('t,args(e))$

istvar(e):=block([ar], ar:args(e), if (ar=['t]) then return(true),
   if (length(ar)>1 or atom(ar[1]) ) then return(false),
   is((op(ar[1])="+") and member('t,args(ar[1]))))$

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

/* System definition */
systdef( [pars]
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

/* Hk submodules */
hk( s ):=block([g,Hi,dx,hk,j],
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
