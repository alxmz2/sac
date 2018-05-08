/**
 * @file operators.mc
 * @author L.A. Marquez-Martinez
 * @date May 20, 2017
 * @brief Miscellaneous operator definitions
 *
 */
/* non-commutative product */
infix("*^",128,127)$ /* binding power to have more precedence than normal product, but less than exponentiation */
"*^"(p1,p2) := block([_pf2,_mp1,_mp2], 
    _pf2: not freeof(del,p2),
    _mp1: matrixp(p1),
    _mp2: matrixp(p2),
    if ( (not freeof(del,p1)) )
       then error("wrong arguments"),
    pol1:coefpow(p1),
    if (_mp1 and _mp2) then oper:"."
                       else oper:"*",
    if _pf2 then return(ratsimp(sum(apply(oper,[pol1[1][i],tshift(p2,pol1[2][i])]),i,1,length(pol1[1]))))
            else return(ratsimp(sum(apply(oper,[pol1[1][i],tshift(p2,pol1[2][i])*_D^pol1[2][i]]),i,1,length(pol1[1]))))
)$

/* wedge product */
 wedge([ar]):= block([l,u,v,nu,nv,wprod], 
  l:length(ar),
  if l<2 then error("error: at least two arguments are required"),
  u:pop(ar),  

  if l=2 then v:pop(ar)
         else v:tree_reduce(wedge,ar),
  if (p_degree(u)=0) then return(u*^v),
  if (p_degree(v)=0) then return(v*^u),
  u:dot_fact(u,0),
  if not(freeof(_D,u)) then u:dot_fact(transpose(u[1])*^u[2],0),
  v:dot_fact(v,0),  
  if not(freeof(_D,v)) then v:dot_fact(transpose(v[1])*^v[2],0),
  nu:length(u[1]),
  nv:length(v[1]),
  wprod:0,
  for i:1 thru nu do
     for j:1 thru nv do
       wprod: wprod
             +u[1][i,1]*v[1][j,1]
             * apply(del, flatten([args(u[2][i,1]),args(v[2][j,1])])),
  return(wprod) 
 )$

/* differential operator */
 _d(f):=block([rv,tmp,tmp2,suma,i],
  if (not freeof(del,f)) then
     (suma:0,
      tmp: dot_fact(f,0),
      tmp2:matrixmap(_d,tmp[1]),
      for i:1 thru length(tmp[1]) do
          suma:suma+wedge(tmp2[i,1],tmp[2][i,1]),
      return(suma)
     )
     else
    (rv:showtvars(f),
     if (rv=[]) then return(0) else
                     return(matrix(gradfnc(f,rv)).transpose(matrix(map(del,rv))))
    )
  )$


/**
* @brief Computes the k-th time-derivative of a function
* @author L.A. Marquez-Martinez
*
* \details Given a system S:\f$\{\dot{x}=f(x_\tau,u_\tau)\f$ and a function
 \f$h(x_\tau,u_\tau^{(i)})\f$, find the time-derivative of h along
 the trajectories of S:
 \f[
  \mbox{d\_dt}(h,S)=\sum_{j=0}^s
    \left(
      \sum_{i=1}^n\frac{\displaystyle \partial h}{\displaystyle \partial x(t-j)}f(t-j)
    + \sum_{k=0}^r\frac{\displaystyle \partial h}{\displaystyle \partial u^{(k)}(t-j)}u^{(k+1)}(t-j)
    \right)
 \f]
*
*
* <b>Usage</b>
* @code
* (%i1) load("sac.mc")$
*
* @endcode
*
* @param f function
* @param S system whose trajectories constrain the state.
* @param k (optional) number of times to derivate. Default is 1.
* @return \f$h^{(k)}\f$, following the trajectories of S.
* @see lie
* @note won't complain with p-forms, but will give wrong results.
*/

/*v function d_dt(function f, system S, int k){}  */

/*v // */ d_dt([ar])
:=block([f,S,k,l,c,d,cdt,ddt,vl,vu,a],
  if length(ar)<2 then error ("expected at least 2 arguments"),
  f:pop(ar),
  S:pop(ar),
  if ar=[] then k:1 else k:pop(ar),
  if k=0 then return(f),
  if k<0 then error("k cannot be negative"),
  if k>1 then f:d_dt(f,S,k-1),
  if matrixp(f) then return(matrixmap(lambda([u],d_dt(u,S)),f)),
  if not freeof(del,f) then
    ( /* case of 1-forms */
      if (p_degree(f)>1)
        then error("not implemented for p-forms")
        else
         ( [c,d]:dot_fact(f,0),
           cdt:matrixmap(lambda([u],d_dt(u,S)),c),
           ddt:matrixmap(lambda([u],_d(d_dt(inpart(u,1),S))),d),
           return(ratsimp(c*^ddt+cdt*^d))
         )
    ),
  vl:showtvars(f),
  l:length(vl),
  vu:[],
  for i:1 thru l do
    (
     a:pop(vl),
     if not (member(tshift(a,-rel_shift(a)),S@statevar)) then push(a,vu)
    ),
  l:length(vu),
  return(ratsimp(lie(f,S)+sum(ratcoef(f,vu[i])*diff(vu[i],t),i,1,l)))
)$
/**
 * @brief Returns the integral of a list of closed 1-forms 
 *
 * This function returns the integral form of its argument, which can be a closed 1-form or 
 * a list of closed 1-forms.
 *
 * <b>Usage</b>
 * @code
 * 
(%i1) load("sac.mc")$
(%i2) L:[del(u(t)),x[1](t-1)*del(x[1](t-1))]$
(%i3) antider(L);
                                       2
                                      x (t - 1)
                                       1
(%o3)                          [u(t), ---------]
                                          2
 *
 * @endcode
 * 
 * If one of the arguments is not a closed 1-form (even if it is integrable), then it throws an error.
 *
 * @code
(%i4) antider(x[1](t)*del(x[2](t-1)));

argument is not a [list of ] closed 1-form[s]
@endcode
 *
 * @param L closed 1-form or list of closed 1-forms.
 * @return   Integral of the argument(s).
 * @note The integration is done using the routine @b potential, from the "vect" package.
 * @see _d
 */
/*v list */ antider(
/*v list */ L ):=block([F,vlist,lv,c,d],
if not(is_closed(L)) then error("argument is not a [list of ] closed 1-form[s]"),
if listp(L) 
   then return(maplist(antider,L))
   else (
       if (L=0) then return(0),
       [c,d]:dot_fact(L,0),
       vlist:showtvars([c,matrixmap(args,d)]),
       lv:length(vlist),
       F:subst(makelist(vlist[i]=concat(x,i),i,1,lv),flatten(args(c))),
       scalefactors(makelist(concat(x,i),i,1,lv)),
       return(subst(makelist(concat(x,i)=vlist[i],i,1,lv),potential(F)))
       )
)$


