/**
 * @file operators.mc
 * @author L.A. Marquez-Martinez
 * @date May 20, 2017
 * @brief Miscellaneous operator definitions
 *
 */

 /**
  * @brief  Wedge product
  * @author L.A. Marquez-Martinez
  *
  * Defines the wedge product \f$\Lambda:\mathcal{E}^p\times\mathcal{E}^q\to\mathcal{E}^{p+q}\f$.
  *
  * <b>Usage</b>
  * 
  * @code
  * (%i1) load("sac.mc")$
  * (%i2) wedge(del(x[1](t-1)), del(x[1](t-2),x[2](t)) );
  * (%o2)                - del(x (t - 2), x (t - 1), x (t))
  *                             1          1          2
  * @endcode
  * @param u p-form
  * @param v q-form
  * @return   wedge product \f$u\wedge v\f$
  * @note \f$d(x)\wedge d(y)\f$ is written as \f$d(x,y)\f$.
  * @bugs It will produce invalid results if the coefficients contain _D.
  * @todo create routine that will expand \f$p(\delta)w(t)\f$ as \f$\sum_i p_iw(t-i)\f$.
  */
/*v (p1+p2+..+ps)-form wedge(p1-form w1, p2-form w2,..., ps-form ws){}  */
/*v // */ wedge([ar]):= block([l,u,v,nu,nv,wprod], /* hack for doxygen */
  l:length(ar),
  if l<2 then error("error: at least two arguments are required"),
  u:pop(ar),  

  if l=2 then v:pop(ar)
         else v:tree_reduce(wedge,ar),
  if ( (p_degree(u)=0) or (p_degree(v)=0) ) then return(u*v),
  u:dot_fact(u,0),
  v:dot_fact(v,0),  
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

 /**
  * @fn infix("*^")
  * @brief  Operator *^
  * @author L.A. Marquez-Martinez
  *
  * Defines the non-commutative operator *^. This allows to multiply polynomials in
  * \f$\mathcal{K}[\delta)\f$.
  *
  * <b>Usage</b>
  * p1 *^ p2
  * @code
  * (%i1) load("sac.mc")$
  * (%i2) x(t-1)*_D *^ x(t-2);
  * (%o2)              x(t - 3) x(t - 1) _D
  * @endcode
  * @param p1 polynomial with scalar or matrix coefficients.
  * @param p2 polynomial with scalar or matrix coefficients, or p-form.
  * @return   non-commutative product
  * @note \f$\delta\f$ is written as _D
  */
/*v infix("*^") := p1 *^ p2 {} */  /* this is a hack for Doxygen */
infix("*^",128,127)$ /* binding power to have more precedence than normal product, but less than exponentiation */
/*v // */ "*^"(p1,p2) := block([_pf2,_mp1,_mp2], /* this one too */
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

/**
 * @brief Protects a symbol so it cannot be assigned.
 *
 * Some symbols are reserved for the software.  Assigning them would lead
 * to weird and hard-to-track bugs.  This command avoids this problem by
 * reserving the word.  So far we reserve the words <tt>t, del, true, false</tt>.
 * @param s symbol
 * @returns protected symbol s
 * @see unprotect
 */
/*v protect ( symbol  s ):=block([] */
/*v )$ */
/**
 * @brief Unprotects a symbol.
 *
 * Removes the protection of the argument s.
 * @param s protected symbol
 */
/*v unprotect( symbol s ):=block([] */
/*v )$ */

/**
 * @brief Shifts its argument in time
 *
 * If only one argument is given, shifts a function in time by one unit.
 * If a second argument s is given, it shifts the first argument by s units of time.
 *
 *<b>Usage</b>
 * @code
 * (%i2) tshift(x(t-2),4);
 * (%o2)    x(t - 6)
 * (%i3) tshift([matrix([x(t-1)],[u(t)]),x[3](t-3)]);
 *                            [ x(t - 2) ]
 * (%o6)                     [[          ], x (t - 4)]
 *                            [ u(t - 1) ]   3
 * @endcode
 * @param f any valid function, polynomial, matrix, p-form, or a list of these elements
 * @param s (optional) positive integer
 * @returns the same argument, with all functions of t shifted by 1 unit, or by s units if a second argument is given.
 * @see protect
 */
/*v expr tshift(expr f, integer s):= block([], */
/*v )$ */

/**
 * @brief Computes the differential form of a function.
 * @author L.A. Marquez-Martinez
 *
 * Given \f$f(z_\tau)\f$, this routine computes df:
 \f[
 df = \sum_{i=1}^n \sum_{j=0}^s \frac{\displaystyle \partial f}{\partial z_i(t-j)} dz_i(t-j)
 \f]
 *
 * The partial derivatives are taken against the variables which explicitely depend on \f$t\f$.
 *
 *
 * <b>Usage</b>
 * @code
 (%i1) load("sac.mc")$
 (%i2) _d(x[1](t-2)*u(t)+x[1](t));
 (%o2)       del(x (t)) + x (t - 2) del(u(t)) + u(t) del(x (t - 2))
                  1        1                              1
 * @endcode
 *
 * If \f$ f\f$ is a p-form, then it returns its differential, which is a (p+1)-form.
 *
 * @code
(%i3)  dlist:[x[1](t),(x[1](t)+x[2](t))*del(x[1](t)),(sin(u(t-1))+x[2](t-1)*u(t)^2)*del(x[2](t),u(t-1))];
(%o3) [x (t), (x (t) + x (t)) del(x (t)), 
        1       2       1          1

                                    2
                      - (x (t - 1) u (t) + sin(u(t - 1))) del(u(t - 1), x (t))]
                          2                                              2
(%i4) maplist(_d,dlist);
(%o4) [del(x (t)), - del(x (t), x (t)), 
            1             1      2

(- 2 x (t - 1) u(t) del(u(t - 1), x (t), u(t)))
      2                            2

    2
 - u (t) del(x (t - 1), u(t - 1), x (t))]
              2                    2

 * @endcode
 * @param f function
 * @return df
 */
/*v one_form */  _d(
/*v function */ f):=block([rv,tmp,tmp2,suma,i],
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


