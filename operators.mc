
/**
 * @file operators.mc
 * @author A. Garate-Garcia and L.A. Marquez-Martinez
 * @date May 20, 2017
 * @brief Miscellaneous operator definitions
 *
 */
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
/*v infix("*^") := p1 *^ p2 {} */  /* this is a trick for Doxygen */
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
 df = \sum_{i=1}^n \sum_{j=0}^s \frac{\displaystyle \partial f}{\partial z_i(t-j)}\delta^j dz_i
 \f]
 *
 *
 * <b>Usage</b>
 * @code
 (%i1) load("sac.mc")$
 (%i2) _d(x[1](t-2)*u(t)+x[1](t));
 (%o2)       del(x (t)) + x (t - 2) del(u(t)) + u(t) del(x (t - 2))
                  1        1                              1
 (%i3) dot_fact(%);
                   [                   2     ]  [ del(u(t))  ]
 (%o3)            [[ x (t - 2)  u(t) _D  + 1 ], [            ]]
                   [  1                      ]  [ del(x (t)) ]
                                                [      1     ]
 * @endcode
 *
 * @param f function
 * @return df
 * @todo support for p-forms
 */
/*v one_form */  _d(
/*v function */ f):=block([rv],
  if (not freeof(del,f)) then error("p-forms not yet supported"),
  rv:showratvars(f),
  return(matrix(grad(f,rv)).transpose(matrix(map(del,rv)))))$


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

/*v // */ d_dt([args])
:=block([f,S,i,k,l,c,d,cdt,ddt,vl,vu],
  if length(args)<2 then error ("expected at least 2 arguments"),
  f:pop(args),
  S:pop(args),
  if args=[] then k:1 else k:pop(args),
  if k=0 then return(f),
  if k<0 then error("k cannot be negative"),
  if k>1 then f:d_dt(f,S,k-1),
  if matrixp(f) then return(matrixmap(lambda([u],d_dt(u,S)),S@g)),
  if not freeof(del,f) then
    ( /* case of 1-forms */
      [c,d]:dot_fact(f),
      cdt:map(lambda([u],d_dt(u,S)),c),
      ddt:matrixmap(lambda([u],_d(d_dt(inpart(u,1),S))),d),
      return(ratsimp(c*^ddt+cdt*^d))
    ),
  vl:showratvars(f),
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
