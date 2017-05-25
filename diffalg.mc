/**
 * @file diffalg.mc
 * @author A. Marquez, A. Garate, and R. Cuesta
 * @date May 20, 2017
 * @brief Differential algebra routines
 *
 * This file contains the functions related to differential algebra for
 * working with p-forms (\f$p\in{Z}^+\f$).
 */
/**
 * @brief Computes the gradient of a function
 * @author L.A. Marquez-Martinez
 *
 * Given a list of variables [\f$v=[v_1,\ldots,v_s]\f$], it will return the partial
 * derivative of function \f$f(\cdot)\f$ with respect to them:
 \f[
  \mbox{grad}(f,v)=\left[\frac{\displaystyle \partial f}{\displaystyle\partial v_1 },\ldots,\frac{\displaystyle \partial f}{\displaystyle\partial v_s }\right]
 \f]
 *
 * <b>Usage</b>
 * @code
(%i1) load("sac.mc")$
(%i2) grad(((x[2](t-1))^2),[x[1](t-1),x[2](t-1)]);
(%o2)                [0, 2 x (t - 1)]
                            2
(%i3) grad(matrix([x[2](t-1)^2],[x[1](t-1)*x[2](t)]),[x[1](t-1),x[2](t-1)]);
                            [   0    2 x (t - 1) ]
                            [           2        ]
(%o3)                       [                    ]
                            [ x (t)       0      ]
                            [  2                 ]
 *
 * @endcode
 * @param f scalar or vector function with entries in \f$\mathcal{K}\f$.
 * @param v list of s variables \f$v=[v_1,\ldots,v_s]\f$
 * @return List with the partial derivatives of \f$f\f$ with respect to the elements
   of the list \f$x=[var_1,\ldots,var_s]\f$ :
 *  \f[
    grad(f,x)=\left[\frac{\displaystyle \partial f}{\displaystyle \partial var_1},\ldots,
    \frac{\displaystyle \partial f}{\displaystyle \partial var_s}\right]
  \f]
 * @note If the list coincides with all the state variables, the result is the gradient of the function.
 */
/*v list */ grad(
/*v function */ f,
/*v list */ v) := block([],
    if matrixp(f)
      then (
            if length(transpose(f))>1 then error("only scalar or vector functions")
            else apply('matrix,makelist(grad(f[i,1],v),i,1,length(f)))
            )
      else map(lambda([u],diff(f,u)), v)
)$

/**
 * @brief Factorizes a 1-form.
 * @author L.A. Marquez-Martinez
 *
 * Given a 1-form \f$\omega=\sum_{i=1}^s a_i\,dz_i=a dz\f$ returns the matrices
 * \f$ a = [a_1,\ldots,a_s]\f$ and \f$dz=[d(z_1),\ldots,d(z_s)]\f$.
 *
 * <b>Usage</b>
 * @code
(%i1) load("sac.mc")$
(%i2) _d(x[1](t-3)*u(t-2));
(%o2)          x (t - 3) del(u(t - 2)) + u(t - 2) del(x (t - 3))
                1                                      1
(%i3) dot_fact(%);
                [             2             3 ]  [ del(u(t))  ]
(%o3)          [[ x (t - 3) _D   u(t - 2) _D  ], [            ]]
                [  1                          ]  [ del(x (t)) ]
                                                 [      1     ]
 *
 * @endcode
 *
 * @param w 1-form
 * @return list of matrices [a,dz] such that w = a dz
 * @see _d
 */
/*v matrix_list */ dot_fact(
/*v 1-form */ w):=block([rv,e,dl,pw,cf,tmp],
  if matrixp(w) or listp(w)
    then error("wrong type of argument: 1-form expected"),
  rv:showratvars(w),
  dl:[],
  for e in rv do if not freeof(del,e) then push(e,dl),
  cf:map(lambda([u],ratcoef(w,u)),dl),
  if (dl=[]) then return(w),
  pw:map(rel_shift,dl),
  cf:map(lambda([i,j],_D^j*i),cf,pw),
  dl:map(lambda([i,j],tshift(i,-j)),dl,pw),
  tmp:unique(dl),
  if length(tmp) # length(dl) /* in case we have dx(t-i) and dx(t-j) */
    then return(dot_fact(matrix(cf).transpose(matrix(dl))))
    else return([matrix(cf),transpose(matrix(dl))])
  )$

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
   * (%i1) load("sac.mc")$
   *
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
 * @param k (optional) number of times to derivate
 * @return \f$h^{(k)}\f$, following the trajectories of S.
 * @see lie
 * @note won't complain with p-forms, but will give wrong results.
 */

 /*v sys d_dt(function f, system S, int k){}  */

 /*v // */ d_dt([args])
   :=block([f,S,k,l,c,d,cdt,ddt],
    if length(args)<2 then error ("expected at least 2 arguments"),
    f:pop(args),
    S:pop(args),
    if args=[] then k:1 else k:pop(args),
    if k>1 then f:d_dt(f,S,k-1),
    if not freeof(del,f) then
      ( /* case of 1-forms */
        [c,d]:dot_fact(f),
        cdt:map(lambda([u],d_dt(u,S)),c),
        ddt:matrixmap(lambda([u],_d(d_dt(inpart(u,1),S))),d),
        return(c.ddt+cdt.d)
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
