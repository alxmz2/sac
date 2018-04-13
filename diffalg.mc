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
 * @brief Returns the degree p of a p-form.
 * @author L.A. Marquez-Martinez
 *
 * Given a p-form \f$\omega\in\mathcal{E}^p\f$, it returns the
 * integer \f$p\f$.
 *
 * <b>Usage</b>
 * @code
 * (%i1) load("sac.mc")$
 * (%i2) p_degree(del(x[1](t),x[3](t-1)));
 * (%o2)            2
 * @endcode
 *
 * @param v p-form
 * @return  p
 * @todo
 *
 */
/*v int    */ p_degree(
/*v p-form */ v):=block([d1],
  d1:transpose(dot_fact(v)[2])[1],
  if freeof(del,d1)
     then return(0)
     else (d1:unique(maplist(length,d1)),
           if length(d1)>1 then error("wrong argument: mixed p-forms")
                           else return(d1[1])
          )
)$

/**
 * @brief Computes the gradient of a function
 * @author L.A. Marquez-Martinez
 *
 * Given a list of variables [\f$v=[v_1,\ldots,v_s]\f$], it will return the partial
 * derivative of function \f$f(\cdot)\f$ with respect to them:
 \f[
  \mbox{gradfnc}(f,v)=\left[\frac{\displaystyle \partial f}{\displaystyle\partial v_1 },\ldots,\frac{\displaystyle \partial f}{\displaystyle\partial v_s }\right]
 \f]
 *
 * <b>Usage</b>
 * @code
(%i1) load("sac.mc")$
(%i2) gradfnc(((x[2](t-1))^2),[x[1](t-1),x[2](t-1)]);
(%o2)                [0, 2 x (t - 1)]
                            2
(%i3) gradfnc(matrix([x[2](t-1)^2],[x[1](t-1)*x[2](t)]),[x[1](t-1),x[2](t-1)]);
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
    gradfnc(f,x)=\left[\frac{\displaystyle \partial f}{\displaystyle \partial var_1},\ldots,
    \frac{\displaystyle \partial f}{\displaystyle \partial var_s}\right]
  \f]
 * @note If the list coincides with all the state variables, the result is the gradient of the function.
 */
/*v list */ gradfnc(
/*v function */ f,
/*v list */ v) := block([],
    if matrixp(f)
      then (
            if length(transpose(f))>1 then error("only scalar or vector functions")
            else apply('matrix,makelist(gradfnc(f[i,1],v),i,1,length(f)))
            )
      else map(lambda([u],diff(f,u)), v)
)$

/**
* @brief Gradient using _D operator
* @author L.A. Marquez-Martinez
*
* Computes gradient of a function \f$\sum_{i=0}^{\tau_M}\partial f/\partial v(t-i)\delta^i\f$.
*
*
* <b>Usage</b>
* @code
(%i1) load("sac.mc")$
(%i2) ncgrad(x[2](t-2)*x[1](t)+x[1](t-2),[x[1](t),x[2](t)]);
                            2                      2
(%o2)                    [_D  + x (t - 2), x (t) _D ]
                                 2          1

*
* @endcode
*
* @param f expression
* @param v variable
* @return matrix with gradient
*/
/*v matrix     */ ncgrad(
/*v expression */ f,
/*v var        */ v
    ):=block([],
     tmax:maxd(f),
     for i:0 thru tmax do (
        return(sum(gradfnc(f,tshift(v,i))*_D^i,i,0,tmax))
     )
)$

/**
 * @brief Factorizes a p-form.
 * @author L.A. Marquez-Martinez
 *
 * Given a p-form \f$\omega=\sum_{i=1}^s a_i\,dz_i=a dz\f$ returns the matrices
 * \f$ a = [a_1,\ldots,a_s]\f$ and \f$dz=[d(z_1),\ldots,d(z_s)]\f$.
 *
 * If the argument is a column vector of p-forms \f$v\f$, it returns the matrices
 * \f$ M \f$ and \f$dz=[d(z_1),\ldots,d(z_s)]\f$ such that \f$\omega=M\,dz\f$.
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
 * @param w p-form
 * @param flag If 0, then factorization does not include _D operator.
 * @return list of matrices [a,dz] such that w = a dz
 * @see _d
 */
/*v matrix_list */ dot_fact(
/*v p-form  */ [w]):=block([dfact,rv,e,dl,pw,cf,tmp,flag],
  if (last(w)#0) then dfact:1,
  w:pop(w),
  flag:false,
  if matrixp(w) then
    ( if ( length(transpose(w))>1 or apply("or",map(lambda([u],freeof(del,u)),args(w))))
        then error("wrong type of argument: 1-form or column vector of 1-forms expected")
        else flag:true
    ),
  if listp(w) then error("wrong type of argument: 1-form or column vector of 1-forms expected"),
  rv:showtvars(w),
  /* this replaces _d(a+b) by _d(a)+_d(b) but not _d(a,b) */
  for e in rv do
     if not freeof(del,e) then
        if (length(args(e))=1) then w:ratsubst(_d(inpart(e,1)),e,w),
  /* select all differential terms */
  dl:sort(sublist(showtvars(w),lambda([u], is(not freeof(del,u))))),
  cf:map(lambda([u],ratcoef(w,u)),dl),
  if (dl=[]) then return([w,matrix([1])]),
  if  (dfact=1) then (
  pw:map(rel_shift,dl),
  cf:map(lambda([i,j],_D^j*i),cf,pw),
  dl:map(lambda([i,j],tshift(i,-j)),dl,pw)),
  tmp:unique(dl),
  if length(tmp) # length(dl) /* in case we have dx(t-i) and dx(t-j) */
    then return(dot_fact(transpose(matrix(cf)).transpose(matrix(dl))))
    else if flag=true
            then (
              tmp:pop(cf),
              for i in cf do tmp:addcol(tmp,i),
              return([tmp,transpose(matrix(dl))])
              )
            else return([transpose(matrix(cf)),transpose(matrix(dl))])
  )$
