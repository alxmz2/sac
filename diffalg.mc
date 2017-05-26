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
