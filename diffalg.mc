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
 * Given a list of variables , it will return the partial
 * derivative of function \f$f(\cdot)\f$ with respect to them.
 *
 * <b>Usage</b>
 * @code
 * (%2) grad(((x[2](t-1))^2),[x[1](t-1),x[2](t-1)]);
 * (%o2)                [0, 2 x (t - 1)]
 *                             2
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
 * @warning In the present context for time-delay systems, the gradient of a function is not defined.
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
