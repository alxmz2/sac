/**
 * @file utils-lisp.mc
 * @author Alejandro Marquez
 * @date May 18, 2017
 * @brief Lisp functions
 *
 * This file Several utilities written in LISP defined in utils.lisp
 */
/**
 * @brief tshift([args])
 *
 * Shifts a function in time.
 *
 *<b>Usage</b>
 * @code
 * (%i3) tshift(x(t));
 * (%o3)    x(t - 1)
 * (%i4) tshift(x(t)*u(t-2),4);
 * (%o4)    u(t - 6) x(t - 4)
 * @endcode
 * @param _name Descripción chida del first parameter of the function.
 * @param param2 The second one, which follows @p param1.
 * @return   Explicación de lo que regresa la función
 *  \f[
    |I_2|=\left| \int_{0}^T \psi(t)
             \left\{
                u(a,t)-
                \int_{\gamma(t)}^a
                \frac{d\theta}{k(\theta,t)}
                \int_{a}^\theta c(\xi)u_t(\xi,t)\,d\xi
             \right\} dt
          \right|
  \f]
 * @see Box_The_Last_One
 * @see http://website/
 * @note Something to note.
 * @warning Warning.
 */
tshift([_args]):=
  block([_largs,_sol],
  modedeclare(__largs,number),
  return(fix_derivative_shift(_sol))
)$
