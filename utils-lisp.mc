/**
 * @file utils-lisp.mc
 * @author Alejandro Marquez
 * @date May 18, 2017
 * @brief Lisp functions
 *
 * This file Several utilities written in LISP defined in utils.lisp
 */
/**
 * @brief Protects a symbol so it cannot be assigned to anything.
 *
 * Some symbols are reserved for the software.  Assigning them would lead
 * to weird and hard-to-track bugs.  This command avoids this problem by
 * reserving the word.  So far we reserve the words <tt>t, del, true, false</tt>.
 * @param s symbol
 * @returns protected symbol s
 * @see unprotect
 */
protect (
/*v symbolp */ x
):=block([]
)$
/**
 * @brief Unprotects a symbol.
 *
 * Removes the protection of the argument x.
 * @param x protected symbol
 */
unprotect(
/*v symbolp */ x
):=block([]
)$
/**
 * @brief tshift([args])
 *
 * Shifts a function in time.
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
 * @param [args] any valid function, polynomial, matrix, p-form, or a list of these elements
 * @param s (optional) positive integer
 * @returns the same argument, with all functions of t shifted by 1 unit, or by s units if a second argument is given.
 */
tshift([_args],s):=
  block([_largs,_sol],
  modedeclare(__largs,number),
  return(fix_derivative_shift(_sol))
)$
