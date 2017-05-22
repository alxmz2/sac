
/**
 * @file operators.mc
 * @author Araceli Garate and Alejandro Marquez
 * @date May 20, 2017
 * @brief Miscellaneous operator definitions
 *
 */
 /**
  * @fn infix("*^")
  * @brief  Operator *^
  *
  * Defines the non-commutative operator *^. This allows to multiply polynomials in
  * \f$\mathcal{K}[\delta)\f$.
  *
  * <b>Usage</b>
  * p1 *^ p2
  * @code
  * (%i2) x(t-1)*_D *^ x(t-2);
  * (%o2)              x(t - 3) x(t - 1) _D
  * @endcode
  * @param p1 polynomial with scalar or matrix coefficients.
  * @param p2 polynomial with scalar or matrix coefficients, or p-form.
  * @return   non-commutative product
  * @note \f$\delta\f$ is written as _D
  * @warning Bugs found! does not work with polynomial *^ p-forms.
  */
  infix("*^")$
p1 *^p2 := block([_pf2,_mp1,_mp2,_tmp,_argsp2,_m,_k,_lpol2,_j,_lrm,_lcm,_pol2,_hp1],
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
 * @brief Protects a symbol so it cannot be assigned to anything.
 *
 * Some symbols are reserved for the software.  Assigning them would lead
 * to weird and hard-to-track bugs.  This command avoids this problem by
 * reserving the word.  So far we reserve the words <tt>t, del, true, false</tt>.
 * @param s symbol
 * @returns protected symbol s
 * @see unprotect
 */
/*v protect ( symbolp  x ):=block([] */
/*v )$ */
/**
 * @brief Unprotects a symbol.
 *
 * Removes the protection of the argument x.
 * @param x protected symbol
 */
/*v unprotect( symbolp x ):=block([] */
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
 * @param [args] any valid function, polynomial, matrix, p-form, or a list of these elements
 * @param s (optional) positive integer
 * @returns the same argument, with all functions of t shifted by 1 unit, or by s units if a second argument is given.
 */
/*v tshift(var [_args], integer s):= block([_largs,_sol], */
/*v )$ */
