/**
 * @file ncalg.mc
 * @author Araceli Garate and Alejandro Marquez
 * @date May 20, 2017
 * @brief Definitions for non-commutative algebra.
 *
 */
 /**
  * @brief Decomposes a polynomial in monomials
  *
  * Given a polynomial p, it returns a list of the nonzero coefficients
  * \f$\mathcal{K}[\delta)\f$.
  *
  * <b>Usage</b>
  * @code
  *
  * (%i2) coefpow(_D^3+2);
  * (%o2)                        [[2, 1], [0, 3]]
  * (%i3) coefpow(matrix([_D^3+3],[_D]));
  *                         [ 3 ]  [ 0 ]  [ 1 ]
  * (%o159)               [[[   ], [   ], [   ]], [0, 1, 3]]
  *                         [ 0 ]  [ 1 ]  [ 0 ]
  *
  * @endcode
  * @param p Polynomial with scalar or matrix coefficients.
  * @return [c,p], where c is a list of coefficients of _D, and p the list of the corresponding powers
  */
 coefpow(
 /*v polynomial */ p) := block([hp1,cero,lista,l,i,pot],
   hp1 : hipow(p,_D),
   if matrixp(p) then cero:zeromatrix(length(p),length(transpose(p))) else cero:0,
   if operatorp(p,"+") then lista:args(p) /* even -a-b returns (+) */
   else lista:[p],
   l:[],
   pot:[],
   for i:0 thru hp1 do
       ( ci:ratcoef(p,_D,i),
         if ci#cero then (
            l:append(l,[ci]),
            pot:append(pot,[i])
            )
        ),
   return( [l,pot] )
  )$
