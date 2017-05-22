/**
 * @file ncalg.mc
 * @author Araceli Garate and Alejandro Marquez
 * @date May 20, 2017
 * @brief Definitions for non-commutative algebra.
 *
 */
/**
  *
  * @brief Decomposes a polynomial
  *
  * Given a polynomial \f$p[\delta)=\sum_i p_i\delta^i\f$, it returns a list of the nonzero coefficients,
  * \f$p_i\in\mathcal{K}[\delta)\f$, and another list with the corresponding exponent of \f$p_i\f$
  * in ascending order.
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
  * @param p Polynomial \f$p[\delta)=\sum_i p_i\delta^i\f$ with scalar or matrix coefficients.
  * @return \f$[c,e]\f$, where \f$c=[p_i\mid p_i\neq0]\f$, and \f$e=[i\mid p_i\neq0]\f$, in ascending order.
  */
coefpow(
/*v polynomial */ p) := block([hp1,cero,lista,l,i,c,e],
 hp1 : hipow(p,_D),
 if matrixp(p) then cero:zeromatrix(length(p),length(transpose(p))) else cero:0,
 if operatorp(p,"+") then lista:args(p) /* even -a-b returns (+) */
 else lista:[p],
 c:[],
 e:[],
 for i:0 thru hp1 do
     ( ci:ratcoef(p,_D,i),
       if ci#cero then (
          c:append(c,[ci]),
          e:append(e,[i])
          )
      ),
 return( [c,e] )
)$
/**
 * @brief Computes (left) Ore and Bezout polynomials.
 * @author Araceli Garate
 *
 * Let \f$a,b\in\mathcal{K}[\delta)\f$. The Ore polynomials \f$\alpha[1],\beta[1]\f$
 * and the Bezout polynomials \f$\alpha[2],\beta[2]\f$
 * Example of usage:
 * @code
 * systdef(eq,matrix([x[1](t)*u(t-1)],[u(t-2)]))$
 * eq_F;
 * eq_dF;
 * @endcode
 * @param a polynomial
 * @param b \f$\in\mathcal{K}[\delta)\f$
 * @return  Matrix \f$P[\delta)\f$ such that
 \f[
   P[\delta) \begin{bmatrix} a\\ b\end{bmatrix} = \begin{bmatrix} glcd(a,b)\\ 0\end{bmatrix}
 \f]
 * @note Something to note.
 * @warning Under development...
 */
lorebez(/*v polynomial */ _a,
/*v polynomial */ _b) := block([_ans,_k,_qr,_V,_T],
  _b:rat(_b),
  if (_b=0) then error("Division by 0"),
  _T:ident(2),
  _a:rat(_a),
/*    if ((hipow(_a,_D)=0) and (_a#0)) then  return(matrix([1/_a,0],[-_b/_a,1])),
    if  (hipow(_b,_D)=0) then (
                   return(matrix([0,1/_b],[1,-_a/_b]))), */
   _V:matrix([_a],[_b]),
    while (_V[2,1]#0) do (
           _qr:euclid(_V[1,1],_V[2,1]),
           _T:matrix([0,1],[1,-_qr[1,1]])*^_T,
           _V[1,1]:_V[2,1],
           _V[2,1]:_qr[1,2]
           ),
    _V[1,1]:rat(_V[1,1],showratvars(_V[1,1])),
    _k:ratcoef(_V[1,1],_D,hipow(_V[1,1],_D)),
    /*_T[1]:factor(1/(_k)*_T[1]),*/			/*creo que tengo la sol. usando makelist(1/(_k)*_T[1,i],i,1,2)*/
    _T[1]:factor(makelist(1/(_k)*_T[1,i],i,1,2)),
    return(rat(map(factor,_T)))
    )	$
