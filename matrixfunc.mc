/**
 * @file matrixfunc.mc
 * @author Alejandro Marquez
 * @date May 18, 2017
 * @brief General matrix functions.
 *
 * Functions for matrices with elements in \f$\mathcal{K}[\delta)\f$
 */
/**
 * @brief Computes (left) Ore and Bezout polynomials.
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
