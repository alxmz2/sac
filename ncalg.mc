/**
 * @file ncalg.mc
 * @author A. Marquez, A. Garate
 * @date May 20, 2017
 * @brief Definitions for non-commutative algebra.
 *
 */
/**
  * @brief Decomposes a polynomial
  * @author A. Marquez
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
  * @return [c,e] List where \f$c=[p_i\mid p_i\neq0]\f$, and \f$e=[i\mid p_i\neq0]\f$, in ascending order.
  */
/*v list */ coefpow(
/*v polynomial */ p) := block([hp1,cero,lista,l,i,c,e],
 hp1 : hipow(p,_D),
 if matrixp(p)
    then cero:zeromatrix(length(p),length(transpose(p)))
    else cero:0,
 if (symbolp(p) or numberp(p))
    then lista:[p]
    elseif operatorp(p,"+")
        then lista:args(p)
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
 * @brief Euclid's division.
 * @author A. Garate-Garcia and L.A. Marquez-Martinez
 *
 * Given two poynomials a, b\f$\in\mathcal{K}[\delta)\f$, perform the Euclid's division
 * to find q, r\f$\in\mathcal{K}[\delta)\f$ such that  a=q b+r, where the polynomial
 * degree of pol.d(r) is strictly less than pol.d(b).
 *
 * <b>Usage</b>
 * @code
 * (%i1) load("sac.mc")$
 * (%i2) a:x(t)*_D^2+2$
 * (%i3) b:u(t)*_D-1$
 * (%i4) d:euclid(a,b);
 *             [ x(t) _D        x(t)       2 u(t) u(t - 1) + x(t) ]
 * (%o4)/R/    [ -------- + -------------  ---------------------- ]
 *             [ u(t - 1)   u(t - 1) u(t)      u(t) u(t - 1)      ]
 * (%i5) fullratsimp(d[1][1]*^b+d[1][2]);
 *                                        2
 * (%o5)                           x(t) _D  + 2
 * @endcode
 *
 * @param a polynomial \f$\in\mathcal{K}[\delta)\f$
 * @param b polynomial \f$\in\mathcal{K}[\delta)\f$
 * @return M polynomial matrix M=[q,r] such that a = q b + r, and deg(r)<deg(b).
 */
/*v matrix */ euclid(
                    /*v polynomial */ a,
                    /*v polynomial */ b) := block([_bm,_l,_p,_q,_r],
    if (b=0) then error("division by 0"),
    _q:0,
    _r:rat(a,_D),
    b:rat(b,_D),
    _pdb:hipow(b,_D),
    _bm:ratcoef(b,_D,_pdb),
    _pdr:hipow(_r,_D),
    _l:_pdr-_pdb,
    while ( (_l>=0) and (_r # 0)) do(
      _p:(ratcoef(_r,_D,_pdr)/tshift(_bm,_l))*^(_D^_l),
      _r:_r-_p*^b,
      _q:_q+_p,
      _pdr:hipow(_r,_D),
      _l:_pdr-_pdb
      ),
     return(matrix([_q,_r]))
   )$
/**
  * @brief Swaps rows of a matrix
  *
  * Description of  @c this which is equivalent to <tt>this</tt>
  * We can also include text verbatim,
  * @verbatim like this@endverbatim
  *
  * <b>Usage</b>
  * @code
  *
  * @endcode
  *
  * if it is just one word, then you can just do @b this.
  * @param M Descripción chida del first parameter of the function.
  * @param r1 The second one, which follows @p param1.
  * @param r2 The second one, which follows @p param1.
  * @return   Explicación de lo que regresa la función
  *  \f[
    x^2
     \f]
  * @see Box_The_Last_One
  * @see http://website/
  * @note Something to note.
  * @warning Warning.
  */
/*v matrix swaprow( matrix M, int r1, int r2 ) {} */

/**
 * @brief Computes (left) Ore and Bezout polynomials.
 * @author A. Garate
 *
 * Let \f$a,b\in\mathcal{K}[\delta)\f$. We call \f$\alpha,\beta\f$
 * Ore polynomials if they satisfy the left-Ore condition:
 \f[ \alpha\,a + \beta \, b = 0 \f]
 * and we call them Bezout polynomials if they satisfy
 \f[ \alpha\,a + \beta \, b = gcld(a,b) \f]
 * where glcd(a,b) stands for greatest left common divisor of (a,b).
 *
 * <b>Usage</b>
 * @code
 * (%i1) load("sac.mc")$
 * (%i2) a:_D^2+1$
 * (%i3) b:x(t)$
 * (%i4) lorebez(a,b);
                       [             1            ]
                       [ 0          ----          ]
                       [            x(t)          ]
 *  (%o4)/R/           [                          ]
                       [                   2      ]
                       [      x(t - 2) + _D  x(t) ]
                       [ 1  - ------------------- ]
                       [         x(t - 2) x(t)    ]
 * (%i5) lorebez(a,b)*^matrix([a],[b]);
 *                                 [ 1 ]
 * (%o5)                           [   ]
 *                                 [ 0 ]
 * @endcode
 * <b>Hierarchy</b>
 *
 * @param a polynomial \f$\in\mathcal{K}[\delta)\f$
 * @param b polynomial \f$\in\mathcal{K}[\delta)\f$
 * @return  Matrix \f$P[\delta)\f$ such that
 \f[
   P[\delta) \begin{bmatrix} a\\ b\end{bmatrix} = \begin{bmatrix} glcd(a,b)\\ 0\end{bmatrix}
 \f]
 *
 */
/*v matrix */ lorebez(
                      /*v polynomial */ a,
                      /*v polynomial */ b) := block([_ans,_k,_qr,_V,_T],
  b:rat(b),
  if (b=0) then error("Division by 0"),
  _T:ident(2),
  a:rat(a),
  _V:matrix([a],[b]),
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
