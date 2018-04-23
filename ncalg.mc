/**
 * @file ncalg.mc
 * @author A. Marquez, A. Garate
 * @date May 20, 2017
 * @brief Definitions for non-commutative algebra.
 *
 */
/**
  * @fn coefpow 
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
/*v polynomial */ p) := block([hp1,cero,c,e],
 hp1 : hipow(p,_D),
 if matrixp(p)
    then cero:zeromatrix(length(p),length(transpose(p)))
    else cero:0,
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
 * to find q, r\f$\in\mathcal{K}[\delta)\f$ such that  a=qb+r, where the polynomial
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
  _T[1]:factor(makelist(1/(_k)*_T[1,i],i,1,2)),
  return(rat(map(factor,_T)))
  )$

/**
 * @brief Finds the positions where a specific element occurs within a matrix
 * @author L.A. Marquez-Martinez
 *
 * Given a matrix @c M, an element @c e, and an index @c idx, returns a list of all the pairs \f$[i,j]\f$ such that \f$M[i,j]=k\f$, and
 * \f$i,j \geq r\f$.  If no index is given, it will be taken as idx=1.
 *
 *
 * <b>Usage</b>
 * @code
 (%i1) load("sac.mc")$
 (%i2) M:genmatrix(lambda([i,j],(2*i-j)),4,10);
                [ 1  0  - 1  - 2  - 3  - 4  - 5  - 6  - 7  - 8 ]
                [                                              ]
                [ 3  2   1    0   - 1  - 2  - 3  - 4  - 5  - 6 ]
 (%o2)          [                                              ]
                [ 5  4   3    2    1    0   - 1  - 2  - 3  - 4 ]
                [                                              ]
                [ 7  6   5    4    3    2    1    0   - 1  - 2 ]
 (%i3) find_el(M,-3);
 (%o3)                    [[1, 5], [2, 7], [3, 9]]
 (%i4) find_el(M,1,2);
 (%o4)                    [[2, 3], [3, 5], [4, 7]]
 (%i5) find_el(M,-3,2,6);
 (%o5)                        [[2, 7], [3, 9]]
 *
 * @endcode
 *
 * @param M matrix
 * @param e any valid expression
 * @param idr (optional) initial row
 * @param idc (optional) initial column
 * @return list of pairs \f$[i,j]\f$ satisfying \f$M[i,j]=e, i,j\geq idx\f$.
 */
/*v list find_el(matrix M, expr e, int idr, int idc) */
/* // */ find_el([pars]
      ):=block([M,e,idx,n,m,L],
   M:pop(pars),
   if not(matrixp(M)) then error("first argument must be a matrix"),
   e:pop(pars),
   if pars=[] then idr:1 else idr:pop(pars),
   if pars=[] then idc:1 else idc:pop(pars),
   [n,m]:matrix_size(M),
   L:[],
   for j:idc thru m do
      for i:idr thru n do if (M[i,j]=e) then push([i,j],L),
   reverse(L)
)$

/**
 * @brief Computes the upper triangular form
 * @author L.A. Marquez-Martinez
 *
 * Returns a structure with 3 elements: P, S, and Q, such that for the given matrix M,
   P M Q = S
  with S in upper-triangular form, and the elements of the main diagonal are normalized.
 *
 * <b>Usage</b>
 * @code
 * (%i1) load("sac.mc")$
 *
 * @endcode
 *
 * @param M matrix \f$\in\mathcal{K}[\delta)\f$
 * @return list with the Smith form of \f$M,\ P,\ Q,\  P^{-1},\ Q^{-1}\f$.
 *
 * @warning en desarrollo, posibles bugs.
 * @todo optimize algorithm,
 */
/*v matrix_list */ nctriangularize(
/*v matrix      */ M
                  ):=block([tmp,PQ,Mp,l,L,Ll,ans,m,n,p,limit],

   /* creates structure(P,S,Q) */
   [n,m]:matrix_size(M),
   ans:new(PSQ),
   ans@S:copy(M),
   ans@P:ident(n),
   ans@Q:ident(m),
   limit:min(m,n),
   for i:1 thru limit do (  /* iterate following the diagonal  */
     /* finds powers of _D. Infinity for 0  */

     /* finds nonzero polynomial of lowest hipow */
     L:locate_matrix_entry(ans@S,i,i,n,m,lambda([u],if u=0 then inf else hipow(u,_D)),'min),
     if L=false then L:[i,i],
     ans:psqswap(ans,[i,i],L),
     for j:i+1 thru n do (
        if ans@S[j,i] # 0 then (
          PQ:lorebez(ans@S[i,i],ans@S[j,i]),
          tmp:PQ*^apply(matrix,[ans@P[i],ans@P[j]]),
          ans@P[i]:tmp[1],
          ans@P[j]:tmp[2],
          tmp:PQ*^apply(matrix,[ans@S[i],ans@S[j]]),
          ans@S[i]:tmp[1],
          ans@S[j]:tmp[2]
        )
     )
   ), /* next i */
   if ((ans@S[limit,limit] # 0) and (hipow(ans@S[limit,limit],_D)=0)) then (
      ans@P[limit]:ans@P[limit]/ans@S[limit,limit],
      ans@S[limit]:ans@S[limit]/ans@S[limit,limit]),
   ans@S:expand(ans@S),
   return(ans)
)$

/**
* @brief Swap matrices of a PSQ structure
* @author L.A. Marquez-Martinez
*
* Given a PSQ structure, and 2 lists of indexes [r1,c1] and [r2,c2], swap rows and columns
of the elements P, S, Q to have P', S', and Q' such that
S'= P' S Q', where S' is obtained by swapping rows r1 and r2, and columns c1 and c2.
*
* <b>Usage</b>
* @code
* (%i1) load("sac.mc")$
*
* @endcode
*
* @param s Structure PSQ
* @param l1 list [r1,c1], the first row and column pair.
* @param l2 list [r2,c2], the second row and column pair.
* @return psq_struct with swapped rows and columns.
* @see
* @note
* @warning
*/
/*v psq_struct psqswap( psq_struct s, list l1, list l2) {} */
/*v // */ psqswap([pars]
            ) :=block( [ans,r1,r2,c1,c2],
    ans:pop(pars),
    [r1,c1]:pop(pars),
    [r2,c2]:pop(pars),
    if [r1,c1]=[r2,c2] then return(ans),
    ans@P:rowswap(ans@P,r1,r2),
    ans@S:rowswap(ans@S,r1,r2),
    ans@S:columnswap(ans@S,c1,c2),
    ans@Q:columnswap(ans@Q,c1,c2),
    return(ans)
  )$

  /**
  * @brief Computes the inverse matrix
  * @author A. Garate-Garcia, R. Cuesta-Garcia, and L.A. Marquez-Martinez
  *
  * This routine computes the inverse of a matrix with entries in \f$\mathcal{K}[\delta)\f$
    if it exists. Otherwise, signals an error.
  *
  *
  * <b>Usage</b>
  * @code
(%i1) load("sac.mc")$
(%i2) M:matrix([1+_D,-_D],[_D,1-_D])*^matrix([x[2](t),u(t)],[1,-u(t-2)]);
      [ (x (t - 1) - 1) _D + x (t)    (u(t - 3) + u(t - 1)) _D + u(t)   ]
      [   2                   2                                         ]
(%o2) [                                                                 ]
      [   (x (t - 1) - 1) _D + 1    (u(t - 3) + u(t - 1)) _D - u(t - 2) ]
      [     2                                                           ]
(%i3) ncinverse(M);
      [   (u(t) + u(t - 2)) _D - u(t - 2)  u(t) + (u(t) + u(t - 2)) _D ]
      [ - -------------------------------  --------------------------- ]
      [        u(t) + u(t - 2) x (t)          u(t) + u(t - 2) x (t)    ]
      [                         2                              2       ]
(%o3) [                                                                ]
      [        1 + (x (t) - 1) _D             x (t) + (x (t) - 1) _D   ]
      [              2                         2        2              ]
      [       ---------------------         - ----------------------   ]
      [       u(t) + u(t - 2) x (t)           u(t) + u(t - 2) x (t)    ]
      [                        2                               2       ]
(%i4) ncinverse(M)*^M;
                                   [ 1  0 ]
(%o4)                              [      ]
                                   [ 0  1 ]
(%i5)

  * @endcode
  *
  * @param M matrix with entries in \f$\mathcal{K}[\delta)\f$
  * @return Inverse of M, if it exists.
  * @see nctriangularize
  */
/*v matrix */  ncinverse(
/*v matrix */ M
    ):= block([tf,m,n],
  [m,n]:matrix_size(M),
  if not matrixp(M) then error("argument must be a matrix"),
  if m # n then error("cannot invert a non-square matrix"),
  tf:nctriangularize(M),
  for i:1 thru n do if tf@S[i,i] # 1 then error("non unimodular matrix"),
  genmatrix(lambda([i,j],if j>i then -tf@S[i,j] else tf@S[i,j]),n,n)*^tf@P*^tf@Q
  )$

  /**
  * @brief
  * @author L.A. Marquez-Martinez
  *
  * Returns the row rank over \f$\mathcal{K}[\delta)\f$ of M.
  *
  *
  * <b>Usage</b>
  * @code
(%i1) load("sac.mc")$
(%i2) ncrow_rank(matrix([_D^2,1],[_D,1+_D],[_D,1-_D]));
(%o2)                                 2
(%i3) M:matrix([u(t),u(t-1),_D],[u(t-1)*_D,u(t-2)*_D,_D^2]);
                       [    u(t)       u(t - 1)    _D  ]
(%o3)                  [                               ]
                       [                             2 ]
                       [ u(t - 1) _D  u(t - 2) _D  _D  ]
(%i4) ncrow_rank(M);
(%o4)                                 1
  * @endcode
  *
  * @param M matrix
  * @return row rank over \f$\mathcal{K}[\delta)\f$ of M
  */
/*v int */ ncrow_rank(
/*v matrix */ M
           ):=block([ans,rnk],
    ans:nctriangularize(M),
    rnk:apply(min,matrix_size(M)),
    for i:1 thru rnk do (
       if ans@S[i,i]=0 then return(rnk:i-1)
    ),
    return(rnk)
)$
/**
* @brief
* @author L.A. Marquez-Martinez
*
* Returns a basis for the left kernel of a matrix with entries in \f$\mathcal{K}[\delta)\f$.
*
*
* <b>Usage</b>
* @code
(%i1) load("sac.mc")$
(%i2) M:matrix([_D],[u(t)],[1+_D]);
                                  [   _D   ]
                                  [        ]
(%o2)                             [  u(t)  ]
                                  [        ]
                                  [ _D + 1 ]
(%i3) left_kernel(M);
                        [              _D            ]
                        [ 1       - --------       0 ]
                        [           u(t - 1)         ]
(%o3)                   [                            ]
                        [      u(t - 1) + u(t) _D    ]
                        [ 0  - ------------------  1 ]
                        [        u(t - 1) u(t)       ]
* @endcode
*
* @param M matrix
* @return base of the left-kernel of M.
*/
/*v matrix */ left_kernel(
/*v matrix */ M
             ):=block([ans,rnk,n,P],
   ans:nctriangularize(M),
   rnk:apply(min,matrix_size(M)),
   P:args(ans@P),
   for i:1 thru rnk do (
      if ans@S[i,i]=0 then return(),
      pop(P)
    ),
   if P=[] then return(zeromatrix(1,rnk))
           else return(apply(matrix,P))
)$
