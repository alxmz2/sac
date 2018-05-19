/*
 * file ncalg.mc
 * author A. Marquez, A. Garate
 * date May 20, 2017
 * brief Definitions for non-commutative algebra.
*/

/* Decomposes a polynomial in coefficents and exponents of _D */
coefpow( p ) := block([hp1,cero,c,e],
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

/* Euclid's division. */
euclid(  a, b) := block([_bm,_l,_p,_q,_r],
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

/* Computes (left) Ore and Bezout polynomials. */
lorebez( a, b ) := block([_ans,_k,_qr,_V,_T],
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

/* Finds the positions where a specific element occurs within a matrix */
find_el([pars]) := block([M,e,idx,n,m,L],
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

/* Computes the upper triangular form  */
nctriangularize( M ):=block([tmp,PQ,Mp,l,L,Ll,ans,m,n,p,limit],
   [n,m]:matrix_size(M),
   ans:new(PSQ),    /* creates structure(P,S,Q) */
   ans@S:copy(M),
   ans@P:ident(n),
   ans@Q:ident(m),
   limit:min(m,n),
   for i:1 thru limit do (
     /* iterate following the diagonal  */
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

/* Swap matrices of a PSQ structure */
psqswap([pars]) :=block( [ans,r1,r2,c1,c2],
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

/* @brief Computes the inverse matrix */
ncinverse( M ):= block([tf,m,n],
  [m,n]:matrix_size(M),
  if not matrixp(M) then error("argument must be a matrix"),
  if m # n then error("cannot invert a non-square matrix"),
  tf:nctriangularize(M),
  for i:1 thru n do if tf@S[i,i] # 1 then error("non unimodular matrix"),
  genmatrix(lambda([i,j],if j>i then -tf@S[i,j] else tf@S[i,j]),n,n)*^tf@P*^tf@Q
  )$

/* computes the row rank over K[_D) */
ncrow_rank( M ):=block([ans,rnk],
    ans:nctriangularize(M),
    rnk:apply(min,matrix_size(M)),
    for i:1 thru rnk do (
       if ans@S[i,i]=0 then return(rnk:i-1)
    ),
    return(rnk)
)$

/* Returns a basis for the left kernel of a matrix */
left_kernel( M ):=block([ans,rnk,n,P],
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
