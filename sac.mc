/* Symbolic Analysis and Control package (SAC)  */
/* 2018 © Copyright CICESE, All rights reserved */
/* Licensed under GPL License v3.0 */
defstruct(sys (affine,f,dF,g,fg,h,n,m,p,statevar,controlvar,outputvar,taumax,hk))$
defstruct(PSQ (P,S,Q))$
declare(del,antisymmetric)$
if (get(vect,version)=false) then load("vect")$
load("utils.lisp")$
protect(t)$
protect(del)$
protect(_D)$
texput(_D,"\\delta")$

/* non-commutative product */
infix("*^",128,127)$ /* binding power to have more precedence than normal product, but less than exponentiation */
"*^"(p1,p2) := block([pf2,mp1,mp2],
    pf2: not freeof(del,p2),
    mp1: matrixp(p1),
    mp2: matrixp(p2),
    if not freeof(del,p1)
       then error("wrong arguments: left operand cannot be a p-form"),
    pol1:coefpow(p1),
    if (mp1 and mp2) then oper:"."
                     else oper:"*",
    if pf2 then return(ratsimp(sum(apply(oper,[pol1[1][i],tshift(p2,pol1[2][i])]),i,1,length(pol1[1]))))
            else return(ratsimp(sum(apply(oper,[pol1[1][i],tshift(p2,pol1[2][i])*_D^pol1[2][i]]),i,1,length(pol1[1]))))
 )$

/* Integrates [a list of] 1-forms */
antider( [L] ) := block([F,vlist,lv,c,d],
  L:flatten(L),
  if not(isclosed(L)) then error("argument is not a [list of ] closed 1-form[s]"),
  if length(L)>1
    then return(maplist(antider,L))
    else (
         L:pop(L), /* converts [dw] to dw */
         if (L=0) then return(0),
         [c,d]:dotfact(L,0),
         vlist:showtvars([c,matrixmap(args,d)]),
         lv:length(vlist),
         F:subst(makelist(vlist[i]=concat(x,i),i,1,lv),flatten(args(c))),
         scalefactors(makelist(concat(x,i),i,1,lv)),
         return(subst(makelist(concat(x,i)=vlist[i],i,1,lv),potential(F)))
         )
  )$

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
       if ci#cero
         then (
            c:append(c,[ci]),
            e:append(e,[i]) )
      ),
  return( [c,e] )
 )$

/* differential operator */
_d(f):=block([rv,tmp,tmp2,suma,i],
 if (not freeof(del,f))
    then( suma:0,
          tmp: dotfact(f,0),
          tmp2:matrixmap(_d,tmp[1]),
          for i:1 thru length(tmp[1]) do
             suma:suma+wedge(tmp2[i,1],tmp[2][i,1]),
          return(suma))
    else (rv:showtvars(f),
          if (rv=[])
             then return(0)
             else return(matrix(gradfnc(f,rv)).transpose(matrix(map(del,rv)))) )
 )$

/* time-derivative along a vector field */
ddt([ar]) := block([f,S,k,l,c,d,cdt,ddt,vl,vu,a],
   if length(ar)<2 then error ("expected at least 2 arguments"),
   f:pop(ar),
   S:pop(ar),
   if ar=[] then k:1 else k:pop(ar),
   if k=0 then return(f),
   if k<0 then error("k cannot be negative"),
   if k>1 then f:ddt(f,S,k-1),
   if matrixp(f) then return(matrixmap(lambda([u],ddt(u,S)),f)),
   if not freeof(del,f) then
     ( /* case of 1-forms */
     if (pdegree(f)>1)
        then error("not implemented for p-forms")
        else
           ( [c,d]:dotfact(f,0),
             cdt:matrixmap(lambda([u],ddt(u,S)),c),
             ddt:matrixmap(lambda([u],_d(ddt(inpart(u,1),S))),d),
             return(ratsimp(c*^ddt+cdt*^d))
           )
     ),
   vl:showtvars(f),
   l:length(vl),
   vu:[],
   for i:1 thru l do (
     a:pop(vl),
     if not (member(tshift(a,-relshift(a)),S@statevar)) then push(a,vu) ),
   l:length(vu),
   return(ratsimp(lie(f,S)+sum(ratcoef(f,vu[i])*diff(vu[i],t),i,1,l)))
   )$

/* Factorizes a p-form. */
dotfact( [w] ):=block([dfact,rv,e,dl,pw,cf,tmp,flag],
  if (last(w)#0) then dfact:1,
  w:pop(w),
  flag:false,
  if matrixp(w)
     then ( if ( length(transpose(w))>1 or apply("or",map(lambda([u],((u#[0]) and freeof(del,u))),args(w))))
              then error("wrong type of argument: 1-form or column vector of 1-forms expected")
              else flag:true ),
  if listp(w) then error("wrong type of argument: 1-form or column vector of 1-forms expected"),
  rv:showtvars(w),
  /* this replaces _d(a+b) by _d(a)+_d(b) but not _d(a,b) */
  for e in rv do
     if not freeof(del,e)
        then if (length(args(e))=1)
                then w:ratsubst(_d(inpart(e,1)),e,w),
  /* select all differential terms */
  dl:sort(sublist(showtvars(w),lambda([u], is(not freeof(del,u))))),
  cf:map(lambda([u],ratcoef(w,u)),dl),
  if (dl=[]) then return([w,matrix([1])]),
  if (dfact=1) then (
      pw:map(lambda([u],max(relshift(u),0)),dl),
      cf:map(lambda([i,j],_D^j*i),cf,pw),
      dl:map(lambda([i,j],tshift(i,-j)),dl,pw)),
  tmp:unique(dl),
  if length(tmp) # length(dl) /* in case we have dx(t-i) and dx(t-j) */
     then return(dotfact(transpose(matrix(cf)).transpose(matrix(dl))))
     else if flag=true
             then (
                  tmp:pop(cf),
                  for i in cf do tmp:addcol(tmp,i),
                  return([tmp,transpose(matrix(dl))])
                  )
             else return([transpose(matrix(cf)),transpose(matrix(dl))])
  )$

/* Euclid's division. */
euclid( a, b) := block([bm,l,p,pdb,pdr,q,r],
  if (b=0) then error("division by 0"),
  q:0,
  r:rat(a,_D),
  b:rat(b,_D),
  pdb:hipow(b,_D),
  bm:ratcoef(b,_D,pdb),
  pdr:hipow(r,_D),
  l:pdr-pdb,
  while ( (l>=0) and (r # 0)) do(
    p:(ratcoef(r,_D,pdr)/tshift(bm,l))*^(_D^l),
    r:r-p*^b,
    q:q+p,
    pdr:hipow(r,_D),
    l:pdr-pdb ),
  return(matrix([q,r]))
 )$

/* Finds the positions where a specific element occurs within a matrix */
findel([pars]) := block([M,e,idx,n,m,L],
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

/* find max index of a subscripted variable s in expression f */
findmaxidx(f,s):=block([l,maxidx:minf],
   l:showalltvars(f),
   for e in l do(
       if subvarp(op(e))
         then (if op(op(e))=s then maxidx:max(maxidx,part(op(e),1)))
         else (if op(e)=s then maxidx:max(maxidx,1))),
   return(maxidx)
   )$

/* Computes the gradient of a function */
gradfnc( f, v) := block([],
  if matrixp(f)
    then (
          if length(transpose(f))>1 then error("only scalar or vector functions")
          else apply('matrix,makelist(gradfnc(f[i,1],v),i,1,length(f))) )
    else map(lambda([u],diff(f,u)), v)
  )$

/* Hk submodules */
hk( s ):=block([g,Hi,dx,hk,j],
  g:copy(s@g),
  gi:copy(g),
  dx:transpose(matrix(map(del,s@statevar))),
  hk:[],
  dimhk:s@m,
  for i:1 thru s@n do (
     if i=s@n then j:inf else j:i+1,
     Hi:leftkernel(g)*^dx,
     if Hi = 0 then (hk:append(hk,[[inf,0]]),return()),
     if matrixp(Hi) then hk:append(hk,[[j,flatten(args(matrixmap(num,Hi)))]])
                    else hk:append(hk,[[j,[num(Hi)]]]),
     gi:s@dF*^gi-ddt(gi,s),
     g:addcol(g,gi)),
  s@hk:hk
  )$

/* Tests a system for generic observability. */
isobservable( S ):= is(ncrowrank(ncgrad(apply(addrow,makelist(ddt(S@h,S,i),i,0,S@n)),S@statevar))=S@n)$

/*Checks the strong accessibility condition for a given system. */
isaccessible( S  ):= if not listp(S@hk) then hk(S) else is(last(s@hk)[2]=0)$

/* Checks whether a 1-form or list of 1-forms is integrable or not. */
isintegrable( [L] ) := block([ww,dw,dww],
  L:flatten(L),
  if not(unique(maplist(pdegree,L))=[1]) then error("argument must be a 1-form or list of 1-forms"),
  if (length(L)>1)
     then(
           ww:lreduce(wedge,L),
           dw:maplist(_d,L),
           dww:unique(maplist(lambda([u],wedge(ww,u)),dw)),
           return(is(dww=[0]))
          )
     else  return(is(wedge(_d(L[1]),L[1])=0))
  )$

/* Checks whether a 1-form or list of 1-forms are closed or not. */
isclosed( [L] ) := block([],
  L:flatten(L),
  if (length(L)>1)
     then return(is(unique(maplist(isclosed,L))=[true]))
     else ( L:first(L),
           return (is((L=0) or ( (pdegree(L)=1) and (_d(L))=0) ) ))
  )$

/* Returns a basis for the left kernel of a matrix */
leftkernel( M ):=block([ans,rnk,n,P],
   ans:nctriangularize(M),
   rnk:apply(min,matrix_size(M)),
   P:args(ans@P),
   for i:1 thru rnk do (
      if ans@S[i,i]=0 then return(),
      pop(P) ),
   if P=[] then return(zeromatrix(1,rnk))
           else return(apply(matrix,P))
   )$

/* Lie derivative */
lie([pars]) := block([i,l,k,h,S,p],
   l:length(pars),
   if (l<2) then error("expects at least 2 arguments"),
   h:pop(pars),
   S:pop(pars),
   if pars=[] then k:1 else k:pop(pars),
   if k=0 then return(h),
   if k<0 then error("k cannot be negative"),
   if k>1 then h:lie(h,S,k-1),
   if matrixp(h)
      then
        h:matrixmap(lambda([u],lie(u,S)),h)
      else (
        if (l=3) and (k>1) then h:lie(h,S,k-1),
        p:maxd(h),
        return(sum(matrix(gradfnc(h,tshift(S@statevar,i))).tshift(S@fg,i),i,min(0,p),max(0,p)))  )
   )$

/* Computes (left) Ore and Bezout polynomials. */
lorebez( a, b ) := block([ans,k,qr,V,T],
  b:rat(b),
  if (b=0) then error("Division by 0"),
  T:ident(2),
  a:rat(a),
  V:matrix([a],[b]),
  while (V[2,1]#0) do (
         qr:euclid(V[1,1],V[2,1]),
         T:matrix([0,1],[1,-qr[1,1]])*^T,
         V[1,1]:V[2,1],
         V[2,1]:qr[1,2]
         ),
  V[1,1]:rat(V[1,1],showratvars(V[1,1])),
  k:ratcoef(V[1,1],_D,hipow(V[1,1],_D)),
  T[1]:factor(makelist(1/(k)*T[1,i],i,1,2)),
  return(rat(map(factor,T)))
  )$

/* Finds maximal delay in one expression. */
maxd( f ) := -apply(min,flatten(subst([t=0],maplist(args,showalltvars(f))))) $

/* Gradient using _D operator */
ncgrad( f,  v ):= sum(gradfnc(f,tshift(v,i))*_D^i,i,0,maxd(f))$

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
ncrowrank( M ):=block([ans,rnk],
   ans:nctriangularize(M),
   rnk:apply(min,matrix_size(M)),
   for i:1 thru rnk do (
   if ans@S[i,i]=0 then return(rnk:i-1)),
   return(rnk)
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
        ans@S[j]:tmp[2] ) )
    ), /* next i */
  if ((ans@S[limit,limit] # 0) and (hipow(ans@S[limit,limit],_D)=0)) then (
    ans@P[limit]:ans@P[limit]/ans@S[limit,limit],
    ans@S[limit]:ans@S[limit]/ans@S[limit,limit]),
  ans@S:expand(ans@S),
  return(ans)
  )$

/* Returns the degree of a p-form. */
pdegree( v ):=block([d1],
  d1:transpose(dotfact(v)[2])[1],
  if freeof(del,d1)
    then return(0)
    else (d1:unique(maplist(length,d1)),
  if length(d1)>1 then error("wrong argument: mixed p-forms")
                  else return(d1[1]))
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

/* Finds the relative shift in one expression. */
relshift( f ) := apply(min, map(maxd,showalltvars(f)))$

/* Auxiliary functions for showtvars and showalltvars */
makealltvarslist(e):=if not atom(e)
     then ((if matchvar(e) then push(e,mylist)
                           else map(makealltvarslist,args(e))))$
matchvar(e):=not atom(e)
  and (istvar(e) or subvarp(op(e)))  /* istvar replaces (diff(args(e),t)=[1])  */
  and not freeof('t,args(e))$
istvar(e):=block([ar], ar:args(e), if (ar=['t]) then return(true),
  if (length(ar)>1 or atom(ar[1]) ) then return(false),
  is((op(ar[1])="+") and member('t,args(ar[1]))))$
maketvarslist(e):=if not atom(e)
   then (if op(e) = nounify(diff)
            then push(e,mylist)
            else if op(e) = nounify(del)
                    then ( if istvar(args(e)) or length(args(e))>1
                              then push(e,mylist)
                              else mylist:flatten([showratvars(subst(_d(args(e)[1]),e,e)),mylist]) )
                    else (if matchvar(e)
                       then push(e,mylist)
                       else map(maketvarslist,args(e))))$

/* This function lists all time-dependent variables. */
showalltvars(e) := block ([mylist : [] ], makealltvarslist(e), unique(mylist))$

/* This function is similar to showalltvars, but keeping del() or diff() operators. */
showtvars(e) := block ([mylist : [] ], maketvarslist(e), unique(mylist))$

/* System definition */
systdef( [pars]
  ) := block([vars,name,tmp,varlist,s,ss],
  name:new(sys),           /* create new system structure        */
  if length(pars)=1
     then (                /* if only one parameter is given,    */
          eq:pop(pars),    /* it is f or [f,h]                   */
          vars:[x,u,y]     /* set default variables              */
     )
     else (
          eq:pop(pars),    /* f, or [f,h]                        */
          vars:pop(pars)   /* [state, control, output]           */
     ),
  if (listp(eq)) then (
     if length(eq) # 2 then error("two arguments are expected"), /* [f,h] */
     name@fg:pop(eq),
     name@h:pop(eq),
     name@p:length(name@h),
    if not matrixp(name@h) then name@h:matrix([name@h])
    )
   else (
    name@fg:eq,
    name@h:0,
    name@p:0
    ),
  name@n:length(name@fg),
  tmp:pop(vars),
  name@statevar:makelist(tmp[i](t),i,1,name@n),
  tmp:pop(vars), /* control var name */
  if vars # []
    then
      (
       name@outputvar:pop(vars),
       if name@h #0 then name@outputvar:makelist(name@outputvar[i](t),i,1,length(name@h))
      ),
  /* find maximal delay */
  name@taumax:maxd(name@fg),
  /* substitute u(t-i) by u[1](t-i) */
  name@fg:subst(makelist(tmp(t-i)=tmp[1](t-i),i,0,name@taumax),name@fg),
  /* find size of control input */
  name@m:findmaxidx(name@fg,tmp),
  /* compute dF */
  name@dF:sum(gradfnc(name@fg,tshift(name@statevar,i))*_D^i,i,0,name@taumax),
  /* compute g(_D) */
  if (name@m >0 ) then (
     name@controlvar:makelist(tmp[i](t),i,1,name@m),
     name@g:sum(gradfnc(name@fg,tshift(name@controlvar,i))*_D^i,i,0,name@taumax)
  ) else name@g:zeromatrix(name@n,1),
  /* system is not affine if g depends on u */
  name@affine: is(findmaxidx(name@g,tmp)<0),
  if name@affine then
      name@f:subst(map(lambda([u],u=0),flatten(makelist(tshift(name@controlvar,i),i,0,name@taumax))),name@fg),
  /* compute hk's  */
  hk(name),
  return(name)
  )$

/* wedge product */
wedge([ar]):= block([l,u,v,nu,nv,wprod],
  ar:flatten(ar),
  l:length(ar),
  if l<2 then error("error: at least two arguments are required"),
  u:pop(ar),
  if l=2 then v:pop(ar)
         else v:treereduce(wedge,ar),
  if (pdegree(u)=0) then return(u*^v),
  if (pdegree(v)=0) then return(v*^u),
  u:dotfact(u,0),
  if not(freeof(_D,u)) then u:dotfact(transpose(u[1])*^u[2],0),
  v:dotfact(v,0),
  if not(freeof(_D,v)) then v:dotfact(transpose(v[1])*^v[2],0),
  nu:length(u[1]),
  nv:length(v[1]),
  wprod:0,
  for i:1 thru nu do
     for j:1 thru nv do
       wprod: wprod
             +u[1][i,1]*v[1][j,1]
             * apply(del, flatten([args(u[2][i,1]),args(v[2][j,1])])),
  return(wprod)
  )$
