/**
 * @file diffalg.mc
 * @author A. Marquez, A. Garate, and R. Cuesta
 * @date May 20, 2017
 * @brief Differential algebra routines
 *
 * This file contains the functions related to differential algebra for
 * working with p-forms (\f$p\in{Z}^+\f$).
 */

/* Returns the degree of a p-form. */
p_degree( v ):=block([d1],
  d1:transpose(dot_fact(v)[2])[1],
  if freeof(del,d1)
     then return(0)
     else (d1:unique(maplist(length,d1)),
           if length(d1)>1 then error("wrong argument: mixed p-forms")
                           else return(d1[1])
          )
)$

/* Computes the gradient of a function */
gradfnc( f, v) := block([],
    if matrixp(f)
      then (
            if length(transpose(f))>1 then error("only scalar or vector functions")
            else apply('matrix,makelist(gradfnc(f[i,1],v),i,1,length(f)))
            )
      else map(lambda([u],diff(f,u)), v)
)$

/* Gradient using _D operator */
ncgrad( f,  v ):=block([],
     tmax:maxd(f),
     for i:0 thru tmax do (
        return(sum(gradfnc(f,tshift(v,i))*_D^i,i,0,tmax))
     )
)$

/* Factorizes a p-form. */
dot_fact( [w] ):=block([dfact,rv,e,dl,pw,cf,tmp,flag],
  if (last(w)#0) then dfact:1,
  w:pop(w),
  flag:false,
  if matrixp(w) then
    ( if ( length(transpose(w))>1 or apply("or",map(lambda([u],freeof(del,u)),args(w))))
        then error("wrong type of argument: 1-form or column vector of 1-forms expected")
        else flag:true
    ),
  if listp(w) then error("wrong type of argument: 1-form or column vector of 1-forms expected"),
  rv:showtvars(w),
  /* this replaces _d(a+b) by _d(a)+_d(b) but not _d(a,b) */
  for e in rv do
     if not freeof(del,e) then
        if (length(args(e))=1) then w:ratsubst(_d(inpart(e,1)),e,w),
  /* select all differential terms */
  dl:sort(sublist(showtvars(w),lambda([u], is(not freeof(del,u))))),
  cf:map(lambda([u],ratcoef(w,u)),dl),
  if (dl=[]) then return([w,matrix([1])]),
  if  (dfact=1) then (
     pw:map(rel_shift,dl),
     cf:map(lambda([i,j],_D^j*i),cf,pw),
     dl:map(lambda([i,j],tshift(i,-j)),dl,pw)
     ),
  tmp:unique(dl),
  if length(tmp) # length(dl) /* in case we have dx(t-i) and dx(t-j) */
    then return(dot_fact(transpose(matrix(cf)).transpose(matrix(dl))))
    else if flag=true
            then (
              tmp:pop(cf),
              for i in cf do tmp:addcol(tmp,i),
              return([tmp,transpose(matrix(dl))])
              )
            else return([transpose(matrix(cf)),transpose(matrix(dl))])
  )$
