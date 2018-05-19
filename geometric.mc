/**
 * @file geometric.mc
 * @author A. Marquez, E. García
 * @date May 20, 2017
 * @brief Definitions for geometric tools for  non-commutative algebra.
 *
 */
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
        if (l=3) and (k>1)
        then h:lie(h,S,k-1),
        p:maxd(h),
        return(sum(matrix(gradfnc(h,tshift(S@statevar,i))).tshift(S@fg,i),i,min(0,p),max(0,p)))
      )
)$
