/**
 * @file analysis.mc
 * @author L.A. Marquez-Martinez A. Garate-Garcia, R. Cuesta-Garcia.
 * @date May 2017
 * @brief Routines for analysis of NLTDS
 *
 */

/* Finds maximal delay in one expression. */
 maxd( f ) := block([],
   return(-apply(min,flatten(subst([t=0],maplist(args,showalltvars(f))))))
)$

/* Finds the relative shift in one expression. */
 rel_shift( f ) := apply(min, map(maxd,showalltvars(f)))$

/* Tests a system for generic observability. */
 is_observable( S ):=block([],
    is(ncrow_rank(ncgrad(apply(addrow,makelist(d_dt(S@h,S,i),i,0,S@n)),S@statevar))=S@n)
)$

/*Checks the strong accessibility condition for a given system. */
is_accessible( s  ):=block([],
   if not listp(s@hk) then hk(s),
   is(last(s@hk)[2]=0)
)$

/* Checks whether a 1-form or list of 1-forms is integrable or not. */
is_integrable( [L] ) := block([ww,dw,dww],
L:flatten(L),
if not(unique(maplist(p_degree,L))=[1]) then error("argument must be a 1-form or list of 1-forms"),
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
is_closed( [L] ) := block([],
L:flatten(L),
if (length(L)>1)
   then return(is(unique(maplist(is_closed,L))=[true]))
   else (L:first(L),
         if ((L=0) or ( (p_degree(L)=1) and (_d(L))=0) )
           then return(true)
           else return(false)
        )
)$
