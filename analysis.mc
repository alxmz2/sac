/**
 * @file analysis.mc
 * @author A. Garate-Garcia, R. Cuesta-Garcia, and L.A. Marquez-Martinez
 * @date May 2017
 * @brief Routines for analysis of NLTDS
 *
 *
 */

/**
 * @brief Finds maximal delay in one expression.
 * @author L.A. Marquez-Martinez
 *
 * Given a function, matrix, or p-form, it finds the maximum delay found in any
 * time-dependant variable.
 *
 * <b>Usage</b>
 * @code
  (%i1) load("sac.mc")$
  (%i2) maxd(x[3](t-1)*u(t-4));
  (%o2)                        4
  (%i3) maxd(
 * @endcode
 *
 * @param f function, matrix, or p-form
 * @return int maximum delay found in f.
 * @warning not tested for p>1.
 */
/*v int */ maxd(
/*v var */      f ) := block([_fn,_mdf,vlist],
   vlist:flatten(maplist(lambda([u],if atom(u) then 0 else if (u[1]='t) then u else  args(u)),showtvars(f))),
   push(0,vlist),
   _mdf:apply(min,subst([t=0],flatten(maplist(lambda([u],if atom(u) then 0 else if (u[1]='t) then u else  args(u)),vlist))))
)$
/**
 * @brief Finds the relative shift in one expression.
 * @author L.A. Marquez-Martinez
 *
 * the relative shift of a function \f$f(z_tau)\f$ is defined as the maximal
 * forward time shift such that the resulting function is still causal.
 * Mathematically, \f[
 \mbox{rel\_shift}(f(z_\tau)) = f(t)= max\{k\in\mathbb{Z}^+\ \mid\
 d(f(t+k)\in span_{\mathcal{K}[\delta)}\{dz\})
 \f]
 *
 *
 * <b>Usage</b>
 * @code
 * (%i1) load("sac.mc")$
 * (%i2) rel_shift(x[3](t-1)*u(t-4));
 * (%o2)                                 1
 * @endcode
 *
 * @param f function, matrix, or p-form
 * @return int minimum delay found in f.
 */
/*v int */ rel_shift(
/*v var */      f ) := block([vlist,dlist],
   vlist:showtvars(f),
   dlist:sublist(vlist,lambda([u],not (freeof(del,u)))),
   if not(dlist=[]) then vlist:flatten([vlist,map(args,dlist)]),
   apply(min, map(maxd,vlist))
)$

/**
 * @brief Tests a system for generic observability.
 *
 * Given a system \f$\Sigma\f$ with output \f$y\f$, it checks the generic
 * observability condition given by:
 *
 * \f[
    rank_{\mathcal{K[\delta)}}\frac{\partial [y\ \dot y \cdots y^{(n-1)}]}{\partial x}=n
   \f]
 *
 * <b>Usage</b>
 * @code
 * (%i1) load("sac.mc")$
 * (%i2) f:matrix([x[2](t-2)*u[1](t)],[u[2](t-3)])$
 * (%i3) h:x[1](t-1)$
 * (%i4) S:systdef([f,h],[x,u,y])$
 * (%i5) is_observable(S);
 * (%o5)                                true
 * @endcode
 *
 * @param s system
 * @return True if the observability condition is satisfied, false otherwise.
 * @see systdef
 * @warning not tested for time-delay systems
 */
/*v boolean */ is_observable(
/*v system  */ s ):=block([],
    is(ncrow_rank(ncgrad(apply(addrow,makelist(d_dt(s@h,s,i),i,0,2)),s@statevar))=s@n)
)$

/**
 * @brief Checks the strong accessibility condition for a given system.
 *
 * Tests if the system \f$\Sigma\f$ satisfies the strong accessibility condition:
 *
 * \f[
    \mathcal{H}_\infty = 0
   \f]
 * <b>Usage</b>
 * @code
 * (%i1) load("sac.mc")$
 * (%i2) f:matrix([x[2](t-2)*u[1](t)],[u[2](t-3)])$
 * (%i3) S:systdef(f)$
 * (%i4) is_accessible(S);
 * (%o4)                                true
 * @endcode
 *
 * @param s System
 * @return  True if the given system is accessible, false otherwise.
 * @see systdef
 * @see is_observable
 */
/*v boolean */ is_accessible(
/*v system  */ s  ):=block([],
   if not listp(s@hk) then hk(s),
   is(last(s@hk)[2]=0)
)$

/**
 * @brief Checks whether a 1-form or list of 1-forms is integrable or not.
 *
 * This routine checks if a 1-form or list of 1-forms is integrable using
 * the Frobenious theorem.
 *
 * <b>Usage</b>
 * @code
 * (%i1) load("sac.mc")$
 *
 * @endcode
 *
 * @param L 1-form or list of 1-forms.
 * @return   True if \f$span_\mathcal{K}\{l\}\f$ is integrable.
 * @todo The integrability condition is still to be implemented for time-delay systems.
 * @warning  It is not valid for time-delay systems yet.
 */
/*v boolean */ is_integrable(
/*v list    */ L ) := block([ww,dw,dww],
if not(listp(L)) then L:[L],
if not(unique(maplist(p_degree,L))=[1]) then error("argument must be a 1-form or list of 1-forms"),
if (length(L)>1)
   then(
         ww:lreduce("~^",L),
         dw:maplist(_d,L),
         dww:unique(maplist(lambda([u],ww~^u),dw)),
         return(is(dww=[0]))
        )
   else  return(is((_d(L[1])~^L[1])=0))
)$

/**
 * @brief Checks whether a 1-form or list of 1-forms are closed or not.
 *
 * This routine checks if a 1-form or list of 1-forms are closed, that is,
 *\f[ d(\omega_i)=0,
  \f]
 *
 * <b>Usage</b>
 * @code
 * (%i1) load("sac.mc")$
 *
 * @endcode
 *
 * @param L 1-form or list of 1-forms \f$[\omega_1,\ldots,\omega_s]\f$.
 * @return   True if \f$d\omega_i=0,\ \forall\,i=1,\ldots,s\f$.
 * @note Any argument which is not a 1-form, will return false.
 */
/*v boolean */ is_closed(
/*v list    */ L ) := block([],
if listp(L)
   then return(maplist(is_closed,L))
   else (
         if ( (p_degree(L)=1) and (_d(L))=0)
           then return(true)
           else return(false)
        )
)$
