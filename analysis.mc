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
 * @author A. Garate-Garcia, R. Cuesta-Garcia, and L.A. Marquez-Martinez
 *
 *
 *
 *
 * <b>Usage</b>
 * @code
 * (%i1) load("sac.mc")$
 * (%i2) maxd(x[3](t-1)*u(t-4));
 * (%o2)                        4
 * @endcode
 *
 * @param f function, matrix, or p-form
 * @return int maximum delay found in f.
 * @warning not tested for p>1.
 */
/*v int */ maxd(
/*v var */      f ) := block([_fn,_mdf,_md,_lenfnc,_fl],
    _fl:inflag,
    inflag:true,
    if atom(f) then(
        return(0)
    ),
    _mdf:-inf,
    _lenfnc:length(f),
    if(nterms(f)=1 and _lenfnc>1 and inpart(f,2)=t)
       then return(max(_mdf,maxd(inpart(f,1)))),
    for _m:1 thru _lenfnc do(
        _fn:inpart(f,_m),
        if (_fn=t)and(_m=2) then(
            _md:-inpart(f,1),
            _mdf:max(_mdf,_md)
        )else(
            _mdf:max(_mdf,maxd(_fn))
        )
    ),
    inflag:_fl,
    _mdf
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
 * @return int maximum delay found in f.
 * @warning not tested for p>1.
 */
/*v int */ rel_shift(
/*v var */      f ) := apply(min, map(maxd,showratvars(f))
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
 * if it is just one word, then you can just do @b this.
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
