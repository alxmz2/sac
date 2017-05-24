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
 * (%i21) maxd(x[3](t-1)*u(t-4));
 * (%o21)                                 4
 * @endcode
 *
 * @param f function, matrix, or p-form
 * @return int maximum delay found in f.
 * @warning not tested for p>1.
 */
/*v int */ maxd(
 /*v var */ f ) := block([_fn,_mdf,_md,_lenfnc,_fl],
    _fl:inflag,
    inflag:true,
    if atom(f) then(
        return(0)
    ),
    _mdf:0,
    _lenfnc:length(f),
    if(nterms(f)=1 and _lenfnc>1 and inpart(f,2)=t) then(
        flag1:1,
        return(max(_mdf,maxd(inpart(f,1))))
    ),
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
);