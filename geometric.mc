/**
 * @file geometric.mc
 * @author A. Marquez, E. García
 * @date May 20, 2017
 * @brief Definitions for geometric tools for  non-commutative algebra.
 *
 */
 /**
  * @brief
  * @author L.A. Marquez-Martinez
  *
  * Computes the Lie derivative
  *
  *
  * <b>Usage</b>
  * @code
  (%i1) load("sac.mc")$
  (%i2) f:matrix([x[2](t-2)*u(t)],[x[1](t)]);
                                [ x (t - 2) u(t) ]
                                [  2             ]
  (%o2)                         [                ]
                                [     x (t)      ]
                                [      1         ]
  (%i3) s:systdef(f,[x,u])$
  (%i4) lie(x[2](t-3)*_D,s,2);
  (%o5)                        x (t - 5) u(t - 3) _D
                                2
  *
  * @endcode
  *
  * @param h polynomial in \f$\mathcal{K}[\delta)\f$.
  * @param S Refernce system
  * @param k (optional) number of times to derivate.
  * @return k-th time-derivative of h following the trajectory of S
  */
/*v polynomial lie(polynomial h, system S, int k):= { (*/
/*v // */ lie([pars]) := block([l,k,h,S,p],
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
        return(sum(matrix(gradient(h,tshift(S@statevar,i))).tshift(S@fg,i),i,0,p))
      )
  )$
