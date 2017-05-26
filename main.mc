/*! \mainpage Symbolic Analysis and Control package
 *
 * This project aims to provide symbolic tools for people working with nonlinear time-delay systems.
 * It implements routines for:
 * -  manipulating polynomials (scalar or matrices) in \f$\mathcal{K}[\delta)\f$.
 *
 * -  analysis of nonlinear time-delay systems.
 *
 * The \ref jerarquia shows all the provided routines and their relationships.
 *
 */

/*! \page jerarquia General hierarchy diagram
 *
 * Functions and their relationships.
 *
 \dot
 digraph sac {
     rankdir=LR;
     tshift [URL = "\ref tshift"];
     grad [ URL = "\ref grad"];
     d [ label="_d" URL = "\ref _d"];
     ncProd [ label="*^" URL = "\ref infix(\"*^\")"];
     coefpow [URL = "\ref coefpow"];
     lorebez [ URL="\ref lorebez" ];
     euclid  [label="Euclid" URL="\ref euclid"];
     SwapRow [URL = "\ref swaprow"];
     ddt [label="d/dt" URL="\ref d_dt"];
     dot_fact [URL="\ref dot_fact"];
     rel_shift [URL="\ref rel_shift"];
     Lie [label="Lie" URL="\ref lie"];
     maxd [URL="\ref maxd"];
     systdef [URL="\ref systdef"];
     edge [arrowhead=open];
     ncProd -> {coefpow,tshift} ;
     protect;
     unprotect;
     d -> grad;
     grad -> grad;
     euclid -> tshift;
     ddt -> {dot_fact,lie,rel_shift,tshift,ddt};
     dot_fact -> {dot_fact,rel_shift,tshift,ncProd};
     lorebez -> euclid;
     Lie -> {grad,Lie,maxd,tshift};
     maxd -> maxd;
     rel_shift -> maxd;
     systdef -> {tshift,grad,maxd}
 }
 \enddot
 *
 */

/*! \page ayc Analysis and control routines
 *
 * Hk
 */
