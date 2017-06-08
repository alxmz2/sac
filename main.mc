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
     ncProd [ label="*^" URL = "\ref infix"];
     coefpow [URL = "\ref coefpow"];
     ddt [label="d/dt" URL="\ref d_dt"];
     lorebez [ URL="\ref lorebez" ];
     euclid  [label="Euclid" URL="\ref euclid"];
     dot_fact [URL="\ref dot_fact"];
     rel_shift [URL="\ref rel_shift"];
     Lie [label="Lie" URL="\ref lie"];
     maxd [URL="\ref maxd"];
     systdef [URL="\ref systdef"];
     nctriangularize [URL="\ref nctriangularize"];
     ncinverse [URL="\ref ncinverse"];
     find_max_idx [URL="\ref find_max_idx"];
     edge [arrowhead=open];
     ncProd -> {coefpow,tshift} ;
     protect;
     unprotect;
     d -> grad;
     grad -> grad;
     ddt -> {dot_fact,Lie,rel_shift,tshift,ddt};
     dot_fact -> {dot_fact,rel_shift,tshift,ncProd};
     euclid -> tshift;
     lorebez -> euclid;
     Lie -> {grad,Lie,maxd,tshift};
     maxd -> maxd;
     rel_shift -> maxd;
     systdef -> {tshift,grad,maxd,find_max_idx};
     ncinverse -> {ncProd,nctriangularize};
     nctriangularize -> {lorebez,psqswap,infix}
 }
 \enddot
 *
 */

/*! \page ayc Analysis and control routines
 *
 * Hk
 */
