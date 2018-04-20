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
     antider [URL="\ref antider"];
     coefpow [URL = "\ref coefpow"];
     d [ label="_d" URL = "\ref _d"];
     ddt [label="d/dt" URL="\ref d_dt"];
     dot_fact [URL="\ref dot_fact"];
     euclid  [label="Euclid" URL="\ref euclid"];
     find_el [URL="\ref find_el"];
     find_max_idx [URL="\ref find_max_idx"];
     gradfnc [ URL = "\ref gradfnc"];
     hk [URL="\ref hk"];
     isaccessible [URL="\ref is_accessible"];
     isclosed [URL="\ref is_closed"];
     isintegrable [URL="\ref is_integrable"];
     isobservable [URL="\ref is_observable"];
     left_kernel  [URL="\ref left_kernel"];
     Lie [label="Lie" URL="\ref lie"];
     lorebez [URL="\ref lorebez"];
     maxd [URL="\ref maxd"];
     ncgrad [URL="\ref ncgrad"];
     ncinverse [URL="\ref ncinverse"];
     ncProd [ label="*^" URL = "\ref infix"];
     ncRowRank [URL="\ref ncrow_rank"];
     nctriangularize [URL="\ref nctriangularize"];
     pdeg [URL="\ref p_degree"];
     psqswap [URL="\ref psqswap"];
     rel_shift [URL="\ref rel_shift"];
     systdef [URL="\ref systdef"];
     showalltvars [URL="\ref showalltvars"];
     showtvars [URL="\ref showtvars"];
     tshift [URL = "\ref tshift"];
     wedge [URL="\ref infix"];

     edge [arrowhead=open];
     antider -> {antider};
     coefpow -> {ncProd};
     d -> {ddt,isintegrable};
     ddt -> {ddt,isobservable};
     dot_fact -> {antider,d,ddt,ncProd,pdeg,wedge};
     euclid -> {lorebez};
     find_max_idx -> {systdef};
     gradfnc -> {d,dot_fact,gradfnc,Lie,ncgrad,systdef};
     hk -> {isaccessible};
     isclosed -> {antider,isclosed};
     Lie -> {Lie};
     lorebez -> {nctriangularize};
     maxd -> {Lie,maxd,ncgrad,rel_shift,systdef};
     ncgrad -> {isobservable};
     ncRowRank -> {isobservable};
     nctriangularize -> {left_kernel,ncinverse,ncRowRank};
     ncProd -> {ddt,dot_fact,euclid,lorebez,ncgrad,ncinverse,nctriangularize};
     pdeg -> {ddt,isclosed,isintegrable};
     psqswap -> {nctriangularize};
     rel_shift -> {dot_fact};  
     showalltvars -> {maxd};   
     showtvars -> {antider,d,dot_fact,rel_shift,showalltvars,showtvars};
     tshift -> {ddt,dot_fact,euclid,Lie,ncgrad,ncProd,systdef};
     wedge -> {isintegrable}
 }
 \enddot
 *
 */

/*! \page ayc Analysis and control routines
 *
 * In progress....
 */
