
/**
 * @file ncalg.mc
 * @author Alejandro Marquez
 * @date May 20, 2017
 * @brief Definitions for non-commutative algebra.
 *
 */

 /**
  * @fn infix("*^")
  * @brief  Operator *^
  *
  * Defines the non-commutative operator *^. This allows to multiply polynomials in
  * \f$\mathcal{K}[\delta)\f$.
  *
  * <b>Usage</b>
  * p1 *^ p2
  * @code
  * (%i2) x(t-1)*_D *^ x(t-2);
  * (%o2)              x(t - 3) x(t - 1) _D
  * @endcode
  * @param p1 polynomial with scalar or matrix coefficients.
  * @param p2 polynomial with scalar or matrix coefficients, or p-form.
  * @return   non-commutative product
  * @warning Bugs found! does not work with polynomial *^ p-forms.
  */
  infix("*^")$ /*v infix("*^"):= */
  p1 *^ p2 := block([_tmp,_argsp2,_m,_k,_lpol2,_j,_lrm,_lcm,_pol2,_hp1],
      _tmp:0,
      _hp1:hipow(p1,_D),
      if _hp1=0 then (
          if (matrixp(p1) and matrixp(p2)) then(
              _tmp:p1.p2
          )else(
              _tmp:p1*p2
          ),
          return(ratsimp(_tmp))
      ),
      /* print("big problem"),*/
      for i:0 thru hipow(p1,_D) do(
          if (i#0 and freeof(del,p2)=true) then(
              if matrixp(p2)=false then(
                  p2:(tshift(p2))*_D
              )else(
                  if i=1 then(
                      _pol2:(tshift(p2))*_D
                  )else(
                      _pol2:(tshift(_pol2))*_D
                  )
              )
          ),
          /*-------------------------------------------------------------*/
          if (i#0 and freeof(del,p2)=false) then(
              if matrixp(p2)=false then(
                  _m:inflag,
                  inflag:false,
                  _argsp2:args(p2),
                  _lpol2:length(_argsp2),
                  inflag:_m,
                  if _lpol2>1 then(
                      p2:0,
                      for _k:1 thru _lpol2 do(
                          if freeof(del,_argsp2[_k])=true then(
                              _argsp2[_k]:(tshift(_argsp2[_k]))*_D
                          )else(
                              _argsp2[_k]:tshift(_argsp2[_k])
                          ),
                          p2:_argsp2[_k]+p2
                      )
                  )else(
                      p2:tshift(p2)
                  )
              )else(
                  _lrm:length(p2),_lcm:length(transpose(p2)),
                  if i=1 then(
                      _pol2:zeromatrix(_lrm,_lcm)
                  ),
                  for _k:1 thru _lcm do(
                      for _j:1 thru _lrm do(
                          if freeof(del,p2[_j,_k])=false then(
                              if i=1 then(
                                  _pol2[_j,_k]:tshift(p2[_j,_k])
                              )else(
                                  _pol2[_j,_k]:tshift(_pol2[_j,_k])
                              )
                          )else(
                              if i=1 then(
                                  _pol2[_j,_k]:tshift(p2[_j,_k])*_D
                              )else(
                                  _pol2[_j,_k]:tshift(_pol2[_j,_k])*_D
                              )
                          )
                      )
                  )
              )
          ),
          /*---------------------------------------------------------------*/
          if (matrixp(p1) and matrixp(p2)) then(
              if i#0 then(
                  _tmp:expand(ratsimp(_tmp+ratcoeff(p1,_D,i)._pol2))
              )else(
                  _tmp:expand(ratsimp(_tmp+ratcoeff(p1,_D,i).p2))
              )
          )else(
              if matrixp(p2) and i#0 then(
                  _tmp:expand(ratsimp(_tmp+ratcoeff(p1,_D,i)*_pol2))
              )else(
                  _tmp:expand(ratsimp(_tmp+ratcoeff(p1,_D,i)*p2))
              )
          )
      ),
      _tmp
  )$
