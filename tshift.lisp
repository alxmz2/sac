
(defun $replace-t (olista)
  (let (lista m mm)
    (setq lista (reverse olista))
    (setq m (second lista))
    (if (and (listp m) (member '$T m))
        (setq m '$T))     ;return if not a time derivative
    (setq mm (car lista))
    (setq lista (rest (rest lista)))
    (setq lista (reverse (append (list mm m ) lista)))
   )
)

(defun $fix_derivative_shift (expr)
    (cond ((listp expr)
        (if (equal (car expr) '(%DERIVATIVE SIMP))
            ($replace-t expr)
            (mapcar #'$fix_derivative_shift expr)))
        (t expr))
)

;; funcion traducida por maxima 
;; funcion original:
;; tshift([_args]):=
;;   block([_largs,_sol],
;;    modedeclare(__largs,number),
;;    _largs:length(_args),
;;    if (_largs=1) then (_sol:subst(t=t-1,_args[1]))
;;                   else
;;          (_sol:subst(t=t-_args[2],_args[1])),
;;   return(fix_derivative_shift(_sol)))$
(PROGN
  (DEFPROP $TSHIFT T TRANSLATED)
  (ADD2LNC '$TSHIFT $PROPS)
  (DEFMTRFUN ($TSHIFT $ANY MDEFINE T NIL) ($_ARGS)
      (DECLARE (SPECIAL $_ARGS))
      ((LAMBDA ($_LARGS $_SOL)
         (DECLARE (SPECIAL $_SOL $_LARGS))
         (PROG ()
           (SETQ $_LARGS (MFUNCTION-CALL $LENGTH $_ARGS))
           (COND
             ((LIKE $_LARGS 1)
              (SETQ $_SOL
                    (SIMPLIFY
                        (MFUNCTION-CALL $SUBSTITUTE
                            (SIMPLIFY
                                (LIST '(MEQUAL) (TRD-MSYMEVAL $T '$T)
                                      (ADD* (TRD-MSYMEVAL $T '$T) -1)))
                            (MAREF $_ARGS 1)))))
             (T (SETQ $_SOL
                      (SIMPLIFY
                          (MFUNCTION-CALL $SUBSTITUTE
                              (SIMPLIFY
                                  (LIST '(MEQUAL) (TRD-MSYMEVAL $T '$T)
                                        (ADD* (TRD-MSYMEVAL $T '$T)
                                         (*MMINUS (MAREF $_ARGS 2)))))
                              (MAREF $_ARGS 1))))))
           (RETURN
             (SIMPLIFY (MFUNCTION-CALL $FIX_DERIVATIVE_SHIFT $_SOL)))))
       '$_LARGS '$_SOL)))
