
(defun $protect (x)
    (if (symbolp x)
        (putprop x 'neverset 'assign)
        (merror "Argument to protect must be a symbol")
    )
)

(defun $unprotect (x)
    (if (symbolp x)
        (remprop x 'assign)
        (merror "Argument to unprotect must be a symbol")
    )
)
