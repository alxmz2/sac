defstruct(sys (f,g,fg,h,statevar,controlvar))$
load("utils.lisp")$
load("operators.mc")$
load("ncalg.mc")$
load("diffalg.mc")$
load("system-utils.mc")$
protect(t)$
texput(_D,"\\delta")$
_wvarsta:x$
_wvarcon:u$
_wpert:q$
