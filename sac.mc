defstruct(sys (affine,f,dF,g,fg,h,n,m,p,statevar,controlvar,outputvar,taumax,hk))$
defstruct(PSQ (P,S,Q))$
load("lrats")$
load("utils.lisp")$
load("operators.mc")$
load("ncalg.mc")$
load("diffalg.mc")$
load("geometric.mc")$
load("analysis.mc")$
load("system-utils.mc")$
protect(t)$
texput(_D,"\\delta")$
