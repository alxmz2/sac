file_search_maxima: append (file_search_maxima,
    ["/home/lmarquez/Documents/s/software/sac/###.mc", "/home/lmarquez/Documents/s/software/sac/###.lisp"])$
defstruct(sys (affine,f,dF,g,fg,h,n,m,p,statevar,controlvar,outputvar,taumax,hk))$
defstruct(PSQ (P,S,Q))$
declare(del,antisymmetric)$
load("lrats")$
if (get(vect,version)=false) then load("vect")$
load("utils.lisp")$
load("operators.mc")$
load("ncalg.mc")$
load("diffalg.mc")$
load("geometric.mc")$
load("analysis.mc")$
load("system-utils.mc")$
protect(t)$
protect(del)$
texput(_D,"\\delta")$
