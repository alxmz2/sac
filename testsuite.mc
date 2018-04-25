load("sac.mc")$

testfnc(name,testcond):=if testcond then print(name,": passed") else error(name,": failed")$

/* define some test data */
flist:[1,a,x[1](t),u(t),x[1](t)*u(t-1),x[1](t-1)^2,x[2](t+3)*diff(u[1](t-2),t)]$
dlist:[x[1](t),(x[1](t)+x[2](t))*del(x[1](t)),x[2](t-1)*u(t)^2*del(x[2](t))]$
w:sin(x[1](t-1))*u(t-3)*q*x[2]*del(diff(u(t-2),t,2),diff(u(t-8),t))$

/* utils */

testfnc(showtvars, 
  showtvars(w)=[u(t-3),x[1](t-1),del('diff(u(t-8),t,1),'diff(u(t-2),t,2))])$

testfnc(showalltvars,
  showalltvars(w)=[u(t-3),x[1](t-1),u(t-2),u(t-8)])$


/* operators */

df1:maplist(_d,flist)$
df:flatten([df1,maplist(_d,dlist)])$

testfnc(_d, 
   df=[0,0,del(x[1](t)),del(u(t)),u(t-1)*del(x[1](t))+x[1](t)*del(u(t-1)),2*x[1](t-1)*del(x[1](t-1)),
       x[2](t+3)*del('diff(u[1](t-2),t,1))
       +'diff(u[1](t-2),t,1)*del(x[2](t+3)),del(x[1](t)),-del(x[1](t),x[2](t)),
       u(t)^2*del(x[2](t-1),x[2](t))-2*x[2](t-1)*u(t)*del(x[2](t),u(t))])$

pop(flist)$pop(flist)$pop(df1)$pop(df1)$
testfnc(antider,flist=antider(df1))$

testfnc("antider (detect not closed argument)",(errcatch(antider(u(t)*del(u(t-1))))=[]))$

