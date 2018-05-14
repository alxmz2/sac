load("sac.mc")$

testfnc(name,testcond):=if testcond then print(name,": passed") else error(name,": failed")$

/* define some test data */
flist:[1,a,x[1](t),u(t),x[1](t)*u(t-1),x[1](t-1)^2,x[2](t+3)*diff(u[1](t-2),t)]$
dlist:[x[1](t),(x[1](t)+x[2](t))*del(x[1](t)),x[2](t-1)*u(t)^2*del(x[2](t))]$
w:sin(x[1](t-1))*u(t-3)*q*x[2]*del(diff(u(t-2),t,2),diff(u(t-8),t))$

/* utils */

testfnc(tshift,
 flatten([tshift(flist),tshift(dlist,2)])=[1,a,x[1](t-1),u(t-1),u(t-2)*x[1](t-1),x[1](t-2)^2,x[2](t+2)*'diff(u[1](t-3),t,1),x[1](t-2),(x[2](t-2)+x[1](t-2))*del(x[1](t-2)),x[2](t-3)*u(t-2)^2*del(x[2](t-2))])$

testfnc(showtvars,
  showtvars(w)=[u(t-3),x[1](t-1),del('diff(u(t-8),t,1),'diff(u(t-2),t,2))])$

testfnc(showalltvars,
  showalltvars(w)=[u(t-8),u(t-3),u(t-2),x[1](t-1)])$

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

testfnc(is_closed,
        not(is_closed([u(t-1)*del(u(t-1)),u(t)*del(u(t-1))]))
        and is_closed(df1)
        )$

list:[1,a,_D,x(t)*_D,del(x(t)),del(u(t),x(t-1))]$
wedges:[]$
for i in list do for j in reverse(list) do push(wedge(i,j),wedges)$
wedgefn:[-del(x(t-1),u(t)),-a*del(x(t-1),u(t)),-del(x(t-2),u(t-1)),
        -x(t)*del(x(t-2),u(t-1)),-del(x(t-1),u(t),x(t)),0,del(x(t)),
        a*del(x(t)),del(x(t-1)),x(t)*del(x(t-1)),0,-del(x(t-1),u(t),x(t)),
        x(t)*_D,a*x(t)*_D,x(t)*_D^2,x(t-1)*x(t)*_D^2,x(t)*del(x(t-1)),
        -x(t)*del(x(t-2),u(t-1)),_D,a*_D,_D^2,x(t-1)*_D^2,del(x(t-1)),
        -del(x(t-2),u(t-1)),a,a^2,a*_D,a*x(t)*_D,a*del(x(t)),
        -a*del(x(t-1),u(t)),1,a,_D,x(t)*_D,del(x(t)),-del(x(t-1),u(t))]$

testfnc(wedge, wedges=wedgefn)$
