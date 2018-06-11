testfnc(name,testcond):=if testcond then print(name,": passed") else error(name,": failed")$

/* define some test data */
flist:[1,a,x[1](t),u(t),x[1](t)*u(t-1),x[1](t-1)^2,x[2](t+3)*diff(u[1](t-2),t)]$
dlist:[x[1](t),(x[1](t)+x[2](t))*del(x[1](t)),x[2](t-1)*u(t)^2*del(x[2](t))]$
w:sin(x[1](t-1))*u(t-3)*q*x[2]*del(diff(u(t-2),t,2),diff(u(t-8),t))$

testfnc("Noncommutative product (*^)",
   (u(t)*_D*^matrix(flist)=matrix([u(t)*_D,a*u(t)*_D,x[1](t-1)*u(t)*_D,u(t-1)*u(t)*_D,
           u(t-2)*x[1](t-1)*u(t)*_D,x[1](t-2)^2*u(t)*_D,
           u(t)*x[2](t+2)*'diff(u[1](t-3),t,1)*_D])) and
          (u(t)*_D*^matrix(dlist)=matrix([x[1](t-1)*u(t),(x[2](t-1)+x[1](t-1))*u(t)*del(x[1](t-1)),
          x[2](t-2)*u(t-1)^2*u(t)*del(x[2](t-1))]))
          )$

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

testfnc(coefpow,
  coefpow(matrix([_D,1,_D^3+_D^2],_D^3*^[x[5](t-1),_D*^u(t),4]))=
       [[matrix([0,1,0],[0,0,0]),matrix([1,0,0],[0,0,0]),
         matrix([0,0,1],[0,0,0]),matrix([0,0,1],[x[5](t-4),0,4]),
         matrix([0,0,0],[0,u(t-4),0])],[0,1,2,3,4]]
   )$

testfnc(dotfact,
   dotfact(matrix([_d(u(t)*u(t-1))],
                  [_d(x(t-3)*x(t-2))]))
         =[matrix([u(t)*_D+u(t-1), 0],
                  [0,              x(t-2)*_D^3+x(t-3)*_D^2]),
        matrix([del(u(t))],[del(x(t))])] )$

testfnc(euclid,
  unique(map(lambda([u],fullratsimp(euclid(u[1],u[2])-u[3])),[
     [_D^7+1, u(t)*_D^5-x(t), matrix([_D^2/u(t-2),(u(t-2)+x(t-2)*_D^2)/u(t-2)])],
     [_D^2-_D,u(t),matrix([-(_D*u(t-2)-u(t-1)*_D^2)/(u(t-1)*u(t-2)),0])] ] ))
     =[matrix([0,0])])$
     
S:systdef(matrix([x[2](t-1)],[u[1](t)]))$
testfnc(ddt,
  ddt(x[1](t-1)^2,S,4)=2*x[1](t-1)*'diff(u[1](t-2),t,2)
                       +8*x[2](t-2)*'diff(u[1](t-2),t,1)+6*u[1](t-2)^2)$

push(a,flist)$ push(1,flist)$
list:[1,a,_D,x(t)*_D,del(x(t)),del(u(t),x(t-1))]$
wedges:[]$
for i in list do for j in reverse(list) do push(wedge(i,j),wedges)$

testfnc(isclosed,
        not(isclosed([u(t-1)*del(u(t-1)),u(t)*del(u(t-1))]))
        and isclosed(df1)
        )$

M:matrix([_D],[u(t-1)*_D+x(t)],[x(t-1)]);

testfnc(matchvar,
    map(matchvar,[1,a,x[3],u(t),v[3](t),t^2,(t+1),(t-3)^4,x[4](t+3)])=[false, false, false, true, true, false, false, false, true])$

testfnc(ncrowrank,  ncrowrank(M*^transpose(M))=1)$

testfnc(showalltvars,
showalltvars(w)=[u(t-8),u(t-3),u(t-2),x[1](t-1)])$

testfnc(showtvars,
  showtvars(w)=[u(t-3),x[1](t-1),del('diff(u(t-8),t,1),'diff(u(t-2),t,2))])$

testfnc(tshift,
   flatten([tshift(flist),tshift(dlist,2)])=[1,a,x[1](t-1),u(t-1),u(t-2)*x[1](t-1),x[1](t-2)^2,x[2](t+2)*'diff(u[1](t-3),t,1),x[1](t-2),(x[2](t-2)+x[1](t-2))*del(x[1](t-2)),x[2](t-3)*u(t-2)^2*del(x[2](t-2))])$

wedgefn:[-del(x(t-1),u(t)),-a*del(x(t-1),u(t)),-del(x(t-2),u(t-1)),
        -x(t)*del(x(t-2),u(t-1)),-del(x(t-1),u(t),x(t)),0,del(x(t)),
        a*del(x(t)),del(x(t-1)),x(t)*del(x(t-1)),0,-del(x(t-1),u(t),x(t)),
        x(t)*_D,a*x(t)*_D,x(t)*_D^2,x(t-1)*x(t)*_D^2,x(t)*del(x(t-1)),
        -x(t)*del(x(t-2),u(t-1)),_D,a*_D,_D^2,x(t-1)*_D^2,del(x(t-1)),
        -del(x(t-2),u(t-1)),a,a^2,a*_D,a*x(t)*_D,a*del(x(t)),
        -a*del(x(t-1),u(t)),1,a,_D,x(t)*_D,del(x(t)),-del(x(t-1),u(t))]$

testfnc(wedge, wedges=wedgefn)$
