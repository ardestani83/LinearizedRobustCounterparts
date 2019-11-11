function AppB_MultiProdAssembly_Ex

m=2; n=3;
c = [25 3]';
q = [380 800 1200]';
l = [0 0 0]';
s = [4 1]';
A = [9 9; 0 5; 9 4]';
u1 = [425 805 1240]';
u2 = [140 310]';
u3 = u1+A'*u2-(q-l+A'*s);
dbar = [9000 10000 8000]';
dhat = [8000 8000 8000]';
M = 100000;
Gamma = 2;
u1_ = q-l+A'*s;
u2_ = max((A>0).*(ones(2,1)*u1_')./(A+eps),[],2);
d = q-l-A'*s;
B = [eye(size(A,2)); A; -eye(size(A,2))];
u = getMaxU(B,d,max([u1;u2;u3]));
u1 = u(1:n); u2 = u(n+[1:m]); u3 = u(n+m+[1:n]);

%aarc solution
rome_begin; 
h_ = rome_model('AARC'); % Create Rome Model
newvar x(m,1);
newvar deltan(n) uncertain;
newvar y(n,deltan') linearrule; 
d = dbar+diag(dhat)*(-deltan);
rome_maximize(-c'*x + (q-l)'*y+s'*(x-A*y));
rome_constraint(y<=d);
rome_constraint(A*y<=x);
rome_constraint(deltan<=1);
rome_constraint(deltan>=0);
rome_constraint(sum(deltan)<= Gamma);
rome_constraint(x>=0);
rome_constraint(y>=0);
rome_constraint(x<=M);
h_.solve;
aarc.x = h_.eval(x);
aarc.y = h_.eval(y);
aarc.fval = h_.objective;
rome_end;
aarc.fval_true = evalWorstCase(aarc.x,n,m,c,q,l,s,A,dbar,dhat);



%mlrc solution
rome_begin; 
h_ = rome_model('MLRC'); % Create Rome Model
newvar x(m,1);
newvar deltan(n) uncertain;
newvar y(n,deltan') linearrule; 
newvar z1(n,deltan') linearrule; 
newvar z2(m,deltan') linearrule; 
newvar z3(n,deltan') linearrule; 
d = dbar-diag(dhat)*deltan;
%rome_maximize(-c'*x + (q-l)'*y+s'*(x-A*y)-u1'*z1-u2'*z2);
rome_maximize(-c'*x + (q-l)'*y+s'*(x-A*y)-u1'*z1-u2'*z2-u3'*z3);
rome_constraint(y<=d+z1);
rome_constraint(A*y<=x+z2);
rome_constraint(deltan<=1);
rome_constraint(deltan>=0);
rome_constraint(sum(deltan)<= Gamma);
rome_constraint(x>=0);
rome_constraint(y+z3>=0);
rome_constraint(x<=M);
rome_constraint(z1>=0);
rome_constraint(z2>=0);
rome_constraint(z3>=0);
h_.solve;
mlrc.x = h_.eval(x);
mlrc.y = h_.eval(y);
mlrc.fval = h_.objective;
rome_end;
mlrc.fval_true = evalWorstCase(mlrc.x,n,m,c,q,l,s,A,dbar,dhat);

%exact solution
deltas = [0 0 0; 1 1 0; 1 0 1; 0 1 1; 0 0 1; 1 0 0; 0 1 0]';
rome_begin; 
h_ = rome_model('Optimal'); % Create Rome Model
newvar x(m,1);
newvar y(n,size(deltas,2));
newvar t(1);
ds = dbar*ones(1,size(deltas,2))-diag(dhat)*deltas;
xs = x*ones(1,size(y,2));
rome_maximize(-c'*x + t);
rome_constraint(y(:)<=ds(:));
rome_constraint(t<=(q-l)'*y+s'*(xs-A*y));
rome_constraint(A*y<=xs);
rome_constraint(x>=0);
rome_constraint(y>=0);
rome_constraint(x<=M);
h_.solve;
exact.x = h_.eval(x);
exact.y = h_.eval(y);
exact.fval = h_.objective;
rome_end;
exact.fval_true = evalWorstCase(exact.x,n,m,c,q,l,s,A,dbar,dhat);

display('The vector of penalties is as follows:');
u1 = u(1:3)'
u2 = u(4:5)'
u3 = u(6:end)'

display('Number of parts of each type are as follows: AARC | MLRC | Exact');
[aarc.x mlrc.x exact.x]

display('Optimal bound obtained from each model is as follows: AARC | MLRC | Exact');
[aarc.fval mlrc.fval exact.fval]

display('Worst-case profit achieved by solution of each model is as follows: AARC | MLRC | Exact');
[aarc.fval_true mlrc.fval_true exact.fval_true]


return

function [fval,ys] = evalWorstCase(x,n,m,c,q,l,s,A,dbar,dhat)

deltas = [0 0 0; 1 1 0; 1 0 1; 0 1 1; 0 0 1; 1 0 0; 0 1 0]';

for k=1:size(deltas,2)
rome_begin; 
h_ = rome_model('MLRC'); % Create Rome Model
newvar y(n,1); 
d = dbar-diag(dhat)*deltas(:,k);
rome_maximize(-c'*x + (q-l)'*y+s'*(x-A*y));
rome_constraint(y<=d);
rome_constraint(A*y<=x);
rome_constraint(y>=0);
h_.solve;
ys(:,k) = h_.eval(y);
fvals(k)= h_.objective;
rome_end;
end;
fval = min(fvals);
return