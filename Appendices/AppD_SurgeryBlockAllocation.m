function AppD_SurgeryBlockAllocation


n = 3;
m=2;

cv = 1000;
cf = 390000;
dmin = [0; 240; 320];
dmax = [160; 352; 512];
T = 480;
Gamma = 2;

for numberORs=1:2 
%Denton model
rome_begin; 
h = rome_model('Denton'); 
newvar x(m,1) binary;
newvar y(n,m) binary; 
newvar alpha1;
newvar kau(n,m);
newvar gamma1(m,1);
rome_minimize(cf*sum(x)+Gamma*alpha1+sum(gamma1));
rome_constraint(y<=ones(n,1)*x');
rome_constraint(sum(y,2)==1);
rome_constraint(alpha1+kau >= cv*((dmax-dmin)*ones(1,m)).*y);
for k=1:m
    rome_constraint(sum(kau(:,k),1)<= cv*(T*x(k)-dmin'*y(:,k))+gamma1(k));
end;
rome_constraint(alpha1>=0);
rome_constraint(kau>=0);
rome_constraint(gamma1>=0);
rome_constraint(sum(x)==numberORs);
h.solve;
DentonProb.fval(numberORs,1) = h.objective;
rome_end;

[fval, x, y] = solveSurgeryBlockAllocationFast(n,m,cf,cv,T,dmin,dmax,Gamma,numberORs,'LPRC');
aarc.fval(numberORs,1) = fval;

[fval, x, y] = solveSurgeryBlockAllocationFast(n,m,cf,cv,T,dmin,dmax,Gamma,numberORs,'exact');
exact.fval(numberORs,1) = fval;
end;

display(sprintf('Open one OR: RORA=%d,  AARC=%d,   Exact=%d',DentonProb.fval(1),aarc.fval(1),exact.fval(1)));
display(sprintf('Open two ORs: RORA=%d,  AARC=%d,   Exact=%d',DentonProb.fval(2),aarc.fval(2),exact.fval(2)));
return

function [fval, x, y] = solveSurgeryBlockAllocationFast(n,m,cf,cv,T,dmin,dmax,Gamma,numberORs,approxMode)
if numberORs == 1
    xs = [1 0; 0 1];
    ys(:,:,1) = [1 0; 1 0; 1 0];
    ys(:,:,2) = [0 1; 0 1; 0 1];
elseif numberORs == 2
    xs = ones(6,1)*[1 1];
    ys(:,:,1) = [0 1; 0 1; 1 0];
    ys(:,:,2) = [0 1; 1 0; 0 1];
    ys(:,:,3) = [0 1; 1 0; 1 0];
    ys(:,:,4) = [1 0; 0 1; 0 1];
    ys(:,:,5) = [1 0; 0 1; 1 0];
    ys(:,:,6) = [1 0; 1 0; 0 1];
end;


for k=1:size(xs,1)
    x = xs(k,:)'; y = ys(:,:,k);
if strcmp(approxMode,'exact')
        [wc,dx,dy] = getWorstCase(x,y,n,m,cf,cv,T,dmin,dmax,Gamma);
elseif strcmp(approxMode,'LPRC')
        [wc,dx,dy] = getLPRCApprox(x,y,n,m,cf,cv,T,dmin,dmax,Gamma);
end;
    fvals(k) = wc;
end;
[fval,index] = min(fvals);
x = xs(min(index),:)'; y = ys(:,:,min(index));
return


function [wc,dx,dy] = getWorstCase(x,y,n,m,cf,cv,T,dmin,dmax,Gamma)
dhat = dmax-dmin;
rome_begin; 
h = rome_model('Exact'); 
newvar z(m) binary;
newvar Delta(n,m);
newvar delta(n);
rome_maximize(cf*sum(x)+cv*(dmin'*y*z-T*x'*z+dhat'*(y.*Delta)*ones(m,1)));
rome_constraint(sum(delta) <= Gamma);
rome_constraint(delta<= 1);
rome_constraint(delta>= 0);
rome_constraint(Delta>= 0);
rome_constraint(Delta<= ones(n,1)*z');
rome_constraint(Delta<= delta*ones(1,m));
rome_constraint(sum(Delta,1)'<= Gamma*z);
h.solve;
prob.z = h.eval(z);
prob.delta = h.eval(delta);
prob.Delta = h.eval(Delta);
prob.fval = h.objective;
rome_end;
z = prob.z; Delta = prob.Delta; delta = prob.delta;
wc = prob.fval;
dx = cf*ones(m,1)-cv*T*z;
dy = cv*dmin*z'+(dhat*ones(1,m)).*Delta;
dy = dy(:);

return

function [wc,dx,dy] = getLPRCApprox(x,y,n,m,cf,cv,T,dmin,dmax,Gamma)

dhat = dmax-dmin;

Delta=sdpvar(n,m,'full');
z = sdpvar(m,1);
delta = sdpvar(n,1);
objective = cf*sum(x)+cv*(dmin'*y*z-T*x'*z+dhat'*(y.*Delta)*ones(m,1));
constraints = [];
constraints = [constraints, delta >= 0];
constraints = [constraints, delta <= 1];
constraints = [constraints, sum(delta) <= Gamma];
constraints = [constraints, z >= 0];
constraints = [constraints, z <= 1];
constraints = [constraints, vec(Delta)>=0];
constraints = [constraints, vec(Delta)<=vec(ones(n,1)*z')];
constraints = [constraints, vec(Delta) <= vec(delta*ones(1,m))];
constraints = [constraints, ones(1,n)*Delta<=Gamma*z'];
constraints = [constraints, sum(delta)-sum(Delta,1)'<= Gamma*(1-z)];
constraints = [constraints, 1-delta*ones(1,m)-ones(n,1)*z'+Delta>=0];
options=sdpsettings('verbose',0,'solver', 'cplex');
va=optimize(constraints,-objective, options);
tau = value(objective);
Delta=value(Delta);
z = value(z);
delta = value(delta);

wc = cf*sum(x)+cv*(dmin'*y*z-T*x'*z+dhat'*(y.*Delta)*ones(m,1));
dx = cf*ones(m,1)-cv*T*z;
dy = cv*dmin*z'+(dhat*ones(1,m)).*Delta;
dy = dy(:);

return



