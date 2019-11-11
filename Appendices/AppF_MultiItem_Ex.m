function AppF_MultiItem_Ex
n=3;
r=[80 80 80]'; c=[70 50 20]'; s=[20 15 10]'; p=[60 60 50]';
Gamma=2; dbar=[80 80 60]'; dhat=[60 60 40]';
h=r-s; cx = c-r;
mixingMat = 0.5*[1 1 0; 0 1 1; 1 0 1]; 

aarc = solveAARC(n,Gamma,h,p,dbar,dhat,mixingMat,cx);
aarc.fval_true = findExactWorstCase(aarc.x,n,cx,h,p,Gamma,dbar,dhat,mixingMat);

[fval,x] = solveExactWorstCase(n,cx,h,p,Gamma,dbar,dhat,mixingMat);
exact.fval = fval; exact.fval_true = exact.fval;

sdp_small = getORSDPmodel(n,Gamma,dhat,dbar,mixingMat,cx,h,p);
sdp_small.fval_true = findExactWorstCase(sdp_small.x,n,cx,h,p,Gamma,dbar,dhat,mixingMat);

sdp = getSDPmodel(n,Gamma,dhat,dbar,mixingMat,cx,h,p);
sdp.fval_true = findExactWorstCase(sdp.x,n,cx,h,p,Gamma,dbar,dhat,mixingMat);

%performances bound in terms of profit
display(sprintf('Optimal bound on worst-case profit:  AARC=%d,  SDP-A&D=%d,  SDP-LRC=%d,  Exact=%d',-aarc.fval,-sdp_small.fval,-sdp.fval,-exact.fval));

%actual worst-case performance in terms of profit
display(sprintf('Worst-case profit of solution:  AARC=%d,  SDP-A&D=%d,  SDP-LRC=%d,  Exact=%d',-aarc.fval_true,-sdp_small.fval_true,-sdp.fval_true,-exact.fval_true));
return

function ans = getMatrixOfAllNsequencesKchoice(N,K)
if N==1, ans = [1:K]'; return; end;
ans = [];
for k=1:K, tmp = getMatrixOfAllNsequencesKchoice(N-1,K); ans = [ans; [k*ones(size(tmp,1),1)  tmp]];
end;
return

function fval = findExactWorstCase(x,n,cx,h,p,Gamma,dbar,dhat,mixingMat)
fval = solveExactWorstCase(n,cx,h,p,Gamma,dbar,dhat,mixingMat,x);
return


function [fval,x] = solveExactWorstCase(n,cx,h,p,Gamma,dbar,dhat,mixingMat,xfixed)
if ~exist('xfixed'), xfixed = []; end;

tmpComb = getMatrixOfAllNsequencesKchoice(n,2);
rome_begin; 
h_ = rome_model('exact');
if isempty(xfixed), newvar x(n);
else, x = xfixed; end;
newvar t(1,1);
newvar deltap(n) uncertain;
newvar deltan(n) uncertain;
d = dbar+diag(dhat)*mixingMat*(deltap-deltan);
rome_minimize(cx'*x + t);
for k=1:size(tmpComb,1)
    rome_constraint(t>=sum(diag((tmpComb(k,:)'==1).*h)*(x-d)+diag((tmpComb(k,:)'==2).*p)*(d-x)));
end;
rome_constraint(deltan+deltap<=1);
rome_constraint(deltan>=0);
rome_constraint(deltap>=0);
rome_constraint(sum(deltan)+sum(deltap)== Gamma);
h_.solve;
fval = h_.objective;
if isempty(xfixed), x = h_.eval(x); 
else, x = xfixed; end;
rome_end;
return

function aarc = solveAARC(n,Gamma,h,p,dbar,dhat,mixingMat,cx)

%aarc solution
    rome_begin;
    h_ = rome_model('AARC'); 
    newvar x(n,1);
    newvar deltap(n) uncertain;
    newvar deltan(n) uncertain;
    newvar y(n,[deltap' deltan']) linearrule;
    d = dbar+diag(dhat)*mixingMat*(deltap-deltan);
    rome_minimize(cx'*x + sum(y));
    rome_constraint(y>=diag(h)*(x-d));
    rome_constraint(y>=diag(p)*(d-x));
    rome_constraint(deltan+deltap<=1);
    rome_constraint(deltan>=0);
    rome_constraint(deltap>=0);
    rome_constraint(sum(deltan)+sum(deltap)== Gamma);
    rome_constraint(x>=0);
    t1 = cputime;
    h_.solve; aarc.time=cputime-t1;
    aarc.x = h_.eval(x);
    aarc.y = h_.eval(y);
    aarc.fval = h_.objective;
    rome_end;

return


function sdp_small = getORSDPmodel(n,Gamma,dhat,dbar,mixingMat,cx,h,p)
%OR paper
c2 = [];
c2(:,:,1) = -(diag(h.*dhat)*mixingMat)'; %
c2(:,:,2) = (diag(p.*dhat)*mixingMat)';
dz = [-h.*dbar p.*dbar];

    
    z = sdpvar(n,2);
   deltap = sdpvar(n,1);
   deltan = sdpvar(n,1);
   Deltap = sdpvar(n,n,2,'full'); 
   Deltan = sdpvar(n,n,2,'full');
   Lambdas = {}; for k=1:n, Lambdas{k} = sdpvar(2,2,'symmetric'); end;
   Lambdap = sdpvar(n,n,'symmetric');
   Lambdan = sdpvar(n,n,'symmetric');
   f = c2(:)'*(Deltap(:)-Deltan(:))+dz(:)'*z(:);
   C = [];
   C = [C, deltap >= 0, ...
   deltan >= 0, ...
   deltap+deltan<=1, ...
   sum(deltap+deltan)==Gamma, ...
   sum(z,2) == 1, ...
   sum(Deltap,3)==deltap*ones(1,n), ...
   sum(Deltan,3)==deltan*ones(1,n), ...
   vec(Deltap) >= 0; vec(Deltan) >= 0];
   for k=1:n, C = [C, vec(Deltap(k,:,:)+Deltan(k,:,:)) <= z(:)]; end;
   C = [C, vec(sum(Deltap+Deltan,1))==Gamma*z(:)];
   for k=1:n
   C = [C, [Lambdas{k} reshape(Deltap(:,k,:),n,2)' z(k,:)'; reshape(Deltap(:,k,:),n,2) Lambdap deltap; z(k,:) deltap' 1] >= 0, ...
   [Lambdas{k} reshape(Deltan(:,k,:),n,2)' z(k,:)'; reshape(Deltan(:,k,:),n,2) Lambdan deltan; z(k,:) deltan' 1] >= 0, ...
   Lambdas{k} == diag(z(k,:))];
   end;
   C = [C, diag(Lambdap)<=deltap; diag(Lambdan)<=deltan,...
       vec(Lambdap)>=0, vec(Lambdan)>=0,...
       vec(z)>=0, vec(z)<=1];
    C = [C, cx+h.*z(:,1)-p.*z(:,2)>=0];
    x_index = length(C);
   output = optimize(C, -f,sdpsettings('solver','mosek','verbose',0));
sdp_small.time=output.solvertime;
sdp_small.x=dual(C(x_index));
sdp_small.fval= value(c2(:)'*(Deltap(:)-Deltan(:))+dz(:)'*z(:));

return

function sdp = getSDPmodel(n,Gamma,dhat,dbar,mixingMat,cx,h,p)
P = [eye(n) eye(n); -eye(2*n); ones(1,2*n); -ones(1,2*n)]; q = [ones(n,1); zeros(2*n,1); Gamma; -Gamma];
A = [diag(h); -diag(p)]; a = [-diag(h)*dbar; diag(p)*dbar];
B = [-eye(n); -eye(n)];
Psi = [diag(h.*dhat)*mixingMat -diag(h.*dhat)*mixingMat; -diag(p.*dhat)*mixingMat diag(p.*dhat)*mixingMat];
PsiT = Psi;
c = -cx; d = -ones(n,1);

%optimize in x through dualizing outer minimization
zeta = sdpvar(2*n,1);
mylambda = sdpvar(2*n,1);
Delta = sdpvar(2*n,2*n,'full');
Lambda = sdpvar(2*n,2*n,'symmetric');
Xi = sdpvar(2*n,2*n,'symmetric');
f = trace(Psi*Delta) - a'*mylambda;
C = [];
C = [C, B'*mylambda == d, ...
    P*zeta <= q, ...
    -mylambda <=0, ...
    vec(Delta*B) == vec(zeta*d'), ...
    vec(P*Delta) <= vec(q*mylambda'), ... 
   [Lambda Delta' mylambda; Delta Xi zeta; mylambda' zeta' 1] >= 0, ...
   vec(Lambda*B) == vec(mylambda*d'), ...
   vec(P*Xi*P' + q*q') >= vec(P*zeta*q' + q*zeta'*P'), ...
   vec(Lambda) >= 0];
%   dual variable x;
C = [C, c <= A'*mylambda];
x_index = length(C);
output = optimize(C,f,sdpsettings('solver','mosek','verbose',0));
sdp.time=output.solvertime;
sdp.x=dual(C(x_index));
sdp.fval= value(-(trace(Psi*Delta) - a'*mylambda));

return

