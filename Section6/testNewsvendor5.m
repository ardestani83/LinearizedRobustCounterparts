function [RESULT,TIME,avgCorrel]=testNewsvendor5(n,Gamma,correlLevel)
if ~exist('n'), n = 40; end;
if ~exist('Gamma'), Gamma = 0.5*n; end;
if ~exist('correlLevel'), correlLevel = 0.1; end;


while 1==1
    r=100*ones(n,1);
    c=ceil(60*rand(n,1));
    s=ceil(c.*rand(n,1));
    p=ceil(c.*rand(n,1));
    dbar=100*ones(n,1);
    dhat=ceil(dbar.*rand(n,1));
    
    r=100*rand(n,1);
    c=r.*rand(n,1);
    s=c.*rand(n,1);
    p=c.*rand(n,1);
    dbar=100*rand(n,1);
    dhat=dbar.*rand(n,1);
    
    
    h=r-s; cx = c-r;
    
    if isinf(correlLevel)
        %IF CORRELLEVAL=inf then no correlation
        mixingMat = eye(n); avgCorrel = 0;
    else
        %IF CORRELLEVAL<inf than correlation is inversely proportional to CORRELLEVAL
        correlTmp = vineBeta(n, correlLevel);
        avgCorrel = sum(sum(abs(correlTmp-eye(n))))/(n^2-n);
        sigmas = dhat;
        covarTmp = (sigmas*sigmas').*correlTmp;
        sqCovar = real(sqrtm(covarTmp));
        mixingMat = sqCovar./(diag(sqCovar)*ones(1,n));
        mixingMat=mixingMat./(sum(abs(mixingMat),2)*ones(1,n));
        
    end;
    
    
    %aarc solution under Gamma=n
    aarc_test = solveAARC(n,n,h,p,dbar,dhat,mixingMat,cx);
    if aarc_test.fval <= -1, break; end;%find a problem where worst-case profit is above 1$
    
end;

aarc = solveAARC(n,Gamma,h,p,dbar,dhat,mixingMat,cx);

if 1==1
    [opt_obj,opt_x,timesUp,solTime]=Newsvendor_ColumnConstGen(n,Gamma,cx,h,p,dbar,dhat, mixingMat,1e3);
    exact.time = solTime;
    if timesUp, exact.time = inf; end;
    exact.x = opt_x;
    exact.fval = opt_obj;
    exact.timesUp = timesUp;
end;

sdp_small = getORSDPmodel(n,Gamma,dhat,dbar,mixingMat,cx,h,p);

sdp = getSDPmodel(n,Gamma,dhat,dbar,mixingMat,cx,h,p);

%performances in profit
RESULT= -[exact.fval aarc.fval sdp_small.fval sdp.fval];
TIME=[exact.time aarc.time sdp_small.time sdp.time];
return

function aarc = solveAARC(n,Gamma,h,p,dbar,dhat,mixingMat,cx)
P = [eye(n) eye(n); -eye(2*n); ones(1,2*n); -ones(1,2*n)]; q = [ones(n,1); zeros(2*n,1); Gamma; -Gamma];
A = [diag(h); -diag(p)]; a = [-diag(h)*dbar; diag(p)*dbar];
B = [-eye(n); -eye(n)];
Psi = [diag(h.*dhat)*mixingMat -diag(h.*dhat)*mixingMat; -diag(p.*dhat)*mixingMat diag(p.*dhat)*mixingMat];
PsiT = Psi;
c = -cx; d = -ones(n,1);
u = ones(2*n,1);

x = sdpvar(n,1,'full');
Y = sdpvar(size(B,2),size(P,2),'full');
y = sdpvar(size(B,2),1,'full');
S = sdpvar(size(P,1),size(A,1),'full');
s = sdpvar(size(P,1),1,'full');
h = c'*x+d'*y-q'*s;
C = [];
C = [vec(P'*S)==vec(Y'*B'-Psi'),...
    a+A*x+B*y+S'*q<=0, ...
    vec(P'*s) == vec(- Y'*d), ...
    s>=0, ...
    vec(S)>=0, ...
    x>=0];
output = optimize(C,-h,sdpsettings('solver','+cplex','verbose',0));
%optimize(C,-h,sdpsettings('solver','mosek','verbose',0,'mosek.MSK_DPAR_OPTIMIZER_MAX_TIME',10))
aarc.time=output.solvertime;
aarc.x=value(x);
aarc.fval= -value(h);
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
%u = ones(2*n,1); %u is redundant for this model

%use yalmip to call mosek
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
%    mylambda <= u, ...
    vec(Delta*B) == vec(zeta*d'), ...
    vec(P*Delta) <= vec(q*mylambda'), ... 
%    vec(P*Delta) >= vec(q*mylambda'-(q-P*zeta)*u'), ...
   [Lambda Delta' mylambda; Delta Xi zeta; mylambda' zeta' 1] >= 0, ...
   vec(Lambda*B) == vec(mylambda*d'), ...
   vec(P*Xi*P' + q*q') >= vec(P*zeta*q' + q*zeta'*P'), ...
   vec(Lambda) >= 0];%, ...
%   vec(Lambda) <= vec(u*mylambda'), ...
%   vec(Lambda+u*u')>=vec(mylambda*u'+u*mylambda')];
C = [C, c <= A'*mylambda];
x_index = length(C);
t1 = cputime;
output = optimize(C,f,sdpsettings('solver','mosek','verbose',0));%,'mosek.MSK_DPAR_OPTIMIZER_MAX_TIME',10))    
sdp.time=output.solvertime;
sdp.x=dual(C(x_index));
sdp.fval= value(-(trace(Psi*Delta) - a'*mylambda));
return

function [opt_obj,opt_x,timesUp,solTime]=Newsvendor_ColumnConstGen(n,Gamma,cx,h,p,dbar,dhat, mixingMat,big_M,x_fixed)
tolerance = 1e-3;
if n<=20, maxTime = 300;
elseif n<=40, maxTime = 3600;
elseif n<=100, maxTime = 6*3600;
else, maxTime = 24*3600;    
end;
if ~exist('x_fixed'), x_fixed = []; end;

zbars = zeros(n,1);
time_Master =[]; time_Slave = [];
solTime = 0; timesUp = 1==0;
opt_obj = inf;
t1 = cputime;
bounds = [];
while 1==1
    
    Kprime = size(zbars,2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Optimize using subset of vertices as uncertainty set
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ds = dbar*ones(1,Kprime)+diag(dhat)*mixingMat*zbars;
    
        x = sdpvar(n,1);
        y = sdpvar(n,Kprime,'full');
        t = sdpvar(1,1);
        C = [];
        for k=1:Kprime
            C = [C, t>=cx'*x+sum(y(:,k)), ...
                y(:,k)>=diag(h)*(x-ds(:,k)),...
                y(:,k)>=diag(p)*(ds(:,k)-x)];
        end;
        C = [C, x>=0];
        if ~isempty(x_fixed)
            C = [C, x == x_fixed];
        end;
        diagnose=optimize(C,t,sdpsettings('solver','+cplex','verbose',0));
        x = value(x);
        cost_lowerbound=value(t);
        solTime_ = diagnose.solvertime;
    
    solTime = solTime+solTime_;
    time_Master(end+1) = solTime_;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Identify worst-case z
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [cost_wc,z,timesUp,solTime_] = solveWorstCaseCostBudgetSet(x,n,Gamma,cx,h,p,dbar,dhat, mixingMat,big_M,maxTime-solTime);
    solTime = solTime+solTime_;
    time_Slave(end+1) = solTime_;
    if timesUp, break; end;
    
    %test if wc z is in zbar
    if min(sqrt(sum((zbars-z*ones(1,size(zbars,2))).^2)))<sqrt(eps)
        cost_wc = cost_lowerbound;
    end;
    bounds(:,end+1) = [cost_wc; cost_lowerbound];
    
    
    if cost_wc-cost_lowerbound<tolerance %small enough
        opt_obj = cost_wc;
        opt_x = x;
        timesUp = 1==0;
        break;
    elseif solTime>maxTime
        opt_obj = cost_wc;
        opt_x = x;
        timesUp = 1==1;
        break;
    end;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Add new z to set of vertices
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    zbars = [zbars z];
end
return

function [profit_wc,z,timesUp,solTime] = solveWorstCaseCostBudgetSet(x,n,Gamma,cx,h,p,dbar,dhat, mixingMat,big_M,timeLeft)
Big_Mplus = p.*(dbar+dhat-x)-h.*(x-dbar-dhat);
Big_Mminus = h.*(x-dbar+dhat)-p.*(dbar-dhat-x);


deltap = sdpvar(n,1);
deltan = sdpvar(n,1);
y= sdpvar(n,1);
zplus = binvar(n,1);
zminus = binvar(n,1);
C = [];
C = [C, y>=diag(h)*(x-(dbar+diag(dhat)*mixingMat*(deltap-deltan)))];
C = [C, y>=diag(p)*((dbar+diag(dhat)*mixingMat*(deltap-deltan))-x)];
C = [C, zplus+zminus == 1];
C = [C, y-diag(h)*(x-(dbar+diag(dhat)*mixingMat*(deltap-deltan)))<=Big_Mplus.*(1-zplus)];
C = [C, y-diag(p)*((dbar+diag(dhat)*mixingMat*(deltap-deltan))-x)<=Big_Mminus.*(1-zminus)];
C = [C, deltan+deltap<=1];
C = [C, deltan>=0, deltap>=0];
C = [C, sum(deltan)+sum(deltap)== Gamma];
h = cx'*x + sum(y);
diagnostic = optimize(C,-h,sdpsettings('solver','+cplex','verbose',0,'cplex.timelimit',1.15*timeLeft));
timesUp = diagnostic.problem == 3;
solTime = diagnostic.solvertime;
z = value(deltap-deltan);
profit_wc=value(h);

return
