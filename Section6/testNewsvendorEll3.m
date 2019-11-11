function [RESULT,TIME,avgCorrel]=testNewsvendorEll(n,Gamma,correlLevel,algNames)
if ~exist('n'), n = 10; end;
if ~exist('Gamma'), Gamma = 2.5; end;
if ~exist('correlLevel'), correlLevel = 1; end;
if ~exist('algNames'), algNames = {'AARC'}; end;

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
        sigmas = dhat;%/norminv(0.95);
        covarTmp = (sigmas*sigmas').*correlTmp;
        sqCovar = real(sqrtm(covarTmp));
        mixingMat = sqCovar./(diag(sqCovar)*ones(1,n));
        mixingMat=mixingMat./(sum(abs(mixingMat),2)*ones(1,n));
    end;
    
    %aarc solution under Gamma=n
    aarc_test = solveAARC(n,n,h,p,dbar,dhat,mixingMat,cx);
    if aarc_test.fval <= -1, break; end;%find a problem where worst-case profit is above 1$
end;

for algK=1:length(algNames)
    if strcmp(algNames{algK},'Exact')
        [opt_obj,opt_x,timesUp,solTime]=Newsvendor_ColumnConstGen_Ell(n,Gamma,cx,h,p,dbar,dhat, mixingMat,1e3);
        exact.time = solTime;
        if timesUp, exact.time = inf; end;
        exact.x = opt_x;
        exact.fval = opt_obj;
        exact.timesUp = timesUp;
        algResult = exact;
    elseif strcmp(algNames{algK},'AARC')
        algResult = solveAARC(n,Gamma,h,p,dbar,dhat,mixingMat,cx);
    elseif strcmp(algNames{algK},'SDP-A&D')
        algResult = getORSDPmodel(n,Gamma,dhat,dbar,mixingMat,cx,h,p,1==1);
    elseif strcmp(algNames{algK},'SDP-LRC2')
        algResult = getORSDPmodel(n,Gamma,dhat,dbar,mixingMat,cx,h,p,1==0);
    elseif strcmp(algNames{algK},'SDP-LRC')
        algResult = getSDPmodel(n,Gamma,dhat,dbar,mixingMat,cx,h,p);
    end;
    %performances in profit
    RESULT(1,algK) = -algResult.fval;
    TIME(1,algK) = algResult.time;
end;
return

function aarc = solveAARC(n,Gamma,h,p,dbar,dhat,mixingMat,cx)
P = [eye(n) eye(n); -eye(2*n)]; q = [ones(n,1); zeros(2*n,1)];
A = [diag(h); -diag(p)]; a = [-diag(h)*dbar; diag(p)*dbar];
B = [-eye(n); -eye(n)];
Psi = [diag(h.*dhat)*mixingMat -diag(h.*dhat)*mixingMat; -diag(p.*dhat)*mixingMat diag(p.*dhat)*mixingMat];
PsiT = Psi;
c = -cx; d = -ones(n,1);
%u = ones(2*n,1); %this is redundant for newsvendor problem "complete recourse"

diffZetas = [eye(2*n)];

%use yalmip to call mosek
%optimize in x through dualizing outer minimization
zeta = sdpvar(2*n,1);
mylambda = sdpvar(2*n,1);
Delta = sdpvar(2*n,2*n,'full');
Lambda = sdpvar(2*n,2*n,'symmetric');
Xi = sdpvar(2*n,2*n,'symmetric');
%   dual variable x;
f = trace(Psi*Delta) - a'*mylambda;
C = [];
C = [C, B'*mylambda == d, ...
    P*zeta <= q, ...
    -mylambda <=0, ...
    %    mylambda <= u, ...
    vec(Delta*B) == vec(zeta*d'), ...
    vec(P*Delta) <= vec(q*mylambda'), ...
    %    vec(P*Delta) >= vec(q*mylambda'-(q-P*zeta)*u'), ...
    vec(Lambda*B) == vec(mylambda*d'), ...
    vec(P*Xi*P' + q*q') >= vec(P*zeta*q' + q*zeta'*P'), ...
    vec(Lambda) >= 0, ...
    %   vec(Lambda) <= vec(u*mylambda'), ...
    %   vec(Lambda+u*u')>=vec(mylambda*u'+u*mylambda'), ...
    norm(diffZetas*zeta)<=Gamma];
for k=1:size(Delta,2)
    C = [C, norm(diffZetas*Delta(:,k))<=Gamma*mylambda(k)];
end;
for k=1:size(P,1)
    C = [C, norm(diffZetas*(q(k)*zeta-Xi*P(k,:)'))<=Gamma*(q(k)-P(k,:)*zeta)];
end;
%   x: c <= A'*mylambda;
C = [C, c <= A'*mylambda];
x_index = length(C);
output = optimize(C,f,sdpsettings('solver','cplex+','verbose',0));%,'mosek.MSK_DPAR_OPTIMIZER_MAX_TIME',10))
aarc.time=output.solvertime;
aarc.x=dual(C(x_index));
aarc.fval= value(-(trace(Psi*Delta) - a'*mylambda));
return


function sdp_small = getORSDPmodel(n,Gamma,dhat,dbar,mixingMat,cx,h,p,isNaive)
%Inspired by OR paper
P = [eye(n) eye(n); -eye(2*n)]; q = [ones(n,1); zeros(2*n,1)];
A = [diag(h); -diag(p)]; a = [-diag(h)*dbar; diag(p)*dbar];
B = [-eye(n); -eye(n)];
Psi = [diag(h.*dhat)*mixingMat -diag(h.*dhat)*mixingMat; -diag(p.*dhat)*mixingMat diag(p.*dhat)*mixingMat];
PsiT = Psi;
c = -cx; d = -ones(n,1);
%u = ones(2*n,1);

indexes = [];
for k=1:n
    indexes(:,end+1) = [k,k+n,2*n+1:2*n+n,4*n+1];
    indexes(:,end+1) = [k,k+n,3*n+1:3*n+n,4*n+1];
end;

diffZetas = [eye(2*n)];

%use yalmip to call mosek
%optimize in x through dualizing outer minimization
zeta = sdpvar(2*n,1);
mylambda = sdpvar(2*n,1);
Delta = sdpvar(2*n,2*n,'full');
Lambda = sdpvar(2*n,2*n,'symmetric');
Xi = sdpvar(2*n,2*n,'symmetric');
Z = sdpvar(4*n+1,4*n+1,'symmetric');
%   dual variable x;
f = trace(Psi*Delta) - a'*mylambda;
C = [];
C = [C, B'*mylambda == d, ... %(6b)
    P*zeta <= q, ... %(6c)
    -mylambda <=0, ... %(6d)
    %    mylambda <= u, ...
    vec(Delta*B) == vec(zeta*d'), ... %(8b)
    vec(P*Delta) <= vec(q*mylambda'), ... %(8c)
    %    vec(P*Delta) >= vec(q*mylambda'-(q-P*zeta)*u'), ...
    Z == [Lambda Delta' mylambda; Delta Xi zeta; mylambda' zeta' 1], ...
    %Z >= 0;
    vec(Lambda*B) == vec(mylambda*d'), ... %(8d)
    vec(P*Xi*P' + q*q') >= vec(P*zeta*q' + q*zeta'*P'), ... %(8g)
    vec(Lambda) >= 0, ... %(8e)
    %   vec(Lambda) <= vec(u*mylambda'), ...
    %   vec(Lambda+u*u')>=vec(mylambda*u'+u*mylambda'), ...
    norm(diffZetas*zeta)<=Gamma]; %(35c)
for k=1:size(indexes,2)
    C = [C, Z(indexes(:,k),indexes(:,k))>=0];
end;
if ~isNaive
    for k=1:size(Delta,2)
        C = [C, norm(diffZetas*Delta(:,k))<=Gamma*mylambda(k)]; %(35f)
    end;
    for k=1:size(P,1)
        C = [C, norm(diffZetas*(q(k)*zeta-Xi*P(k,:)'))<=Gamma*(q(k)-P(k,:)*zeta)]; %(extra for linear constraints)
    end;
end;
%   x: c <= A'*mylambda;
C = [C, c <= A'*mylambda];
x_index = length(C);
output = optimize(C,f,sdpsettings('solver','mosek','verbose',0));%,'mosek.MSK_DPAR_OPTIMIZER_MAX_TIME',10))
sdp_small.time=output.solvertime;
sdp_small.x=dual(C(x_index));
sdp_small.fval= value(-(trace(Psi*Delta) - a'*mylambda));
return

function sdp = getSDPmodel(n,Gamma,dhat,dbar,mixingMat,cx,h,p)
P = [eye(n) eye(n); -eye(2*n)]; q = [ones(n,1); zeros(2*n,1)];
A = [diag(h); -diag(p)]; a = [-diag(h)*dbar; diag(p)*dbar];
B = [-eye(n); -eye(n)];
Psi = [diag(h.*dhat)*mixingMat -diag(h.*dhat)*mixingMat; -diag(p.*dhat)*mixingMat diag(p.*dhat)*mixingMat];
PsiT = Psi;
c = -cx; d = -ones(n,1);
%u = ones(2*n,1); %redundant

diffZetas = [eye(2*n)];

%use yalmip to call mosek
%optimize in x through dualizing outer minimization
zeta = sdpvar(2*n,1);
mylambda = sdpvar(2*n,1);
Delta = sdpvar(2*n,2*n,'full');
Lambda = sdpvar(2*n,2*n,'symmetric');
Xi = sdpvar(2*n,2*n,'symmetric');
%   dual variable x;
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
    vec(Lambda) >= 0, ...
    %   vec(Lambda) <= vec(u*mylambda'), ...
    %   vec(Lambda+u*u')>=vec(mylambda*u'+u*mylambda'), ...
    norm(diffZetas*zeta)<=Gamma];
for k=1:size(Delta,2)
    C = [C, norm(diffZetas*Delta(:,k))<=Gamma*mylambda(k)];
end;
for k=1:size(P,1)
    C = [C, norm(diffZetas*(q(k)*zeta-Xi*P(k,:)'))<=Gamma*(q(k)-P(k,:)*zeta)];
end;

%   x: c <= A'*mylambda;
C = [C, c <= A'*mylambda];
x_index = length(C);
output = optimize(C,f,sdpsettings('solver','mosek','verbose',0));%,'mosek.MSK_DPAR_OPTIMIZER_MAX_TIME',10))
sdp.time=output.solvertime;
sdp.x=dual(C(x_index));
sdp.fval= value(-(trace(Psi*Delta) - a'*mylambda));
return

function [opt_obj,opt_x,timesUp,solTime]=Newsvendor_ColumnConstGen_Ell(n,Gamma,cx,h,p,dbar,dhat, mixingMat,big_M,x_fixed)
tolerance = 1e-3;
if n<=20, maxTime = 600;
elseif n<=40, maxTime = 900;
elseif n<=100, maxTime = 3600;
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
    time_Master(end+1) = cputime-t1;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Identify worst-case z
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [cost_wc,z,timesUp,solTime_] = solveWorstCaseCostBudgetSet(x,n,Gamma,cx,h,p,dbar,dhat, mixingMat,big_M,maxTime-solTime);
    solTime = solTime+solTime_;
    time_Slave(end+1) = solTime_;
    if timesUp
        opt_obj = cost_wc;
        opt_x = x;
        timesUp = 1==1;
        break;
    end;
    
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
C = [C, norm([deltap; deltan]) <= Gamma];
h = cx'*x + sum(y);
diagnostic = optimize(C,-h,sdpsettings('solver','+cplex','verbose',0,'cplex.timelimit',1.15*timeLeft));
timesUp = diagnostic.problem == 3;
solTime = diagnostic.solvertime;
z = value(deltap-deltan);
profit_wc=value(h);
return
