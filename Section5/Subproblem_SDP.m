function [tau,theta,lambda,lambda_delta,delta,psi,theta_delta]=Subproblem_SDP(m,n,eta,D_bar,D_hat,x,Gamma,u1,u2,Method)
x=max(x,zeros(m,1));
cx = 0;
d = eta(:);
P = [-eye(n); eye(n); ones(1,n)];
q = [zeros(n,1); ones(n,1); Gamma];
Ax = [-D_bar; -x; zeros(n*m,1)];
B = [];
for k=1:n
    tmp = ones(m,1)*([1:n]==k);
    B = [B; vec(tmp)'];
end
for k=1:m
    tmp = ([1:m]==k)'*ones(1,n);
    B = [B; vec(tmp)'];
end
for k=1:n*m
    B = [B; -vec([1:n*m]==k)'];
end
Psi = [-diag(D_hat); zeros(m,n); zeros(m*n,n)];
u = [u1; u2];
Q = [eye(n+m) zeros(n+m,n*m)];


n_const = size(B,1);

    %SDP-LRC
zeta=sdpvar(n,1,'full');
lambda=sdpvar(n_const,1,'full');
Delta=sdpvar(n,n_const,'full');
Lambda=sdpvar(n_const,n_const);
Xi=sdpvar(n,n);
objective = cx+vec(Psi')'*Delta(:)-Ax'*lambda;
constraints = [];
%6b-6d
constraints = [constraints, B'*lambda == d, P*zeta<=q, lambda>=0];
%8b - 8f
constraints = [constraints, vec(Delta*B) == vec(zeta*d'),
    vec(P*Delta) <= vec(q*lambda'),
    vec(Lambda*B) == vec(lambda*d'),
    vec(Lambda) >= 0,
    vec(P*Xi*P'+q*q') >= vec(P*zeta*q'+q*zeta'*P')];
%19
%constraints = [constraints, lambda<= u];
constraints = [constraints, Q*lambda<= u];
%21
%constraints = [constraints, vec(P*Delta) >= vec(q*lambda'-(q-P*zeta)*u')];
constraints = [constraints, vec(P*Delta*Q') >= vec(q*lambda'*Q'-(q-P*zeta)*u')];
%34b 34d
%constraints=[constraints, Lambda(:)<=vec(u*lambda'),
%    vec(Lambda+u*u') >= vec(lambda*u'+u*lambda')];
constraints=[constraints, vec(Q*Lambda)<=vec(u*lambda'),
    vec(Q*Lambda*Q'+u*u') >= vec(Q*lambda*u'+u*lambda'*Q')];
constraints=[constraints, [Lambda Delta' lambda; Delta Xi zeta; lambda' zeta' 1] >= 0];
options=sdpsettings('verbose',0,'solver', 'mosek');
va=optimize(constraints,objective, options);
tau = value(objective);
Delta_ = value(Delta);
zeta_ = value(zeta);
lambda_ = value(lambda);

%vec(Psi')'*Delta(:)-Ax'*lambda
%Z0'*theta_Pool+D_bar'*Lambda_Pool-D_hat'*lambda_delta_Pool>=tau
theta = lambda_(n+1:n+m);
lambda = lambda_(1:n);
lambda_delta = diag(Delta_(1:n,1:n));
delta = zeta_;
psi=ones(m,1)*lambda'+theta*ones(1,n)-eta;
theta_delta = [];
return

