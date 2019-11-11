function [tau,Z0,v]=Masterproblem(J,I,M,c,eta,D_bar,D_hat,K,theta_Pool,Lambda_Pool,lambda_delta_Pool,count,delta)

v= binvar(I,1,'full');
Z0= sdpvar(I,1,'full');
Y= sdpvar(I,J,'full');
tau= sdpvar;


Constraints = Z0 >= 0;
Constraints = [Constraints,Y >= 0];

Objective= -(tau-sum(c.*Z0+K.*v));

Constraints = [Constraints, (Z0'*theta_Pool+D_bar'*Lambda_Pool-D_hat'*lambda_delta_Pool>=tau*ones(1,count+1))];
Constraints = [Constraints, (Z0<=M*v)];

Constraints = [Constraints,(tau<=eta(:)'*Y(:))];
Constraints = [Constraints,(sum(Y)'<=D_bar-D_hat.*delta)];
Constraints = [Constraints,(sum(Y,2)<=Z0)];


ops = sdpsettings('solver','cplex');
ops = sdpsettings(ops,'verbose',0);
ops = sdpsettings(ops,'solver','cplex','cplex.timelimit',172800);
optimize(Constraints,Objective,ops);

tau = value (tau );
Z0 = value (Z0 );
v = value (v);

return