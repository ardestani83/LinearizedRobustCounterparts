function [beta,gama,Delta,delta]=CCG_recourse(I,J,eta,D_bar,D_hat,Gamma,Z0,u1)
Z0=max(Z0,zeros(I,1));
beta= sdpvar(J,1,'full');
Delta= sdpvar(J,1,'full');
gama= sdpvar(I,1,'full');
delta= binvar(J,1,'full');

Constraints = beta >= 0;
Constraints = [Constraints,Delta >= 0];
Constraints = [Constraints,gama >= 0];

Objective= (sum(D_bar.*beta)-sum(D_hat.*Delta)+sum(gama.*Z0));

Constraints = [Constraints, (gama*ones(1,J)+(beta*ones(1,I))'>=eta)];
Constraints = [Constraints, (sum(delta)<=Gamma)];
Constraints = [Constraints, (delta<=1)];
Constraints = [Constraints, (Delta<=u1.*delta)];
Constraints = [Constraints, (Delta<=beta)];


ops = sdpsettings('solver','cplex');
ops = sdpsettings(ops,'verbose',0);
ops = sdpsettings(ops,'solver','cplex','cplex.timelimit',172800);
optimize(Constraints,Objective,ops);

beta=value(beta);
gama=value(gama);
Delta=value(Delta);
delta=value(delta);

return