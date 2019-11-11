function [U_ELAARC,Z_ELAARC,v_ELAARC,X_ELAARC,W_ELAARC]=ELAARC(I,J,eta,c,K,D_bar,D_hat,M,Gamma,u1)

X= sdpvar(I,J,J);
W= sdpvar(I,J,'full');
v= binvar(I,1);
Z0= sdpvar(I,1);
q= sdpvar;
s= sdpvar(I,1);
alfa= sdpvar(J,1);
B= sdpvar(I,J,'full');
r= sdpvar(J,1);
t= sdpvar(I,J,'full');
Beta= sdpvar(J,J,'full');
A= sdpvar(I,J,J,'full');
Ln= sdpvar(J,1);

Constraints = Z0 >= 0;
Constraints = [Constraints,q >= 0];
Constraints = [Constraints,s >= 0];
Constraints = [Constraints,alfa >= 0];
Constraints = [Constraints,B >= 0];
Constraints = [Constraints,r >= 0];
Constraints = [Constraints,t >= 0];
Constraints = [Constraints,Beta >= 0];
Constraints = [Constraints,A >= 0];
Constraints = [Constraints,Ln >= 0];

Objective =-(sum(sum(eta.*W,2))-sum(c.*Z0+K.*v)-Gamma*q-sum(r));
for k=1:J
    Constraints = [Constraints, (q+r(k)>=u1(k).*Ln(k)-sum(sum(squeeze(X(:,:,k).*eta))))];
end
Constraints = [Constraints, (s*Gamma+sum(t,2)+sum(W,2)<=Z0)];
Constraints = [Constraints, (s*ones(1,J)+t>=squeeze(sum(X,2)))];
Constraints = [Constraints, (alfa*Gamma+sum(Beta,2)+sum(W)'<=D_bar)];
Constraints = [Constraints, (alfa+(Beta.*eye(J))*ones(J,1)>=(squeeze(sum(X)).*eye(J))*ones(J,1)+D_hat-Ln)];
Constraints = [Constraints, (alfa*ones(1,J)+Beta- Beta.*eye(J)>=squeeze(sum(X))-squeeze(sum(X)).*eye(J))];
Constraints = [Constraints, (-Gamma*B-sum(A,3)+W>=0)];
for k=1:J
    Constraints = [Constraints, (B+A(:,:,k)>=-X(:,:,k))];
end
Constraints = [Constraints, (Z0<=M*v)];


ops = sdpsettings('solver','cplex');
ops = sdpsettings(ops,'verbose',0);
ops = sdpsettings(ops,'solver','cplex','cplex.timelimit',172800);
optimize(Constraints,Objective,ops);

U_ELAARC = -value(Objective);
v_ELAARC = value(v);
Z_ELAARC = value(Z0);
X_ELAARC = value(X);
W_ELAARC = value(W);
return