function [U_MLRC,Z_MLRC,v_MLRC,X_MLRC,W_MLRC]=MLRC(I,J,eta,c,K,D_bar,D_hat,M,Gamma,u1,u2)

X= sdpvar(I,J,J);
W= sdpvar(I,J,'full');

X1= sdpvar(J,J,'full');
W1= sdpvar(J,1);
X2= sdpvar(I,J,'full');
W2= sdpvar(I,1);
v= binvar(I,1);


AA1= sdpvar(J,J,'full');
BB1= sdpvar(J,1);
AA2= sdpvar(I,J,'full');
BB2= sdpvar(I,1);
Z0= sdpvar(I,1);
q= sdpvar;
s= sdpvar(I,1);
alphaa= sdpvar(J,1);
BB= sdpvar(I,J,'full');
r= sdpvar(J,1);
t= sdpvar(I,J,'full');
betaa= sdpvar(J,J,'full');
AA= sdpvar(I,J,J,'full');

Constraints = Z0 >= 0;
Constraints = [Constraints,AA1 >= 0];
Constraints = [Constraints,BB1 >= 0];
Constraints = [Constraints,AA2 >= 0];
Constraints = [Constraints,BB2 >= 0];
Constraints = [Constraints,Z0 >= 0];
Constraints = [Constraints,q >= 0];
Constraints = [Constraints,s >= 0];
Constraints = [Constraints,alphaa >= 0];
Constraints = [Constraints,BB >= 0];
Constraints = [Constraints,r >= 0];
Constraints = [Constraints,t >= 0];
Constraints = [Constraints,betaa >= 0];
Constraints = [Constraints,AA >= 0];

Objective= -(sum(sum(eta.*W,2))-sum(c.*Z0+K.*v)-Gamma*q-sum(r)-u1'*W1-u2'*W2);
for k=1:J
    Constraints = [Constraints, (q+r(k)>=-sum(sum(squeeze(X(:,:,k).*eta)))+u1'*X1(:,k)+u2'*X2(:,k))];
end
Constraints = [Constraints, (s*Gamma+sum(t,2)+sum(W,2)-W2<=Z0)];
Constraints = [Constraints, (s*ones(1,J)+t>=squeeze(sum(X,2))-X2)];
Constraints = [Constraints, (alphaa*Gamma+sum(betaa,2)+sum(W)'-W1<=D_bar)];
Constraints = [Constraints, (alphaa+(betaa.*eye(J))*ones(J,1)>=(squeeze(sum(X)).*eye(J))*ones(J,1)+D_hat-(X1.*eye(J))*ones(J,1))];
Constraints = [Constraints, (alphaa*ones(1,J)+betaa- betaa.*eye(J)>=squeeze(sum(X))-squeeze(sum(X)).*eye(J)-X1+X1.*eye(J))];
Constraints = [Constraints, (-Gamma*BB-sum(AA,3)+W>=0)];
for k=1:J
    Constraints = [Constraints, (BB+AA(:,:,k)>=-X(:,:,k))];
end
Constraints = [Constraints, (-Gamma*BB1-sum(AA1,2)+W1>=0)];
Constraints = [Constraints, (BB1*ones(1,J)+AA1>=-X1)];
Constraints = [Constraints, (-Gamma*BB2-sum(AA2,2)+W2>=0)];
Constraints = [Constraints, (BB2*ones(1,J)+AA2>=-X2)];
Constraints = [Constraints, (Z0<=M*v)];


ops = sdpsettings('solver','cplex');
ops = sdpsettings(ops,'verbose',0);
ops = sdpsettings(ops,'solver','cplex','cplex.timelimit',172800);
optimize(Constraints,Objective,ops);

U_MLRC= -value(Objective);
v_MLRC= value(v);
Z_MLRC=value(Z0 );
X_MLRC=value(X);
W_MLRC=value(W);
return