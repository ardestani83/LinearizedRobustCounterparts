function [Z0_master,v_master,tau_master]=CCG_Masterproblem(I,J,eta,c,K,D_bar,D_hat,M,delta_pool,count)
rome_begin;
L = rome_model;

Z0= sdpvar(I,1,'full');
Y= sdpvar(I,J,count+1,'full');
v= binvar(I,1,'full');
tau= sdpvar;
Objective= -(-sum(c.*Z0+K.*v)+tau);
Constraints = Z0 >= 0;
Constraints = [Constraints,Y >= 0];
for l=0:count
    Constraints = [Constraints,(tau<=sum(sum(eta.*Y(:,:,l+1))))];
    for j=1:J
        Constraints = [Constraints,(sum(Y(:,j,l+1))<=D_bar(j)-D_hat(j)*delta_pool(j,l+1))];
    end
    for i=1:I
        Constraints = [Constraints,(sum(Y(i,:,l+1))<=Z0(i))];
    end
end

Constraints = [Constraints, (Z0<=M*v)];
ops = sdpsettings('solver','cplex');
ops = sdpsettings(ops,'verbose',0);
ops = sdpsettings(ops,'solver','cplex','cplex.timelimit',172800);
optimize(Constraints,Objective,ops);

Z0_master = value(Z0 );
v_master=value(v );
tau_master=value(tau );


return