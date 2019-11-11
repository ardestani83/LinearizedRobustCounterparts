function [tau,theta,lambda,lambda_delta,delta,psi,theta_delta]=Subproblem(I,J,eta,D_bar,D_hat,Z0,Gamma,Method,u1,u2)
Z0=0.5*Z0.*(abs(sign(Z0))+sign(Z0));
delta= sdpvar(J,1,'full');
delta_delta= sdpvar(J,J,'full') ;

lambda= sdpvar(J,1,'full') ;
theta= sdpvar(I,1,'full') ;
lambda_delta= sdpvar(J,J,'full') ;
theta_delta= sdpvar(I,J,'full') ;


tau= sdpvar;
%%
Objective=(tau);

Constraints = lambda >= 0;
Constraints = [Constraints,theta >= 0];
Constraints = [Constraints,vec(lambda_delta) >= 0];
Constraints = [Constraints,vec(theta_delta) >= 0];


Constraints = [Constraints,sum(Z0.*theta)+sum(D_bar.*lambda)-sum(D_hat.*(sum(eye(J).*lambda_delta,2)))<=tau];

Constraints = [Constraints,vec(ones(I,1)*lambda'+theta*ones(1,J))>=vec(eta)];

for k=1:J
    Constraints = [Constraints,vec(ones(I,1)*(lambda_delta(:,k))'+theta_delta(:,k)*ones(1,J))>=vec(eta*delta(k))];
end

Constraints = [Constraints,sum(delta)<=Gamma];
Constraints = [Constraints,delta<=1];
Constraints = [Constraints,delta>=0];

Constraints = [Constraints,sum(lambda_delta,2)<=Gamma*lambda];
Constraints = [Constraints,vec(lambda_delta)<=vec(lambda*ones(1,J))];

Constraints = [Constraints,sum(theta_delta,2)<=Gamma*theta];
Constraints = [Constraints,vec(theta_delta)<=vec(theta*ones(1,J))];

Constraints = [Constraints,vec(ones(I,1)*(sum(lambda_delta,2))'+sum(theta_delta,2)*ones(1,J)-eta...
    *sum(delta))<=vec(Gamma*(ones(I,1)*lambda'+theta*ones(1,J)-eta))];



for k=1:J
    Constraints = [Constraints,vec(ones(I,1)*(lambda_delta(:,k))'+theta_delta(:,k)*ones(1,J)-...
        eta*delta(k))<=vec(ones(I,1)*lambda'+theta*ones(1,J)-eta)];
end

if Method=="E"
    Constraints = [Constraints,vec(lambda_delta)<=vec(u1*delta')];
elseif Method=="M"
    
    lambda_theta= sdpvar(J,I,'full') ;
    lambda_lambda= sdpvar(J,J,'full') ;
    theta_theta= sdpvar(I,I,'full') ;
    
    Z= sdpvar(2*J+I+1,2*J+I+1);
    
    Constraints = [Constraints,[lambda;theta]<=[u1;u2]];
    
    Constraints = [Constraints,vec([lambda_delta;theta_delta])<=vec([u1;u2]*delta')];
    
    Constraints = [Constraints,vec(ones(J,1)*(u1-lambda)'-delta*u1'+lambda_delta')>=0];
    Constraints = [Constraints,vec(ones(J,1)*(u2-theta)'-delta*u2'+theta_delta')>=0];
    
    Constraints = [Constraints,sum(lambda_delta,2)+(Gamma-sum(delta))*u1-Gamma*lambda>=0];
    Constraints = [Constraints,sum(theta_delta,2)+(Gamma-sum(delta))*u2-Gamma*theta>=0];
    
end

ops = sdpsettings('solver','cplex');
ops = sdpsettings(ops,'verbose',0);
ops = sdpsettings(ops,'solver','cplex','cplex.timelimit',172800);
optimize(Constraints,Objective,ops);



tau = value (tau );
theta = value (theta );
lambda = value (lambda );
lambda_delta = value (lambda_delta );
delta = value (delta );
theta_delta = value (theta_delta );

psi=ones(I,1)*lambda'+theta*ones(1,J)-eta;
psi=psi(:);
lambda_delta=sum(eye(J).*lambda_delta,2);
return

