function [UB,Z_opt,v_opt]=Decomposition_algorithm(I,J,Gamma,eta,c,K,D_bar,D_hat,M,u1,u2,Method)
LB=-10^8;
UB=10^8;
count=0;
theta_Pool=[];Lambda_Pool=[];lambda_delta_Pool=[];
while (UB-LB)/abs(LB)>=10^-7
    if count==0
        Z0=zeros(I,1);
        if Method=='S'
            [~,theta,lambda,lambda_delta,delta]=Subproblem_SDP(Lo,N,eta,D_bar,D_hat,Z0,Gamma,u1,u2,'S');
        else
            [~,theta,lambda,lambda_delta,delta]=Subproblem(I,J,eta,D_bar,D_hat,Z0,Gamma,Method,u1,u2);
        end
    end
    theta_Pool=[theta_Pool theta];Lambda_Pool=[Lambda_Pool lambda];lambda_delta_Pool=[lambda_delta_Pool lambda_delta];
    
    [tau,Z0,v]=Masterproblem(J,I,M,c,eta,D_bar,D_hat,K,theta_Pool,Lambda_Pool,lambda_delta_Pool,count,delta);
    UB=tau-sum(c.*Z0+K.*v);
    
    if Method=='S'
        [tau,theta,lambda,lambda_delta,delta]=Subproblem_SDP(Lo,N,eta,D_bar,D_hat,Z0,Gamma,u1,u2,'S');
    else
        [tau,theta,lambda,lambda_delta,delta]=Subproblem(I,J,eta,D_bar,D_hat,Z0,Gamma,Method,u1,u2);
    end
    LB=max(LB,tau-sum(c.*Z0+K.*v));
    count=count+1;
end

v_opt=v;
Z_opt=Z0;
return