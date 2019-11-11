function [UB,Z_opt,v_opt]=CCG(I,J,eta,c,K,D_bar,D_hat,M,Gamma,u1)
LB=-10^10;
UB=inf;
count=0;
delta_pool=[];
while (UB-LB)/abs(LB)>10^-6
    if count==0
        delta=zeros(J,1);
    end
    delta_pool=[delta_pool delta];
    [Z0,v,tau]=CCG_Masterproblem(I,J,eta,c,K,D_bar,D_hat,M,delta_pool,count);
    UB=tau-c'*Z0-K'*v;
    
    [beta,gama,Delta,delta]=CCG_recourse(I,J,eta,D_bar,D_hat,Gamma,Z0,u1);

    LB=max(LB,D_bar'*beta-D_hat'*Delta+gama'*Z0-c'*Z0-K'*v);
    count=count+1;
end
v_opt=v;
Z_opt=Z0;
return
