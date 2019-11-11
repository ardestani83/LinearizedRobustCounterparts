clear
clc
rng(113)

I=5;
J=10;
Total_iter=100;
Exact=zeros(Total_iter,3,3);
AARC_Obj=zeros(Total_iter,3,3);
ELAARC_Obj=zeros(Total_iter,3,3);
MLRC_Obj=zeros(Total_iter,3,3);
SDP_LRC_Obj=zeros(Total_iter,3,3);
SDP_LRC_WC=zeros(Total_iter,3,3);
MLRC_WC=zeros(Total_iter,3,3);
ELAARC_WC=zeros(Total_iter,3,3);
AARC_WC=zeros(Total_iter,3,3);
for iter=1:Total_iter
    c=0.6*ones(I,1);
    K=10000*ones(I,1);
    eta=4+round(300*rand(I,J))/100;
    D_bar=20000*ones(J,1);
    u1=max(eta)';
    u2=max(eta,[],2);
    for iter_Gamma=1:3
        for iter_epsilon=1:3
            Gamma=0.1*(2*iter_Gamma+3)*J;
            epsilon=0.6+0.1*iter_epsilon*ones(J,1);
            D_hat=D_bar.*epsilon;
            M=sum(D_bar.*(1+epsilon));
            [true_WC,Z_opt,v_opt]=CCG(I,J,eta,c,K,D_bar,D_hat,M,Gamma,u1);
            Exact(iter,iter_epsilon,iter_Gamma)=true_WC;
            
            [U_AARC,Z_AARC,v_AARC]=Decomposition_algorithm(I,J,Gamma,eta,c,K,D_bar,D_hat,M,u1,u2,'A');
            AARC_Obj(iter,iter_epsilon,iter_Gamma)=U_AARC;
            
            [U_ELAARC,Z_ELAARC,v_ELAARC]=Decomposition_algorithm(I,J,Gamma,eta,c,K,D_bar,D_hat,M,u1,u2,'E');
            ELAARC_Obj(iter,iter_epsilon,iter_Gamma)=U_ELAARC;
            
            [U_MLRC,Z_MLRC,v_MLRC]=Decomposition_algorithm(I,J,Gamma,eta,c,K,D_bar,D_hat,M,u1,u2,'M');
            MLRC_Obj(iter,iter_epsilon,iter_Gamma)=U_MLRC;
            
            [U_SDP_LRC,Z_SDP_LRC,v_SDP_LRC]=Decomposition_algorithm(I,J,Gamma,eta,c,K,D_bar,D_hat,M,u1,u2,'S');
            SDP_LRC_Obj(iter,iter_epsilon,iter_Gamma)=U_SDP_LRC;
            
            [beta,gama,Delta]=CCG_recourse(I,J,eta,D_bar,D_hat,Gamma,Z_SDP_LRC,u1);
            SDP_LRC_WC(iter,iter_epsilon,iter_Gamma)=D_bar'*beta-D_hat'*Delta+gama'*Z_SDP_LRC-c'*Z_SDP_LRC-K'*v_SDP_LRC;
            
            [beta,gama,Delta]=CCG_recourse(I,J,eta,D_bar,D_hat,Gamma,Z_MLRC,u1);
            MLRC_WC(iter,iter_epsilon,iter_Gamma)=D_bar'*beta-D_hat'*Delta+gama'*Z_MLRC-c'*Z_MLRC-K'*v_MLRC;
            
            [beta,gama,Delta]=CCG_recourse(I,J,eta,D_bar,D_hat,Gamma,Z_ELAARC,u1);
            ELAARC_WC(iter,iter_epsilon,iter_Gamma)=D_bar'*beta-D_hat'*Delta+gama'*Z_ELAARC-c'*Z_ELAARC-K'*v_ELAARC;
            
            [beta,gama,Delta]=CCG_recourse(I,J,eta,D_bar,D_hat,Gamma,Z_AARC,u1);
            AARC_WC(iter,iter_epsilon,iter_Gamma)=D_bar'*beta-D_hat'*Delta+gama'*Z_AARC-c'*Z_AARC-K'*v_AARC;
        end
    end
end

Gap_AARC=(1-AARC_WC./Exact)*100;
Gap_ELAARC=(1-ELAARC_WC./Exact)*100;
Gap_MLRC=(1-MLRC_WC./Exact)*100;
Gap_SDP_LRC=(1-SDP_LRC_WC./Exact)*100;

Table6=[vec(mean(Gap_AARC)),vec(mean(Gap_ELAARC)),vec(mean(Gap_MLRC)),vec(mean(Gap_SDP_LRC)),...
    vec(max(Gap_AARC)),vec(max(Gap_ELAARC)),vec(max(Gap_MLRC)),vec(max(Gap_SDP_LRC))];
Table6=[Table6;mean(Table6)];
save('resultsTable6.mat','AARC_Obj','ELAARC_Obj','MLRC_Obj','SDP_LRC_Obj',...
    'SDP_LRC_WC','MLRC_WC','ELAARC_WC','AARC_WC','Gap_AARC','Gap_ELAARC',...
    'Gap_MLRC','Table6');
filename = 'Tables.xlsx';
xlswrite(filename,Table6,'Table6')
