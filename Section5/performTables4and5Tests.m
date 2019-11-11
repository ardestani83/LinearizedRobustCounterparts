clear
clc
rng(113)

I=10;
J=20;
Total_iter=100;
Exact=zeros(Total_iter,7,5);
AARC_Obj=zeros(Total_iter,7,5);
ELAARC_Obj=zeros(Total_iter,7,5);
MLRC_Obj=zeros(Total_iter,7,5);
MLRC_WC=zeros(Total_iter,7,5);
ELAARC_WC=zeros(Total_iter,7,5);
AARC_WC=zeros(Total_iter,7,5);
for iter=1:Total_iter
    c=0.6*ones(I,1);
    K=10000*ones(I,1);
    eta=4+round(300*rand(I,J))/100;
    D_bar=20000*ones(J,1);
    u1=max(eta)';
    u2=max(eta,[],2);
    for iter_Gamma=1:5
        for iter_epsilon=1:7
            Gamma=0.1*(2*iter_Gamma-1)*J;
            epsilon=0.2+0.1*iter_epsilon*ones(J,1);
            D_hat=D_bar.*epsilon;
            M=sum(D_bar.*(1+epsilon));
            [true_WC,Z_opt,v_opt]=CCG(I,J,eta,c,K,D_bar,D_hat,M,Gamma,u1);
            Exact(iter,iter_epsilon,iter_Gamma)=true_WC;
            
            [U_AARC,Z_AARC,v_AARC]=AARC(I,J,eta,c,K,D_bar,D_hat,M,Gamma);
            AARC_Obj(iter,iter_epsilon,iter_Gamma)=U_AARC;
            
            [U_ELAARC,Z_ELAARC,v_ELAARC]=ELAARC(I,J,eta,c,K,D_bar,D_hat,M,Gamma,u1);
            ELAARC_Obj(iter,iter_epsilon,iter_Gamma)=U_ELAARC;
            
            [U_MLRC,Z_MLRC,v_MLRC]=MLRC(I,J,eta,c,K,D_bar,D_hat,M,Gamma,u1,u2);
            MLRC_Obj(iter,iter_epsilon,iter_Gamma)=U_MLRC;
            
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

tmp1=sort(Gap_AARC);
tmp2=sort(Gap_ELAARC);
tmp3=sort(Gap_MLRC);

Table4=[vec(mean(Gap_AARC)),vec(mean(Gap_ELAARC)),vec(mean(Gap_MLRC)),...
    vec(tmp1(0.9*Total_iter,:,:)),vec(tmp2(0.9*Total_iter,:,:)),vec(tmp3(0.9*Total_iter,:,:)),...
    vec(max(Gap_AARC)),vec(max(Gap_ELAARC)),vec(max(Gap_MLRC))];

Table5=[[sum(Gap_AARC(:)<=0.01)/35,sum(Gap_ELAARC(:)<=0.01)/35,sum(Gap_MLRC(:)<=0.01)/35];...
    [sum(Gap_AARC(:)<=0.1)/35,sum(Gap_ELAARC(:)<=0.1)/35,sum(Gap_MLRC(:)<=0.1)/35];...
    [sum(Gap_AARC(:)<=1)/35,sum(Gap_ELAARC(:)<=1)/35,sum(Gap_MLRC(:)<=1)/35];...
    [sum(Gap_AARC(:)<=10)/35,sum(Gap_ELAARC(:)<=10)/35,sum(Gap_MLRC(:)<=10)/35];...
    mean(Table4(:,1:3));...
    [max(Gap_AARC(:)),max(Gap_ELAARC(:)),max(Gap_MLRC(:))]];
save('resultsTables4and5.mat','AARC_Obj','ELAARC_Obj','MLRC_Obj',...
    'MLRC_WC','ELAARC_WC','AARC_WC','Gap_AARC','Gap_ELAARC',...
    'Gap_MLRC','Table4','Table5');
filename = 'Tables.xlsx';
xlswrite(filename,Table4,'Table4')
xlswrite(filename,Table5,'Table5')
