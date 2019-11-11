clear
clc
rng(113)

I=10;
J=20;
Total_iter=100;
Exact=zeros(Total_iter,9);
AARC_Obj=zeros(Total_iter,9);
ELAARC_Obj=zeros(Total_iter,9);
MLRC_Obj=zeros(Total_iter,9);
MLRC_WC=zeros(Total_iter,9);
ELAARC_WC=zeros(Total_iter,9);
AARC_WC=zeros(Total_iter,9);
for iter=1:Total_iter
    c=0.6*ones(I,1);
    K=10000*ones(I,1);
    eta=4+round(300*rand(I,J))/100;
    D_bar=20000*ones(J,1);
    epsilon=0.3+0.6*rand(J,1);
    D_hat=D_bar.*epsilon;
    M=sum(D_bar.*(1+epsilon));
    u1=max(eta)';
    u2=max(eta,[],2);
    for j=1:9
        Gamma=0.1*j*J;
        [true_WC,Z_opt,v_opt]=CCG(I,J,eta,c,K,D_bar,D_hat,M,Gamma,u1);
        Exact(iter,j)=true_WC;
        
        [U_AARC,Z_AARC,v_AARC]=AARC(I,J,eta,c,K,D_bar,D_hat,M,Gamma);
        AARC_Obj(iter,j)=U_AARC;
        
        [U_ELAARC,Z_ELAARC,v_ELAARC]=ELAARC(I,J,eta,c,K,D_bar,D_hat,M,Gamma,u1);
        ELAARC_Obj(iter,j)=U_ELAARC;
        
        [U_MLRC,Z_MLRC,v_MLRC]=MLRC(I,J,eta,c,K,D_bar,D_hat,M,Gamma,u1,u2);
        MLRC_Obj(iter,j)=U_MLRC;
        
        [beta,gama,Delta]=CCG_recourse(I,J,eta,D_bar,D_hat,Gamma,Z_MLRC,u1);
        MLRC_WC(iter,j)=D_bar'*beta-D_hat'*Delta+gama'*Z_MLRC-c'*Z_MLRC-K'*v_MLRC;
        
        [beta,gama,Delta]=CCG_recourse(I,J,eta,D_bar,D_hat,Gamma,Z_ELAARC,u1);
        ELAARC_WC(iter,j)=D_bar'*beta-D_hat'*Delta+gama'*Z_ELAARC-c'*Z_ELAARC-K'*v_ELAARC;
        
        [beta,gama,Delta]=CCG_recourse(I,J,eta,D_bar,D_hat,Gamma,Z_AARC,u1);
        AARC_WC(iter,j)=D_bar'*beta-D_hat'*Delta+gama'*Z_AARC-c'*Z_AARC-K'*v_AARC;
    end
end

Gap_AARC=(1-AARC_WC./Exact)'*100;
Gap_ELAARC=(1-ELAARC_WC./Exact)'*100;
Gap_MLRC=(1-MLRC_WC./Exact)'*100;

tmp1=sort(Gap_AARC,2);
tmp2=sort(Gap_ELAARC,2);
tmp3=sort(Gap_MLRC,2);

Table2=[mean(Gap_AARC,2),mean(Gap_ELAARC,2),mean(Gap_MLRC,2),...
    tmp1(:,0.9*Total_iter),tmp2(:,0.9*Total_iter),tmp3(:,0.9*Total_iter),...
    max(Gap_AARC,[],2),max(Gap_ELAARC,[],2),max(Gap_MLRC,[],2);...
    zeros(1,9)];

Table3=[[(sum(Gap_AARC(:)<=0.01)+100)/10,(sum(Gap_ELAARC(:)<=0.01)+100)/10,(sum(Gap_MLRC(:)<=0.01)+100)/10];...
    [(sum(Gap_AARC(:)<=0.1)+100)/10,(sum(Gap_ELAARC(:)<=0.1)+100)/10,(sum(Gap_MLRC(:)<=0.1)+100)/10];...
    [(sum(Gap_AARC(:)<=1)+100)/10,(sum(Gap_ELAARC(:)<=1)+100)/10,(sum(Gap_MLRC(:)<=1)+100)/10];...
    [(sum(Gap_AARC(:)<=10)+100)/10,(sum(Gap_ELAARC(:)<=10)+100)/10,(sum(Gap_MLRC(:)<=10)+100)/10];...
    mean(Table2(:,1:3));...
    [max(Gap_AARC(:)),max(Gap_ELAARC(:)),max(Gap_MLRC(:))]];
save('resultsTables2and3.mat','AARC_Obj','ELAARC_Obj','MLRC_Obj',...
    'MLRC_WC','ELAARC_WC','AARC_WC','Gap_AARC','Gap_ELAARC',...
    'Gap_MLRC','Table2','Table3');
filename = 'Tables.xlsx';
xlswrite(filename,Table2,'Table2')
xlswrite(filename,Table3,'Table3')
