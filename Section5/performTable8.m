clear
clc

Facility_number=[2,2,5];
Customers_number=[2,4,10];
Time=zeros(6,3,7);
for instance=1:3
    S=load('InputDataTables7and8.mat');
    I=Facility_number(instance);
    J=Customers_number(instance);
    Gamma=0.2*J;
    c=S.c(1:I);
    K=S.K(1:I);
    eta=S.eta(1:I,1:J);
    D_bar=S.D_bar(1:J);
    u1=max(eta)';
    u2=max(eta,[],2);
    for epsilon_iter=1:5
        
        epsilon=0.1*(2*epsilon_iter-1)*ones(J,1);
        D_hat=D_bar.*epsilon;
        M=sum(D_bar+D_hat);
        
        tic;
        AARC(I,J,eta,c,K,D_bar,D_hat,M,Gamma);
        Time(epsilon_iter,instance,1)=toc;
        if Time(epsilon_iter,instance,1)>=172800, Time(epsilon_iter,instance,1)=nan; end
        
        tic
        ELAARC(I,J,eta,c,K,D_bar,D_hat,M,Gamma,u1);
        Time(epsilon_iter,instance,2)=toc;
        if Time(epsilon_iter,instance,2)>=172800, Time(epsilon_iter,instance,2)=nan; end
        
        tic
        MLRC(I,J,eta,c,K,D_bar,D_hat,M,Gamma,u1,u2);
        Time(epsilon_iter,instance,3)=toc;
        if Time(epsilon_iter,instance,3)>=172800, Time(epsilon_iter,instance,3)=nan; end
        
        tic
        Decomposition_algorithm(I,J,Gamma,eta,c,K,D_bar,D_hat,M,u1,u2,'A');
        Time(epsilon_iter,instance,4)=toc;
        if Time(epsilon_iter,instance,4)>=172800, Time(epsilon_iter,instance,4)=nan; end
        
        tic
        Decomposition_algorithm(I,J,Gamma,eta,c,K,D_bar,D_hat,M,u1,u2,'E');
        Time(epsilon_iter,instance,5)=toc;
        if Time(epsilon_iter,instance,5)>=172800, Time(epsilon_iter,instance,5)=nan; end
        
        tic
        Decomposition_algorithm(I,J,Gamma,eta,c,K,D_bar,D_hat,M,u1,u2,'M');
        Time(epsilon_iter,instance,6)=toc;
        if Time(epsilon_iter,instance,6)>=172800, Time(epsilon_iter,instance,4)=nan; end
        
        tic;
        CCG(I,J,eta,c,K,D_bar,D_hat,M,Gamma,u1);
        Time(epsilon_iter,instance,7)=toc;
        if Time(epsilon_iter,instance,7)>=172800, Time(epsilon_iter,instance,7)=nan; end
        
    end
end
Time(6,:,:)=mean(Time(1:5,:,:));
Table8=[vec(Time(:,:,1)),vec(Time(:,:,2)),vec(Time(:,:,3)),...
    vec(Time(:,:,4)),vec(Time(:,:,5)),vec(Time(:,:,6)),vec(Time(:,:,7))];
save('resultsTables8.mat','Table8');
filename = 'Tables.xlsx';
xlswrite(filename,Table8,'Table8')