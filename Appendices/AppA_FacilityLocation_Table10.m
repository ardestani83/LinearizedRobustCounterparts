function AppA_FacilityLocation_Table10

m=2; n=3;
dbar = 20000*ones(n,1);
dhat = 18000*ones(n,1);
c = 0.6*ones(m,1);
K = 100000*ones(m,1);
eta = [5.9, 5.6, 4.9; 5.6 5.9 4.9];
Gamma = 2;
M = 100000;
u1 = max(eta,[],1)';
u2 = max(eta,[],2);
u3 = u2*ones(1,n)+ones(m,1)*u1'-eta;

deltas = [0 0 0;  0 0 1; 1 0 0; 0 1 0; 1 1 0; 1 0 1; 0 1 1]';


bigM = 2*max(eta(:));
nI = m; nJ = n;
B1_ = [];
for k=1:nJ
    tmp = zeros(nI,nJ); tmp(:,k)=1;
    B1_ = [B1_; tmp(:)']; 
end;
B2_ = []; 
for k=1:nI
    tmp = zeros(nI,nJ); tmp(k,:)=1;
    B2_ = [B2_; tmp(:)']; 
end;
B = [B1_; B2_; -eye(nI*nJ)];
d = eta(:);
u = getMaxU(B,d,bigM);
u1 = u(1:n)
u2 = u(n+[1:m])
u3 = reshape(u(n+m+[1:m*n]),[m n])

%aarc solution
rome_begin; 
h_ = rome_model('AARC'); % Create Rome Model
newvar x(m,1);
newvar v(m,1) binary;
newvar deltan(n) uncertain;
newvar y(m,n,deltan') linearrule; 
d = dbar+diag(dhat)*(-deltan);
rome_maximize(-c'*x -K'*v + eta(:)'*y(:));
rome_constraint(deltan<=1);
rome_constraint(deltan>=0);
rome_constraint(sum(deltan)<= Gamma);
rome_constraint(y(:)>=0);
rome_constraint(sum(y,1)'<=d);
rome_constraint(sum(y,2)<=x);
rome_constraint(x>=0);
rome_constraint(x<=M*v);
h_.solve;
aarc.x = h_.eval(x);
aarc.v = h_.eval(v);
aarc.fval = h_.objective;
rome_end;
aarc.fval_true = evalWorstCase(aarc.x,aarc.v,n,m,c,K,eta,dbar,dhat,deltas);

%aarc with forced one warehouse
rome_begin; 
h_ = rome_model('AARC One Warehouse'); % Create Rome Model
newvar x(m,1);
newvar v(m,1) binary;
newvar deltan(n) uncertain;
newvar y(m,n,deltan') linearrule; 
d = dbar+diag(dhat)*(-deltan);
rome_maximize(-c'*x -K'*v + eta(:)'*y(:));
rome_constraint(deltan<=1);
rome_constraint(deltan>=0);
rome_constraint(sum(deltan)<= Gamma);
rome_constraint(y(:)>=0);
rome_constraint(sum(y,1)'<=d);
rome_constraint(sum(y,2)<=x);
rome_constraint(x>=0);
rome_constraint(x<=M*v);
rome_constraint(v(1)==1);
h_.solve;
aarc2.x = h_.eval(x);
aarc2.fval = h_.objective;
rome_end;

rome_begin; 
h_ = rome_model('AARC One Warehouse Pareto optimal'); % Create Rome Model
newvar x(m,1);
newvar v(m,1) binary;
newvar deltan(n) uncertain;
newvar y(m,n,deltan') linearrule; 
d = dbar+diag(dhat)*(-deltan);
rome_constraint(deltan.mean == mean(deltas,2));
rome_maximize(mean(-c'*x -K'*v + eta(:)'*y(:)));
rome_constraint(-c'*x -K'*v + eta(:)'*y(:)>=aarc2.fval);
rome_constraint(deltan<=1);
rome_constraint(deltan>=0);
rome_constraint(sum(deltan)<= Gamma);
rome_constraint(y(:)>=0);
rome_constraint(sum(y,1)'<=d);
rome_constraint(sum(y,2)<=x);
rome_constraint(x>=0);
rome_constraint(x<=M*v);
rome_constraint(v(1)==1);
h_.solve;
aarc2.x = h_.eval(x);
aarc2.meanfval = h_.objective;
tmp = h_.eval(y); tmp = tmp.linearpart;
aarc2.linearpolicy = tmp;
aarc2.y = []; for k=1:size(deltas,2), aarc2.y(:,:,k) = tmp(:,:,1); for kk=1:size(deltas,1), aarc2.y(:,:,k)=aarc2.y(:,:,k)+tmp(:,:,kk+1)*deltas(kk,k); end; end;
rome_end;

%exact solution
rome_begin; 
h_ = rome_model('Exact'); % Create Rome Model
newvar x(m,1);
newvar v(m,1) binary;
newvar y(m,n,size(deltas,2)); 
newvar t(1);
ds = dbar*ones(1,size(deltas,2))-diag(dhat)*deltas;
rome_maximize(-c'*x -K'*v + t);
rome_constraint(y(:)>=0);
for k=1:size(deltas,2)
    rome_constraint(t<= eta(:)'*reshape(y(:,:,k),[m*n 1]));
    rome_constraint(sum(y(:,:,k),1)'<=ds(:,k));
    rome_constraint(sum(y(:,:,k),2)<=x);
end;
rome_constraint(x>=0);
rome_constraint(x<=M*v);
h_.solve;
exact.x = h_.eval(x);
exact.v = h_.eval(v);
exact.fval = h_.objective;
exact.y = h_.eval(y);
exact.fvalnom = h_.eval(-c'*x -K'*v + eta(:)'*reshape(y(:,:,1),[m*n 1]));
rome_end;
exact.fval_true = evalWorstCase(exact.x,exact.v,n,m,c,K,eta,dbar,dhat,deltas);

rome_begin; 
h_ = rome_model('Exact Pareto Optimal'); % Create Rome Model
newvar x(m,1);
newvar v(m,1) binary;
newvar y(m,n,size(deltas,2)); 
newvar ts(size(deltas,2));
ds = dbar*ones(1,size(deltas,2))-diag(dhat)*deltas;
rome_maximize((1/length(ts))*sum(-c'*x -K'*v + ts));
rome_constraint(-c'*x -K'*v + ts>=exact.fval);
rome_constraint(y(:)>=0);
for k=1:size(deltas,2)
    rome_constraint(ts(k)<= eta(:)'*reshape(y(:,:,k),[m*n 1]));
    rome_constraint(sum(y(:,:,k),1)'<=ds(:,k));
    rome_constraint(sum(y(:,:,k),2)<=x);
end;
rome_constraint(x>=0);
rome_constraint(x<=M*v);
rome_constraint(v(1)==1);
h_.solve;
exact.x = h_.eval(x);
exact.meanfval = h_.objective;
exact.y = h_.eval(y);
exact.fvalnom = h_.eval(-c'*x -K'*v + eta(:)'*reshape(y(:,:,1),[m*n 1]));
rome_end;

%elaarc solution
rome_begin; 
h_ = rome_model('ELAARC'); % Create Rome Model
newvar x(m,1);
newvar v(m,1) binary;
newvar deltan(n) uncertain;
newvar y(m,n,deltan') linearrule; 
newvar z1(n,n); 
d = dbar+diag(dhat)*(-deltan);
rome_maximize(-c'*x -K'*v + eta(:)'*y(:)-u1'*z1*deltan);
rome_constraint(deltan<=1);
rome_constraint(deltan>=0);
rome_constraint(sum(deltan)<= Gamma);
rome_constraint(y(:)>=0);
rome_constraint(sum(y,1)'<=d+z1*deltan);
rome_constraint(sum(y,2)<=x);
rome_constraint(x>=0);
rome_constraint(x<=M*v);
rome_constraint(z1*deltan>=0);
h_.solve;
elaarc.x = h_.eval(x);
elaarc.v = h_.eval(v);
elaarc.fval = h_.objective;
rome_end;
elaarc.fval_true = evalWorstCase(elaarc.x,elaarc.v,n,m,c,K,eta,dbar,dhat,deltas);

%mlrc solution
rome_begin; 
h_ = rome_model('MLRC'); % Create Rome Model
newvar x(m,1);
newvar v(m,1) binary;
newvar deltan(n) uncertain;
newvar y(m,n,deltan') linearrule; 
newvar z1(n,deltan') linearrule; 
newvar z2(m,deltan') linearrule; 
newvar z3(m,n,deltan') linearrule; 
d = dbar+diag(dhat)*(-deltan);
rome_maximize(-c'*x -K'*v + eta(:)'*y(:)-u1'*z1-u2'*z2-u3(:)'*z3(:));
rome_constraint(deltan<=1);
rome_constraint(deltan>=0);
rome_constraint(sum(deltan)<= Gamma);
rome_constraint(y(:)>=0-z3(:));
rome_constraint(sum(y,1)'<=d+z1(:));
rome_constraint(sum(y,2)<=x+z2(:));
rome_constraint(x>=0);
rome_constraint(x<=M*v);
rome_constraint(z1>=0);
rome_constraint(z2>=0);
rome_constraint(z3(:)>=0);
h_.solve;
mlrc.x = h_.eval(x);
mlrc.v = h_.eval(v);
mlrc.fval = h_.objective;
rome_end;
mlrc.fval_true = evalWorstCase(mlrc.x,mlrc.v,n,m,c,K,eta,dbar,dhat,deltas);

rome_begin; 
h_ = rome_model('MLRC under nominal demand'); % Create Rome Model
newvar y(m,n); 
d = dbar;
rome_maximize(-c'*x -K'*v + eta(:)'*y(:));
rome_constraint(y(:)>=0);
rome_constraint(sum(y,1)'<=d);
rome_constraint(sum(y,2)<=x);
h_.solve;
mlrc.nom_fval = h_.objective;
rome_end;

display(sprintf('Exact solution achieves a worst-case of %d and a nominal profit of %d, with capacities loc#1:%d, loc#2:%d',exact.fval,exact.fvalnom,exact.x(1),exact.x(2)));
display(sprintf('AARC finds an optimal worst-case profit of %d, with capacities #1:%d, #2:%d',aarc.fval,aarc.x(1),aarc.x(2)))

display(sprintf('If force one to be open then AARC uses capacity %d and estimates worst-case profit as %d',max(aarc2.x),aarc2.fval));
display(sprintf('AARC policy becomes:'));
index = find(max(aarc2.x)==aarc2.x);
display(sprintf('y_11(delta):=%d+%d delta_1+ %d delta_2 + %d delta_3',aarc2.linearpolicy(index,1,1),aarc2.linearpolicy(index,1,2),aarc2.linearpolicy(index,1,3),aarc2.linearpolicy(index,1,4)));
display(sprintf('y_11(delta):=%d+%d delta_1+ %d delta_2 + %d delta_3',aarc2.linearpolicy(index,2,1),aarc2.linearpolicy(index,2,2),aarc2.linearpolicy(index,2,3),aarc2.linearpolicy(index,2,4)));
display(sprintf('y_11(delta):=%d+%d delta_1+ %d delta_2 + %d delta_3',aarc2.linearpolicy(index,3,1),aarc2.linearpolicy(index,3,2),aarc2.linearpolicy(index,3,3),aarc2.linearpolicy(index,3,4)));

tmp = [deltas' reshape(exact.y(1,:,:),[size(exact.y,2) size(exact.y,3)])';
    deltas' reshape(aarc2.y(1,:,:),[size(aarc2.y,2) size(aarc2.y,3)])'];
csvwrite('AppA_facilityLocY.csv',tmp);
dlmwrite('AppA_facilityLocY.csv',tmp,'precision',8);

display(sprintf('In comparison, ELAARC uses capacity loc#1:%d loc#2:%d and estimates worst-case profit as %d',elaarc.x(1),elaarc.x(2),elaarc.fval));
display(sprintf('In comparison, MLRC uses capacity loc#1:%d loc#2:%d and estimates worst-case profit as %d',mlrc.x(1),mlrc.x(2),mlrc.fval));

return

function fval = evalWorstCase(x,v,n,m,c,K,eta,dbar,dhat,deltas)

rome_begin; 
h_ = rome_model('Optimal'); % Create Rome Model
newvar y(m,n,size(deltas,2)); 
newvar t(1);
ds = dbar*ones(1,size(deltas,2))-diag(dhat)*deltas;
rome_maximize(-c'*x -K'*v + t);
rome_constraint(y(:)>=0);
for k=1:size(deltas,2)
    rome_constraint(t<= eta(:)'*reshape(y(:,:,k),[m*n 1]));
    rome_constraint(sum(y(:,:,k),1)'<=ds(:,k));
    rome_constraint(sum(y(:,:,k),2)<=x);
end;
h_.solve;
fval = h_.objective;
rome_end;
return

