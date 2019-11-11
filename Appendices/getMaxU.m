function u = getMaxU(B,d,M)
if ~exist('M'), M = 5000; end;
[nlambda,ny] = size(B);

M = M+1;

%try to obtain a bound without employing MILP
for k=1:nlambda
    rome_begin;
    h_ = rome_model('maxU_relax'); % Create Rome Model
    newvar lambda(nlambda);
    rome_maximize(lambda(k));
    rome_constraint(B'*lambda == d);
    rome_constraint(lambda>=0);
    h_.solve;
    fval = h_.objective;
    rome_end;
    Ms(k,1) = min(M,fval);
end;

%check that B is full rank
if ~testLinearIndepRows(B')
    error('B prime matrix should have linearly independent rows');
end;

%identify max u using MILP & test for linear independence progressively
nCuts = [];
for k=1:nlambda
    
    cutsA = zeros(0,nlambda); cutsb = zeros(0,1);
    while 1==1
        rome_begin;
        h_ = rome_model('maxU'); % Create Rome Model
        newvar lambda(nlambda);
        newvar v(nlambda) binary;
        rome_maximize(lambda(k));
        rome_constraint(B'*lambda == d);
        rome_constraint(lambda>=0);
        rome_constraint(lambda<=Ms.*(1-v));
        rome_constraint(sum(v)==nlambda-ny);
        rome_constraint(cutsA*v<=cutsb);
        h_.solve;
        lambda_ = h_.eval(lambda);
        v_ = h_.eval(v);
        fval = h_.objective;
        rome_end;
        
        if isempty(cutsb), u_(k,1) = fval; end;
        tmp = eye(nlambda);
        Btmp = [B';tmp(find(v_==1),:)];
        
        if testLinearIndepRows(Btmp), break; end;
        %cutsA(end+1,:)=ones(1,nlambda); cutsA(end,find(v_==0))=-1;
        %cutsb(end+1,1)=sum(v_)-1;
        cutsA(end+1,:)=ones(1,nlambda); cutsA(end,find(v_==0))=0;
        cutsb(end+1,1)=sum(v_)-1;
        
    end;
    nCuts(k)=length(cutsb);
    u(k,1) = fval;
    
    if 1==0
        %do we need to check for linearly independencei in LP?
        M2 = 100;
        rome_begin;
        h_ = rome_model('maxU'); % Create Rome Model
        newvar lambda(nlambda);
        newvar v(nlambda) binary;
        newvar y(nlambda,nlambda);
        rome_maximize(lambda(k));
        rome_constraint(B'*lambda == d);
        rome_constraint(lambda>=0);
        rome_constraint(lambda<=Ms.*(1-v));
        rome_constraint(sum(v)==nlambda-ny);
        rome_constraint(B'*y==0);
        tmp = []; tmp2 = eye(nlambda);
        for k=1:nlambda, tmp(k,:)=vec(tmp2(:,k)*tmp2(:,k)'); end;
        rome_constraint(tmp*y(:) >= v);
        for k=1:nlambda
            tmp = []; indexes = [[1:k-1] [k+1:nlambda]];
            rome_constraint(tmp2(indexes,:)*y(:,k)<=M2*(1-v(indexes)));
            rome_constraint(tmp2(indexes,:)*y(:,k)>=-M2*(1-v(indexes)));
        end;
        h_.solve;
        lambda_ = h_.eval(lambda);
        v_ = h_.eval(v);
        y_ = h_.eval(y);
        fval = h_.objective;
        rome_end;
        u(k,1) = fval;
    end;
    
    
end;

if max(u)==M, error('Big M does not seem large enough'); end;

return

function linInd = testLinearIndepRows(A)
if size(A,1)==size(A,2)
    linInd = det(A) ~= 0;
elseif 1==1
    linInd = rank(A)==size(A,1);
else
    linInd = 1==1;
    for k=1:size(A,1)
rome_begin;
        h_ = rome_model('linInd'); % Create Rome Model
        newvar y(size(A,1),1);
        newvar s(1);
        rome_minimize(s);
        rome_constraint(A'*y <= s);
        rome_constraint(A'*y >= -s);
        rome_constraint(y(k) == 1);
        h_.solve;
        y_ = h_.eval(y);
        fval = h_.objective;
        rome_end;    
        if fval<=sqrt(eps), linInd = 1==0; break; end;
    end;
end;
return

function M = getM(B,d)
%based on corollary 3.2d of Schrijver's Theory of Linear and Integer
%Programming
%this theory is useless because it returns infinity most of the time
nn = size(B,1);
[num,denum] = rat([B' d; eye(nn) zeros(nn,1)]);
sizeBd = 1+ceil(log(abs(num)+1)/log(2))+ceil(log(denum+1)/log(2));
tmp = sum(sizeBd,2);
tmp = sort(tmp(:),1,'descend');
psi = nn^2+sum(tmp(1:nn));
M = 2^(4*nn*psi-2)-1;
return