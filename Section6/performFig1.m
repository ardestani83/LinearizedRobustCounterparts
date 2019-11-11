n=10;
correlLevels = [logspace(-2,2,9) inf];

allSigmas = []; allCorrels = [];
for kk=1:length(correlLevels)
    for k=1:1e3
        correlLevel = correlLevels(kk);
        if isinf(correlLevel)
            mixingMat = eye(n);
        else
            
            dbar=100*rand(n,1);
            dhat=dbar.*rand(n,1);
            sigmas = dhat;
            
            correlTmp = vineBeta(n, correlLevel);
            avgCorrel = sum(sum(abs(correlTmp-eye(n))))/(n^2-n);
            sigmas = dhat;
            covarTmp = (sigmas*sigmas').*correlTmp;
            sqCovar = real(sqrtm(covarTmp));
            mixingMat = sqCovar./(diag(sqCovar)*ones(1,n));
            mixingMat=mixingMat./(sum(abs(mixingMat),2)*ones(1,n));
            
            
        end;
        tmp = diag(dhat)*mixingMat; covarTmp2 = tmp*tmp';
        [sigmas2, correlTmp2] = cov2corr(covarTmp2);
        avgCorrel2 = sum(sum(abs(correlTmp2-eye(n))))/(n^2-n);
        
        sigmas2 = dhat.*sum(abs(mixingMat),2);
        
        
        allSigmas = [allSigmas; sigmas sigmas2];
        allCorrels = [allCorrels; avgCorrel avgCorrel2];
    end;
end;

figure; cdfplot(allCorrels(:,2))
