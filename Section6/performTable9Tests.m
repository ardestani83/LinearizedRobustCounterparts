function performTable9Tests(testId,resultDir)
%% this function generates the results for Table 9
%% it must be called with testId = 1,2,...,11
%% resultDir should be the name of the directory in which results 
%% should be saved

if testId==1
    correlK=0;
    ns = [10,25,50];
else
    correlK=testId-1;
    ns = [100];
end;
    
Gammas = [0.1,0.3,0.5,0.7,0.9];
CorrelLevels = [logspace(-2,2,9) inf];
if correlK>0, CorrelLevels = CorrelLevels(correlK); end;
iterN = 1;
algNames = {'AARC', 'SDP-A&D', 'SDP-LRC2'};
if ~isdir(resultDir), mkdir(resultDir), end;
save_name = [resultDir sprintf('/resultsNewsvendorCorrel_Ell_n_%d_%d_it%d_%d',ns(1),ns(end),iterN,correlK)]

Obj=[];Time=[];
Obj1=[];Time1=[];
%before we start warm up the CPU
testNewsvendorEll3(5);
t1=cputime
for k=1:length(ns)
    n=ns(k);
    for kk=1:length(CorrelLevels)
        correlLevel=CorrelLevels(kk)
    for kkk=1:length(Gammas)
        Gamma=Gammas(kkk)*sqrt(n);
        Obj=[];Time=[];avgCorrel=[];Improv=[];
        for iter=1:iterN
            seedHere = rng;
            rng(seedHere);
            [RESULT,TIME,AVGCORREL]=testNewsvendorEll3(n,Gamma,correlLevel,algNames);
            IMPROV = (RESULT-RESULT(1))./RESULT(1)*100;
            if isinf(max(TIME)) 
                save debug; 
            end;
            Time=[Time;TIME];
            Obj=[Obj;RESULT];
            Improv = [Improv; IMPROV];
            avgCorrel = [avgCorrel;AVGCORREL];
        end
        Improvs(:,kkk,kk,k) = mean(Improv,1);
        Times(:,kkk,kk,k) = mean(Time,1);
        avgCorrels(kkk,kk,k) = mean(avgCorrel,1);
        
        Obj_data(:,:,kkk,kk,k) = Obj;
        Improvs_data(:,:,kkk,kk,k) = Improv;
        Times_data(:,:,kkk,kk,k) = Time;
        display(sprintf('Done with Gamma %d',Gamma));
        end;
        display(sprintf('done with correl %d',kk));
        save(save_name);
    end
    display(sprintf('done with %d after %d sec',n,cputime-t1));
end
display(sprintf('Done with this Test %d for Table 9',testId));
