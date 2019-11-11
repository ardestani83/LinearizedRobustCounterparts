function performFigure4_Tests(testId,resultDir)
if ~exist('resultDir'), resultDir = '.'; end;
correlK = testId;

ns = 10;
Gammas = [0 0.1,0.3,0.5,0.7,0.9 1];
CorrelLevels = [logspace(-2,2,9) inf];
if correlK>0, CorrelLevels = CorrelLevels(correlK); end;
iterN = 20; 
algNames = {'Exact', 'AARC', 'SDP-A&D', 'SDP-LRC2', 'SDP-LRC'};
save_name = [resultDir sprintf('/resultsNewsvendorCorrel_Ell_n%d_it%d_%d',ns(1),iterN,correlK)]

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
