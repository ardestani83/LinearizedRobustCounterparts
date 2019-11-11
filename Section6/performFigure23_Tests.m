function performFigure23_Tests(testId,resultDir)
if testId == 1
    CorrelLevelsIndexes = [1:9];
else
    CorrelLevelsIndexes = [10];
end;

Obj=[];Time=[];
Obj1=[];Time1=[];

ns = [10];
Gammas = [0 0.1,0.3,0.5,0.7,0.9 1];
CorrelLevels = [logspace(-2,2,9) inf];
CorrelLevels = CorrelLevels(CorrelLevelsIndexes);
iterN = 100;

save_name = [resultDir sprintf('/resultsNewsvendor_Poly_n%d_correl%d_%d',ns(1),CorrelLevelsIndexes(1),CorrelLevelsIndexes(end))]


%before we start warm up the CPU
testNewsvendor5(5);

for k=1:length(ns)
    n=ns(k);
    for kk=1:length(CorrelLevels)
        correlLevel=CorrelLevels(kk);
    for kkk=1:length(Gammas)
        Gamma=Gammas(kkk)*n;
        Obj=[];Time=[];avgCorrel=[];Improv=[];
        for iter=1:iterN
            seedHere = rng;
            rng(seedHere);
            [RESULT,TIME,AVGCORREL]=testNewsvendor5(n,Gamma,correlLevel);
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
        save(save_name)
    end;
    display(sprintf('done with %d',n));
end
display(sprintf('Done with this Test %d for Figure 2 and 3',testId));