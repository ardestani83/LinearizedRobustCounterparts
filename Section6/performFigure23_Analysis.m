function performFigure23_Analysis(resultDir)

%Polyhedral 
files = {[resultDir '/resultsNewsvendor_Poly_n10_correl1_9'],
    [resultDir '/resultsNewsvendor_Poly_n10_correl10_10']};
alg_names = {'C&CG','AARC','SDP-A&D/SDP-LRC2','SDP-LRC'};
colors = {'k','b','r','g'};
gamma_x_label = 'Size of budgeted set (\Gamma)';
maxGamma_f = @(n) n;

[Improvs,Times,avgCorrels,Obj_data,Improvs_data,Times_data,CorrelLevels,Gammas,n] = mergeMatResults(files);
algN = length(alg_names);
Gammas = Gammas*maxGamma_f(n);

%move C&CG at end of list
Improvs_data_ = Improvs_data(:,[2:end 1],:,:);
Times_data_ = Times_data(:,[2:end 1],:,:);
alg_names = alg_names([2:end 1]);
colors = colors([2:end 1]);

%plot performance vs correlation
figure;
subplot(1,2,1);plotPerformGivenX(-Improvs_data_,4,mean(avgCorrels(:,1:size(Improvs,3)),1),colors,1==0,1==1,alg_names,'Correlation levels','Optimality gap (in %)')
subplot(1,2,2);plotPerformGivenX(Times_data_,4,mean(avgCorrels(:,1:size(Improvs,3)),1),colors,1==1,1==0,alg_names,'Correlation levels','Solution Time (in sec)')

%plot performance vs Gamma
figure;
subplot(1,2,1);plotPerformGivenX(-Improvs_data_,3,Gammas,colors,1==0,1==1,alg_names,gamma_x_label,'Optimality gap (in %)')
subplot(1,2,2);plotPerformGivenX(Times_data_,3,Gammas,colors,1==1,1==0,alg_names,gamma_x_label,'Solution Time (in sec)')