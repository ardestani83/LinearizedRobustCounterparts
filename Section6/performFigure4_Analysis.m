function performFigure4_Analysis(resultDir)

fileNs = [1:10];%[2 3 4 5 6 7 9];
files = {}; for k=1:length(fileNs), files{k} = [resultDir sprintf('/resultsNewsvendorCorrel_Ell_n10_it20_%d.mat',fileNs(k))]; end;
alg_names = {'C&CG','AARC','SDP-A&D','SDP-LRC2','SDP-LRC'};
colors = {'k','b','r','m','g'};
gamma_x_label = 'Size of ellipsoidal set (\gamma)';
maxGamma_f = @(n) sqrt(n);

[Improvs,Times,avgCorrels,Obj_data,Improvs_data,Times_data,CorrelLevels,Gammas,n] = mergeMatResults(files);
Gammas = Gammas*maxGamma_f(n);

%move C&CG at end of list
Improvs_data_ = Improvs_data(:,[2:end 1],:,:);
Times_data_ = Times_data(:,[2:end 1],:,:);
alg_names = alg_names([2:end 1]);
colors = colors([2:end 1]);

%plot performance vs Gamma
figure;
subplot(1,2,1);plotPerformGivenX(-Improvs_data_(:,:,:,:),3,Gammas,colors(:),1==0,1==1,alg_names(:),gamma_x_label,'Optimality gap (in %)')
subplot(1,2,2);plotPerformGivenX(Times_data_,3,Gammas,colors,1==1,1==0,alg_names,gamma_x_label,'Solution Time (in sec)')