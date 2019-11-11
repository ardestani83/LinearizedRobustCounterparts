function [Improvs_,Times_,avgCorrels_,Obj_data_,Improvs_data_,Times_data_,CorrelLevels_,Gammas,ns] = mergeMatResults(files)
Improvs_ = zeros(0,0,0,0);
Times_ = zeros(0,0,0,0);
avgCorrels_ = zeros(0,0);
Obj_data_ = zeros(0,0,0,0);
Improvs_data_ = zeros(0,0,0,0);
Times_data_ = zeros(0,0,0,0);
CorrelLevels_ = zeros(0,0);


for k=1:length(files)
    load(files{k});
    L1=size(Improvs_,3);
    L2=size(Improvs,3);
    Improvs_(:,:,L1+1:L1+L2,:) = Improvs;
    Times_(:,:,L1+1:L1+L2,:) = Times;
    avgCorrels_(:,L1+1:L1+L2) = avgCorrels;
    CorrelLevels_(L1+1:L1+L2) = CorrelLevels(1:L2);
        
    Obj_data_(:,:,:,L1+1:L1+L2) = Obj_data;
    Improvs_data_(:,:,:,L1+1:L1+L2) = Improvs_data;
    Times_data_(:,:,:,L1+1:L1+L2) = Times_data;
end;
return