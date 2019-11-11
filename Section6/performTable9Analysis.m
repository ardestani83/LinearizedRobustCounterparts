function performTable9Analysis(resultDir)

load([resultDir '/resultsNewsvendorCorrel_Ell_n_10_50_it1_0.mat'])
TIMES_table = reshape(mean(Times,3),[size(Times,1) size(Times,2) size(Times,4)]);
IMPROVS_table = reshape(mean(Improvs,3),[size(Improvs,1) size(Improvs,2) size(Improvs,4)]);

Times_100 = []; Improvs_100 = [];
for k=1:10
    load([resultDir sprintf('/resultsNewsvendorCorrel_Ell_n_100_100_it1_%d.mat',k)],'Times','Improvs');
    Times_100(:,:,k,:) = Times(:,:,1,:);
    Improvs_100(:,:,k,:) = Improvs(:,:,1,:);
end;
TIMES_table(:,:,end+1) = reshape(mean(Times_100,3),[size(Times_100,1) size(Times_100,2) size(Times_100,4)]);
IMPROVS_table(:,:,end+1) = reshape(mean(Improvs_100,3),[size(Improvs_100,1) size(Improvs_100,2) size(Improvs_100,4)]);

tmp = [reshape(TIMES_table,[size(TIMES_table,1) prod(size(TIMES_table(1,:,:)))])' reshape(IMPROVS_table(2:3,:,:),[2 prod(size(TIMES_table(1,:,:)))])']
csvwrite([resultDir '/NewsVendorTable.csv'],tmp)
%You can then import it in NewsVendorTable.xlsx and copy paste the latex code in the sheet.
return