function plotPerformGivenX(Data,X_index,xs,colors,isSemiLog,printPerc,legend_names,x_label,y_label)
if ~exist('isSemiLog'), isSemiLog = 1==0; end;
if ~exist('printPerc'), isSemiLog = 1==1; end;
percs = [0.1, 0.9];
%figure;
Data2 = permute(Data,[2 X_index setdiff([1:4],[2 X_index])]);
Data2 = reshape(Data2,[size(Data,2) size(Data,X_index) prod(size(Data))/size(Data,2)/size(Data,X_index)]);
Data2 = permute(Data2,[3 2 1]);
for k=1:size(Data2,3)
    ys = Data2(:,:,k);
    %mus(1,:)=median(ys,1);
    mus(1,:)=mean(ys,1);
    tmp = sort(ys,1);
    deltas(1,:) = tmp(min(size(tmp,1),ceil(percs(2)*end)),:)-mus;
    deltas(2,:) = mus-tmp(max(1,floor(percs(1)*end)),:);
    PlotTmp(k,:,1) = mus;
    PlotTmp(k,:,2) = mus+deltas(1,:);
    if isSemiLog
    semilogy(xs,mus,colors{k});     hold on;%semilogy(xs,mus+deltas(1,:),['--' colors{k}]);
    else
    plot(xs,mus,colors{k});     hold on;%plot(xs,mus+deltas(1,:),['--' colors{k}]);
    end;
end;
legend(legend_names);
if printPerc
for k=1:size(Data2,3)
    if isSemiLog
    semilogy(xs,PlotTmp(k,:,2),['--' colors{k}]);     hold on;%semilogy(xs,mus+deltas(1,:),['--' colors{k}]);
    else
    plot(xs,PlotTmp(k,:,2),['--' colors{k}]);     hold on;%plot(xs,mus+deltas(1,:),['--' colors{k}]);
    end;
end;
end;
xlabel(x_label); ylabel(y_label);

return