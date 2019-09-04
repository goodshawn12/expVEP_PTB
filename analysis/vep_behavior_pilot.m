function vep_behavior(respMat)

%% ------------------------------------------------------------------------
%                       Behavioral Results
% -------------------------------------------------------------------------

contrast_list = [2,4,8];
stat_percept = zeros(5,3,20);

index_mat = ones(5,3);
for it = 1:length(respMat.stimuli)
   mode_index = respMat.stimuli(it);
   cont_index = find(respMat.contrast(it)==contrast_list);
   
   stat_percept(mode_index,cont_index,index_mat(mode_index,cont_index)) = respMat.rate_percept(it);
   index_mat(mode_index,cont_index) = index_mat(mode_index,cont_index) + 1;
end

mean_percept = mean(stat_percept,3);
std_percept = std(stat_percept,[],3);

figure
hold on
ha = bar(1:5,mean_percept);
% For each set of bars, find the centers of the bars, and write error bars
pause(0.1); %pause allows the figure to be created
for ib = 1:numel(ha)
    %XData property is the tick labels/group centers; XOffset is the offset
    %of each distinct group
    xData = ha(ib).XData+ha(ib).XOffset;
    errorbar(xData,mean_percept(:,ib),std_percept(:,ib),'k.')
end
set(gca,'XTick',[1 2 3 4 5],'XTickLabel',{'FMC', 'MSEQ', 'SSVEP', 'FMC-IMG', 'STATIC'},'FontSize',12);
xlabel('Stimuli Types'); ylabel('Perceived Intensity Rating');
