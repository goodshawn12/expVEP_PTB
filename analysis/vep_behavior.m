function vep_behavior(respMat)

%% ------------------------------------------------------------------------
%                       Behavioral Results
% -------------------------------------------------------------------------

contrast_list = unique(respMat.contrast); 
stimuli_list = unique(respMat.stimuli);
nStim = length(stimuli_list);
nCont = length(contrast_list);
nTr = length(respMat.stimuli);
nTrPerCond = nTr/nStim/nCont;
stat_percept = zeros(nStim,nCont,nTrPerCond);

index_mat = ones(nStim,nCont);
for it = 1:nTr
   mode_index = respMat.stimuli(it);
   cont_index = find(respMat.contrast(it)==contrast_list);
   
   stat_percept(mode_index,cont_index,index_mat(mode_index,cont_index)) = respMat.rate_percept(it);
   index_mat(mode_index,cont_index) = index_mat(mode_index,cont_index) + 1;
end

mean_percept = mean(stat_percept,3);
std_percept = std(stat_percept,[],3);

figure
hold on
ha = bar(1:nStim,mean_percept);
% For each set of bars, find the centers of the bars, and write error bars
pause(0.1); %pause allows the figure to be created
for ib = 1:numel(ha)
    %XData property is the tick labels/group centers; XOffset is the offset
    %of each distinct group
    xData = ha(ib).XData+ha(ib).XOffset;
    errorbar(xData,mean_percept(:,ib),std_percept(:,ib),'k.')
end
set(gca,'XTick',[1 2 3 4 5],'XTickLabel',{'FMC', 'MSEQ', 'SSVEP', 'FMC-IMG'},'FontSize',12);
xlabel('Stimuli Types'); ylabel('Perceived Intensity Rating');
legend('2','8','16');
