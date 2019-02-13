
% addpath EEGLAB
cd('C:\Users\shawn\Desktop\NSF_BIC\expVEP_PTB')

% load EEG from XDF file
EEG = pop_loadxdf('C:\Users\shawn\Desktop\NSF_BIC\expVEP_PTB\data\VEPEXP_013119_Pilot.xdf' , 'streamtype', 'EEG', 'exclude_markerstreams', {});
EEG = pop_select( EEG,'nochannel',{'Trig1'});
pop_saveset(EEG);


%% ------------------------------------------------------------------------
%                       EEG Preprocessing
% -------------------------------------------------------------------------
filename = 'expVEP_013119_pilot.set';
filepath = 'C:\Users\shawn\Desktop\NSF_BIC\expVEP_PTB\data\';
EEG = pop_loadset([filepath filename]);

% separate EEG and EOG channels
EOG = pop_select(EEG,'channel',{'EX3','EX4','EX5','EX6','EX7','EX8'});
EEG = pop_select(EEG,'nochannel',{'EX3','EX4','EX5','EX6','EX7','EX8'});

% high pass filtering
cutoffTh = 1.0; % Hz
EEG = pop_eegfilt( EEG, cutoffTh, 0, [], 0, 0, 0, 'fir1', 0);

% apply notch filter to remove 60Hz line noise and its 2nd, 3rd, and 4th harmonics
lpFilt = designfilt('lowpassiir', 'FilterOrder', 8, ...
    'PassbandFrequency', 55, 'PassbandRipple', 0.2,...
    'SampleRate', EEG.srate);
d1 = designfilt('bandstopiir','FilterOrder',2, ...
    'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
    'DesignMethod','butter','SampleRate',EEG.srate);
d2 = designfilt('bandstopiir','FilterOrder',2, ...
                'HalfPowerFrequency1',119,'HalfPowerFrequency2',121, ...
                'DesignMethod','butter','SampleRate',EEG.srate);
d3 = designfilt('bandstopiir','FilterOrder',2, ...
                'HalfPowerFrequency1',179,'HalfPowerFrequency2',181, ...
                'DesignMethod','butter','SampleRate',EEG.srate);
d4 = designfilt('bandstopiir','FilterOrder',2, ...
                'HalfPowerFrequency1',239,'HalfPowerFrequency2',241, ...
                'DesignMethod','butter','SampleRate',EEG.srate);

tempData = zeros(EEG.nbchan,EEG.pnts);
for ch_id = 1:EEG.nbchan 
    tempData(ch_id,:) = filtfilt(d1, double(EEG.data(ch_id,:)));
%     tempData(ch_id,:) = filtfilt(lpFilt, tempData(ch_id,:));
    tempData(ch_id,:) = filtfilt(d2, tempData(ch_id,:));
    tempData(ch_id,:) = filtfilt(d3, tempData(ch_id,:));
    tempData(ch_id,:) = filtfilt(d4, tempData(ch_id,:));
    fprintf('.');
end

% visualize filter and filtered signals
% fvtool(d3,'Fs',EEG.srate)
[praw,fraw] = periodogram(EEG.data(1,:),[],[],EEG.srate);
[pbutt,fbutt] = periodogram(tempData(1,:),[],[],EEG.srate);
figure, plot(fraw, 20*log10(abs(praw)), fbutt,20*log10(abs(pbutt)),'--')

% store pre-processed EEG
EEG.data = tempData;
EEG = eeg_checkset(EEG);
pop_saveset(EEG);

% identify bad channels using clean_rawdata
% remove artifact using ASR

%% ------------------------------------------------------------------------
%     EEG Analysis - Average Impulse Response to Rising / Falling Edges
% -------------------------------------------------------------------------
filename = 'expVEP_013119_EEG.set';
filepath = 'C:\Users\shawn\Desktop\NSF_BIC\expVEP_PTB\data\';
EEG = pop_loadset([filepath filename]);

load('data\respMat_013119_pilot.mat')
code = respMat{1};

% extract code sequence and rising / falling edges
rise_fmc = cell(1,2);
fall_fmc = cell(1,2);
flat_fmc = cell(1,2);
for loc_it = 1:2
    rise_fmc{loc_it} = find((code.code_fmc{loc_it} - [0 code.code_fmc{loc_it}(1:end-1)]) == 1);
    fall_fmc{loc_it} = find((code.code_fmc{loc_it} - [0 code.code_fmc{loc_it}(1:end-1)]) == -1);
    flat_fmc{loc_it} = find((code.code_fmc{loc_it} - [0 code.code_fmc{loc_it}(1:end-1)]) == 0);
end

% expression for decoding event markers
expression = 'Stim(?<stim>\d+)_CT(?<cont>\d+)_Loc(?<loc>\d+)';
stimList = {'00','01','02','10'};
contList = {'4','8','64'};
locList = {'1','2'};
epochedEEG = cell(length(code.mode),3,2,10);
epochedIndex = zeros(length(code.mode),3,2);
epochWindow = [-0.2, 3];   % sec
for event_it = 1:length(EEG.event)
    tokenNames = regexp(EEG.event(event_it).type,expression,'names');
    if ~isempty(tokenNames)
        stimIndex = find(strcmp(tokenNames.stim,stimList));
        contIndex = find(strcmp(tokenNames.cont,contList));
        locIndex = find(strcmp(tokenNames.loc,locList));
        epochedIndex(stimIndex,contIndex,locIndex) = epochedIndex(stimIndex,contIndex,locIndex) + 1;
        
        dataIndex = floor(EEG.event(event_it).latency) + floor(epochWindow(1)*EEG.srate) : ...
            floor(EEG.event(event_it).latency) + floor(epochWindow(2)*EEG.srate);
        epochedEEG{stimIndex,contIndex,locIndex,epochedIndex(stimIndex,contIndex,locIndex)} = EEG.data(:,dataIndex);
    end
end

%% select channel and average over trials
STIM = 1; 
CONT = 3;
CH = [12,18:20]; % [5,9:15,17:21,24:26]; %[5,10:14,18:20]; %  %[12,17:21,24:26]; 
REFRESH = 60;

% extract and visualize the impulse response to rising edges
% lenSegment = 0.5;   % sec
segWindow = round([-200, 400] * EEG.srate / 1000);    % samples
data_seg_rise = []; % zeros(numSegment, round(lenSegment*EEG.srate));
data_seg_fall = []; % zeros(numSegment, round(lenSegment*EEG.srate));
data_seg_flat = [];

for trial_id = 1:10
    for LOC = 1:2
        data = mean(epochedEEG{STIM,CONT,LOC,trial_id}(CH,:),1);
        
        % extract the impulse response to rising edges
        for it = 1:length(rise_fmc{LOC})
            startIndex = round((rise_fmc{LOC}(it)-1) / REFRESH * EEG.srate) - segWindow(1);
            if startIndex+segWindow(2) <= length(data)
                data_seg_rise = [data_seg_rise; ...
                    data( startIndex+1+segWindow(1) : startIndex+segWindow(2)) ];
            end
        end
    
        % extract the impulse response to falling edges
        for it = 1:length(fall_fmc{LOC})
            startIndex = round((fall_fmc{LOC}(it)-1) / REFRESH * EEG.srate) - segWindow(1);
            if startIndex+segWindow(2) <= length(data)
                data_seg_fall = [data_seg_fall; ...
                    data( startIndex+1+segWindow(1) : startIndex+segWindow(2)) ];
            end
        end
        
        % extract the impulse response to no edges
        for it = 1:length(flat_fmc{LOC})
            startIndex = round((flat_fmc{LOC}(it)-1) / REFRESH * EEG.srate) - segWindow(1);
            if startIndex+segWindow(2) <= length(data)
                data_seg_flat = [data_seg_flat; ...
                    data( startIndex+1+segWindow(1) : startIndex+segWindow(2)) ];
            end
        end
    end
end
avg_data_seg_rise = mean(data_seg_rise,1);
avg_data_seg_fall = mean(data_seg_fall,1);
avg_data_seg_flat = mean(data_seg_flat,1);
avg_data_seg = mean([data_seg_rise; -data_seg_fall],1);

% plot the average impulse response
figure, subplot(2,1,1); hold on, 
plot((segWindow(1)+1:segWindow(2))./EEG.srate*1000,avg_data_seg_rise,'b')
plot((segWindow(1)+1:segWindow(2))./EEG.srate*1000,avg_data_seg_fall,'r')
xlim([-200 400]); legend('rising edge','falling edge'); title('Avg Impulse Response (IR)')
subplot(2,1,2);
plot((segWindow(1)+1:segWindow(2))./EEG.srate*1000, avg_data_seg, 'k')
xlabel('Time (msec)'); ylabel('Amplitude (\muV)'); xlim([-200 400]); title('IR(rising)-IR(falling)')
% plot((segWindow(1)+1:segWindow(2))./EEG.srate*1000,avg_data_seg_flat,'k')

% plot the power spectra 
[prise,frise] = periodogram(avg_data_seg_rise,[],[],EEG.srate);
[pfall,ffall] = periodogram(avg_data_seg_fall,[],[],EEG.srate);
[pdiff,fdiff] = periodogram(avg_data_seg,[],[],EEG.srate);

figure, subplot(2,1,1); hold on, 
plot(frise,20*log10(abs(prise)),'b')
plot(ffall,20*log10(abs(pfall)),'r')
legend('rising edge','falling edge'); title('Avg Impulse Response (IR)')
subplot(2,1,2);
plot(fdiff,20*log10(abs(pdiff)),'k')
xlabel('Frequency (Hz)'); ylabel('Spectral Power (dB)'); title('IR(rising)-IR(falling)')


%% Template matching correlation over time
testCH = 12;

% extract and visualize the impulse response to rising edges
template_rising = avg_data_seg_rise;
template_falling = avg_data_seg_fall;
template_flat = avg_data_seg_flat;
template_diff = avg_data_seg;

len_template = length(avg_data_seg_rise);
len_data = size(epochedEEG{1,1,1,1},2);
corr_rising = zeros(1,len_data-len_template+1);
corr_falling= zeros(1,len_data-len_template+1);
corr_flat = zeros(1,len_data-len_template+1);
corr_diff = zeros(1,len_data-len_template+1);
for trial_id = 1:10
    for LOC = 1:2
        data = mean(epochedEEG{STIM,CONT,LOC,trial_id}(testCH,:),1);
        
        for it = 1:len_data-len_template+1
            corr_rising(it) = corr(data(it:it+len_template-1)',template_rising');
            corr_falling(it) = corr(data(it:it+len_template-1)',template_falling');
            corr_flat(it) = corr(data(it:it+len_template-1)',template_flat');
        end
        
    end
end

figure, hold on,
plot(corr_rising,'b');
plot(corr_falling,'r');
plot(corr_flat,'k');

numCode = floor((len_data-len_template+1)/EEG.srate*REFRESH);
decodedSeqShift = cell(1,floor(EEG.srate/REFRESH));
meanCorr = zeros(1,floor(EEG.srate/REFRESH));
for shift = 1:floor(EEG.srate/REFRESH)
    codeInTime = round((1:numCode)/REFRESH*EEG.srate)-shift+1;
    corrTemp = [corr_rising(codeInTime);corr_flat(codeInTime);corr_falling(codeInTime)];
    [maxCorr,decodedSeqShift{shift}] = max(corrTemp,[],1);
    meanCorr(shift) = mean(maxCorr);
end
[~,opt_shift] = max(meanCorr);


decodedSeq = decodedSeqShift{opt_shift};
figure, plot(decodedSeq);

codeSeq = zeros(1,length(decodedSeq));
for it = 1:length(decodedSeq)
    if decodedSeq(it) == 1  % rising edge
        codeSeq(it) = 1;
    elseif decodedSeq(it) == 3  % falling edge
        codeSeq(it) = 0;
    elseif decodedSeq(it) == 2 && it > 1 % code: 2 no change
        codeSeq(it) = codeSeq(it-1);
    else
        codeSeq(it) = 0;
    end
end

corrSlide = zeros(1,length(code.code_fmc{2})-length(decodedSeq));
for it = 1:length(code.code_fmc{2})-length(decodedSeq)
    corrSlide(it) = corr(code.code_fmc{2}(it:it+length(decodedSeq)-1)', codeSeq');
end

[~,opt_align] = max(corrSlide);

actualCode = code.code_fmc{2}(opt_align:opt_align+length(decodedSeq)-1);
errorSeq = actualCode ~= codeSeq;
errorRate = mean(errorSeq);

%{
[peaks_rising,locs_rising] = findpeaks(corr_rising,'MinPeakHeight',0,'MinPeakDistance',8);
[peaks_falling,locs_falling] = findpeaks(corr_falling,'MinPeakHeight',0,'MinPeakDistance',8);
[peaks_flat,locs_flat] = findpeaks(corr_flat,'MinPeakHeight',0,'MinPeakDistance',8);

figure, hold on,
stem(locs_rising,peaks_rising,'b');
stem(locs_falling,peaks_falling,'r');
stem(locs_flat,peaks_flat,'k');
legend('Rising','Falling','Flat'); xlabel('Samples'); ylabel('Correlation');

[sortedLoc, sortedIndex] = sort([locs_rising, locs_falling]);
decoder = [ones(1,length(locs_rising)),zeros(1,length(locs_falling))];
decoder = decoder(sortedIndex);
figure, stem(sortedLoc,decoder);

tmp = mod(sortedLoc, EEG.srate/REFRESH);
figure, hist(tmp,20);
%}



%% ------------------------------------------------------------------------ 
%       Performing the TRCA-based VEP detection algorithm
% ------------------------------------------------------------------------- 

fprintf('Results of the ensemble TRCA-based method.\n');

% Preparing data
STIM = 1;
CONT = 1;
CH = [5,10:14,18:20];
NTRIAL = 10;
is_ensemble = 1;

data_len_list = [0.05:0.05:3.2];
data_offset = 0; % linspace(0,floor((epochWindow(2)-epochWindow(1)-data_len)*EEG.srate), 25);

var_list = data_len_list;
mean_acc = zeros(4,length(var_list));
mu_ci = zeros(4,length(var_list),2);
for stim_i = 1:4
    
    STIM = stim_i;
    
    for var_i = 1:length(var_list)
        
        data_len = data_len_list(var_i);
        data_range = data_offset + (1:1:floor(data_len*EEG.srate));
        NSAMP = length(data_range);
        
        trial_eeg = zeros(2,length(CH),NSAMP,NTRIAL);
        for loc_i = 1:2
            for tr_i = 1:10
                trial_eeg(loc_i,:,:,tr_i) = epochedEEG{STIM,CONT,loc_i,tr_i}(CH,data_range);
            end
        end
        
        % Leave-one-trial-out cross validation classification accuracy
        labels = [1, 2];
        for cv_i = 1:1:NTRIAL
            
            % Training stage
            traindata = trial_eeg;
            traindata(:, :, :, cv_i) = [];
            model = train_trca(traindata);
            
            % Test stage
            testdata = squeeze(trial_eeg(:, :, :, cv_i));
            estimated = test_trca(testdata, model, is_ensemble);
            
            % Evaluation
            is_correct = (estimated==labels);
            accs(cv_i) = mean(is_correct)*100;
            %         fprintf('Trial %d: Accuracy = %2.2f%%\n', cv_i, accs(cv_i));
            
        end % loocv_i
        
        % Summarize
        alpha_ci = 0.05;
        ci = 1-alpha_ci;
        [mu, ~, muci, ~] = normfit(accs, alpha_ci);
        fprintf('Mean accuracy = %2.2f %% (%2d%% CI: %2.2f - %2.2f %%)\n',...
            mu, ci, muci(1), muci(2));
        
        mean_acc(stim_i,var_i) = mu;
        mu_ci(stim_i,var_i,:) = muci;
    end
end
% figure, plot(data_offset_list/EEG.srate + epochWindow(1), mean_acc,'linewidth',2);
% xlabel('Window Delay (sec)'); ylabel('Cross validation accuracy (''%)');
% set(gca,'fontsize',12)

figure, plot(data_len_list, mean_acc,'linewidth',2);
xlabel('Trial length (sec)'); ylabel('Cross validation accuracy (''%)');
set(gca,'fontsize',12)
ylim([0 100]); 
% legend('FMC-idp','FMC-opp','MSEQ');
legend('FMC-idp','FMC-opp','FMC-chm','MSEQ');
% legend('Cont 4','Cont 8','Cont 64');

%% ------------------------------------------------------------------------
%                   Behavioral Results
% -------------------------------------------------------------------------

% load behavioral results
load('data\respMat_013119_pilot.mat')
code = respMat{1};
resp = respMat{2};

contrast_list = [4,8,64];

stat_percept = zeros(4,3,20);
stat_comfort = zeros(4,3,20);
index_mat = ones(4,3);
for it = 1:length(resp.stimuli)
   mode_index = 3 * resp.stimuli{it}(1) + resp.stimuli{it}(2) + 1;
   cont_index = find(resp.contrast(it)==contrast_list);
   
   stat_percept(mode_index,cont_index,index_mat(mode_index,cont_index)) = resp.rate_percept(it);
   stat_comfort(mode_index,cont_index,index_mat(mode_index,cont_index)) = resp.rate_comfort(it);
   index_mat(mode_index,cont_index) = index_mat(mode_index,cont_index) + 1;
end

mean_percept = mean(stat_percept,3);
std_percept = std(stat_percept,[],3);

mean_comfort = mean(stat_comfort,3);
std_comfort = std(stat_comfort,[],3);

figure
hold on
ha = bar(1:4,mean_percept);
% For each set of bars, find the centers of the bars, and write error bars
pause(0.1); %pause allows the figure to be created
for ib = 1:numel(ha)
    %XData property is the tick labels/group centers; XOffset is the offset
    %of each distinct group
    xData = ha(ib).XData+ha(ib).XOffset;
    errorbar(xData,mean_percept(:,ib),std_percept(:,ib),'k.')
end
set(gca,'XTick',[1 2 3 4],'XTickLabel',{'Normal', 'Opposite', 'Chrome', 'Mseq'},'FontSize',12);
xlabel('Stimuli Types'); ylabel('Perceived Intensity Rating');

figure
hold on
hb = bar(1:4,mean_comfort);
% For each set of bars, find the centers of the bars, and write error bars
pause(0.1); %pause allows the figure to be created
for ib = 1:numel(hb)
    %XData property is the tick labels/group centers; XOffset is the offset
    %of each distinct group
    xData = hb(ib).XData+hb(ib).XOffset;
    errorbar(xData,mean_comfort(:,ib),std_comfort(:,ib),'k.')
end
ylim([2 6]);
set(gca,'XTick',[1 2 3 4],'XTickLabel',{'Normal', 'Opposite', 'Chrome', 'Mseq'},'FontSize',12);
xlabel('Stimuli Types'); ylabel('Comfortability Rating');