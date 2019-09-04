%% addpath 
eeglab

%% load EEG and marker streams from XDF file
[filename, pathname] = uigetfile('*.xdf');

% load and integrate EEG and marker streams
EEG = pop_loadxdf([pathname, filename], 'streamtype', 'EEG', 'exclude_markerstreams', {});
EEG = pop_select( EEG,'nochannel',{'Trig1'});

% select EOG and EEG channels
EOG = pop_select( EEG,'channel',35:37); % EX3-EX5
EEG = pop_select( EEG,'channel',1:34);  % A1-A32, EX1-EX2
disp({EEG.chanlocs.labels});
pop_saveset(EEG,'filename',filename,'filepath',pathname);

%% Preprocessing
% high pass filtering using FIR
filter_hp_cutoff = 1;
EEG = pop_eegfiltnew(EEG,[],filter_hp_cutoff,[],1,0,0); % why set revfilt = 1 (invert filter)?

% low pass filtering using FIR
filter_lp_cutoff = 50;
EEG = pop_eegfiltnew(EEG,[],filter_lp_cutoff,[],0,0,0);
% 
% % remove bad channels
% EEG = clean_flatlines(EEG);
% EEG = clean_channels_nolocs(EEG);
% 
% % remove line noise using cleanLine
% EEG = pop_cleanline(EEG, 'bandwidth',2,'chanlist', 1:EEG.nbchan,'computepower',1,'linefreqs',...
%     [60 120] ,'normSpectrum',0,'p',0.01,'pad',2,'plotfigures',0,'scanforlines',1,'sigtype',...
%     'Channels','tau',100,'verb',0,'winsize',4,'winstep',1, 'ComputeSpectralPower','False');
% 
% % repair artifacts using ASR with mild threshold
% asr_stdcutoff = 20;
% EEG = clean_asr(EEG,asr_stdcutoff);

% store pre-processed EEG
EEG = eeg_checkset(EEG);


%% extract EEG epoch
epochWindow = [-0.2, 3.2];   % define time window for EEG epoch (in sec)
numTrial = 10;
expression = 'Stim(?<stim>\d+)_CT(?<cont>\d+)_Loc(?<loc>\d+)_Img(?<img>\d+)';  % decode event markers

start_event_id = 1; % find(strcmp({EEG.event.type},'Exp1'));

stimList = {'1','2','3','4','5'};   % stimulation pattern: 00: fmc-VEP (independent) 01: fmc-VEP (inverted) 02: fmc-VEP (chrome) 10: mseq
contList = {'2','8','16'};      % contrast level: +- 4, 8, 64
locList = {'1','2','3','4'};    % location: 1: top left screen, 2: top right screen, 3: bottom left, 4: bottom right

epochedEEG = cell(length(stimList),length(contList),length(locList),numTrial);
epochedIndex = zeros(length(stimList),length(contList),length(locList));
for event_it = start_event_id:length(EEG.event)
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
% save('epochedEEG_s1.mat','epochedEEG');
save([filename '.mat'],'epochedEEG');


% Archive - deal with wrong event labels
% expression = 'Stim(?<stim>\d+)_CT(?<cont>\d+)_Loc(?<loc>\d+)_Img';  % decode event markers
% stimIndex = find(strcmp(tokenNames.stim(1),stimList));
% contIndex = find(strcmp(tokenNames.stim(2:end),contList));
% locIndex = find(strcmp(tokenNames.cont,locList));
