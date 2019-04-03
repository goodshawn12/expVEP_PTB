%% load EEG from XDF file
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
