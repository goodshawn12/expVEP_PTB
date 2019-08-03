function recordData = main(config)
%% Gaze VEP Experiment: Pilot

%% Define parameters
stim.BASELINE = config.BASELINE;

stim.FIXATION_TIME = 1; % sec
stim.FIXATION_LOCATION = 25; % percent of width of the screen from left to right
stim.FIXATION_LENGTH = 30;
stim.FIXATION_WIDTH = 5;
stim.FIXATION_LEFT = 25;     % percent of width of the screen from left to right
stim.FIXATION_RIGHT = 75;    % percent of width of the screen from left to right
BG_COLOR = [stim.BASELINE, stim.BASELINE, stim.BASELINE];   % background color

% parameters for generating codes
CODE_LENGTH = floor(config.STIM_LEN * config.REFRESH);
SEED_1 = 10;
SEED_2 = 10^2;
SEED_3 = 10^3;
SEED_4 = 10^4;
NUM_REPEAT = config.NUM_REPEAT_DEMO;    % for showing stimuli only

% parameters for drawing text
BLACK = [0,0,0];
text.TEXT_FONT = 'Arial';
text.FONT_SIZE = 32;

%% ------------------------------------------------------------------------
%              Generate code sequence
%--------------------------------------------------------------------------
% [TODO] (optional) implement a circular shift for FMC, MSEQ codes

% code length for fmc: (2^nBits - 1) * 2
nBits_fmc = ceil(log2(config.STIM_LEN * config.REFRESH / 2 + 1));   % ensure the code does not repeat for given stimulation length
code1_fmc = gen_code_fmc(nBits_fmc,CODE_LENGTH,SEED_1);
code2_fmc = gen_code_fmc(nBits_fmc,CODE_LENGTH,SEED_2);
code3_fmc = gen_code_fmc(nBits_fmc,CODE_LENGTH,SEED_3);
code4_fmc = gen_code_fmc(nBits_fmc,CODE_LENGTH,SEED_4);

% code length for m-sequence: (2^nBits - 1)
nBits_mseq = ceil(log2(config.STIM_LEN * config.REFRESH + 1));   % ensure the code does not repeat for given stimulation length
code1_mseq = gen_msequence(nBits_mseq,CODE_LENGTH,SEED_1);
code2_mseq = gen_msequence(nBits_mseq,CODE_LENGTH,SEED_2);
code3_mseq = gen_msequence(nBits_mseq,CODE_LENGTH,SEED_3);
code4_mseq = gen_msequence(nBits_mseq,CODE_LENGTH,SEED_4);

% default 30Hz with phase difference
code1_ssvep = gen_code_ssvep(9, CODE_LENGTH, config.REFRESH,0);  
code2_ssvep = gen_code_ssvep(10, CODE_LENGTH, config.REFRESH,1);    % delay 1 frame
code3_ssvep = gen_code_ssvep(11, CODE_LENGTH, config.REFRESH,2);  
code4_ssvep = gen_code_ssvep(12, CODE_LENGTH, config.REFRESH,3);    % delay 1 frame


%% ------------------------------------------------------------------------
%              Initialize communication ports (LSL)
%--------------------------------------------------------------------------
if config.ENABLE_LSL
    % Instantiate LSL
    lslObj = lsl_loadlib();
    
    % LabStreamingLayer (LSL) setting for sending event marker
    info = lsl_streaminfo(lslObj, 'MyMarkerStream', 'Markers', 1, 0, 'cf_string', 'gazevep');
    outlet = lsl_outlet(info);
end

%% ------------------------------------------------------------------------
%              Setup PsychoToolbox and Prepare Stim. Screen
%--------------------------------------------------------------------------
fprintf('exp_gazevep: Creating windows.\n');

% Setup PTB with some default values
PsychDefaultSetup(2);

% For demo purpose - tolerate inaccurate timing (change to 0 for exp)
Screen('Preference', 'SkipSyncTests', config.DEMO);
if ~config.DEMO, HideCursor; end

% prepare keyboard information
KbName('UnifyKeyNames');
escKey = KbName('ESCAPE');
oneKey = KbName('1!');
twoKey = KbName('2@');
threeKey = KbName('3#');
fourKey = KbName('4$');
fiveKey = KbName('5%');
sixKey = KbName('6^');
sevenKey = KbName('7&');
eightKey = KbName('8*');
nineKey = KbName('9(');

% if multiple screens opened, select the last one
screenNumber = max(Screen('Screens'));

% obtain the window size and open a window
if ~config.FULLSCREEN
    WINDOW_WIDTH = config.WIN_WIDTH;
    WINDOW_HEIGHT = config.WIN_HEIGHT;
    [window, windowRect] = Screen('OpenWindow', screenNumber, BG_COLOR, [0, 0, WINDOW_WIDTH, WINDOW_HEIGHT], [], 2);
else
    [WINDOW_WIDTH, WINDOW_HEIGHT] = Screen('WindowSize',screenNumber);
    [window, windowRect] = Screen('OpenWindow', screenNumber, BG_COLOR, [], [], 2);
end

[centX, centY] = RectCenter(windowRect);

% Flip to clear
Screen('Flip', window);

% obtain refresh rate of the monitor
refreshRateHz = 1/Screen('GetFlipInterval', window);
fprintf('ONLINE-BCI: Refresh rate of on-screen window = %d [Hz].\n',refreshRateHz);


%% ------------------------------------------------------------------------
%              Prepare off-screen windows for experiment
%--------------------------------------------------------------------------

prepareMsg = sprintf('Now loading stimulation screens...');
prepareScreen = sys_prepInstructionScreen(window, prepareMsg, BG_COLOR, ...
    BLACK, text.TEXT_FONT, text.FONT_SIZE, centX, centY);
Screen('CopyWindow', prepareScreen, window);
Screen('Flip', window);

startMsg_exp0 = cell(1,3);
startMsg_exp0{1} = 'Calibration - Blinking Task';
startMsg_exp0{2} = 'Blink every 2 sec (hint on screen)';
startMsg_exp0{3} = 'Press Any Key To Begin';
startExp0Screen = sys_prepInstructionScreen(window, startMsg_exp0, BG_COLOR, ...
    BLACK, text.TEXT_FONT, text.FONT_SIZE, centX, centY);

exp0_2 = {'2'};
exp0_Screen_2 = sys_prepInstructionScreen(window, exp0_2, BG_COLOR, ...
    BLACK, text.TEXT_FONT, text.FONT_SIZE, centX, centY);

exp0_1 = {'1'};
exp0_Screen_1 = sys_prepInstructionScreen(window, exp0_1, BG_COLOR, ...
    BLACK, text.TEXT_FONT, text.FONT_SIZE, centX, centY);

exp0_blink = {'Blink'};
exp0_Screen_blink = sys_prepInstructionScreen(window, exp0_blink, BG_COLOR, ...
    BLACK, text.TEXT_FONT, text.FONT_SIZE, centX, centY);



startMsg_exp1 = cell(1,2);
startMsg_exp1{1} = 'Experiment - Task 1';
startMsg_exp1{2} = 'Press Any Key To Begin';
startExp1Screen = sys_prepInstructionScreen(window, startMsg_exp1, BG_COLOR, ...
    BLACK, text.TEXT_FONT, text.FONT_SIZE, centX, centY);

endMsg_exp1 = sprintf('End of Task 1. Press any key to continue...');
endExp1Screen = sys_prepInstructionScreen(window, endMsg_exp1, BG_COLOR, ...
    BLACK, text.TEXT_FONT, text.FONT_SIZE, centX, centY);


% Questionanire screen
rateMsg_1 = cell(1,9);
rateMsg_1{1} = 'Rate your feeling toward the stimuli';
rateMsg_1{2} = ' ';
rateMsg_1{3} = '1: Not Perceptible                           ';      % 0
rateMsg_1{4} = '2: Barely Perceptible / Not Annoying';        % 1
rateMsg_1{5} = '3: Perceptible / Slightly Annoying   ';   % 3
rateMsg_1{6} = '4: Quite Strong Stimuli / Annoying   ';           % 10
rateMsg_1{7} = '5: Strong Stimuli / Very Annoying     ';        % 30
rateMsg_1{8} = '6: Very Strong Stimuli / Uncomfortable';        % 30
rateMsg_1{9} = '7: Extremly Strong Stimuli / Very Uncomfortable';        % 30
rateScreen_1 = sys_prepInstructionScreen_align(window, rateMsg_1, BG_COLOR, ...
    BLACK, text.TEXT_FONT, text.FONT_SIZE, centX, centY,1);

demoMsg = cell(1,2);
demoMsg{1} = 'Demo Experimental Stimuli';
demoMsg{2} = 'Press Any Key To Begin';
demoScreen = sys_prepInstructionScreen(window, demoMsg, BG_COLOR, ...
    BLACK, text.TEXT_FONT, text.FONT_SIZE, centX, centY);


%% ------------------------------------------------------------------------
%              Generate Stimulation Textures (stimulus screen)
%--------------------------------------------------------------------------

% generate image frames
[frame_image, frame_contrast, fixation_mask, frame_base] = gen_frame_image( ...
    WINDOW_HEIGHT,WINDOW_WIDTH,config.NUM_IMG,stim.BASELINE,stim.FIXATION_LENGTH,stim.FIXATION_WIDTH);

% generate final textures for gray and image frames
image_textures = gen_image_textures(window, frame_image, frame_contrast, fixation_mask, config.CONTRAST_LIST, config.NUM_IMG);
[gray_textures, baseFixTexture] = gen_gray_textures(window, frame_base, frame_contrast, fixation_mask, config.CONTRAST_LIST);

% % prepare base texture with / without fixation cross at different locations
% [locTexture, noCrossTextures, baseFixTexture, baseTexture] = gen_base_textures(window, stim, config, WINDOW_HEIGHT, WINDOW_WIDTH);
% 
% % prepare grey binary background with different contrast levels and fixation locations
% [contrastTextures, chromeTextures] = gen_contrast_textures(window, stim, config, WINDOW_HEIGHT, WINDOW_WIDTH);

% image and text
% imageTexture = gen_image_textures(image_screen, window, stim, config, WINDOW_HEIGHT, WINDOW_WIDTH);



%% ----------------------------------------------------------------------
%             Start the Main Loop for Experiment
%----------------------------------------------------------------------

Priority(1);
recordData = cell(1,3);

% store stimuli info to output
out.mode = config.MODE;
out.code_fmc = {code1_fmc, code2_fmc, code3_fmc, code4_fmc};
out.code_mseq = {code1_mseq, code2_mseq, code3_mseq, code4_mseq};
out.code_ssvep = {code1_ssvep, code2_ssvep, code3_ssvep, code4_ssvep};
recordData{1} = out;
response = [];

%% ----------------------------------------------------------------------
%             Calibration - blinking task for event synchronization
%----------------------------------------------------------------------

if config.RUN_CALIB
    disp('Enter Calibration: Blinking task');
    
    % Event marker
    if config.ENABLE_LSL, outlet.push_sample({'Exp0'}); end
    
    % present first screen
    Screen('CopyWindow', startExp0Screen, window);
    Screen('Flip', window);
    KbStrokeWait;
    
    for frame_id = 1:30        
        % prep - 2
        Screen('CopyWindow', exp0_Screen_2, window);
        Screen('Flip', window);
        WaitSecs(0.6667);
        
        % prep - 1
        Screen('CopyWindow', exp0_Screen_1, window);
        Screen('Flip', window);
        WaitSecs(0.6667);
        
        % Event marker
        if config.ENABLE_LSL, outlet.push_sample({'Blink'}); end
        
        % blink
        Screen('CopyWindow', exp0_Screen_blink, window);
        Screen('Flip', window);
        WaitSecs(0.6667);
        
        % break the loop if press ESC key
        [~, ~, keyCode] = KbCheck;
        if keyCode(escKey)
            Screen('CloseAll');
            return;
        end
    end
end



%% ----------------------------------------------------------------------
%             Experimental loop - Task 1: contrast levels
%----------------------------------------------------------------------

if config.RUN_EXP_1
    
    disp('Enter Experiment: visual perception task');
    
    % Event marker
    if config.ENABLE_LSL, outlet.push_sample({'Exp1'}); end
    
    nModes = length(config.MODE);

    % obtain trial order, randomized across types of stimuli, contrast levels, and left vs. right fixation
    nTrial = nModes * length(config.CONTRAST_LIST) * 4 * config.NUM_TRIAL_CONTRAST;    % 4: locations
    trial_order = randperm(nTrial);
    
    % present first screen
    Screen('CopyWindow', startExp1Screen, window);
    Screen('Flip', window);
    KbStrokeWait;
    
    respMat = [];
    respMat.stimuli = zeros(1,nTrial);
    respMat.contrast = zeros(1,nTrial);
    respMat.location = zeros(1,nTrial);
    respMat.image = zeros(1,nTrial);
    respMat.rate_percept = zeros(1,nTrial);
    respMat.rt = zeros(1,nTrial);
    
    for trial = 1:nTrial
        
        image_id = mod(trial_order(trial)-1, config.NUM_TRIAL_CONTRAST) + 1;
        loc = mod(ceil(trial_order(trial)/config.NUM_TRIAL_CONTRAST)-1, 4) + 1;
        cont_level = mod(ceil(trial_order(trial)/config.NUM_TRIAL_CONTRAST/4)-1, length(config.CONTRAST_LIST)) + 1;
        mode = ceil(trial_order(trial)/config.NUM_TRIAL_CONTRAST/4/length(config.CONTRAST_LIST));
                
        % Event marker
        if config.ENABLE_LSL, outlet.push_sample({sprintf('Cross_%d',loc)}); end
        
        % present cross for eye fixation
        Screen('DrawTexture', window, baseFixTexture{loc}); % [TODO]
        Screen('Flip',window);
        
        % wait time
        WaitSecs(stim.FIXATION_TIME);
        
        % Event marker
        if config.ENABLE_LSL
            outlet.push_sample({sprintf('Stim%d_CT%d_Loc%d_Img%d', ...
                mode, config.CONTRAST_LIST(cont_level), loc, image_id)});
        end
        
        % present stimuli according to the code sequence and stimulation mode
        for it_code = 1:CODE_LENGTH
            if mode == 1     % fmc
                Screen('DrawTexture',window,gray_textures{cont_level,loc, ...
                    code1_fmc(it_code)+code2_fmc(it_code)*2+code3_fmc(it_code)*4+code4_fmc(it_code)*8+1 });
            elseif mode == 2     % msequence
                Screen('DrawTexture',window,gray_textures{cont_level,loc, ...
                    code1_mseq(it_code)+code2_mseq(it_code)*2+code3_mseq(it_code)*4+code4_mseq(it_code)*8+1 });
            elseif mode == 3     % ssvep
                Screen('DrawTexture',window,gray_textures{cont_level,loc, ...
                    code1_ssvep(it_code)+code2_ssvep(it_code)*2+code3_ssvep(it_code)*4+code4_ssvep(it_code)*8+1 });
            elseif mode == 4    % fmc_image
                Screen('DrawTexture',window,image_textures{cont_level,image_id, ...
                    code1_fmc(it_code)+code2_fmc(it_code)*2+code3_fmc(it_code)*4+code4_fmc(it_code)*8+1 });
            end
            
            % Flip and present the stimuli
            Screen('Flip', window);
            
            % break the loop if press ESC key
            [~, ~, keyCode] = KbCheck;
            if keyCode(escKey)
                Screen('CloseAll');
                return;
            end
        end
        
        % Event marker
        if config.ENABLE_LSL, outlet.push_sample({'Rate'}); end
        
        % present questionnaire
        Screen('CopyWindow', rateScreen_1, window);
        Screen('Flip', window);
        
        % Check the keyboard
        restTic = tic;
        response = 0;
        while true
            [~, ~, keyCode] = KbCheck;
            if keyCode(escKey)
                ShowCursor;
                Screen('CloseAll');
                return;
            end
            if keyCode(oneKey), response = 1; break; end
            if keyCode(twoKey), response = 2; break; end
            if keyCode(threeKey), response = 3; break; end
            if keyCode(fourKey), response = 4; break; end
            if keyCode(fiveKey), response = 5; break; end
            if keyCode(sixKey), response = 6; break; end
            if keyCode(sevenKey), response = 7; break; end
        end
                
        % Record the trial data
        respMat.stimuli(trial) = mode;
        respMat.contrast(trial) = config.CONTRAST_LIST(cont_level);
        respMat.location(trial) = loc;
        respMat.image(trial) = image_id;
        respMat.rate_percept(trial) = response;
        respMat.rt(trial) = toc(restTic);
        
        fprintf('Trial: %d, Mode: %d, Contrast: %d, Location: %d, Image: %d, Percept: %d, RT: %f\n', ...
            trial, mode, config.CONTRAST_LIST(cont_level), loc, image_id, response, respMat.rt(trial));
        
        if mod(trial,config.BREAK_INTERVAL) == 0 && trial ~= nTrial
            recordData{2} = respMat; % store behavioral data for each block (in case interrupted)
            
            if mod(floor(trial/config.BREAK_INTERVAL),config.SESS_INTERNAL) == 0
                % present rest texture - long break
                restMsg_long = cell(1,4);
                restMsg_long{1} = sprintf('Session %d / %d Completed!', floor(trial/config.BREAK_INTERVAL/config.SESS_INTERNAL), floor(nTrial/config.BREAK_INTERVAL/config.SESS_INTERNAL));
                restMsg_long{2} = sprintf('Block %d / %d Completed!', floor(trial/config.BREAK_INTERVAL), floor(nTrial/config.BREAK_INTERVAL));
                restMsg_long{3} = 'You deserve a longer break.';
                restMsg_long{4} = 'Press any key to continue the experiment...';
                restScreen = sys_prepInstructionScreen(window, restMsg_long, BG_COLOR, ...
                    BLACK, text.TEXT_FONT, text.FONT_SIZE, centX, centY);
                
                Screen('CopyWindow', restScreen, window);
                Screen('Flip', window);
                KbStrokeWait;
            else
                % present rest texture - short break
                restMsg_exp1 = cell(1,3);
                restMsg_exp1{1} = sprintf('Block %d / %d Completed!', floor(trial/config.BREAK_INTERVAL), floor(nTrial/config.BREAK_INTERVAL));
                restMsg_exp1{2} = 'Take a rest!';
                restMsg_exp1{3} ='Press any key to continue the experiment...';
                restScreen = sys_prepInstructionScreen(window, restMsg_exp1, BG_COLOR, ...
                    BLACK, text.TEXT_FONT, text.FONT_SIZE, centX, centY);
                
                Screen('CopyWindow', restScreen, window);
                Screen('Flip', window);
                KbStrokeWait;
            end
        end
    end
    recordData{2} = respMat; % store behavioral data - final
    
    % end of experiment screen
    Screen('CopyWindow', endExp1Screen, window);
    Screen('Flip', window);
    KbStrokeWait;
    
end


%% ----------------------------------------------------------------------
%             Only Present Stimuli
%----------------------------------------------------------------------
if config.RUN_DEMO % only show stimuli
    
    disp('Present stimuli only');
    mode = config.MODE_DEMO;

    % present first screen
    Screen('CopyWindow', demoScreen, window);
    Screen('Flip', window);
    KbStrokeWait;
    
    for it_rep = 1:NUM_REPEAT
        for it_code = 1:CODE_LENGTH
            if mode == 1     % fmc
                Screen('DrawTexture',window,gray_textures{config.DEMO_CONTRAST_LEVEL,config.DEMO_LOCATION, ...
                    code1_fmc(it_code)+code2_fmc(it_code)*2+code3_fmc(it_code)*4+code4_fmc(it_code)*8+1 });
            elseif mode == 2     % msequence
                Screen('DrawTexture',window,gray_textures{config.DEMO_CONTRAST_LEVEL,config.DEMO_LOCATION, ...
                    code1_mseq(it_code)+code2_mseq(it_code)*2+code3_mseq(it_code)*4+code4_mseq(it_code)*8+1 });
            elseif mode == 3     % ssvep
                Screen('DrawTexture',window,gray_textures{config.DEMO_CONTRAST_LEVEL,config.DEMO_LOCATION, ...
                    code1_ssvep(it_code)+code2_ssvep(it_code)*2+code3_ssvep(it_code)*4+code4_ssvep(it_code)*8+1 });
            elseif mode == 4    % fmc_image
                Screen('DrawTexture',window,image_textures{config.DEMO_CONTRAST_LEVEL,config.DEMO_IMAGE, ...
                    code1_fmc(it_code)+code2_fmc(it_code)*2+code3_fmc(it_code)*4+code4_fmc(it_code)*8+1 });
            end
            
            % Flip and present the stimuli
            Screen('Flip', window);
            
            % break the loop if press ESC key
            [~, ~, keyCode] = KbCheck;
            if keyCode(escKey)
                Screen('CloseAll');
                return;
            end
        end
    end
end


%% closing
fprintf('ONLINE-BCI: Online BCI has done.\n');
Priority(0);
ShowCursor;
Screen('CloseAll')
sca

end


function code = gen_code_ssvep(freq, code_length, refresh, delay)

EPSILON = 1e-12;

index = 0:1:code_length-1;
index = index + delay;  % delay in frame
tmp = sin(2*pi*freq*(index/refresh)+EPSILON);
code = (tmp>=0);

end

function code = gen_code_fmc(num_bits,code_length,seed)

% generate m-sequence
len_mseq = 2^num_bits-1;
mseq = gen_msequence(num_bits,len_mseq,seed);

% initialize code length 2 times the length of m-sequence
fmc = boolean(zeros(1,2*len_mseq));

% stretch the time sequence and perform frequency modulation (flip every other bits)
fmc(2*(1:len_mseq)-1) = mseq;      % odd indices
fmc(2*(1:len_mseq)) = not(mseq);   % even indices (flipped)

% truncate the code according to the given code length
code = fmc(1:code_length);

end

function code = gen_msequence(num_bits,code_length,seed)

% generate maximum length sequence (m-sequence) - pseudorandom binary sequence
rng(seed)
state = rand(1,num_bits) > 0.5;

% length-m registers produces 2^m - 1
mseq = boolean(zeros(1, 2^num_bits - 1));

% linear feedback shift register operation for generating m-sequence
% (can be improved using bit shift)
for it = 1:length(mseq)
    old_state = state;
    mseq(it) = state(end);
    state(2:end) = old_state(1:end-1);
    if ismember(num_bits, [2,3,4,6,7,15,22])
        state(1) = xor(old_state(num_bits), old_state(num_bits-1));
    elseif ismember(num_bits, [5,11,21])
        state(1) = xor(old_state(num_bits), old_state(num_bits-2));
    elseif ismember(num_bits, [10,17])
        state(1) = xor(old_state(num_bits), old_state(num_bits-3));
    elseif ismember(num_bits, [9])
        state(1) = xor(old_state(num_bits), old_state(num_bits-4));
    elseif ismember(num_bits, [8])
        state(1) = xor(xor(xor(old_state(8), old_state(6)),old_state(5)),old_state(4));
    end
end

% truncate the code according to the given code length
code = mseq(1:code_length);

end


function [image_textures] = gen_image_textures(window, frame_image, frame_contrast, fixation_mask, CONTRAST_LIST, NUM_IMG)

image_textures = cell(length(CONTRAST_LIST),NUM_IMG,16);

[USERVIEW, ~] = memory;

if ((USERVIEW.MemAvailableAllArrays - USERVIEW.MemUsedMATLAB) * 0.5) < ...
        size(frame_image{1},1) * size(frame_image{1},2) * 3 * 16 * NUM_IMG * length(CONTRAST_LIST)
    error('Insufficient Memory for generating stimulation frames!')
else
    for cont_it = 1:length(CONTRAST_LIST)
        contrast = CONTRAST_LIST(cont_it);
        for img_it = 1:NUM_IMG
            for code_it = 1:16
                image_textures{cont_it,img_it,code_it} = Screen('MakeTexture',window, ...
                    uint8( frame_image{img_it} + contrast * frame_contrast{code_it}) );
            end
        end
    end
end

end


function [gray_textures, baseFixTexture] = gen_gray_textures(window, frame_base, frame_contrast, fixation_mask, CONTRAST_LIST)

baseFixTexture = cell(1,4);
gray_textures = cell(length(CONTRAST_LIST),4,16);

for cont_it = 1:length(CONTRAST_LIST)
    contrast = CONTRAST_LIST(cont_it);
    for loc_it = 1:4
        for code_it = 1:16
            gray_textures{cont_it,loc_it,code_it} = Screen('MakeTexture',window, ...
                uint8(fixation_mask{loc_it} .* (frame_base + contrast * frame_contrast{code_it})) );
        end
    end
end

for loc_it = 1:4
    baseFixTexture{loc_it} = Screen('MakeTexture',window, uint8(fixation_mask{loc_it} .* frame_base) );
end

end