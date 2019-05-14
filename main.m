function recordData = main(config)
%% Gaze VEP Experiment: Pilot

%% Define parameters
stim.CONTRAST = config.CONTRAST;  % BASELINE-CONTRAST and BASELINE+CONTRAST
stim.BASELINE = config.BASELINE;

stim.FIXATION_TIME = 0.8; % sec
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
NUM_REPEAT = config.NUM_REPEAT_DEMO;    % for showing stimuli only

% parameters for drawing text
BLACK = [0,0,0];
text.TEXT_FONT = 'Arial';
text.FONT_SIZE = 36;

% Questionnaire:
config.QUESTION_2 = 0;  % display second question

%% read images
image_screen = imread(config.filename_image);

% image_screen = [];
% if config.RUN_DEMO
%     if config.MODE_DEMO(1) == 3   % fmc_single_image
%         image_screen = imread(config.filename_image);
%     elseif config.MODE_DEMO(1) == 4   % fmc_single_text
%         image_screen = imread(config.filename_text);
%     end
% end

%% ------------------------------------------------------------------------
%              Generate code sequence
%--------------------------------------------------------------------------

% code length for fmc: (2^nBits - 1) * 2
nBits_fmc = ceil(log2(config.STIM_LEN * config.REFRESH / 2 + 1));   % ensure the code does not repeat for given stimulation length
code1_fmc = gen_code_fmc(nBits_fmc,CODE_LENGTH,SEED_1);
code2_fmc = gen_code_fmc(nBits_fmc,CODE_LENGTH,SEED_2);

% code length for fmc_slow: (2^nBits - 1) * 4
nBits_slow = ceil(log2(config.STIM_LEN * config.REFRESH / 4 + 1));   % ensure the code does not repeat for given stimulation length
code1_slow = gen_code_fmc(nBits_slow,floor(CODE_LENGTH/2),SEED_1);
code2_slow = gen_code_fmc(nBits_slow,floor(CODE_LENGTH/2),SEED_2);

% code length for m-sequence: (2^nBits - 1)
nBits_mseq = ceil(log2(config.STIM_LEN * config.REFRESH + 1));   % ensure the code does not repeat for given stimulation length
code1_mseq = gen_msequence(nBits_mseq,CODE_LENGTH,SEED_1);
code2_mseq = gen_msequence(nBits_mseq,CODE_LENGTH,SEED_2);

% default 30Hz with phase difference
code1_ssvep = gen_code_ssvep(30, CODE_LENGTH, config.REFRESH,0);  
code2_ssvep = gen_code_ssvep(30, CODE_LENGTH, config.REFRESH,1);    % delay 1 frame

% make sure two codes are not identical
if ~any(code1_fmc~=code2_fmc) || ~any(code1_slow~=code2_slow) || ~any(code1_mseq~=code2_mseq)
    error('Two codes are identical. Select different random seeds.');
end


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

% if ~isempty(image_screen)
%     [WINDOW_HEIGHT, WINDOW_WIDTH, ~] = size(image_screen);
%     [window, windowRect] = Screen('OpenWindow', screenNumber, BG_COLOR, [0, 0, WINDOW_WIDTH, WINDOW_HEIGHT], [], 2);
% end

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

startMsg_exp1 = cell(1,2);
startMsg_exp1{1} = 'Experiment - Task 1';
startMsg_exp1{2} = 'Press Any Key To Begin';
startExp1Screen = sys_prepInstructionScreen(window, startMsg_exp1, BG_COLOR, ...
    BLACK, text.TEXT_FONT, text.FONT_SIZE, centX, centY);

endMsg_exp1 = sprintf('End of Task 1. Press any key to continue...');
endExp1Screen = sys_prepInstructionScreen(window, endMsg_exp1, BG_COLOR, ...
    BLACK, text.TEXT_FONT, text.FONT_SIZE, centX, centY);

startMsg_exp2 = cell(1,2);
startMsg_exp2{1} = 'Experiment - Task 2';
startMsg_exp2{2} = 'Press Any Key To Begin';
startExp2Screen = sys_prepInstructionScreen(window, startMsg_exp2, BG_COLOR, ...
    BLACK, text.TEXT_FONT, text.FONT_SIZE, centX, centY);

endMsg_exp2 = sprintf('End of Task 2. Press any key to continue...');
endExp2Screen = sys_prepInstructionScreen(window, endMsg_exp2, BG_COLOR, ...
    BLACK, text.TEXT_FONT, text.FONT_SIZE, centX, centY);


% Questionanire screen
rateMsg = cell(1,3);
rateMsg{1} = 'Rate the intensity of the flicker that you perceived';
rateMsg{2} = '   1          2          3          4          5    ';
rateMsg{3} = '(no flicker)                   (very strong flicker)';
rateScreen = sys_prepInstructionScreen(window, rateMsg, BG_COLOR, ...
    BLACK, text.TEXT_FONT, text.FONT_SIZE, centX, centY);

rateMsg_1 = cell(1,7);
rateMsg_1{1} = 'Rate your feeling toward the stimuli';
rateMsg_1{2} = ' ';
rateMsg_1{3} = '1: Not Perceptible                           ';      % 0
rateMsg_1{4} = '2: Barely Perceptible / Not Annoying';        % 1
rateMsg_1{5} = '3: Perceptible / Slightly Annoying   ';   % 3
rateMsg_1{6} = '4: Quite Strong Stimuli / Annoying   ';           % 10
rateMsg_1{7} = '5: Strong Stimuli / Very Annoying     ';        % 30
rateScreen_1 = sys_prepInstructionScreen_align(window, rateMsg_1, BG_COLOR, ...
    BLACK, text.TEXT_FONT, text.FONT_SIZE, centX, centY,1);


rateMsg_2 = cell(1,7);
rateMsg_2{1} = 'How long you can look at the stimulus comfortably?';
rateMsg_2{2} = '(Select one option below)';
rateMsg_2{3} = '(1): 0 sec (not tolerable)          ';      % 0
rateMsg_2{4} = '(2): 3 sec (barely tolerable)     ';        % 1
rateMsg_2{5} = '(3): 10 sec (tolerable)                ';   % 3
rateMsg_2{6} = '(4): 30 sec (quite comfortable)';           % 10
rateMsg_2{7} = '(5): 90 sec (comfortable)         ';        % 30
rateScreen_2 = sys_prepInstructionScreen(window, rateMsg_2, BG_COLOR, ...
    BLACK, text.TEXT_FONT, text.FONT_SIZE, centX, centY);

demoMsg = cell(1,2);
demoMsg{1} = 'Demo Experimental Stimuli';
demoMsg{2} = 'Press Any Key To Begin';
demoScreen = sys_prepInstructionScreen(window, demoMsg, BG_COLOR, ...
    BLACK, text.TEXT_FONT, text.FONT_SIZE, centX, centY);


%% ------------------------------------------------------------------------
%              Generate Stimulation Textures (stimulus screen)
%--------------------------------------------------------------------------

% prepare grey binary background with different contrast levels and fixation locations
[contrastTextures, chromeTextures] = gen_contrast_textures(window, stim, config, WINDOW_HEIGHT, WINDOW_WIDTH);

% prepare base texture with / without fixation cross at different locations
[locTexture, noCrossTextures, baseFixTexture, baseTexture] = gen_base_textures(window, stim, config, WINDOW_HEIGHT, WINDOW_WIDTH);

% image and text
imageTexture = gen_image_textures(image_screen, window, stim, config, WINDOW_HEIGHT, WINDOW_WIDTH);

% if ~isempty(image_screen)
%     vepImage{1} = Screen('MakeTexture',window,image_screen-stim.CONTRAST);
%     vepImage{2} = Screen('MakeTexture',window,image_screen+stim.CONTRAST);
%     
%     vepText{1} = Screen('MakeTexture',window,image_screen);
%     vepText{2} = Screen('MakeTexture',window,image_screen+2*stim.CONTRAST);
% end





%% ----------------------------------------------------------------------
%             Start the Main Loop for Experiment
%----------------------------------------------------------------------

Priority(1);
recordData = cell(1,3);

% store stimuli info to output
out.mode = config.MODE;
out.code_fmc = {code1_fmc, code2_fmc};
out.code_slow = {code1_slow, code2_slow};
out.code_mseq = {code1_mseq, code2_mseq};
out.code_ssvep = {code1_ssvep, code2_ssvep};
recordData{1} = out;
response = [];

%% ----------------------------------------------------------------------
%             Experimental loop - Task 1: contrast levels
%----------------------------------------------------------------------

if config.RUN_EXP_1
    
    disp('Enter Experiment: Task 1 - perception test');
    
    % Event marker
    if config.ENABLE_LSL, outlet.push_sample({'Exp1'}); end
    
    nModes = length(config.MODE);
    conditionIndex = [];
    for mode_it = 1:nModes
        for cond_it = 1:length(config.LIST_EACH_MODE{mode_it})
            conditionIndex{end+1} = [mode_it, cond_it];
        end
    end

    % obtain trial order, randomized across types of stimuli, contrast levels, and left vs. right fixation
    nTrial = length(conditionIndex) * config.NUM_TRIAL_CONTRAST * 2;    % 2: left / right locations
    trial_order = repmat( 1:(length(conditionIndex)*2), 1, config.NUM_TRIAL_CONTRAST);
    trial_order = trial_order( randperm(nTrial) );
    loc_order = mod(trial_order-1, 2)+1;
    cond_order = conditionIndex(ceil(trial_order/2));
    
    % present first screen
    Screen('CopyWindow', startExp1Screen, window);
    Screen('Flip', window);
    KbStrokeWait;
    
    respMat = [];
    respMat.stimuli = zeros(1,nTrial);
    respMat.contrast = zeros(1,nTrial);
    respMat.location = zeros(1,nTrial);
    respMat.rate_percept = zeros(1,nTrial);
    respMat.rate_comfort = zeros(1,nTrial);
    respMat.rt = zeros(1,nTrial);
    
    for trial = 1:nTrial
        
        mode = config.MODE{cond_order{trial}(1)}; 
        cont_level = cond_order{trial}(2); 
        loc = loc_order(trial);
        
        % Event marker
        if config.ENABLE_LSL, outlet.push_sample({sprintf('Cross_%d',loc)}); end
        
        % present cross for eye fixation
        Screen('DrawTexture', window, baseFixTexture{loc});
        Screen('Flip',window);
        
        % wait time
        WaitSecs(stim.FIXATION_TIME);
        
        % Event marker
        if config.ENABLE_LSL
            outlet.push_sample({sprintf('Stim%d%d_CT%d_Loc%d', ...
                mode, config.CONTRAST_LIST(cont_level), loc)});
        end
        
        % present stimuli according to the code sequence and stimulation mode
        if mode(2) == 3 % fmc with slow (15 Hz) modulation
            for it_code = 1:length(code1_slow)
                Screen('DrawTexture',window,contrastTextures{ cont_level, code1_slow(it_code)*2+code2_slow(it_code)+1, loc+1 });
                Screen('Flip', window);
                
                Screen('DrawTexture',window,baseFixTexture{loc});
                Screen('Flip', window);
                
                % break the loop if press ESC key
                [~, ~, keyCode] = KbCheck;
                if keyCode(escKey)
                    Screen('CloseAll');
                    return;
                end
            end
        elseif mode(2) == 0 % normal (30Hz modulation), grey, independent codes
            for it_code = 1:CODE_LENGTH
                if mode(1) == 0     % frequency modulated code
                    Screen('DrawTexture',window,contrastTextures{ cont_level, code1_fmc(it_code)*2+code2_fmc(it_code)+1, loc+1 });
                elseif mode(1) == 1     % msequence
                    Screen('DrawTexture',window,contrastTextures{ cont_level, code1_mseq(it_code)*2+code2_mseq(it_code)+1, loc+1 });
                elseif mode(1) == 2     % ssvep
                    Screen('DrawTexture',window,contrastTextures{ cont_level, code1_ssvep(it_code)+2, loc+1 });
%                     Screen('DrawTexture',window,contrastTextures{ cont_level, code1_ssvep(it_code)*2+code2_ssvep(it_code)+1, loc+1 });
                elseif mode(1) == 3     % image
                    Screen('DrawTexture',window,imageTexture{ cont_level, code1_fmc(it_code)*2+code2_fmc(it_code)+1, loc+1 });
                elseif mode(1) == 5     % control (static image)
                    Screen('DrawTexture',window,contrastTextures{ 1, 2, loc+1 });  % always present [0 1] screen with smallest contrast
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
        else
            for it_code = 1:CODE_LENGTH
                if mode(2) == 1 % opposite codes
                    Screen('DrawTexture',window,contrastTextures{ cont_level, code1_fmc(it_code)+2, loc+1 });
                elseif mode(2) == 2 % chromatic modulation
                    Screen('DrawTexture',window,chromeTextures{ cont_level, code1_fmc(it_code)*2+code2_fmc(it_code)+1, loc+1 });
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
        
        % Event marker
        if config.ENABLE_LSL, outlet.push_sample({'Rate'}); end
        
        % present questionnaire
        Screen('CopyWindow', rateScreen_1, window);
        Screen('Flip', window);
        
        % Check the keyboard
        restTic = tic;
        response = 0;
        comfortability = 0;
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
        end
        
        if config.QUESTION_2
            % Event marker
            if config.ENABLE_LSL, outlet.push_sample({'Rate2'}); end
            
            % present questionnaire
            Screen('CopyWindow', rateScreen_2, window);
            Screen('Flip', window);
            
            % listen to keyboard after 0.3 second to avoid accidental button press
            WaitSecs(0.3);
            
            % Check the keyboard
            while true
                [~, ~, keyCode] = KbCheck;
                if keyCode(escKey)
                    ShowCursor;
                    Screen('CloseAll');
                    return;
                end
                if keyCode(oneKey), comfortability = 1; break; end
                if keyCode(twoKey), comfortability = 2; break; end
                if keyCode(threeKey), comfortability = 3; break; end
                if keyCode(fourKey), comfortability = 4; break; end
                if keyCode(fiveKey), comfortability = 5; break; end
            end
        end
        
        % Record the trial data
        respMat.stimuli(trial) = cond_order{trial}(1);
        respMat.contrast(trial) = config.CONTRAST_LIST(cont_level);
        respMat.location(trial) = loc;
        respMat.rate_percept(trial) = response;
        respMat.rate_comfort(trial) = comfortability;
        respMat.rt(trial) = toc(restTic);
        
        fprintf('Trial: %2d, Mode: %d%d, Contrast: %2d, Location: %d, Percept: %d, Comfort: %d, RT: %f\n', ...
            trial, mode, config.CONTRAST_LIST(cont_level), loc, response, comfortability, respMat.rt(trial));
        
        if mod(trial,config.BREAK_INTERVAL) == 0 && trial ~= nTrial
            recordData{2} = respMat; % store behavioral data for each block (in case interrupted)
            
            % presetn rest texture
            restMsg_exp1 = cell(1,2);
            restMsg_exp1{1} = sprintf('Block %d / %d Completed!', floor(trial/config.BREAK_INTERVAL), floor(nTrial/config.BREAK_INTERVAL));
            restMsg_exp1{2} = 'Take a rest! Press any key to continue the experiment...';
            restScreen = sys_prepInstructionScreen(window, restMsg_exp1, BG_COLOR, ...
                BLACK, text.TEXT_FONT, text.FONT_SIZE, centX, centY);
            
            Screen('CopyWindow', restScreen, window);
            Screen('Flip', window);
            KbStrokeWait;
        end
        
    end
    recordData{2} = respMat; % store behavioral data - final
    
    % end of experiment screen
    Screen('CopyWindow', endExp1Screen, window);
    Screen('Flip', window);
    KbStrokeWait;
    
end


%% ----------------------------------------------------------------------
%             Experimental loop - Task 2: fixation locations
%----------------------------------------------------------------------

if config.RUN_EXP_2
    
    disp('Enter Experiment: Task 2 - fixation location');
    
    % Event marker
    if config.ENABLE_LSL, outlet.push_sample({'Exp2'}); end
    
    % obtain randomized trial order for different conditions
    trial_order = repmat( [1:length(config.LOCATION_LIST)], 1, config.NUM_TRIAL_LOCATION);
    trial_order = trial_order( randperm(length(trial_order)) );
    
    % present first screen
    Screen('CopyWindow', startExp2Screen, window);
    Screen('Flip', window);
    KbStrokeWait;
    
    respMat = [];
    respMat.contrast = zeros(1,length(trial_order));
    respMat.location = zeros(1,length(trial_order));
    respMat.rate_percept = zeros(1,length(trial_order));
    respMat.rt = zeros(1,length(trial_order));
    
    for trial = 1:length(trial_order)
        
        loc = trial_order(trial);
        
        % Event marker
        if config.ENABLE_LSL, outlet.push_sample({sprintf('Cross_%d',loc)}); end
        
        % present cross for eye fixation
        Screen('DrawTexture', window, baseFixTexture{2+loc});
        Screen('Flip',window);
        
        % wait one second
        WaitSecs(stim.FIXATION_TIME);
        
        % Event marker
        if config.ENABLE_LSL
            % default to fmc normal independent code modes (to add other modes later)
            outlet.push_sample({sprintf('Stim%d%d_CT%d_Loc%d', config.MODE_LOC, stim.CONTRAST, config.LOCATION_LIST(loc))});
        end
        
        % present stimuli according to the code sequence - normal (30Hz modulation), grey, independent codes
        for it_code = 1:CODE_LENGTH
            if config.MODE_LOC(1) == 0     % frequency modulated code
                Screen('DrawTexture',window,locTexture{ loc, code1_fmc(it_code)*2+code2_fmc(it_code)+1 });
            elseif config.MODE_LOC(1) == 1     % msequence
                Screen('DrawTexture',window,locTexture{ loc, code1_mseq(it_code)*2+code2_mseq(it_code)+1 });
            elseif config.MODE_LOC(1) == 2     % ssvep
                Screen('DrawTexture',window,locTexture{ loc, code1_ssvep(it_code)*2+code2_ssvep(it_code)+1 });
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
        
        % Record the trial data
        respMat.contrast(trial) = stim.CONTRAST;
        respMat.location(trial) = config.LOCATION_LIST(loc);
        
        fprintf('Trial: %2d, Mode: %d%d, Contrast: %2d, Location: %d\n', ...
            trial, config.MODE_LOC, stim.CONTRAST, config.LOCATION_LIST(loc));
        
        % present break slide and wait one second
        Screen('DrawTexture',window,baseTexture);
        Screen('Flip', window);
        WaitSecs(0.5);
        
        if mod(trial,config.BREAK_INTERVAL_LOC) == 0 && trial ~= length(trial_order)
            recordData{3} = respMat; % store behavioral data for each block (in case interrupted)
            
            % presetn rest texture
            restMsg_exp1 = cell(1,2);
            restMsg_exp1{1} = sprintf('Block %d / %d Completed!', floor(trial/config.BREAK_INTERVAL_LOC), floor(length(trial_order)/config.BREAK_INTERVAL_LOC));
            restMsg_exp1{2} = 'Take a rest! Press any key to continue the experiment...';
            restScreen = sys_prepInstructionScreen(window, restMsg_exp1, BG_COLOR, ...
                BLACK, text.TEXT_FONT, text.FONT_SIZE, centX, centY);
            
            Screen('CopyWindow', restScreen, window);
            Screen('Flip', window);
            KbStrokeWait;
        end
        
    end
    recordData{3} = respMat;
    
    % end of experiment screen
    Screen('CopyWindow', endExp2Screen, window);
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
    
    if mode(2) == 3
        for it_rep = 1:NUM_REPEAT
            for it_code = 1:length(code1_slow)
                Screen('DrawTexture',window,noCrossTextures{ code1_slow(it_code)*2+code2_slow(it_code)+1 });
                Screen('Flip', window);
                
                Screen('DrawTexture',window,baseTexture);
                Screen('Flip', window);
                
                % break the loop if press ESC key
                [~, ~, keyCode] = KbCheck;
                if keyCode(escKey)
                    Screen('CloseAll');
                    return;
                end
                
            end
            
        end
    elseif mode(2) == 0 % normal (30Hz modulation), grey, independent codes
        for it_rep = 1:NUM_REPEAT
            for it_code = 1:CODE_LENGTH
                if mode(1) == 0     % frequency modulated code
                    Screen('DrawTexture',window,noCrossTextures{ code1_fmc(it_code)*2+code2_fmc(it_code)+1 });
                elseif mode(1) == 1     % msequence
                    Screen('DrawTexture',window,noCrossTextures{ code1_mseq(it_code)*2+code2_mseq(it_code)+1 });
                elseif mode(1) == 2     % ssvep
                    Screen('DrawTexture',window,noCrossTextures{ code1_ssvep(it_code)*2+code2_ssvep(it_code)+1 });
                elseif mode(1) == 3     % image
                    Screen('DrawTexture',window,vepImage{ code1_fmc(it_code)+1 });
                elseif mode(1) == 4     % text
                    Screen('DrawTexture',window,vepText{ code1_fmc(it_code)+1 });
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
    else
        for it_rep = 1:NUM_REPEAT
            for it_code = 1:CODE_LENGTH
                if mode(2) == 1 % opposite codes
                    Screen('DrawTexture',window,noCrossTextures{ code1_fmc(it_code)+2 });
                elseif mode(2) == 2 % chromatic modulation
                    Screen('DrawTexture',window,chromeTextures{ 2, code1_fmc(it_code)*2+code2_fmc(it_code)+1, 1 });
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

function fixation_mask = gen_fixation_mask(stim, FIXATION_LOCATION, WINDOW_HEIGHT, WINDOW_WIDTH)

fixation_mask = ones(floor(WINDOW_HEIGHT/2)*2, floor(WINDOW_WIDTH/2)*2); % make sure the mask size is a even number to match with stimuli
fixation_mask( floor(WINDOW_HEIGHT/2) + [-floor(stim.FIXATION_WIDTH/2):floor(stim.FIXATION_WIDTH/2)], ...
    floor(FIXATION_LOCATION*WINDOW_WIDTH/100) + [-floor(stim.FIXATION_LENGTH/2):floor(stim.FIXATION_LENGTH/2)]) = 0;
fixation_mask( floor(WINDOW_HEIGHT/2) + [-floor(stim.FIXATION_LENGTH/2):floor(stim.FIXATION_LENGTH/2)] , ...
    floor(FIXATION_LOCATION*WINDOW_WIDTH/100) + [-floor(stim.FIXATION_WIDTH/2):floor(stim.FIXATION_WIDTH/2)] ) = 0;

end

function [contrastTextures, chromeTextures] = gen_contrast_textures(window, stim, config, WINDOW_HEIGHT, WINDOW_WIDTH)

% initialize textures for different contrasts
contrastTextures = cell(length(config.CONTRAST_LIST),4,3);  % 4: number of code combinations, 3: number of fixation locations (no, left, right)
chromeTextures = cell(length(config.CONTRAST_LIST),4,3);

fixation_mask_left = gen_fixation_mask(stim, stim.FIXATION_LEFT, WINDOW_HEIGHT, WINDOW_WIDTH);
fixation_mask_right = gen_fixation_mask(stim, stim.FIXATION_RIGHT, WINDOW_HEIGHT, WINDOW_WIDTH);

for contrast_id = 1:length(config.CONTRAST_LIST)
    
    cont_level = config.CONTRAST_LIST(contrast_id);
    
    % draw the two screen textures needed to flicker - grey stimuli
    plus_half_screen = (stim.BASELINE+cont_level) .* ones(WINDOW_HEIGHT,floor(WINDOW_WIDTH/2));
    minus_half_screen = (stim.BASELINE-cont_level) .* ones(WINDOW_HEIGHT,floor(WINDOW_WIDTH/2));
    
    % draw the two screen textures - chromatically modulated stimuli
    plus_half_screen_c = (stim.BASELINE+cont_level) .* ones(WINDOW_HEIGHT,floor(WINDOW_WIDTH/2),3);
    minus_half_screen_c = (stim.BASELINE-cont_level) .* ones(WINDOW_HEIGHT,floor(WINDOW_WIDTH/2),3);
    plus_half_screen_c(:,:,2) = plus_half_screen_c(:,:,2) - 2*cont_level;     % flip Green contrast
    minus_half_screen_c(:,:,2) = minus_half_screen_c(:,:,2) + 2*cont_level;
    
    % mask for smoothing at boundary
    [smooth_filter_ascend, smooth_filter_decend, smooth_filter_ascend_c, smooth_filter_decend_c] = gen_smooth_mask(config, WINDOW_HEIGHT, WINDOW_WIDTH, cont_level);
    
    % grey stimuli
    % no fixation cross
    contrastTextures{contrast_id, 1, 1} = Screen('MakeTexture',window,uint8( [minus_half_screen, minus_half_screen] ));   % code 0 / code 0
    contrastTextures{contrast_id, 2, 1} = Screen('MakeTexture',window,uint8( [minus_half_screen, plus_half_screen] + smooth_filter_ascend ));    % code 0 / code 1
    contrastTextures{contrast_id, 3, 1} = Screen('MakeTexture',window,uint8( [plus_half_screen, minus_half_screen] + smooth_filter_decend));    % code 1 / code 0
    contrastTextures{contrast_id, 4, 1} = Screen('MakeTexture',window,uint8( [plus_half_screen, plus_half_screen] ));     % code 1 / code 1
    
    % superimpose the fixation cross left
    contrastTextures{contrast_id, 1, 2} = Screen('MakeTexture',window,uint8( fixation_mask_left .* [minus_half_screen, minus_half_screen] ));   % code 0 / code 0
    contrastTextures{contrast_id, 2, 2} = Screen('MakeTexture',window,uint8( fixation_mask_left .* [minus_half_screen, plus_half_screen] + smooth_filter_ascend ));    % code 0 / code 1
    contrastTextures{contrast_id, 3, 2} = Screen('MakeTexture',window,uint8( fixation_mask_left .* [plus_half_screen, minus_half_screen] + smooth_filter_decend ));    % code 1 / code 0
    contrastTextures{contrast_id, 4, 2} = Screen('MakeTexture',window,uint8( fixation_mask_left .* [plus_half_screen, plus_half_screen] ));     % code 1 / code 1
    
    % superimpose the fixation cross right
    contrastTextures{contrast_id, 1, 3} = Screen('MakeTexture',window,uint8( fixation_mask_right .* [minus_half_screen, minus_half_screen] ));   % code 0 / code 0
    contrastTextures{contrast_id, 2, 3} = Screen('MakeTexture',window,uint8( fixation_mask_right .* [minus_half_screen, plus_half_screen] + smooth_filter_ascend ));    % code 0 / code 1
    contrastTextures{contrast_id, 3, 3} = Screen('MakeTexture',window,uint8( fixation_mask_right .* [plus_half_screen, minus_half_screen] + smooth_filter_decend ));    % code 1 / code 0
    contrastTextures{contrast_id, 4, 3} = Screen('MakeTexture',window,uint8( fixation_mask_right .* [plus_half_screen, plus_half_screen] ));     % code 1 / code 1
    
    % chromatic modulated stimuli
    % no fixation cross
    chromeTextures{contrast_id, 1, 1} = Screen('MakeTexture',window,uint8( cat(2,minus_half_screen_c, minus_half_screen_c) ));   % code 0 / code 0
    chromeTextures{contrast_id, 2, 1} = Screen('MakeTexture',window,uint8( cat(2,minus_half_screen_c, plus_half_screen_c) + smooth_filter_ascend_c ));    % code 0 / code 1
    chromeTextures{contrast_id, 3, 1} = Screen('MakeTexture',window,uint8( cat(2,plus_half_screen_c, minus_half_screen_c) + smooth_filter_decend_c ));    % code 1 / code 0
    chromeTextures{contrast_id, 4, 1} = Screen('MakeTexture',window,uint8( cat(2,plus_half_screen_c, plus_half_screen_c) ));     % code 1 / code 1
    
    % superimpose the fixation cross left
    chromeTextures{contrast_id, 1, 2} = Screen('MakeTexture',window,uint8( fixation_mask_left .* cat(2,minus_half_screen_c, minus_half_screen_c) ));   % code 0 / code 0
    chromeTextures{contrast_id, 2, 2} = Screen('MakeTexture',window,uint8( fixation_mask_left .* cat(2,minus_half_screen_c, plus_half_screen_c) + smooth_filter_ascend_c ));    % code 0 / code 1
    chromeTextures{contrast_id, 3, 2} = Screen('MakeTexture',window,uint8( fixation_mask_left .* cat(2,plus_half_screen_c, minus_half_screen_c) + smooth_filter_decend_c ));    % code 1 / code 0
    chromeTextures{contrast_id, 4, 2} = Screen('MakeTexture',window,uint8( fixation_mask_left .* cat(2,plus_half_screen_c, plus_half_screen_c) ));     % code 1 / code 1
    
    % superimpose the fixation cross right
    chromeTextures{contrast_id, 1, 3} = Screen('MakeTexture',window,uint8( fixation_mask_right .* cat(2,minus_half_screen_c, minus_half_screen_c) ));   % code 0 / code 0
    chromeTextures{contrast_id, 2, 3} = Screen('MakeTexture',window,uint8( fixation_mask_right .* cat(2,minus_half_screen_c, plus_half_screen_c) + smooth_filter_ascend_c ));    % code 0 / code 1
    chromeTextures{contrast_id, 3, 3} = Screen('MakeTexture',window,uint8( fixation_mask_right .* cat(2,plus_half_screen_c, minus_half_screen_c) + smooth_filter_decend_c ));    % code 1 / code 0
    chromeTextures{contrast_id, 4, 3} = Screen('MakeTexture',window,uint8( fixation_mask_right .* cat(2,plus_half_screen_c, plus_half_screen_c) ));     % code 1 / code 1
    
end
end


function [locTexture, noCrossTextures, baseFixTexture, baseTexture] = gen_base_textures(window, stim, config, WINDOW_HEIGHT, WINDOW_WIDTH)

% initialize textures for different fixation locations
locTexture = cell(length(config.LOCATION_LIST),4);  % 4: number of code combinations, 3: number of fixation locations
noCrossTextures = cell(1,4);

% draw the two screen textures needed to flicker
plus_half_screen = (stim.BASELINE+stim.CONTRAST) .* ones(WINDOW_HEIGHT,floor(WINDOW_WIDTH/2));
minus_half_screen = (stim.BASELINE-stim.CONTRAST) .* ones(WINDOW_HEIGHT,floor(WINDOW_WIDTH/2));

% generate fixation mask
fixation_mask_left = gen_fixation_mask(stim, stim.FIXATION_LEFT, WINDOW_HEIGHT, WINDOW_WIDTH);
fixation_mask_right = gen_fixation_mask(stim, stim.FIXATION_RIGHT, WINDOW_HEIGHT, WINDOW_WIDTH);

% mask for smoothing at boundary
[smooth_filter_ascend, smooth_filter_decend, ~, ~] = gen_smooth_mask(config, WINDOW_HEIGHT, WINDOW_WIDTH, stim.CONTRAST);

% draw texture without cross
noCrossTextures{1} = Screen('MakeTexture',window,uint8( [minus_half_screen, minus_half_screen] ));   % code 0 / code 0
noCrossTextures{2} = Screen('MakeTexture',window,uint8( [minus_half_screen, plus_half_screen] + smooth_filter_ascend ));    % code 0 / code 1
noCrossTextures{3} = Screen('MakeTexture',window,uint8( [plus_half_screen, minus_half_screen] + smooth_filter_decend));    % code 1 / code 0
noCrossTextures{4} = Screen('MakeTexture',window,uint8( [plus_half_screen, plus_half_screen] ));     % code 1 / code 1

% generate baseline stimuli
baseFixTexture = cell(1, 2+length(config.LOCATION_LIST));
base_screen = stim.BASELINE .* ones(floor(WINDOW_HEIGHT/2)*2,floor(WINDOW_WIDTH/2)*2);
baseTexture = Screen('MakeTexture',window,uint8( base_screen ));   % baseline / baseline
baseFixTexture{1} = Screen('MakeTexture',window,uint8( fixation_mask_left .* base_screen ));
baseFixTexture{2} = Screen('MakeTexture',window,uint8( fixation_mask_right .* base_screen ));

% draw texture with cross at different locations
for loc_id = 1:length(config.LOCATION_LIST)
    
    fixation_loc = config.LOCATION_LIST(loc_id);
    
    fixation_mask = gen_fixation_mask(stim, fixation_loc, WINDOW_HEIGHT, WINDOW_WIDTH);
    
    % superimpose the fixation cross left
    locTexture{loc_id, 1} = Screen('MakeTexture',window,uint8( fixation_mask .* [minus_half_screen, minus_half_screen] ));   % code 0 / code 0
    locTexture{loc_id, 2} = Screen('MakeTexture',window,uint8( fixation_mask .* [minus_half_screen, plus_half_screen] + smooth_filter_ascend ));    % code 0 / code 1
    locTexture{loc_id, 3} = Screen('MakeTexture',window,uint8( fixation_mask .* [plus_half_screen, minus_half_screen] + smooth_filter_decend ));    % code 1 / code 0
    locTexture{loc_id, 4} = Screen('MakeTexture',window,uint8( fixation_mask .* [plus_half_screen, plus_half_screen] ));     % code 1 / code 1
    
    % superimpose the fixation on base screen
    baseFixTexture{2+loc_id} = Screen('MakeTexture',window,uint8( fixation_mask .* base_screen ));
    
end





end



function [smooth_filter_ascend, smooth_filter_decend, smooth_filter_ascend_c, smooth_filter_decend_c] = gen_smooth_mask(config, WINDOW_HEIGHT, WINDOW_WIDTH, cont_level)

    % smooth mask for grey contrast 
    smooth_filter_ascend = zeros(WINDOW_HEIGHT, 2*floor(WINDOW_WIDTH/2));
    smooth_filter_decend = zeros(WINDOW_HEIGHT, 2*floor(WINDOW_WIDTH/2));
    if config.SMOOTH == 1 % linearly smoothed boundary
        % config.SMOOTH_WIDTH should be larger than contrast level
        smooth_width = max(1,floor(config.SMOOTH_WIDTH/2/cont_level));
        for it = 1:cont_level
            index_to_smooth = (it-1)*smooth_width+1 : it*smooth_width;
            smooth_filter_ascend(:,floor(WINDOW_WIDTH/2)-cont_level*smooth_width+index_to_smooth) = it;
            smooth_filter_ascend(:,floor(WINDOW_WIDTH/2)+1+cont_level*smooth_width-index_to_smooth) = -it;
            smooth_filter_decend(:,floor(WINDOW_WIDTH/2)-cont_level*smooth_width+index_to_smooth) = -it;
            smooth_filter_decend(:,floor(WINDOW_WIDTH/2)+1+cont_level*smooth_width-index_to_smooth) = it;
        end
    elseif config.SMOOTH == 2 % probablistic smoothed boundary
        prob_grad = linspace(0,100,config.SMOOTH_WIDTH)/100; % probability gradient
        randWin = rand(WINDOW_HEIGHT, config.SMOOTH_WIDTH);
        for it = 1:floor(config.SMOOTH_WIDTH/2)
            flipIndex = randWin(:,it) < prob_grad(it);
            smooth_filter_ascend(flipIndex, floor(WINDOW_WIDTH/2)-floor(config.SMOOTH_WIDTH/2)+it) = 2*cont_level;
            smooth_filter_decend(flipIndex, floor(WINDOW_WIDTH/2)+floor(config.SMOOTH_WIDTH/2)+1-it) = 2*cont_level;
        end
        for it = 1+floor(config.SMOOTH_WIDTH/2) : config.SMOOTH_WIDTH
            flipIndex = randWin(:,it) < prob_grad(it);
            smooth_filter_ascend(~flipIndex, floor(WINDOW_WIDTH/2)-floor(config.SMOOTH_WIDTH/2)+it) = -2*cont_level;
            smooth_filter_decend(~flipIndex, floor(WINDOW_WIDTH/2)+floor(config.SMOOTH_WIDTH/2)+1-it) = -2*cont_level;
        end
    end
    
    % smooth mask for color contrast
    smooth_filter_ascend_c = zeros(WINDOW_HEIGHT, 2*floor(WINDOW_WIDTH/2),3);
    smooth_filter_decend_c = zeros(WINDOW_HEIGHT, 2*floor(WINDOW_WIDTH/2),3);
    if config.SMOOTH == 1 % linearly smoothed boundary
        % config.SMOOTH_WIDTH should be larger than contrast level
        smooth_width = max(1,floor(config.SMOOTH_WIDTH/2/cont_level));
        for it = 1:cont_level
            index_to_smooth = (it-1)*smooth_width+1 : it*smooth_width;
            smooth_filter_ascend_c(:,floor(WINDOW_WIDTH/2)-cont_level*smooth_width+index_to_smooth,[1 3]) = it;
            smooth_filter_ascend_c(:,floor(WINDOW_WIDTH/2)+1+cont_level*smooth_width-index_to_smooth, [1 3]) = -it;
            smooth_filter_ascend_c(:,floor(WINDOW_WIDTH/2)-cont_level*smooth_width+index_to_smooth,2) = -it;
            smooth_filter_ascend_c(:,floor(WINDOW_WIDTH/2)+1+cont_level*smooth_width-index_to_smooth,2) = it;

            smooth_filter_decend_c(:,floor(WINDOW_WIDTH/2)-cont_level*smooth_width+index_to_smooth, [1 3]) = -it;
            smooth_filter_decend_c(:,floor(WINDOW_WIDTH/2)+1+cont_level*smooth_width-index_to_smooth, [1 3]) = it;
            smooth_filter_decend_c(:,floor(WINDOW_WIDTH/2)-cont_level*smooth_width+index_to_smooth, 2) = it;
            smooth_filter_decend_c(:,floor(WINDOW_WIDTH/2)+1+cont_level*smooth_width-index_to_smooth, 2) = -it;
        end
    elseif config.SMOOTH == 2 % probablistic smoothed boundary
        prob_grad = linspace(0,100,config.SMOOTH_WIDTH)/100; % probability gradient
        randWin = rand(WINDOW_HEIGHT, config.SMOOTH_WIDTH);
        for it = 1:floor(config.SMOOTH_WIDTH/2)
            flipIndex = randWin(:,it) < prob_grad(it);
            smooth_filter_ascend_c(flipIndex, floor(WINDOW_WIDTH/2)-floor(config.SMOOTH_WIDTH/2)+it, [1 3]) = 2*cont_level;
            smooth_filter_decend_c(flipIndex, floor(WINDOW_WIDTH/2)+floor(config.SMOOTH_WIDTH/2)-it, [1 3]) = 2*cont_level;
            smooth_filter_ascend_c(flipIndex, floor(WINDOW_WIDTH/2)-floor(config.SMOOTH_WIDTH/2)+it, 2) = -2*cont_level;
            smooth_filter_decend_c(flipIndex, floor(WINDOW_WIDTH/2)+floor(config.SMOOTH_WIDTH/2)-it, 2) = -2*cont_level;

        end
        for it = 1+floor(config.SMOOTH_WIDTH/2) : config.SMOOTH_WIDTH
            flipIndex = randWin(:,it) < prob_grad(it);
            smooth_filter_ascend_c(~flipIndex, floor(WINDOW_WIDTH/2)-floor(config.SMOOTH_WIDTH/2)+it, [1 3]) = -2*cont_level;
            smooth_filter_decend_c(~flipIndex, floor(WINDOW_WIDTH/2)+floor(config.SMOOTH_WIDTH/2)+1-it, [1 3]) = -2*cont_level;
            smooth_filter_ascend_c(~flipIndex, floor(WINDOW_WIDTH/2)-floor(config.SMOOTH_WIDTH/2)+it, 2) = 2*cont_level;
            smooth_filter_decend_c(~flipIndex, floor(WINDOW_WIDTH/2)+floor(config.SMOOTH_WIDTH/2)+1-it, 2) = 2*cont_level;
        end
    end
        
end


function imageTexture = gen_image_textures(image_screen, window, stim, config, WINDOW_HEIGHT, WINDOW_WIDTH)

% initialize textures for different contrasts
imageTexture = cell(length(config.CONTRAST_LIST),4,3);

fixation_mask_left = gen_fixation_mask(stim, stim.FIXATION_LEFT, WINDOW_HEIGHT, WINDOW_WIDTH);
fixation_mask_right = gen_fixation_mask(stim, stim.FIXATION_RIGHT, WINDOW_HEIGHT, WINDOW_WIDTH);

% resize image if the window size does not match
image_background = zeros(WINDOW_HEIGHT,2*floor(WINDOW_WIDTH/2),3);
[image_height, image_width, ~] = size(image_screen);
if image_height >= WINDOW_HEIGHT && image_width >= WINDOW_WIDTH
    image_background = double(image_screen( ...
        floor(image_height/2)-floor(WINDOW_HEIGHT/2) + (1:WINDOW_HEIGHT), ...
        floor(image_width/2)-floor(WINDOW_WIDTH/2) + (1:WINDOW_WIDTH), :));
elseif image_height < WINDOW_HEIGHT && image_width < WINDOW_WIDTH
    image_background( ...
        floor(WINDOW_HEIGHT/2)-floor(image_height/2) + (1:image_height), ...
        floor(WINDOW_WIDTH/2)-floor(image_width/2) + (1:image_width), :) ...
        = double(image_screen);
else
    error('Please adjust image size');
end


for contrast_id = 1:length(config.CONTRAST_LIST)
    
    cont_level = config.CONTRAST_LIST(contrast_id);
        
    % draw the two screen textures - chromatically modulated stimuli
    plus_left_screen_c = (image_background(:,1:floor(WINDOW_WIDTH/2),:)+cont_level) .* ones(WINDOW_HEIGHT,floor(WINDOW_WIDTH/2),3);
    minus_left_screen_c = (image_background(:,1:floor(WINDOW_WIDTH/2),:)-cont_level) .* ones(WINDOW_HEIGHT,floor(WINDOW_WIDTH/2),3);
    plus_right_screen_c = (image_background(:,floor(WINDOW_WIDTH/2)+1:end,:)+cont_level) .* ones(WINDOW_HEIGHT,floor(WINDOW_WIDTH/2),3);
    minus_right_screen_c = (image_background(:,floor(WINDOW_WIDTH/2)+1:end,:)-cont_level) .* ones(WINDOW_HEIGHT,floor(WINDOW_WIDTH/2),3);
    
    % mask for smoothing at boundary
    [~, ~, smooth_filter_ascend_c, smooth_filter_decend_c] = gen_smooth_mask(config, WINDOW_HEIGHT, WINDOW_WIDTH, cont_level);
    
    % no fixation cross
    imageTexture{contrast_id, 1, 1} = Screen('MakeTexture',window,uint8( cat(2,minus_left_screen_c, minus_right_screen_c) ));   % code 0 / code 0
    imageTexture{contrast_id, 2, 1} = Screen('MakeTexture',window,uint8( cat(2,minus_left_screen_c, plus_right_screen_c) + smooth_filter_ascend_c ));    % code 0 / code 1
    imageTexture{contrast_id, 3, 1} = Screen('MakeTexture',window,uint8( cat(2,plus_left_screen_c, minus_right_screen_c) + smooth_filter_decend_c ));    % code 1 / code 0
    imageTexture{contrast_id, 4, 1} = Screen('MakeTexture',window,uint8( cat(2,plus_left_screen_c, plus_right_screen_c) ));     % code 1 / code 1
    
    % superimpose the fixation cross left
    imageTexture{contrast_id, 1, 2} = Screen('MakeTexture',window,uint8( fixation_mask_left .* cat(2,minus_left_screen_c, minus_right_screen_c) ));   % code 0 / code 0
    imageTexture{contrast_id, 2, 2} = Screen('MakeTexture',window,uint8( fixation_mask_left .* cat(2,minus_left_screen_c, plus_right_screen_c) + smooth_filter_ascend_c ));    % code 0 / code 1
    imageTexture{contrast_id, 3, 2} = Screen('MakeTexture',window,uint8( fixation_mask_left .* cat(2,plus_left_screen_c, minus_right_screen_c) + smooth_filter_decend_c ));    % code 1 / code 0
    imageTexture{contrast_id, 4, 2} = Screen('MakeTexture',window,uint8( fixation_mask_left .* cat(2,plus_left_screen_c, plus_right_screen_c) ));     % code 1 / code 1
    
    % superimpose the fixation cross right
    imageTexture{contrast_id, 1, 3} = Screen('MakeTexture',window,uint8( fixation_mask_right .* cat(2,minus_left_screen_c, minus_right_screen_c) ));   % code 0 / code 0
    imageTexture{contrast_id, 2, 3} = Screen('MakeTexture',window,uint8( fixation_mask_right .* cat(2,minus_left_screen_c, plus_right_screen_c) + smooth_filter_ascend_c ));    % code 0 / code 1
    imageTexture{contrast_id, 3, 3} = Screen('MakeTexture',window,uint8( fixation_mask_right .* cat(2,plus_left_screen_c, minus_right_screen_c) + smooth_filter_decend_c ));    % code 1 / code 0
    imageTexture{contrast_id, 4, 3} = Screen('MakeTexture',window,uint8( fixation_mask_right .* cat(2,plus_left_screen_c, plus_right_screen_c) ));     % code 1 / code 1
    
end

end