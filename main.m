function recordData = main(config)
%% Gaze VEP Experiment: Pilot

%% Define parameters
stim.CONTRAST = config.CONTRAST;  % BASELINE-CONTRAST and BASELINE+CONTRAST
stim.BASELINE = config.BASELINE;

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
NUM_REPEAT = 10;    % for showing stimuli only

% parameters for drawing text
BLACK = [0,0,0];
text.TEXT_FONT = 'Arial';
text.FONT_SIZE = 36;


%% read images
image_screen = [];
if strcmp(config.MODE, 'fmc_single_image')
    image_screen = imread(config.filename_image);
elseif strcmp(config.MODE, 'fmc_single_text')
    image_screen = imread(config.filename_text);
end


%% ------------------------------------------------------------------------
%              Generate code sequence
%--------------------------------------------------------------------------

if strcmp(config.MODE, 'fmc_opposite_slow') || strcmp(config.MODE, 'fmc_binary_slow')
    % code length for fmc_slow: (2^nBits - 1) * 4
    nBits = ceil(log2(config.STIM_LEN * config.REFRESH / 4 + 1));   % ensure the code does not repeat for given stimulation length
    code = gen_code_fmc(nBits,CODE_LENGTH,SEED_1);
    code2 = gen_code_fmc(nBits,CODE_LENGTH,SEED_2);
elseif strcmp(config.MODE, 'fmc_binary') || strcmp(config.MODE, 'fmc_binary_chrome') || ...
        strcmp(config.MODE, 'fmc_opposite') || strcmp(config.MODE, 'fmc_opposite_chrome') || ...
        strcmp(config.MODE, 'fmc_single_image') || strcmp(config.MODE, 'fmc_single_text')
    % code length for fmc: (2^nBits - 1) * 2
    nBits = ceil(log2(config.STIM_LEN * config.REFRESH / 2 + 1));   % ensure the code does not repeat for given stimulation length
    code = gen_code_fmc(nBits,CODE_LENGTH,SEED_1);
    code2 = gen_code_fmc(nBits,CODE_LENGTH,SEED_2);
elseif strcmp(config.MODE, 'mseq_opposite') || strcmp(config.MODE, 'mseq_binary')
    % code length for m-sequence: (2^nBits - 1)
    nBits = ceil(log2(config.STIM_LEN * config.REFRESH + 1));   % ensure the code does not repeat for given stimulation length
    code = gen_msequence(nBits,CODE_LENGTH,SEED_1);
    code2 = gen_msequence(nBits,CODE_LENGTH,SEED_2);
elseif strcmp(config.MODE, 'ssvep')
    % default 10Hz and 12Hz
    code = gen_code_ssvep(10, CODE_LENGTH, config.REFRESH);
    code2 = gen_code_ssvep(12, CODE_LENGTH, config.REFRESH);
else
    error('Mode input is incorrect');
end
if ~any(code~=code2) % make sure two codes are not identical
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
if isempty(image_screen)
    if ~config.FULLSCREEN
        WINDOW_WIDTH = config.WIN_WIDTH;
        WINDOW_HEIGHT = config.WIN_HEIGHT;
        [window, windowRect] = Screen('OpenWindow', screenNumber, BG_COLOR, [0, 0, WINDOW_WIDTH, WINDOW_HEIGHT], [], 2);
    else
        [WINDOW_WIDTH, WINDOW_HEIGHT] = Screen('WindowSize',screenNumber);
        [window, windowRect] = Screen('OpenWindow', screenNumber, BG_COLOR, [], [], 2);
    end
else
    [WINDOW_HEIGHT, WINDOW_WIDTH, ~] = size(image_screen);
    [window, windowRect] = Screen('OpenWindow', screenNumber, BG_COLOR, [0, 0, WINDOW_WIDTH, WINDOW_HEIGHT], [], 2);
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

demoMsg = cell(1,2);
demoMsg{1} = 'Demo Experimental Stimuli';
demoMsg{2} = 'Press Any Key To Begin';
demoScreen = sys_prepInstructionScreen(window, demoMsg, BG_COLOR, ...
    BLACK, text.TEXT_FONT, text.FONT_SIZE, centX, centY);


%% ------------------------------------------------------------------------
%              Generate Stimulation Textures (stimulus screen)
%--------------------------------------------------------------------------

% prepare grey binary background with different contrast levels and fixation locations
[contrastTextures, locTexture, noCrossTextures, baseFixTexture, baseTexture] = gen_textures(window, stim, config, WINDOW_HEIGHT, WINDOW_WIDTH);

if ~isempty(image_screen)
    % image
    vepImage{1} = Screen('MakeTexture',window,image_screen-stim.CONTRAST);
    vepImage{2} = Screen('MakeTexture',window,image_screen+stim.CONTRAST);

    % text
    vepText{1} = Screen('MakeTexture',window,image_screen);
    vepText{2} = Screen('MakeTexture',window,image_screen+2*stim.CONTRAST);
end


% chromatically modulated stimuli
plus_half_screen_c = (stim.BASELINE+stim.CONTRAST) .* ones(WINDOW_HEIGHT,floor(WINDOW_WIDTH/2),3);
minus_half_screen_c = (stim.BASELINE-stim.CONTRAST) .* ones(WINDOW_HEIGHT,floor(WINDOW_WIDTH/2),3);
plus_half_screen_c(:,:,2) = plus_half_screen_c(:,:,2) - 2*stim.CONTRAST;     % flip Green contrast
minus_half_screen_c(:,:,2) = minus_half_screen_c(:,:,2) + 2*stim.CONTRAST;

chromeTextures{1} = Screen('MakeTexture',window,uint8( cat(2,minus_half_screen_c, minus_half_screen_c) ));   % code 0 / code 0
chromeTextures{2} = Screen('MakeTexture',window,uint8( cat(2,minus_half_screen_c, plus_half_screen_c) ));    % code 0 / code 1
chromeTextures{3} = Screen('MakeTexture',window,uint8( cat(2,plus_half_screen_c, minus_half_screen_c) ));    % code 1 / code 0
chromeTextures{4} = Screen('MakeTexture',window,uint8( cat(2,plus_half_screen_c, plus_half_screen_c) ));     % code 1 / code 1



%% ----------------------------------------------------------------------
%             Start the Main Loop for Experiment
%----------------------------------------------------------------------

Priority(1);
recordData = cell(1,2);

%% ----------------------------------------------------------------------
%             Experimental loop - Task 1: contrast levels
%----------------------------------------------------------------------

if config.RUN_EXP_1    
    
    disp('Enter Experiment: Task 1 - contrast levels');
    
    % Event marker
    if config.ENABLE_LSL, outlet.push_sample({'Exp1'}); end

    % obtain randomized trial order for different conditions
    trial_order = repmat( [1:length(config.CONTRAST_LIST)*2], 1, config.NUM_TRIAL_CONTRAST);
    trial_order = trial_order( randperm(length(trial_order)) );
    
    % present first screen
    Screen('CopyWindow', startExp1Screen, window);
    Screen('Flip', window);
    KbStrokeWait;
    
    respMat = [];
    respMat.contrast = zeros(1,length(trial_order));
    respMat.location = zeros(1,length(trial_order));
    respMat.rate_percept = zeros(1,length(trial_order));
    respMat.rt = zeros(1,length(trial_order));
    respMat.mode = config.MODE;
    respMat.code = {code, code2};
    
    for trial = 1:length(trial_order)
        
        cont_level = mod(trial_order(trial)-1, length(config.CONTRAST_LIST)) + 1;
        loc = floor( (trial_order(trial)-1) / length(config.CONTRAST_LIST)) + 1;

        % Event marker
        if config.ENABLE_LSL, outlet.push_sample({sprintf('Cross_%d',loc)}); end

        % present cross for eye fixation
        Screen('DrawTexture', window, baseFixTexture{loc});
        Screen('Flip',window);
        
        % wait one second
        WaitSecs(1);
        
        % Event marker
        if config.ENABLE_LSL
            outlet.push_sample({sprintf('Stim%d%d%d%d_CT%d_Loc%d', config.MODE_STIM, ...
                config.MODE_SIDE, config.MODE_COLOR, config.MODE_NORM, config.CONTRAST_LIST(cont_level), loc)}); 
        end

        % present stimuli according to the code sequence
        for it_code = 1:CODE_LENGTH
            
            % Draw it on the back buffer
            if strcmp(config.MODE, 'fmc_opposite') || strcmp(config.MODE, 'mseq_opposite')
                Screen('DrawTexture',window,contrastTextures{ cont_level, code(it_code)+2, loc+1 });
            elseif strcmp(config.MODE, 'fmc_binary') || strcmp(config.MODE, 'mseq_binary') || strcmp(config.MODE, 'ssvep')
                Screen('DrawTexture',window,contrastTextures{ cont_level, code(it_code)*2+code2(it_code)+1, loc+1 });
            end
            
            % Tell PTB no more drawing commands will be issued until the next flip
            Screen('DrawingFinished', window);
            
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
        Screen('CopyWindow', rateScreen, window);
        Screen('Flip', window);
        
        % Check the keyboard
        restTic = tic;
        while true
            if toc(restTic) > 30, break; end
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

        % Record the trial data
        respMat.contrast(trial) = config.CONTRAST_LIST(cont_level);
        respMat.location(trial) = loc;
        respMat.rate_percept(trial) = response;
        respMat.rt(trial) = toc(restTic);
                  
        fprintf('Trial: %2d, Contrast: %2d, Location: %d, Rating: %d, RT: %f\n', ...
            trial, config.CONTRAST_LIST(cont_level), loc, response, respMat.rt(trial));
        
    end
    recordData{1} = respMat;
    
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
    respMat.mode = config.MODE;
    respMat.code = {code, code2};

    for trial = 1:length(trial_order)
        
        loc = trial_order(trial);
        
        % Event marker
        if config.ENABLE_LSL, outlet.push_sample({sprintf('Cross_%d',loc)}); end

        % present cross for eye fixation
        Screen('DrawTexture', window, baseFixTexture{2+loc});
        Screen('Flip',window);
        
        % wait one second
        WaitSecs(1);
        
        % Event marker
        if config.ENABLE_LSL
            outlet.push_sample({sprintf('Stim%d%d%d%d_CT%d_Loc%d', config.MODE_STIM, ...
                config.MODE_SIDE, config.MODE_COLOR, config.MODE_NORM, stim.CONTRAST, config.LOCATION_LIST(loc))});
        end
        
        % present stimuli according to the code sequence
        for it_code = 1:CODE_LENGTH
            
            % Draw it on the back buffer
            if strcmp(config.MODE, 'fmc_opposite') || strcmp(config.MODE, 'mseq_opposite')
                Screen('DrawTexture',window,locTexture{ loc, code(it_code)+2});
            elseif strcmp(config.MODE, 'fmc_binary') || strcmp(config.MODE, 'mseq_binary') || strcmp(config.MODE, 'ssvep')
                Screen('DrawTexture',window,locTexture{ loc, code(it_code)*2+code2(it_code)+1 });
            end
            
            % Tell PTB no more drawing commands will be issued until the next flip
            Screen('DrawingFinished', window);
            
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
                  
        fprintf('Trial: %2d, Contrast: %2d, Location: %d\n', ...
            trial, stim.CONTRAST, config.LOCATION_LIST(loc));
        
        % present break slide and wait one second
        Screen('DrawTexture',window,baseTexture);
        Screen('Flip', window);
        WaitSecs(1);
        
    end     
    recordData{2} = respMat;

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
    
    % present first screen
    Screen('CopyWindow', demoScreen, window);
    Screen('Flip', window);
    KbStrokeWait;

    if strcmp(config.MODE, 'fmc_opposite_slow') || strcmp(config.MODE, 'fmc_binary_slow')
        for it_rep = 1:NUM_REPEAT
            for it_code = 1:CODE_LENGTH
                
                if strcmp(config.MODE, 'fmc_opposite_slow')
                    Screen('DrawTexture',window,noCrossTextures{ code(it_code)+2 });
                    Screen('Flip', window);
                    
                    Screen('DrawTexture',window,baseTexture);
                    Screen('Flip', window);
                    
                elseif strcmp(config.MODE, 'fmc_binary_slow')
                    Screen('DrawTexture',window,noCrossTextures{ code(it_code)*2+code2(it_code)+1 });
                    Screen('Flip', window);
                    
                    Screen('DrawTexture',window,baseTexture);
                    Screen('Flip', window);
                end
                
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
                % Draw it on the back buffer
                if strcmp(config.MODE, 'fmc_opposite') || strcmp(config.MODE, 'mseq_opposite')
                    Screen('DrawTexture',window,noCrossTextures{ code(it_code)+2 });
                elseif strcmp(config.MODE, 'fmc_binary') || strcmp(config.MODE, 'mseq_binary') || strcmp(config.MODE, 'ssvep')
                    Screen('DrawTexture',window,noCrossTextures{ code(it_code)*2+code2(it_code)+1 });
                elseif strcmp(config.MODE, 'fmc_opposite_chrome')
                    Screen('DrawTexture',window,chromeTextures{ code(it_code)+2 });
                elseif strcmp(config.MODE, 'fmc_binary_chrome')
                    Screen('DrawTexture',window,chromeTextures{ code(it_code)*2+code2(it_code)+1 });
                elseif strcmp(config.MODE, 'fmc_single_image')
                    Screen('DrawTexture',window,vepImage{ code(it_code)+1 });
                elseif strcmp(config.MODE, 'fmc_single_text')
                    Screen('DrawTexture',window,vepText{ code(it_code)+1 });
                end
                
                
                % Tell PTB no more drawing commands will be issued until the next flip
                Screen('DrawingFinished', window);
                
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
sca

end


function code = gen_code_ssvep(freq, code_length, refresh)

EPSILON = 1e-12;

index = 0:1:code_length-1;
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

function [contrastTextures, locTexture, noCrossTextures, baseFixTexture, baseTexture] = gen_textures(window, stim, config, WINDOW_HEIGHT, WINDOW_WIDTH)

% initialize textures for different contrasts
contrastTextures = cell(length(config.CONTRAST_LIST),4,3);  % 4: number of code combinations, 3: number of fixation locations

fixation_mask_left = gen_fixation_mask(stim, stim.FIXATION_LEFT, WINDOW_HEIGHT, WINDOW_WIDTH);
fixation_mask_right = gen_fixation_mask(stim, stim.FIXATION_RIGHT, WINDOW_HEIGHT, WINDOW_WIDTH);

for contrast_id = 1:length(config.CONTRAST_LIST)
    
    cont_level = config.CONTRAST_LIST(contrast_id);
    
    % draw the two screen textures needed to flicker
    plus_half_screen = (stim.BASELINE+cont_level) .* ones(WINDOW_HEIGHT,floor(WINDOW_WIDTH/2));
    minus_half_screen = (stim.BASELINE-cont_level) .* ones(WINDOW_HEIGHT,floor(WINDOW_WIDTH/2));
    
    % inearly smoothed boundary
    smooth_filter_ascend = zeros(WINDOW_HEIGHT, 2*floor(WINDOW_WIDTH/2));
    smooth_filter_decend = zeros(WINDOW_HEIGHT, 2*floor(WINDOW_WIDTH/2));
    if config.SMOOTH
        % config.SMOOTH_WIDTH should be larger than contrast level
        smooth_width = max(1,floor(config.SMOOTH_WIDTH/2/cont_level));
        for it = 1:cont_level
            index_to_smooth = (it-1)*smooth_width+1 : it*smooth_width;
            smooth_filter_ascend(:,floor(WINDOW_WIDTH/2)-cont_level*smooth_width+index_to_smooth) = it;
            smooth_filter_ascend(:,floor(WINDOW_WIDTH/2)+1+cont_level*smooth_width-index_to_smooth) = -it;
            smooth_filter_decend(:,floor(WINDOW_WIDTH/2)-cont_level*smooth_width+index_to_smooth) = -it;
            smooth_filter_decend(:,floor(WINDOW_WIDTH/2)+1+cont_level*smooth_width-index_to_smooth) = it;
        end
    end
    
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
    
end


% initialize textures for different fixation locations
locTexture = cell(length(config.LOCATION_LIST),4);  % 4: number of code combinations, 3: number of fixation locations
noCrossTextures = cell(1,4);

% draw the two screen textures needed to flicker
plus_half_screen = (stim.BASELINE+stim.CONTRAST) .* ones(WINDOW_HEIGHT,floor(WINDOW_WIDTH/2));
minus_half_screen = (stim.BASELINE-stim.CONTRAST) .* ones(WINDOW_HEIGHT,floor(WINDOW_WIDTH/2));

% inearly smoothed boundary
smooth_filter_ascend = zeros(WINDOW_HEIGHT, 2*floor(WINDOW_WIDTH/2));
smooth_filter_decend = zeros(WINDOW_HEIGHT, 2*floor(WINDOW_WIDTH/2));
if config.SMOOTH
    % config.SMOOTH_WIDTH should be larger than contrast level
    smooth_width = max(1,floor(config.SMOOTH_WIDTH/2/stim.CONTRAST));
    for it = 1:stim.CONTRAST
        index_to_smooth = (it-1)*smooth_width+1 : it*smooth_width;
        smooth_filter_ascend(:,floor(WINDOW_WIDTH/2)-stim.CONTRAST*smooth_width+index_to_smooth) = it;
        smooth_filter_ascend(:,floor(WINDOW_WIDTH/2)+1+stim.CONTRAST*smooth_width-index_to_smooth) = -it;
        smooth_filter_decend(:,floor(WINDOW_WIDTH/2)-stim.CONTRAST*smooth_width+index_to_smooth) = -it;
        smooth_filter_decend(:,floor(WINDOW_WIDTH/2)+1+stim.CONTRAST*smooth_width-index_to_smooth) = it;
    end
end

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