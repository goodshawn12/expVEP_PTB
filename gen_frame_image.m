function [frame_image, frame_contrast, fixation_mask, frame_base] = gen_frame_image(H,W,N,cont_base,fix_len,fix_thick)

% % parameter settings
% H = 960;    % window height
% W = 1280;   % window width
% N = 10;     % # trials / images
% cont_base = 128;    % baseline intensity
% fix_len = 30;       % fixation cross - length
% fix_thick = 5;      % fixation cross - thickness


%%%%%%%%%%%%%%%%%%%%%%%
% generate image frames
%%%%%%%%%%%%%%%%%%%%%%%

% load images
raw_image = cell(1,N);
for it = 1:N
    raw_image{it} = imread(sprintf('img_%02d.jpg',it));
end
[imgH, imgW, ~] = size(raw_image{1});

% organize images
frame_base = cont_base * ones(H,W,3);
frame_image = cell(1,N);
loc_h = floor((H/2 - imgH)/2);
loc_w = floor((W/2 - imgW)/2);

for it = 1:N
    frame_image_tmp = frame_base;
    
    frame_image_tmp( loc_h+(1:imgH), loc_w+(1:imgW), : ) = raw_image{it};
    frame_image_tmp( loc_h+(1:imgH), loc_w+(1:imgW)+W/2, : ) = raw_image{it};
    frame_image_tmp( loc_h+(1:imgH)+H/2, loc_w+(1:imgW), : ) = raw_image{it};
    frame_image_tmp( loc_h+(1:imgH)+H/2, loc_w+(1:imgW)+W/2, : ) = raw_image{it};
    frame_image{it} = frame_image_tmp;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate contrast modulation frames
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
frame_contrast = cell(1,16);

for it = 1:16
    tl = rem(it-1,2)+1;
    tr = rem(ceil(it/2)-1,2)+1;
    bl = rem(ceil(it/4)-1,2)+1;
    br = rem(ceil(it/8)-1,2)+1;
    
    clear tmp;
    tmp = zeros(H,W,3);
    tmp( loc_h+(1:imgH), loc_w+(1:imgW), : ) = sign(tl-1.5);
    tmp( loc_h+(1:imgH), loc_w+(1:imgW)+W/2, : ) = sign(tr-1.5);
    tmp( loc_h+(1:imgH)+H/2, loc_w+(1:imgW), : ) = sign(bl-1.5);
    tmp( loc_h+(1:imgH)+H/2, loc_w+(1:imgW)+W/2, : ) = sign(br-1.5);
    frame_contrast{it} = tmp;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate fixation frames
%%%%%%%%%%%%%%%%%%%%%%%%%%%
fixation_mask = cell(1,4);
for it = 1:4
    it1 = rem(it-1,2)+1;
    it2 = rem(ceil(it/2)-1,2)+1;
    fixation_mask{it} = ones(H,W);
    fixation_mask{it}( ...
        floor(H/4) + (it2-1)*H/2 + [-ceil(fix_len/2):floor(fix_len/2)], ...
        floor(W/4) + (it1-1)*W/2 + [-ceil(fix_thick/2):floor(fix_thick/2)]) = 0;
    fixation_mask{it}( ...
        floor(H/4) + (it2-1)*H/2 + [-ceil(fix_thick/2):floor(fix_thick/2)], ...
        floor(W/4) + (it1-1)*W/2 + [-ceil(fix_len/2):floor(fix_len/2)]) = 0;
end


