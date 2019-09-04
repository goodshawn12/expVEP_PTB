function handle = sys_prepInstructionScreen_align(window, msg,...
    bg_color, text_color, text_font, font_size, centX, centY, align_left)

handle = Screen(window, 'OpenOffScreenWindow', bg_color);
Screen(handle, 'TextColor', text_color);
Screen(handle, 'TextFont', text_font);
Screen(handle, 'TextSize', font_size);

if ~iscell(msg)
    bounds = Screen(handle, 'TextBounds', msg);
    Screen('DrawText', handle, msg, ...
        centX-bounds(RectRight)/2, centY-bounds(RectBottom)/2, text_color);
else
    nLine = length(msg);
    xbounds = zeros(1,nLine);
    for it = 1:nLine
        bounds = Screen(handle, 'TextBounds', msg{it});
        xbounds(it) = bounds(RectRight)/2;
    end
    
    max_bounds = max(xbounds);
    for it = 1:nLine
        if ~align_left
            Screen('DrawText', handle, msg{it}, ...
                centX-xbounds(it), centY-(nLine/2-it+1)*bounds(RectBottom), text_color);
        else    % align text to the left of the first line
            Screen('DrawText', handle, msg{it}, ...
                centX-max_bounds, centY-(nLine/2-it+1)*bounds(RectBottom), text_color);
        end
    end
    
end