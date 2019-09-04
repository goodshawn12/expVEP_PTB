function handle = sys_prepInstructionScreen(window, msg,...
    bg_color, text_color, text_font, font_size, centX, centY)

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
    for it = 1:nLine
        bounds = Screen(handle, 'TextBounds', msg{it});
        Screen('DrawText', handle, msg{it}, ...
            centX-bounds(RectRight)/2, centY-(nLine/2-it+1)*bounds(RectBottom), text_color);
    end
end