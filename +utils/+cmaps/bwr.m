function cm = bwr(Resolution)
%cm = BWR(Resolution) Returns the blue-white-red colormap with the desired
%resolution
%
%cm = BWR Returns the blue-white-red colormap with resolution 256

    if nargin < 1
        Resolution = 256;
    end
    
    Colors = [0 0 0.8;    % Blue
              1 1 1;    % White
              0.8 0 0];   % Red
    cm = utils.cmaps.blend(Colors, Resolution);

end

