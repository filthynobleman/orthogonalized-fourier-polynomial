function cmap = blend(Colors, Resolution)
%cmap = BLEND(Colors, Resolution) Generates a colormap with the given
%resolution which passes through all the colors in the given matrix. Notice
%the matrix Colors is n-by-3.
%
%cmap = BLEND(Colors) Uses the default resolution of 256 entries


    % Infer default resolution, if needed
    if nargin < 2
        Resolution = 256;
    end
    
    % Get the number of colors and the number of entries between each color
    NumColors = size(Colors, 1);
    ColorEntries = floor(Resolution / (NumColors - 1));
    
    % Initialize the colormap
    if NumColors == 1
        error("Cannot blend one color. Please, use utils.cmaps.constant for this.");
    end
    cmap = zeros(Resolution, 3);
    
    % Fill the colormap
    for i = 1:NumColors-1
        % Compute the begin and end indices
        Begin = (i - 1) * ColorEntries + 1;
        End = Begin + ColorEntries;
        % Compute the channels
        Channels.R = linspace(Colors(i, 1), Colors(i + 1, 1), ColorEntries)';
        Channels.G = linspace(Colors(i, 2), Colors(i + 1, 2), ColorEntries)';
        Channels.B = linspace(Colors(i, 3), Colors(i + 1, 3), ColorEntries)';
        % Fill the map
        cmap(Begin:End-1, :) = [Channels.R, Channels.G, Channels.B];
    end
    
    
    

end

