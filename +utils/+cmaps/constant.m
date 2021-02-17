function cm = constant(Color, Resolution)
%cm = CONSTANT(Color) Create a colormap with a constant given RGB color.
%The map resolution is 256. Color can be a scalar, in that case it is
%interpreted as a color to assign to all the three channels.
%
%cm = CONSTANT(Color, Resolution) Also specifies the map resolution.


    if size(Color, 1) > 1
        error("A constant colormap can only have one color.");
    end
    if size(Color, 2) == 1
        Color = repmat(Color, 1, 3);
    end
    
    if nargin < 2
        Resolution = 256;
    end
    
    cm = repmat(Color, Resolution, 1);
end

