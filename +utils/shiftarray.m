function shifted = shiftarray(arr, shift)
%SHIFTARRAY Summary of this function goes here
%   Detailed explanation goes here

    if shift == 0
        shifted = arr;
        return;
    end
    
    n = length(arr);
    if shift > 0
        % Shifting down means taking last k elements goes to top
        bottom = reshape(arr((end - shift + 1):end), shift, 1);
        top = reshape(arr(1:(end - shift)), n - shift, 1);
    else
        shift = - shift;
        % Shifting up means first k elements goes to bottom
        bottom = reshape(arr((shift + 1):end), n - shift, 1);
        top = reshape(arr(1:shift), shift, 1);
    end
    shifted = [bottom; top];

end

