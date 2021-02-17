function M = rotate(M, axis, angle, angle_type)
%M = ROTATE(M, axis, angle)
%   Rotate the mesh along the given axis and with the given angle. The
%   angle is assumed to be in degrees.
%
%M = ROTATE(M, axis, angle, angle_type)
%   Rotate the mesh along the given axis and with the given angle. The
%   angle interpretation is 'degree' or 'radians', as specified by
%   angle_type.

    if nargin < 4
        angle_type = 'degree';
    end
    
    if lower(angle_type) == "degree"
        angle = deg2rad(angle);
    elseif lower(angle_type) ~= "radians"
        error("%s is not a valid angle type.", angle_type);
    end
    
    % Create the rotation matrix
    if lower(axis) == "x"
        R = [1, 0, 0;
             0, cos(angle), -sin(angle);
             0, sin(angle), cos(angle)];
    elseif lower(axis) == "y"
        R = [cos(angle), 0, sin(angle);
             0, 1, 0;
             -sin(angle), 0, cos(angle)];
    elseif lower(axis) == "z"
        R = [cos(angle), -sin(angle), 0;
             sin(angle), cos(angle), 0,
             0, 0, 1];
    else
        error("%s is not a valid axis name.", axis);
    end

    M.VERT = M.VERT * R';
    M.X = M.VERT(:, 1);
    M.Y = M.VERT(:, 2);
    M.Z = M.VERT(:, 3);
end

