function frames2gif(varargin)
%FRAMES2GIF Summary of this function goes here
%   Detailed explanation goes here

    args = parse_args(varargin{:});
    check_args_validity(args);
    
    for i = min(args.range):max(args.range)
        frame_name = strcat(args.dir, "/", args.prefix, args.frameName, ...
                            string(i), args.suffix, ".", args.extension);
        img = imread(frame_name, args.extension);
        [im, cmap] = rgb2ind(img, 256);
        if i == min(args.range)
            imwrite(im, cmap, args.output, 'gif', ...
                    'LoopCount', inf, ...
                    'DelayTime', 1 / args.fps);
        else
            imwrite(im, cmap, args.output, 'gif', ...
                    'WriteMode', 'append', ...
                    'DelayTime', 1 / args.fps);
        end
    end
end

function check_args_validity(args)
    if args.dir == ""
        error("No input directory given.");
    elseif args.frameName == ""
        error("No input frame name given.");
    elseif args.range(1) > args.range(2)
        error("The frame range is invalid.");
    elseif args.extension == ""
        error("No frame extension given.");
    elseif args.output == ""
        error("No output filename given.");
    elseif args.fps <= 0 || isinf(args.fps)
        error("FPS and Delay must be strictly positive and finite.");
    end
end

function args = parse_args(varargin)
    args.dir = "";
    args.frameName = "frame";
    args.prefix = "";
    args.suffix = "";
    args.range = [0 0];
    args.extension = "png";
    args.output = "generated.gif";
    args.overwrite = false;
    args.fps = 1;
    
    for i = 1:2:length(varargin)
        if varargin{i} == "InputDir"
            args.dir = string(varargin{i + 1});
        elseif varargin{i} == "FrameName"
            args.frameName = string(varargin{i + 1});
        elseif varargin{i} == "FramePrefix"
            args.prefix = string(varargin{i + 1});
        elseif varargin{i} == "FrameSuffix"
            args.suffix = string(varargin{i + 1});
        elseif varargin{i} == "FrameRange"
            args.range = varargin{i + 1};
        elseif varargin{i} == "FrameExtension"
            args.extension = string(varargin{i + 1});
        elseif varargin{i} == "OutputFile"
            args.output = string(varargin{i + 1});
            if ~endsWith(args.output, ".gif", "IgnoreCase", true)
                args.output = strcat(args.output, ".gif");
            end
        elseif varargin{i} == "FramePerSecond"
            args.fps = varargin{i + 1};
        elseif varargin{i} == "FrameDelay"
            args.fps = 1 / varargin{i + 1};
        elseif varargin{i} == "Overwrite"
            args.overwrite = varargin{i + 1};
        end
    end
    
    if exist(args.output, 'file') && ~args.overwrite
        i = 0;
        while true
            if exist(replace(args.output, ...
                             ".gif", ...
                             strcat("_", ...
                                    string(i), ...
                                    ".gif")), ...
                     'file')
                i = i + 1;
            else
                break;
            end
        end
    end
end