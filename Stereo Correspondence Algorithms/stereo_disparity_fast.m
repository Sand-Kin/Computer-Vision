function [Id] = stereo_disparity_fast(Il, Ir, bbox)
% STEREO_DISPARITY_FAST Fast stereo correspondence algorithm.
%
%  Id = STEREO_DISPARITY_FAST(Il, Ir, bbox) computes a stereo disparity image
%  from left stereo image Il and right stereo image Ir.
%
%  Inputs:
%  -------
%   Il    - Left stereo image, m x n pixels, colour or greyscale.
%   Ir    - Right stereo image, m x n pixels, colour or greyscale.
%   bbox  - Bounding box, relative to left image, top left corner, bottom
%           right corner (inclusive). Height is u, width is v.
%
%  Outputs:
%  --------
%   Id  - Disparity image (map), u x v pixels, greyscale.

% Hints:
%
%  - Loop over each image row, computing the local similarity measure, then
%    aggregate. At the border, you may replicate edge pixels, or just avoid
%    using values outside of the image.
%
%  - You may hard-code any parameters you require (e.g., disparity range) in
%    this function.
%
%  - Use whatever window size you think might be suitable.
%
%  - Don't optimize for runtime, optimize for clarity.

%convert images to doubles first (uint8 causes known bug)
Il = double(Il);
Ir = double(Ir);

%set max disparity and window size
dmax = 63;
winsize = 5;

%create a filter that sums all points within window when applied using convolution
convsum = ones(2*winsize+1, 2*winsize+1);

%housekeeping
sizeL = size(Il);

%iterate through disparities, starting with 0
for d = 0:dmax
    
    if d ~= 0  
        %keep Ir the same size while shifting it over every iteration by dropping a column, and padding with zeros
        Ir = [zeros(sizeL(1),1), Ir(:, 1:end-1)];
    end
    
    %compute absolute differences between images
    AD = abs(Il - Ir);
    %compute SAD using convolution and filter
    kindaSAD = conv2(convsum, AD);
    %trim the excess values- s.t. the sad array is the same size as the images
    SAD = kindaSAD(winsize+1:end-winsize, winsize+1:end-winsize);
    
    %housekeeping
    sadsize = size(SAD);
    
    if d == 0
        %initialize array with smallest sad values and disparity map
        mostSAD = SAD;
        Id = zeros(sadsize(1), sadsize(2));
    else 
        for y = 1:sadsize(1)
            for x = 1:sadsize(2)
                if SAD(y,x) < mostSAD(y,x)
                    %iterate through the SAD for each pixel, and compare it with the lowest SAD for that pixel so far
                    %if a lower SAD is reached, update the lowest SAD and the disparity map
                    mostSAD(y,x) = SAD(y,x);
                    Id(y,x) = d;    
                end
            end
        end
    end
    
end
%crop top (and bottom, and sides)
Id = uint8(Id(bbox(2,1):bbox(2,2), bbox(1,1):bbox(1,2)));
imshow(Id)


end

