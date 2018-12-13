function [Id] = stereo_disparity_best(Il, Ir, bbox)
% STEREO_DISPARITY_BEST Alternative stereo correspondence algorithm.
%
%  Id = STEREO_DISPARITY_BEST(Il, Ir, bbox) computes a stereo disparity image 
%  from left stereo image Il and right stereo image Ir.
%
%  Inputs:
%  -------
%   Il    - Left stereo image, m x n pixels, colour or greyscale.
%   Ir    - Right stereo image, m x n pixels, colour or greyscale.
%   bbox  - Bounding box, relative to left image, top left corner, bottom
%           right corner (inclusive). Width is v.
%
%  Outputs:
%  --------
%   Id  - Disparity image (map), m x v pixels, greyscale.

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

%Simple solution: apply a sharpening filter to both images and use SAD again. *!PRESTO!*
filter = [0 -1 0; -1 5 -1; 0 -1 0]

Il = conv2(filter, Il);
Ir = conv2(filter, Ir);

%the code below is the same as part 1, hence the lack of comments
Il = double(Il);
Ir = double(Ir);

dmax = 63;
winsize = 5;
convsum = ones(2*winsize+1, 2*winsize+1);
sizeL = size(Il);

for d = 0:dmax
    if d ~= 0  
        Ir = [zeros(sizeL(1),1), Ir(:, 1:end-1)];
    end
    AD = abs(Il - Ir);
    kindaSAD = conv2(convsum, AD);
    SAD = kindaSAD(winsize+1:end-winsize, winsize+1:end-winsize);
    sadsize = size(SAD);
    if d == 0
        mostSAD = SAD;
        Id = zeros(sadsize(1), sadsize(2));
    else 
        for y = 1:sadsize(1)
            for x = 1:sadsize(2)
                if SAD(y,x) < mostSAD(y,x)
                    mostSAD(y,x) = SAD(y,x);
                    Id(y,x) = d;
                    
                end
            end
        end
    end
    
end
Id = uint8(Id(bbox(2,1):bbox(2,2), bbox(1,1):bbox(1,2)));
imshow(Id)
  
end
