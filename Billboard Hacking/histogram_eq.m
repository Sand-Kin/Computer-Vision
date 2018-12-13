function [J] = histogram_eq(I)
% HISTOGRAM_EQ Histogram equalization for greyscale image.
%
%   Perform histogram equalization on the 8-bit greyscale intensity image I
%   to produce a contrast-enhanced image J. Full details of the algorithm are
%   provided in the Szeliski text.
%
%   Inputs:
%   -------
%    I  - Single-band (greyscale) intensity image, 8-bit (i.e., uint8).
%
%   Outputs:
%   --------
%    J  - Contrast-enhanced greyscale intensity image, 8-bit (i.e., uint8).

%--- FILL ME IN ---

%necessary variables. N is number of pixels.
sizeI = size(I);
N = sizeI(1)*sizeI(2);

%creating histogram function, h(x)
%(x is the intensity value, element of [0,255])
%(h(x) is the number of pixels with x intensity value)
%NOTE: due to matlab indexing, x values actually range from [1,256]. That's why you'll see I(x,y) + 1.
h = zeros(256, 1);
for x = 1:sizeI(1)
    for y = 1:sizeI(2)
        h(I(x,y)+1) = h(I(x,y)+1) + 1;
    end
end

%creating cumulative distribution function using the formula on page 108 of szeliski 
sumh = 0;
c = zeros(256, 1);
for i = 1:size(h)
    sumh = sumh + h(i);
    c(i) = sumh / N;
end

%assigning intensities to J based on the cumulative distribution
%(basically the inverse of creating a histogram)
for x = 1:sizeI(1)
    for y = 1:sizeI(2)
        J(x,y) = c(I(x,y)+1);
    end
end

%scale J so that its value are in the range [0,255]
J = uint8(J.*255) %J isn't scaled properly


%------------------

end

