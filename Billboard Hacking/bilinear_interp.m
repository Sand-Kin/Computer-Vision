function [b] = bilinear_interp(I, pt)
% BILINEAR_INTERP Performs bilinear interpolation for a given image point.
%
%   Given the (x, y) location of a point in an input image, use the surrounding
%   4 pixels to conmpute the bilinearly-interpolated output pixel intensity.
%
%   Note that images are (usually) integer-valued functions (in 2D), therefore
%   the intensity value you return must be an integer (use round()).
%
%   This function is for a *single* image band only - for RGB images, you will 
%   need to call the function once for each colour channel.
%
%  Inputs:
%  -------
%   I   - Single-band (greyscale) intensity image, 8-bit (i.e., uint8).
%   pt  - 2x1 point in input image (x, y), with subpixel precision.
%
%  Outputs
%  -------
%   b  - Interpolated brightness or intensity value (whole number >= 0).

%--- FILL ME IN ---

%alternate algorithm from https://en.wikipedia.org/wiki/Bilinear_interpolation
%Basically, approximate a function f(x,y) from the surrounding points to determine intensity
%f(x,y) = a0 + a1*x + a2*y + a3*x*y
%AX=B where f(x,y) is B, A is points around desired subpixel, X is weights a0-a3.

%find closest existing points to the given one
x = pt(1);
y = pt(2);
x1 = floor(x);
x2 = ceil(x);
y1 = floor(y);
y2 = ceil(y);

%Create array A with points around desired subpixel
A = [1 , x1, y1, x1*y1;
     1 , x1, y2, x1*y2;
     1 , x2, y1, x2*y1;
     1 , x2, y2, x2*y2]
 
%create array B with f(x,y) values 
B = double([I(y1,x1); I(y2,x1); I(y1,x2); I(y2,x2)])

%Solve AX = B for X to get weights
X = A\B;

%desired point in proper format
point = [1, pt(1), pt(2), pt(1)*pt(2)];

%multiply point by weights for bilinear interpolation
b = round(point*X)

%------------------

end

