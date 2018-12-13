function [pt] = saddle_point(I)
% SADDLE_POINT Locate saddle point in an image patch.
%
%   [pt] = SADDLE_POINT(I) finds the subpixel center of a cross-junction in the
%   image patch I, by fitting a hyperbolic paraboloid to it, and then finding the 
%   critical point of that paraboloid.
%
%   Note that the location of 'p' is relative to (0.5, 0.5) at the upper left
%   corner of the patch, i.e., the pixels are treated as covering an area of
%   one unit square.
%
%   Inputs:
%   -------
%    I  - mxn image patch (grayscale, double or integer class).
%
%   Outputs:
%   --------
%    pt  - 2x1 subpixel location of saddle point in I (x, y coords).

%--- FILL ME IN ---


%PSI is I

%just need to fit the function alphax^2 etc... with linear least squares.. how?

im_size = size(I);
J_square = zeros(6,6);
J_psi = zeros(6,1);

for x = 1:im_size(1)
    for y = 1:im_size(2)

        J = [x^2 x*y y^2 x y 1];
        
        J_square = J_square + (J'*J);
        J_psi = J_psi + (J'*double(I(y,x)));

    
    end
end

invJ = inv(J_square);

theta = invJ * J_psi;

pt = -inv([2*theta(1) theta(2); theta(2) 2*theta(3)]) * [theta(4); theta(5)];



%------------------

end

