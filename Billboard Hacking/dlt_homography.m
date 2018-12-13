function [H, A] = dlt_homography(I1pts, I2pts)
% DLT_HOMOGRAPHY Perspective Homography between two images.
%
%   Given 4 points from 2 separate images, compute the perspective homography
%   (warp) between these points using the DLT algorithm.
%
%   Inputs:
%   -------
%    I1pts  - 2x4 array of points from Image 1 (each column is x, y).
%    I2pts  - 2x4 array of points from Image 2 (1-to-1 correspondence).
%
%   Outputs:
%   --------
%    H  - 3x3 perspective homography (matrix map) between image coordinates.
%    A  - 8x9 DLT matrix used to determine homography.


%Step 1: Construct A matrix
for i = 1:4
    %Construct Ai from given points
    x = I1pts(1,i);
    y = I1pts(2,i);
    u = I2pts(1,i);
    v = I2pts(2,i);
    Ai = [-x, -y, -1, 0, 0, 0, u*x, u*y, u;
          0, 0, 0, -x, -y, -1, v*x, v*y, v];
       
    %Stack Ai on top of eachother to get A
    if (i == 1)
        A = Ai;
    else
        A = [A;Ai]; 
    end    
end  

%Step 2: Use SVD to solve for H. (Ah = 0)
[U,S,V] = svd(A);
%Last singular vector of V is the solution for H
H = reshape(V(:,9),3,3)';

%------------------

end
