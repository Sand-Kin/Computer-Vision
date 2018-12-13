function [Ipts] = cross_junctions(I, boundPoly, Wpts)
% CROSS_JUNCTIONS Find cross-junctions in image with subpixel accuracy.
%
%   [Ipts] = CROSS_JUNCTION(I, boundPoly, Wpts) locates a series of cross- 
%   junction points on a planar calibration target, where the target is
%   bounded in the image by the specified 4-sided polygon. The number of
%   cross-junctions identified should be equal to the number of world points.
%
%   Note also that the world and image points must be in *correspondence*,
%   that is, the first world point should map to the first image point, etc.
%
%   Inputs:
%   -------
%    I          - Image (grayscale, double or integer class).
%    boundPoly  - 2x4 array defining bounding polygon for target (clockwise).
%    Wpts       - 3xn array of world points (in 3D, on calibration target).
%
%   Outputs:
%   --------
%    Ipts  - 2xn array of cross-junctions (x, y), relative to the upper left
%            corner of I.

%--- FILL ME IN ---

%compute homography between normal 640x480 rectange and bounding box
sizeI = size(I);
x = [1, sizeI(2), sizeI(2), 1;
      1, 1, sizeI(1), sizeI(1)];
  
winpix = 15;
  
[H, A] = dlt_homography(x, boundPoly);
newI = I;

%make new image, which is the calibration target, but right side up (not turned in any way
%this is done using project 1 homography techniques

for y = 1:sizeI(1)
    for x = 1:sizeI(2)
           
        orig = [x;y;1];
        warp = H*orig;
            
        %normalize x' (warp)
        warp = warp/warp(3);

        newI(y,x) = bilinear_interp(I, [warp(1);warp(2)]);   
    end 
end

%compute harris points for our new, unturned calibration target
sigma = 2;
threshhold = 3200;
sizeI;
harrispts = harris(newI, 2, 1000, true);
sizeH = size(harrispts); %harris pts are x, y)

Xpts = [0;0];
sizeX = size(Xpts);


for pt = 1:sizeH(2)
    
    %I played with the harris point values such that only corners corresponding
    %to the X corners show up
    %since we have a rectangular image of the calibration target, we're basically 
    %zooming in, and throwing out any points along the edges, so that the X junctions remain
    if (harrispts(1,pt) < 60) | (harrispts(1,pt) > sizeI(2)-60)
        continue
    elseif (harrispts(2,pt) < 40) | (harrispts(2,pt) > sizeI(1)-40)
        continue
    else
        counter = 0;
        %here we go through all of the harrispts we have, and if none are close to it
        %we add it to the list of points. We don't need the best value here, since we'll use a patch around this point
        %This gives us a distinct, 48 harris points, one for each cross junction
        for i = 1:sizeX(2)
            d = norm(Xpts(:,i) - harrispts(:,pt));
            if d < 15
                counter = 1;
                break
            end
        end
        if counter == 0
            Xpts = [Xpts, harrispts(:,pt)];
        end
        sizeX = size(Xpts);     
    end
end

%here we are reordering the points such that the final points are in row major order
Xpts = Xpts(:,2:sizeX(2))';
Xpts = [Xpts(:,2), Xpts(:,1)];
reX = sortrows(Xpts)
reX = [reX(:,2), reX(:,1)];
reX = [sortrows(reX(1:8,:));
       sortrows(reX(9:16,:));
       sortrows(reX(17:24,:));
       sortrows(reX(25:32,:));
       sortrows(reX(33:40,:));
       sortrows(reX(41:48,:))];
reX = [reX(:,2), reX(:,1)]

sizeREX = size(reX);

%finally, we use the homography again to find the approximate cross junction point in the original image,
%round it to the nearest pixel (doesn't need to be exact), and pass a patch around that pixel to the saddle point function
%the output of the saddlepoint function is then normalized and added to the Ipts array.
for i = 1:sizeREX(1)
    orig = [reX(i,2); reX(i,1); 1];
    warp = H*orig;
    warp = warp/warp(3);
    x = round(warp(1));
    y = round(warp(2));
    LOC = saddle_point(I(y-winpix:y+winpix, x-winpix:x+winpix));
    XLOC = [x-winpix+LOC(1); y-winpix+LOC(2)];
    if i == 1
        Ipts = XLOC;
    else
        Ipts = [Ipts, XLOC];
    end
end

end

function [b] = bilinear_interp(I, pt)
x = pt(1);
y = pt(2);
x1 = floor(x);
x2 = ceil(x);
y1 = floor(y);
y2 = ceil(y);
A = [1 , x1, y1, x1*y1;
     1 , x1, y2, x1*y2;
     1 , x2, y1, x2*y1;
     1 , x2, y2, x2*y2];
B = double([I(y1,x1); I(y2,x1); I(y1,x2); I(y2,x2)]);
X = A\B;
point = [1, pt(1), pt(2), pt(1)*pt(2)];
b = round(point*X);
end
function [H, A] = dlt_homography(I1pts, I2pts)
for i = 1:4
    x = I1pts(1,i);
    y = I1pts(2,i);
    u = I2pts(1,i);
    v = I2pts(2,i);
    Ai = [-x, -y, -1, 0, 0, 0, u*x, u*y, u;
          0, 0, 0, -x, -y, -1, v*x, v*y, v];
    if (i == 1)
        A = Ai;
    else
        A = [A;Ai]; 
    end    
end  
[U,S,V] = svd(A);
H = reshape(V(:,9),3,3)';
end
function [pt] = saddle_point(I)
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
end
