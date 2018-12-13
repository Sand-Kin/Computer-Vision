% Billboard hack script file.

% Bounding box in Y & D Square image.
bbox = [404, 490, 404, 490;
         38,  38, 354, 354];

% Point correspondences.
Iyd_pts = [416, 485, 488, 410;
            40,  61, 353, 349];
Ist_pts = [2, 219, 219,   2; 
           2,   2, 410, 410];

Iyd = imread('yonge_dundas_square.jpg');
Ist = imread('uoft_soldiers_tower_dark.png');

% Make a copy for hacking... your script MUST use the Ihack variable name.
Ihack = Iyd; 

% Let's do the histogram equalization first...
hist = histogram_eq(Ist);

% Compute the perspective homography we need...

% In this case, x' in (x' = Hx) is the original image since we're inverse warping
% hence, the order of the points in the dlt_homography call
[H, A] = dlt_homography(Iyd_pts, Ist_pts);


% Main 'for' loop to do the warp and insertion

%proper format of points for inpolygon function
xv = Iyd_pts(1, :);
yv = Iyd_pts(2, :);

for y = 38:354 %for loop constraints from bbox
    for x = 404:490
        %check if point is in polygon
        if (inpolygon(x,y,xv,yv) == 1)
            
            %inverse warp to original image point
            orig = [x;y;1];
            warp = H*orig;
            
            %normalize x' (warp)
            warp(1) = warp(1)/warp(3);
            warp(2) = warp(2)/warp(3);
            
            %if x' (warp) is a pixel, just insert point, no need for bilinear interpolation
            if (warp(1) == round(warp(1))) && (warp(2) == round(warp(2)))
                Ihack(y,x, [1,2,3]) = hist(warp(2),warp(1));
                
            %if x' is a subpixel, computer bilinear interpolation of subpixel    
            else
                
                %assign intensity to the R, G, and B layer (hence the [1,2,3])
                Ihack(y,x, [1,2,3]) = bilinear_interp(hist, [warp(1);warp(2)]);
            end
        end
    end
end
%}        
% this could be vectorized to be faster if needed.

% You may wish to make use of the inpolygon() function (see MATLAB help).

figure;  imshow(Ihack);

%---------- Functions Go Below ----------

% ADD YOUR DLT_HOMOGRAPHY CODE BELOW
function [H, A] = dlt_homography(I1pts, I2pts)

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

end

% ADD YOUR BILINEAR_INTERP CODE BELOW
function [b] = bilinear_interp(I, pt)

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
     1 , x2, y2, x2*y2];
 
%create array B with f(x,y) values 
B = double([I(y1,x1); I(y2,x1); I(y1,x2); I(y2,x2)]);

%Solve AX = B for X to get weights
X = A\B;

%desired point in proper format
point = [1, pt(1), pt(2), pt(1)*pt(2)];

%multiply point by weights for bilinear interpolation
b = round(point*X);

end

% ADD YOUR HISTOGRAM_EQ CODE BELOW
function [J] = histogram_eq(I)

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
J = uint8(J.*255); %J isn't scaled properly

end
