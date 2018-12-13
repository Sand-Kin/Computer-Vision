function [Twc] = pose_estimate_nls(K, Twcg, Ipts, Wpts)
%  POSE_ESTIMATE_NLS Estimate camera pose from 2D-3D correspondences via NLS.
%
%   [Twc] = POSE_ESTIMATE_NLS(Twcg, Ipts, Wpts) performs a nonlinear least squares 
%   optimization procedure to determine the best estimate of the camera pose in 
%   the calibration target frame, given 2D-3D point correspondences.
%
%   Inputs:
%   -------
%    K     - 3x3 camera intrinsic calibration matrix.
%    Twcg  - 4x4 homogenous pose matrix, initial guess for camera pose.
%    Ipts  - 2xn array of cross-junction points (with subpixel accuracy).
%    Wpts  - 3xn array of world points (one-to-one correspondence with Ipts).
%
%   Outputs:
%   --------
%    Twc  - 4x4 homogenous pose matrix, estimate of camera pose in target frame.

% 1. Set initial parameter vector and maximum iterations. Use rpy_from_dcm.
R = Twcg(1:3,1:3);
t = Twcg(1:3,4);
rpy = rpy_from_dcm(R);
param = [t; rpy];
maxIters = 200;


% 2. Allocate space for residuals vector and full Jacobian matrix.
sizeW = size(Wpts);
N = sizeW(2);
xPTS = reshape(Ipts,2*N,1);
WHOMO = [Wpts; ones(1, N)];

% 3. Loop until convergence.
iter = 1;

while true

  % 4. Project each known landmark point into camera image, given 
  %    current pose estimate.
  Tcwg = inv(Twcg);
  xCAM = Tcwg*WHOMO;
  xHOMO = K*xCAM(1:3,:);
  xNEW = xHOMO(1:2,:)./xHOMO(3,:);
  xGUESS = reshape(xNEW,2*N,1);
  
  % 5. Solve system of normal equations for this iteration.
  %compute jacobian
  for i = 1:N
      J(2*i-1:2*i,:) = find_jacobian(K, Twcg, Wpts(:,i));
  end
  
  % 6. Check - converged? Use variable 'diff'.
  %compute difference using residuals
  res = xPTS - xGUESS;
  diff = inv(J'*J)*J'*res;
  if norm(diff) < 1e-12
    break;
  elseif iter == maxIters
    break;
  end
  
  %update best guess for parameter vector for next iteration
  param = param+diff;
  wcg(1:3,1:3) = dcm_from_rpy(param(4:6));
  Twcg(1:3,4) = param(1:3);
  
  
  iter = iter + 1;
end
%update Twc to best guess found (either by hitting max iterations
%or by being under e-12 difference
Twc = Twcg
end

%---------- Functions Go Below ----------

% NOTE: THE CODE BELOW IS UNCOMMENTED FOR SPACE. 
% If you want to see the comments, check the respective sections (i.e. part 1 has saddlepoints commented)
function [J] = find_jacobian(K, Twc, Wpt)
R = Twc(1:3, 1:3);
t = Twc(1:3, 4);
x = Wpt-t;
xCAM = R'*x;
[dRdr, dRdp, dRdy] = dcm_jacob_rpy(R);
dr = dRdr'*x;
dp = dRdp'*x;
dy = dRdy'*x;

Tcw = inv(Twc);
rand1 = Tcw(3,2);
for i = 1:2
    rand2 = Tcw(i,2);
    J(i,1) = -K(i,i) *(Tcw(i,1)/xCAM(3) - xCAM(i)/ xCAM(3)^2* Tcw(3,1));
    J(i,2) = -K(i,i) *(rand2/xCAM(3) - xCAM(i)/ xCAM(3)^2* rand1);
    J(i,3) = -K(i,i) *(Tcw(i,1)/xCAM(3) - xCAM(i)/ xCAM(3)^2* Tcw(3,3));
    J(i,4) = K(i,i) * (dr(i)/xCAM(3) - xCAM(i)/ xCAM(3)^2 *dr(3));
    J(i,5) = K(i,i) * (dp(i)/xCAM(3) - xCAM(i)/ xCAM(3)^2 *dp(3));
    J(i,6) = K(i,i) * (dy(i)/xCAM(3) - xCAM(i)/ xCAM(3)^2 *dy(3));
end
end
function [Ipts] = cross_junctions(I, boundPoly, Wpts) 
sizeI = size(I);
x = [1, sizeI(2), sizeI(2), 1;
      1, 1, sizeI(1), sizeI(1)];
winpix = 15;
[H, A] = dlt_homography(x, boundPoly);
newI = I;
for y = 1:sizeI(1)
    for x = 1:sizeI(2)
        orig = [x;y;1];
        warp = H*orig;
        warp = warp/warp(3);
        newI(y,x) = bilinear_interp(I, [warp(1);warp(2)]);   
    end 
end
sigma = 2;
threshhold = 3200;
sizeI;
harrispts = harris(newI, 2, 1000, true);
sizeH = size(harrispts); %harris pts are x, y)
Xpts = [0;0];
sizeX = size(Xpts);
for pt = 1:sizeH(2)
    if (harrispts(1,pt) < 60) | (harrispts(1,pt) > sizeI(2)-60)
        continue
    elseif (harrispts(2,pt) < 40) | (harrispts(2,pt) > sizeI(1)-40)
        continue
    else
        counter = 0;
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


