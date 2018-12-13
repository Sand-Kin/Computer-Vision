function [J] = find_jacobian(K, Twc, Wpt)
%  FIND_JACOBIAN Determine Jacobian for NLS camera pose optimization.
%
%   [J] = FIND_JACOBIAN(K, Twc, Wpt) computes the Jacobian of an image 
%   plane point with respect to the current camera pose estimate, given
%   a landmark point. The projection model is the simple pinhole model.
%
%   Inputs:
%   -------
%    K    - 3x3 camera intrinsic calibration matrix.
%    Twc  - 4x4 homogenous pose matrix, current guess for pose in world frame.
%    Wpt  - 3x1 world point on calibration target (one of n).
%
%   Outputs:
%   --------
%    J  - 2x6 Jacobian matrix (columns are tx, ty, tz, r, p, q).

%--- FILL ME IN ---

%extract rotation matrix and translation from homogeneous matrix
R = Twc(1:3, 1:3);
t = Twc(1:3, 4);

%compute world point with respect to the camera
x = Wpt-t;
xCAM = R'*x;

%compute partial derivatives with respect to r, p, q(? should be y? I digress)
[dRdr, dRdp, dRdy] = dcm_jacob_rpy(R);
dr = dRdr'*x;
dp = dRdp'*x;
dy = dRdy'*x;

%fill in jacobian using quotient rule
for i = 1:2
    J(i,1) = -K(i,i) *(R(1,i)/xCAM(3) - xCAM(i)/ xCAM(3)^2 * R(1,3));
    J(i,2) = -K(i,i) *(R(2,i)/xCAM(3) - xCAM(i)/ xCAM(3)^2* R(2,3));
    J(i,3) = -K(i,i) *(R(3,i)/xCAM(3) - xCAM(i)/ xCAM(3)^2* R(3,3));
    J(i,4) = K(i,i) * (dr(i)/xCAM(3) - xCAM(i)/ xCAM(3)^2 *dr(3));
    J(i,5) = K(i,i) * (dp(i)/xCAM(3) - xCAM(i)/ xCAM(3)^2 *dp(3));
    J(i,6) = K(i,i) * (dy(i)/xCAM(3) - xCAM(i)/ xCAM(3)^2 *dy(3));
end

end

