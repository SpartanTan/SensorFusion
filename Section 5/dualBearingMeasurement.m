function [hx, Hx] = dualBearingMeasurement(x, s1, s2)
%DUOBEARINGMEASUREMENT calculates the bearings from two sensors, located in 
%s1 and s2, to the position given by the state vector x. Also returns the
%Jacobian of the model at x.
%
%Input:
%   x           [n x 1] State vector, the two first element are 2D position
%   s1          [2 x 1] Sensor position (2D) for sensor 1
%   s2          [2 x 1] Sensor position (2D) for sensor 2
%
%Output:
%   hx          [2 x 1] measurement vector
%   Hx          [2 x n] measurement model Jacobian
%
% NOTE: the measurement model assumes that in the state vector x, the first
% two states are X-position and Y-position.

% Your code here

n = length(x);
xk = x(1);
yk = x(2);

s1_x=s1(1);
s1_y=s1(2);
s2_x=s2(1);
s2_y=s2(2);

atany1 = yk-s1_y;
atanx1 = xk-s1_x;
atany2 = yk-s2_y;
atanx2 = xk-s2_x;

hx = [atan2(atany1, atanx1);
      atan2(atany2, atanx2)];
Hx = [ -atany1/(atanx1^2+atany1^2) atanx1/(atanx1^2+atany1^2) zeros(1,n-2);
       -atany2/(atanx2^2+atany2^2) atanx2/(atanx2^2+atany2^2) zeros(1,n-2)];


end