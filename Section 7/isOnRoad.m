function [u] = isOnRoad(x,y)
% Input:    vectors with x and y positions
%
% Output:   a vector u such that u(i) = 1 if (x(i),y(i)) is on the road
%           and 0 otherwise. 
%


%   Make sure that x and y are column vectors
n   =   length(x);      
x = reshape(x,n,1); 
y = reshape(y,n,1);

%   The number of buildings (including two rectangles in the middle)
m = 9;             

%   To check if any vector is in any building we create
%   matrices of size n x m:
X = x*ones(1,m);
Y = y*ones(1,m);

%   We should check that we are on the map
bounds = ([1+1i 1+9*1i 11+9*1i 11+1i]);

%   And that we are not in any of these houses
house = zeros(m,5);
house(1,:) = ([2+5.2*1i 2+8.3*1i 4+8.3*1i 4+5.2*1i 2+5.2*1i]);%House 1
house(2,:) = ([2+3.7*1i 2+4.4*1i 4+4.4*1i 4+3.7*1i 2+3.7*1i]);%House 2
house(3,:) = ([2+2*1i 2+3.2*1i 4+3.2*1i 4+2*1i 2+2*1i]);%House 3
house(4,:) = ([5+1i 5+2.2*1i 7+2.2*1i 7+1i 5+1i]);%House 4
house(5,:) = ([5+2.8*1i 5+5.5*1i 7+5.5*1i 7+2.8*1i 5+2.8*1i]);%House 5
house(6,:) = ([5+6.2*1i 5+9*1i 7+9*1i 7+6.2*1i 5+6.2*1i]);%House 6
house(7,:) = ([8+4.6*1i 8+8.4*1i 10+8.4*1i 10+4.6*1i 8+4.6*1i]);%House 7
house(8,:) = ([8+2.4*1i 8+4*1i 10+4*1i 10+2.4*1i 8+2.4*1i]);%House 8
house(9,:) = ([8+1.7*1i 8+1.8*1i 10+1.8*1i 10+1.7*1i 8+1.7*1i]);%House 9

%   Let us check if we are in any of the houses:
X1 = X >= ones(n,1)*real(house(:,1))';
X2 = X <= ones(n,1)*real(house(:,3))';
Y1 = Y >= ones(n,1)*imag(house(:,1))';
Y2 = Y <= ones(n,1)*imag(house(:,2))';
XX = X1 & X2;               % Finds houses that match the x-vector
YY = Y1 & Y2;               % Finds houses that match the y-vector
UU = XX & YY;               % Finds houses that match both x and y
u1 = ~any(UU, 2); % Sets u(i)=0 if (x(i),y(i)) is in a house

%   We should also make sure that the vectors are in the map
x3 = x > ones(n,1)*real(bounds(1))';
x4 = x < ones(n,1)*real(bounds(3))';
y3 = y > ones(n,1)*imag(bounds(1))';
y4 = y < ones(n,1)*imag(bounds(2))';

xx = x3 & x4;        %   Checks that the x-coordinates are in the map
yy = y3 & y4;        %   and that the y-coordinates are in the map
u2 = xx & yy;        %   Both must be inside

% Finally, we set the output to zero if (x,y) is either in a building
% or outside the map:
u = u1 & u2;


