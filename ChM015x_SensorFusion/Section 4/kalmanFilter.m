function [X, P] = kalmanFilter(Y, x_0, P_0, A, Q, H, R)
%KALMANFILTER Filters measurements sequence Y using a Kalman filter. 
%
%Input:
%   Y           [m x N] Measurement sequence
%   x_0         [n x 1] Prior mean
%   P_0         [n x n] Prior covariance
%   A           [n x n] State transition matrix
%   Q           [n x n] Process noise covariance
%   H           [m x n] Measurement model matrix
%   R           [m x m] Measurement noise covariance
%
%Output:
%   x           [n x N] Estimated state vector sequence
%   P           [n x n x N] Filter error convariance
%

%% Parameters
N = size(Y,2);

n = length(x_0);
m = size(Y,1);

%% Data allocation
x = zeros(n,N);
P = zeros(n,n,N);

x(:,1) = x_0;
P(:,:,1) = P_0;
for i=2:N+1
    
    %% Prediction
    x(:,i) = A*x(:,i-1);
    P(:,:,i) = A*P(:,:,i-1)*A'+Q;

    %% Update
    S_k = H*P(:,:,i)*H'+R;
    K_k = P(:,:,i)*H'*inv(S_k);
    v_k = Y(:,i-1) - H*x(:,i);

    x(:,i) = x(:,i) + K_k*v_k;
    P(:,:,i) = P(:,:,i) - K_k*S_k*K_k';
end
X=x(:,2:N+1);
P=P(:,:,2:N+1);
end

