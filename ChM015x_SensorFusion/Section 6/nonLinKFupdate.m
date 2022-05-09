function [x, P] = nonLinKFupdate(x, P, y, h, R, type)
%NONLINKFUPDATE calculates mean and covariance of predicted state
%   density using a non-linear Gaussian model.
%
%Input:
%   x           [n x 1] Predicted mean
%   P           [n x n] Predicted covariance
%   y           [m x 1] measurement vector
%   h           Measurement model function handle
%               [hx,Hx]=h(x) 
%               Takes as input x (state), 
%               Returns hx and Hx, measurement model and Jacobian evaluated at x
%               Function must include all model parameters for the particular model, 
%               such as sensor position for some models.
%   R           [m x m] Measurement noise covariance
%   type        String that specifies the type of non-linear filter
%
%Output:
%   x           [n x 1] updated state mean
%   P           [n x n] updated state covariance
%
    n = length(x);
    m = length(y);
    switch type
        case 'EKF'
           [hx, Hx] = h(x);
           S_k = Hx * P * Hx' + R;
           Kk = P * Hx' * inv(S_k);
           
           x = x + Kk * (y - hx);
           P = P - Kk * S_k * Kk';         
           
        case 'UKF'
            [SP, W] = sigmaPoints(x, P, type);
            % y_hat
            y_hat = zeros(m,1);
            for i  = 1: 2*n+1
                y_hat = y_hat + h(SP(:,i))*W(i);
            end
            
            % P_xy
            P_xy = zeros(n,m);
            for i = 1:2*n+1
                P_xy = P_xy + (SP(:,i)-x)*(h(SP(:,i))-y_hat)'*W(i);
            end
            
            % S_k
            S_k = zeros(m,m);
            for i = 1:2*n+1
                S_k = S_k + (h(SP(:,i))-y_hat)*(h(SP(:,i))-y_hat)'*W(i);
            end
            S_k = R +S_k;
            
            % updated state mean
            x = x + P_xy * inv(S_k) * (y - y_hat);
            P = P - P_xy * inv(S_k) * P_xy';
        case 'CKF'
            [SP, W] = sigmaPoints(x, P, type);
            % y_hat
            y_hat = zeros(m,1);
            for i  = 1: 2*n
                y_hat = y_hat + h(SP(:,i))*W(i);
            end
            
            % P_xy
            P_xy = zeros(n,m);
            for i = 1:2*n
                P_xy = P_xy + (SP(:,i)-x)*(h(SP(:,i))-y_hat)'*W(i);
            end
            
            % S_k
            S_k = zeros(m,m);
            for i = 1:2*n
                S_k = S_k + (h(SP(:,i))-y_hat)*(h(SP(:,i))-y_hat)'*W(i);
            end
            S_k = R +S_k;
            
            % updated state mean
            x = x + P_xy * inv(S_k) * (y - y_hat);
            P = P - P_xy * inv(S_k) * P_xy';
                
        otherwise
            error('Incorrect type of non-linear Kalman filter')
    end

end

