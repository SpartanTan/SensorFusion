function [x, P] = nonLinKFprediction(x, P, f, Q, type)
%NONLINKFPREDICTION calculates mean and covariance of predicted state
%   density using a non-linear Gaussian model.
%
%Input:
%   x           [n x 1] Prior mean
%   P           [n x n] Prior covariance
%   f           Motion model function handle
%               [fx,Fx]=f(x) 
%               Takes as input x (state), 
%               Returns fx and Fx, motion model and Jacobian evaluated at x
%               All other model parameters, such as sample time T,
%               must be included in the function
%   Q           [n x n] Process noise covariance
%   type        String that specifies the type of non-linear filter
%
%Output:
%   x           [n x 1] predicted state mean
%   P           [n x n] predicted state covariance
%
    n=length(x);
    switch type
        case 'EKF'
            [fx,Fx] = f(x);
            x = fx; 
            P = Fx*P*Fx'+Q;
            
        case 'UKF'
            [SP, W] = sigmaPoints(x, P, type); % SP: [n x 2n+1] W:[1 x 2n+1]
            % mean
            x = zeros(n ,1);
            for i = 1:2*n+1
                x = x + SP(:,i) * W(i);
            end
            % covariance
            P = zeros(n,n);
            for i = 1:1:2*n+1
                P = P + (f(SP(:,i))-x)*(f(SP(:,i))-x)'*W(i);
            end
            P = Q+P;
            
        case 'CKF'
            [SP, W] = sigmaPoints(x, P, type);
            % mean
            x = zeros(n ,1);
            for i = 1:2*n
                x = x + SP(:,i) * W(i);
            end
            % covariance
            P = zeros(n,n);
            for i = 1:1:2*n
                P = P + (f(SP(:,i))-x)*(f(SP(:,i))-x)'*W(i);
            end
            P = Q+P;
            
        otherwise
            error('Incorrect type of non-linear Kalman filter')
    end

end