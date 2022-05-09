function [SP,W] = sigmaPoints(x, P, type)
% SIGMAPOINTS computes sigma points, either using unscented transform or
% using cubature.
%
%Input:
%   x           [n x 1] Prior mean
%   P           [n x n] Prior covariance
%
%Output:
%   SP          [n x 2n+1] UKF, [n x 2n] CKF. Matrix with sigma points
%   W           [1 x 2n+1] UKF, [1 x 2n] CKF. Vector with sigma point weights 
%
    n=length(x);
    switch type        
        case 'UKF'
            n=length(x);
            P_sqrt=sqrtm(P);
            W_0=1-n/3;
            factor=sqrt(n/(1-W_0));
            % allocate space for SP
            SP = zeros(n, 2*n);
            % assign SP
            for i=1:n
                SP(:,i)=x+factor*P_sqrt(:,i);
                SP(:,i+n)=x-factor*P_sqrt(:,i);
            end
            SP=[x SP]; % add the center point
            W = [W_0 ones(1, 2*n)*(1-W_0)/(2*n)];
            
        case 'CKF'
            P_sqrt=sqrtm(P);
            sqrt_n=sqrt(n);
            % allocate space for SP
            SP=zeros(n,2*n);
            % assign SP
            for i=1:n
                SP(:,i)=x+sqrt_n*P_sqrt(:,i);
                SP(:,i+n)=x-sqrt_n*P_sqrt(:,i);
            end
            W_i=1/(2*n);
            W=ones(1,2*n)*W_i;
            
        otherwise
            error('Incorrect type of sigma point')
    end

end