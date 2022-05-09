function [xfp, Pfp, Xp, Wp] = pfFilter(x_0, P_0, Y, proc_f, proc_Q, meas_h, meas_R, ...
                             N, bResample, plotFunc)
%PFFILTER Filters measurements Y using the SIS or SIR algorithms and a
% state-space model.
%
% Input:
%   x_0         [n x 1] Prior mean
%   P_0         [n x n] Prior covariance
%   Y           [m x K] Measurement sequence to be filtered
%   proc_f      Handle for process function f(x_k-1)
%   proc_Q      [n x n] process noise covariance
%   meas_h      Handle for measurement model function h(x_k)
%   meas_R      [m x m] measurement noise covariance
%   N           Number of particles
%   bResample   boolean false - no resampling, true - resampling
%   plotFunc    Handle for plot function that is called when a filter
%               recursion has finished.
% Output:
%   xfp         [n x K] Posterior means of particle filter
%   Pfp         [n x n x K] Posterior error covariances of particle filter
%   Xp          [n x N x K] Particles for posterior state distribution in times 1:K
%   Wp          [N x K] Non-resampled weights for posterior state x in times 1:K

% Your code here, please. 
% If you want to be a bit fancy, then only store and output the particles if the function
% is called with more than 2 output arguments.
n = size(x_0,1);
K = size(Y,2);

% allocate memory
xfp = zeros(n,K);
Pfp = zeros(n,n,K);
Xp = zeros(n,N,K);
Wp = zeros(N,K);

% sample around the prior mean
if size(x_0,2) == 1
        Xp(:,:,1) = mvnrnd(x_0,P_0,N)';
    else
        Xp(:,:,1) = x_0;
    end

Wp(:,1) = 1/N * ones(N,1);
j = 1:N; % index

for k = 2:K+1
    X_kmin1 = Xp(:,:,k-1); % start from k = 1
    W_kmin1 = Wp(:,k-1)';
    % if resample
    if bResample
        [X_kmin1, W_kmin1, j] = resampl(X_kmin1, W_kmin1);
    end
    % perform a particle filter step with the measurement
    [Xp(:,:,k), Wp(:,k)] = pfFilterStep(X_kmin1, W_kmin1, Y(:,k-1), proc_f, proc_Q, meas_h, meas_R);
    
    % plot particles with function handle
    if ~isempty(plotFunc)
        plotFunc(k-1, Xp(:,:,k), Xp(:,:,k-1), Wp(:,k)', j); 
    end
    % estimate mean and covariance given the particles
    xfp(:,k) = sum(Xp(:,:,k).*Wp(:,k)',2); % sum all the particles in time k
    Pfp(:,:,k) = Wp(:,k)'.*(Xp(:,:,k) - xfp(:,k))*(Xp(:,:,k) - xfp(:,k))'; % not know how to get this
    
end
% remove prior from vector
    xfp = xfp(:,2:end);
    Pfp = Pfp(:,:,2:end);
    Xp  = Xp(:,:,2:end);
    Wp  = Wp(:,2:end);



end

function [Xr, Wr, j] = resampl(X, W)
%RESAMPLE Resample particles and output new particles and weights.
% resampled particles. 
%
%   if old particle vector is x, new particles x_new is computed as x(:,j)
%
% Input:
%   X   [n x N] Particles, each column is a particle.
%   W   [1 x N] Weights, corresponding to the samples
%
% Output:
%   Xr  [n x N] Resampled particles, each corresponding to some particle 
%               from old weights.
%   Wr  [1 x N] New weights for the resampled particles.
%   j   [1 x N] vector of indices refering to vector of old particles

% Your code here!

    N = size(X,2);

    segment = [0 cumsum(W)/sum(W)];
    samples = rand([1 N]);
    j = zeros(1,N);

    for i = 1:N
        j(i) = find(segment <= samples(i), 1, 'last');
    end
    Wr = 1/N*ones(1,N);
    Xr = X(:,j);

end


function [X_k, W_k] = pfFilterStep(X_kmin1, W_kmin1, yk, proc_f, proc_Q, meas_h, meas_R)
% Input:
%   X_kmin1     [n x N] Particles for state x in time k-1
%   W_kmin1     [1 x N] Weights for state x in time k-1
%   y_k         [m x 1] Measurement vector for time k
%   proc_f      Handle for process function f(x_k-1)
%   proc_Q      [n x n] process noise covariance
%   meas_h      Handle for measurement model function h(x_k)
%   meas_R      [m x m] measurement noise covariance
%
% Output:
%   X_k         [n x N] Particles for state x in time k
%   W_k         [1 x N] Weights for state x in time k
    % Copy your code from previous task!
    n = size(X_kmin1,1);
    N = size(X_kmin1,2);

    X_k = zeros(n,N);
    W_k = zeros(1,N);

    % draw N samples p(x_k|x_{k-1})
    X_k = mvnrnd(proc_f(X_kmin1)',proc_Q)';
    likelihood = mvnpdf(yk',meas_h(X_k)',meas_R)'; % 1*m
    W_k = W_kmin1.*likelihood;
    W_k = W_k/sum(W_k);
end