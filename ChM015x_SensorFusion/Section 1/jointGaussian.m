function [mu, Sigma] = jointGaussian(mu_x,mu_r,A,b,sigma2_x,sigma2_r)
%jointGaussian calculates the joint Gaussian density as defined
%in problem 1.3a. 
%
%Input
%   MU_X        Expected value of x
%   SIGMA2_X    Covariance of x
%   SIGMA2_R    Covariance of the noise r
%
%Output
%   MU          Mean of joint density 
%   SIGMA       Covariance of joint density


%Your code here
A = [1 0;A b];
mu_x = [mu_x;0];
Sigma_x = blkdiag(sigma2_x,sigma2_r);
[mu, Sigma] = affineGaussianTransform(mu_x, Sigma_x, A, mu_r);
end