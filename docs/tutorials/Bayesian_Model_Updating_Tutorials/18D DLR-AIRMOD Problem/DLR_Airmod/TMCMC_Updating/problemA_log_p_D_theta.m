function logL = problemA_log_p_D_theta(D, theta)
% Calculation of the log_likelihood for the example in problemA.m
%
% USAGE:
% logL = problemA_log_p_D_theta(D, theta)
%
% INPUTS:
% D     = experimental observations   nobs x dim_x
% theta = epistemic parameters        npar x dim_theta
%
% OUTPUTS:
% logL(i)  = loglikelihood for the set of parameters theta(i,:) and the
%            data D, i = 1, 2, ...npar.        logL = npar x 1

%--------------------------------------------------------------------------
% who                    when         observations
%--------------------------------------------------------------------------
% Diego Andres Alvarez   Jul-24-2013  First algorithm
%--------------------------------------------------------------------------
% Diego Andres Alvarez - daalvarez@unal.edu.co

%%
npar = size(theta,1);  % number of thetas to evaluate
logL = zeros(npar,1);
for i = 1:npar
   logL(i) = sum(log(p_x_theta_pdf(D, theta(i,:))));
   if isinf(logL(i))
      logL(i) = -1e10;
   end
end

return;

%%
function p = p_x_theta_pdf(x, theta_i)
% x       = set of observations         nobs x dim_x
% theta_i = point in epistemic space       1 x dim_theta

%% Simulate M points from the aletatory space (0,1]^4
% NOTE: I am not making this because it is easier to simulate directly 
% the random variables than using the inverse transform method. Also, the
% random variables are independent, so there is no need of a copula
% in (0,1]^4

M = 50000; % number of points employed to create the ksdensity of p_x_theta

%% p1 (beta distribution)
p1_mu  = theta_i(1);
p1_var = theta_i(2);

% syms a b p1_mu p1_var
% sol = solve(p1_mu == a/(a+b), p1_var == a*b/((a+b+1)*((a+b)^2)), a, b);
% factor(sol.a), factor(sol.b)
a = -(p1_mu*(p1_var - p1_mu + p1_mu^2))/p1_var;
b = ((p1_mu - 1)*(p1_var - p1_mu + p1_mu^2))/p1_var;

if not((a>1) && (b>1)) % Unimodality
   p = 0; %zeros(size(x));
   return;
end

p1 = betarnd(a,b,M,1);

%% p2 (interval)
p2 = repmat(theta_i(3),M,1);

%% p3 (uniform random variable)
p3 = rand(M,1); % = unifrnd(0,1,M,1);

%% p4, p5 (bidimensional gaussian distribution)
p4_mu = theta_i(4);  p4_var = theta_i(5);
p5_mu = theta_i(6);  p5_var = theta_i(7);  p4p5_rho = theta_i(8);
c12  = p4p5_rho*sqrt(p4_var*p5_var);
covmat = [p4_var c12; c12 p5_var];
%{
[~, pp] = chol(covmat);
if pp ~= 0
   keyboard
   error('c12 should de symmetric positive semi-definite');
end
%}
p4p5 = mvnrnd([p4_mu p5_mu], covmat, M);

%% Map these simulations through the system
savedState = rng;                   % grab the stream state
output  = p_to_x1([p1 p2 p3 p4p5]);
rng(savedState);

%% Estimate the PDF p_x_theta_pdf(x | theta)
% f = ksdensity(x,xi) specifies the vector xi of values, where the density 
% estimate of the data in x is to be evaluated
p = ksdensity(output, x);  % p(i) = p_x_theta_pdf(x(i,:) | theta)

return;
