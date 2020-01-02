% clear variables; clc; close all

load('metamodel_light','Xnn_no_data');

[~,hostname]=system('hostname'); hostname(end)=[];           % name of the computer
directory = 'workfolder';

%diary(fullfile(directory,[hostname '_AIRMOD.txt'])); % store the output of MATLAB in this file

%% What time and date is it?
datetimenow = datestr(clock);

%% Load the experimental data (as a column vector):
load('SMUdata.mat');      
D = dataCloud{1,1};
% remove failed experiments (all NaN on a line)
idxNaN = all(isnan(D),2);
D(idxNaN,:) = [];
Nexpsamples = size(D,1);
% keep only selected modes 
% selected_modes = 5;
selected_modes = [1:8,10:12,14,19,20];
D = D(:,selected_modes);

%% Defining the prior PDF p(theta)
Xrvset = Xnn_no_data(1).Xinput.Xrvset.(Xnn_no_data(1).Xinput.CnamesRandomVariableSet{1});

fixed_variance = false;

if ~fixed_variance
    std_lb = 1e-4; 
    std_ub = 10*sqrt(nanvar(D.^2))./nanmean(D.^2);
    std_ub([11,13]) = 10*std_ub([11,13]);
    
    p_theta    = @(x) prod([Xrvset.evalpdf('Mxsamples',x(:,1:Xrvset.Nrv)),...
                            unifpdf(x(:,Xrvset.Nrv+(1:length(selected_modes))),std_lb, std_ub)], 2); 

    p_thetarnd = @(N) [get(Xrvset.sample('Nsamples', N),'MsamplesPhysicalSpace'),...
                   unifrnd(std_lb, repmat(std_ub, N, 1))];
else
    p_theta = @(x)Xrvset.evalpdf('Mxsamples',x);
    p_thetarnd = @(N) get(Xrvset.sample('Nsamples', N),'MsamplesPhysicalSpace');
end
%% The loglikelihood of D given theta
log_p_D_theta = @(theta) Airmod_log_p_D_theta(D, theta, Xnn_no_data, selected_modes);

%% The Posterior Distribution function:
target_PDF = @(x) log(p_theta(x)) + log_p_D_theta(x);

%% Bayesian estimation of theta: bayesian model updating using TMCMC
tic
Nsamples = 500;
fprintf('Nsamples TMCMC = %d\n', Nsamples);
samples_ftheta_D = TMCMCsampler(log_p_D_theta, p_theta, p_thetarnd, Nsamples);
save(fullfile(directory,[hostname '_results_AIRMOD.mat']));
toc

%% Bayesian estimation of theta: bayesian model updating using MCMC
%{
% Define the Proposal distribution function for MCMC sampler:
stdev = [108.145683697584, 233.718980920495, 6.33841668350872, 12.6605610403785, 12.5094840585852, 104641.863773591,...
    253268590.022243, 0.0169267842503432, 0.00513742932278792, 0.00365231166885186, 0.000900071656680624, 0.00191637792190504,...
    0.00211445344424288, 1819331.98764385, 976536.954886673, 391813.010840492, 11720500.4290358, 11118394.6551786,...
    0.0307932426341752, 0.0458087146717363, 0.0224966678408943, 0.0132682332934532, 0.00109451068943748, 0.00277050499911286,...
    0.00132551564654006, 0.000785680092944234, 0.000686134111839630, 0.00526755402673266, 0.00485781300669912, 0.00605108290600254,...
    0.135461557612244, 0.00847536777699360];

Tuning_parameter = (10^-4).*stdev.*eye(32);

proppdf = @(x,y) mvnpdf(x,y,Tuning_parameter);
proprnd = @(x) mvnrnd(x, Tuning_parameter);

tic
Nsamples = 500;
NumberOfChain = 1;
Start = p_thetarnd(NumberOfChain);
BurnIn = 1000;
fprintf('Nsamples MCMC = %d\n', Nsamples);
[samples_ftheta_D,accept] = mhsample(Start,Nsamples,'logpdf',target_PDF,'proppdf',proppdf,...
'proprnd',proprnd,'symmetric',1, 'burnin', BurnIn, 'nchain', NumberOfChain);
save(fullfile(directory,[hostname '_results_AIRMOD.mat']));
toc

fprintf('Accept Level = %d\n', accept);
%}

%% Bayesian estimation of theta: bayesian model updating using SMC
%{
stdev = [108.145683697584, 233.718980920495, 6.33841668350872, 12.6605610403785, 12.5094840585852, 104641.863773591,...
    253268590.022243, 0.0169267842503432, 0.00513742932278792, 0.00365231166885186, 0.000900071656680624, 0.00191637792190504,...
    0.00211445344424288, 1819331.98764385, 976536.954886673, 391813.010840492, 11720500.4290358, 11118394.6551786,...
    0.0307932426341752, 0.0458087146717363, 0.0224966678408943, 0.0132682332934532, 0.00109451068943748, 0.00277050499911286,...
    0.00132551564654006, 0.000785680092944234, 0.000686134111839630, 0.00526755402673266, 0.00485781300669912, 0.00605108290600254,...
    0.135461557612244, 0.00847536777699360];

no_iterations = 1;
q_cov = (10^-2).*stdev.*eye(32);

tic
Nsamples = 500;
fprintf('Nsamples SMC = %d\n', Nsamples);

smc = SMCsampler('nsamples',Nsamples,...
    'prior_values',p_thetarnd(Nsamples),'prior',p_theta,'loglikelihood',log_p_D_theta,...
    'no_iterations',no_iterations,'prop_covariance',q_cov);
smc = smc.generate_samples();
samples_ftheta_D = smc.theta;
save(fullfile(directory,[hostname '_results_AIRMOD.mat']));
toc
%}

%% run the model with the posterior distributions of the inputs
Xinput_post = Xnn_no_data(1).Xinput.sample('Nsamples',1);
Xinput_post.Xsamples.MsamplesPhysicalSpace = samples_ftheta_D(:,1:18);
output = zeros(Xinput_post.Nsamples,length(selected_modes));
for imode = 1:length(selected_modes)
    Xout = Xnn_no_data(selected_modes(imode)).apply(Xinput_post);
    output(:,imode)  = Xout.getValues('CSnames',Xnn_no_data(selected_modes(imode)).Coutputnames(1));
end

figure;
plotmatrix(samples_ftheta_D(:,1:18));
title('Plotmatrix of samples');

figure;
plotEigenvaluesScatter(output,D);
title('Eigenvalues Scatter Plot');

diary off;
