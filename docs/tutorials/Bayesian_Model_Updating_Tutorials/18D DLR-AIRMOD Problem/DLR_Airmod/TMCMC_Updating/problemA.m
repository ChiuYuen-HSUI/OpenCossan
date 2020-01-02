clear all; clc; close all

%addpath('/home/daa/matlab_work/stoch_mech/challenge/');
addpath('/home/dalvarez/uq_challenge/challenge/mexa64');

Nexpsamples = 50;
hostname = getenv('HOSTNAME');                       % name of the computer
directory = ['NASA_refined_' num2str(Nexpsamples)]

diary([directory '/' hostname '_problemA_refined_NASA_' num2str(Nexpsamples) '_samples.txt']); % store the output of MATLAB in this file

%% What time and date is it?
datetimenow = datestr(clock);


%% Load the experimental data (as a column vector):
load('x1samples1');         % loads x1sams1 25x1
if Nexpsamples == 25
   D = x1sams1;
elseif Nexpsamples == 50
   load('x1samples2');      % loads x1sams2 25x1
   D = [x1sams1; x1sams2];  % all data
else
	error('The number of samples must be either 25 or 50');
end

%% VARIABLE NAME  LOWER_BOUND  UPPER_BOUND
variables = { ...
   'p1_mu'           3/5          4/5           % original information
   'p1_var'          1/50         1/25          % original information
%  'p1_mu'         3.1/5         3.2/5          % improved interval NASA
%  'p1_var'         0.8/25        0.9/25        % improved interval NASA
   'p2'               0            1
   'p4_mu'           -5            5            % original information
   'p4_var'          1/400         4            % original information
%  'p4_mu'           3.9           4.5          % improved interval NASA
%  'p4_var'        3.5/100         5/100        % improved interval NASA
   'p5_mu'           -5            5
   'p5_var'          1/400         4
   'p4p5_rho'        -1            1            % original information
%  'p4p5_rho'        -0.5         0.5           % improved interval NASA
};
%p1:  3.1/5<=E[p1]<=3.2/5 0.8/25<=v[p1]<=0.9/25
%p4:  3.9<=E[p4]<=4.5, 3.5/100<=v[p4]<=5/100, |rho|<=0.5

%% Defining the prior PDF p(theta)
lb = cell2mat(variables(:,2))';
ub = cell2mat(variables(:,3))';

p_theta    = @(x) problemA_p_theta_pdf(x, lb, ub);
p_thetarnd = @(N) problemA_p_theta_rnd(lb, ub, N);

%% The loglikelihood of D given theta
log_p_D_theta = @(theta) problemA_log_p_D_theta(D, theta);

%% Maximum likelihood estimation of theta (using genetic algorithms)
if matlabpool('size') == 0  % activate all available cores
   matlabpool open;
end;

PopulationSize = 50;
Generations    = 10;
fprintf('PopulationSize = %d, Generations = %d\n', PopulationSize, Generations);
options = gaoptimset(@ga);
options = gaoptimset(options, ...
            'PopulationSize', PopulationSize, ...
            'Generations',    Generations, ...
            'EliteCount',     5, ...
            'MutationFcn',    @mutationadaptfeasible, ...
            'Display',        'iter', ...
            'Vectorized',     'off', ...
            'UseParallel',    'always');
%           'PlotFcns',{@gaplotbestf,@gaplotmaxconstr}, ...

% minus = ga makes a minimization, and we want to maximize
minus_log_p_D_theta = @(theta) -problemA_log_p_D_theta(D, theta);

tic
[thetaML, p_thetaML, exitflag, output, final_pop, scores] = ...
            ga(minus_log_p_D_theta, ...
            8,  ...
            [],[],  ... % A, b,
            [],[],  ... % Aeq, beq,
            lb,ub,  ...
            [],     ... % nonlcon,
            options);
toc


fprintf('The maximum likelihood estimation of theta gives =\n');
for i = 1:8
   fprintf('%8s = %f\n', variables{i,1}, thetaML(i));
end
fprintf('p_thetaML = %f\n', p_thetaML);
save([directory '/' hostname '_results_ga_' num2str(Nexpsamples) '_samples.mat']);


%% Plot final population
opengl software; % to avoid problems with the graphics card controller
for i = 1:8
   figure
   hist(final_pop(:,i), ceil(sqrt(PopulationSize)));
   hold on;
   plot(thetaML(i), 0, 'rx');
   grid on;
   title(variables{i,1},'Interpreter','none');

   % Save the plot in a .eps file (it adds automatically the .eps)
   print('-depsc', [directory '/' hostname '_ga_' num2str(Nexpsamples) '_samples_' variables{i,1} '.eps']); % or '-dpsc2',
end


%% Bayesian estimation of theta: bayesian model updating using TMCMC
tic
Nsamples = 10000;
fprintf('Nsamples TMCMC = %d\n', Nsamples);
samples_ftheta_D = problemA_tmcmc(log_p_D_theta, p_theta, p_thetarnd, Nsamples);
save([directory '/' hostname '_results_' num2str(Nexpsamples) '_samples.mat']);
toc

%% Plot the results of the maximum likelihood and the bayesian estimations
for i = 1:8
   figure
   subplot(2,1,1);
   ksdensity(samples_ftheta_D(:,i), 'support', [lb(i) ub(i)]);
   hold on;
   plot(thetaML(i), 0, 'rx');
   grid on;
   title(variables{i,1},'Interpreter','none');

   subplot(2,1,2);
   hist(samples_ftheta_D(:,i), ceil(sqrt(Nsamples)));
   hold on;
   plot(thetaML(i), 0, 'rx');
   grid on;
   title(variables{i,1},'Interpreter','none');

   % Save the plot in a .eps file (it adds automatically the .eps)
   print('-depsc', [directory '/' hostname '_tmcmc_' num2str(Nexpsamples) '_samples_' variables{i,1} '.eps']); % or '-dpsc2',
end

titulos = { ...
   'E[p_1]'
   'Var[p_1]'
   'p_2'
   'E[p_4]'
   'Var[p_4]'
   'E[p_5]'
   'Var[p_5]'
   '\rho[p_4,p_5]' };

figure
for i = 1:8
   hh = subplot(2,4,i);
   if (i>=5)
      pp = get(hh, 'pos');
      pp(2) = pp(2) + 0.17;
      set(hh, 'pos', pp);
   end
   histnorm(samples_ftheta_D(:,i), ceil(sqrt(Nsamples)), 'plot');
   hold on;
   grid on;
   %title(titulos{i},'Interpreter','none');
   xlabel(titulos{i});
   if (i==1) || (i==5)
      ylabel('normalised bin counts');
   end
end

% Save the plot in a .eps file (it adds automatically the .eps)
%print('-depsc', [directory '/' hostname '_tmcmc_' num2str(Nexpsamples) '_samples_all_variables.eps']); % or '-dpsc2',
saveas(gcf,[directory '/' hostname '_tmcmc_' num2str(Nexpsamples) '_samples_all_variables.fig'],'fig')
saveas(gcf,[directory '/' hostname '_tmcmc_' num2str(Nexpsamples) '_samples_all_variables.pdf'],'eps')

diary off;
