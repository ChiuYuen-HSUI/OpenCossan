%% run FE model with prior samples
Spath = fileparts(which('AIRMOD_Bayesian_part1a_LHS.m'));
OpenCossan.setWorkingPath(fullfile(Spath,'workfolder'));
load(fullfile(Spath,'fullmodel.mat'))

% load prior samples MsamplesPhysicalSpace
load MsamplesPhysicalSpaceExpanded5 

Xrvset = GaussianMixtureRandomVariableSet('Mdataset',MsamplesPhysicalSpace,...
    'Cmembers',{'theta01','theta02','theta03','theta04',...
    'theta05','theta06','theta07','theta08','theta09','theta10','theta11',...
    'theta12','theta13','theta14','theta15','theta16','theta17','theta18'},...
    'Vconstraints',zeros(18,1),'Mcoeff',-1*eye(18),'Lrejection',true);

Xinput=Input('CSmembers',{'Xrvset'},'CXmembers',{Xrvset});
Xinput =Xinput.sample('Nsamples',2000);

Xmodel.Xevaluator.Vconcurrent=[10 0];
 
Xout_FE_prior5 = Xmodel.apply(Xinput);

save Xout_FE_prior5 Xout_FE_prior5

