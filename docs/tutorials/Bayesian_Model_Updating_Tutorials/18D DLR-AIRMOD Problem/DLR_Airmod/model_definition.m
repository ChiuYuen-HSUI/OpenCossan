%% Define AIRMOD FE model
function [Xmodel]=model_definition(Spath,'warp.q')

%% Model definition
%% Connector definition
%{
Identification of the parameters used in the FE with the parameters of the
paper

These paramenters are listed as already updated, I tried to understand what
parameter each line represents by looking at their updated value (order of
magnitude) and their location on the model (id number)

Original file line      | Possible identification

PELAS    160    160.095 | maybe 03 (if they are grouped as in paper)
PELAS    171    171.654 | maybe 04 (if they are grouped as in paper)
PELAS    172    171.234 | maybe 05 (if they are grouped as in paper)
PELAS    1      2.65+6  | maybe 06
PELAS    3      5.65+8  | surely 07
PELAS    4      4+7     | 14, 15, 16
PELAS    5      9.42+6  | 14, 15, 16
PELAS    6      2.43+6  | 14, 15, 16
PELAS    10     7.06+6  | 17, 18
PELAS    11     5.23+6  | 17, 18
MAT1     2      3313.44 | either 01, 02 (01 if they are grouped as in paper)
MAT1     3      4397.10 | either 01, 02 (02 if they are grouped as in paper)
PMASS    103    0.16825 | either 08, 09, 10 (08 if they are grouped as in paper)
PMASS    1041   0.19470 | either 08, 09, 10 (09 if they are grouped as in paper)
PMASS    1042   0.17434 | either 08, 09, 10 (10 if they are grouped as in paper)
PMASS    105    0.03387 | either 11, 12, 13 (11 if they are grouped as in paper)
PMASS    106    7.35-3  | either 11, 12, 13 (12 if they are grouped as in paper)
PMASS    107    7.55-3  | either 11, 12, 13 (13 if they are grouped as in paper)

UPDATE:
Parameters identification confirmed!
%}
Xinj = Injector('Sscanfilepath',Spath,...
    'Sscanfilename','PARAM.DAT.cossan',...
    'Sfile','PARAM.DAT');

Nmodes = 30; % nr. of extracted modes
Xresponse = Response(...
    'Sname','eigenfrequencies',...
    'Clookoutfor',{'         R E A L   E I G E N V A L U E S'},...
    'Nrow',6, 'Ncol',68,...
    'Nrepeat',Nmodes,...
    'Sformat','%12e');

Xext = Extractor('Sfile','AIRMOD103.f06',...
    'CXresponse',{Xresponse});

Xconnector = Connector('Stype','nastran',...
    'Ssolverbinary','/Apps/msc/MSC_Nastran/20131/bin/nast20131',... '/usr/software/Nastran/20122/bin/nast20122',...
    'Sexecmd','%Ssolverbinary %Smaininputfile %Sexeflags',...
    'Sexeflags','scr=yes news=no bat=no old=no',...
    'SerrorFileExtension','f06',...
    'SerrorString','FATAL',...
    'SmainInputPath',Spath,...
    'SmainInputFile','AIRMOD103.bdf',...
    'CSadditionalFiles',{'AIRMOD_SOLID_MODEL_ONLY.BDF'},...
    'LkeepSimulationFiles',false);

Xconnector = Xconnector.add(Xinj);
Xconnector = Xconnector.add(Xext);
%% Post-processing MIO
% This script splits the eigenfrequencies, saved in a dataseries by the
% connector, into individual scalar outputs
CeigfNames = cell(1,Nmodes);
for i = 1:Nmodes
    CeigfNames{i} = ['freq' num2str(i)];
end

Xmio = Mio('Spath',Spath,'Sfile','mio_post.m',...
    'CinputNames',{'eigenfrequencies'},...
    'Coutputnames',CeigfNames,...
    'Lfunction',false,'LioStructure',true,'LioMatrix',false);


%% Evaluator
Xeval = Evaluator('CXmembers',{Xconnector,Xmio},'CSmembers',{'Xconnector','Xmio'},...
    'XjobManagerInterface',JobManagerInterface('Stype','GridEngine'),...
    'CSqueues',{'warp.q',''},'CShostnames',{'iru2',''},...
    'Vconcurrent',[24 0],'LremoteInjectExtract',true);
%% Prior
% Define the prior random variables. The mean values are assigned according
% to table 2.

% We set the lowerbound to 95% under the mean and upperbound to 100% over
% the mean. The magic numbers come from:
[mean_factor,v] = unifstat(0.05,2.0); % normalized bounds
CoV = sqrt(v)/mean_factor;

%% Defining the prior PDF p(theta)
load LatinHCfromSensitivity 
sensitivity_samples=Xout.getValues('Cnames',{'theta01','theta02','theta03', ...
    'theta04','theta05','theta06','theta07','theta08','theta09','theta10',...
    'theta11','theta12','theta13','theta14','theta15','theta16','theta17','theta18'});

OrderOfMagnitude = 10.^floor(log10(sensitivity_samples(1,:)));
sensitivity_samples = sensitivity_samples./repmat(OrderOfMagnitude,size(sensitivity_samples,1),1);

Xrvset = GaussianMixtureRandomVariableSet('Mdataset',sensitivity_samples,...
    'Cmembers',{'theta01_norm','theta02_norm','theta03_norm','theta04_norm',...
    'theta05_norm','theta06_norm','theta07_norm','theta08_norm','theta09_norm','theta10_norm','theta11_norm',...
    'theta12_norm','theta13_norm','theta14_norm','theta15_norm','theta16_norm','theta17_norm','theta18_norm'},...
    'Vconstraints',zeros(18,1),'Mcoeff',-1*eye(18),'Lrejection',true);

% Vperm=randperm(size(sensitivity_samples,1));
% 
% Xrvset2 = GaussianMixtureRandomVariableSet('Mdataset',sensitivity_samples(Vperm(1:50),:),...
%     'Cmembers',{'theta01_norm','theta02_norm','theta03_norm','theta04_norm',...
%     'theta05_norm','theta06_norm','theta07_norm','theta08_norm','theta09_norm','theta10_norm','theta11_norm',...
%     'theta12_norm','theta13_norm','theta14_norm','theta15_norm','theta16_norm','theta17_norm','theta18_norm'},...
%     'Vconstraints',zeros(18,1),'Mcoeff',-1*eye(18),'Lrejection',true);

% Using a variance that is 10 times the variance original in order to
% smooth down the 
Sigma=Xrvset.gmDistribution.Sigma;
Sigma=Sigma*10;

XrvsetLargeVariance = GaussianMixtureRandomVariableSet('Mdataset',sensitivity_samples,'Mcovariance',Sigma,...
    'Cmembers',{'theta01_norm','theta02_norm','theta03_norm','theta04_norm',...
    'theta05_norm','theta06_norm','theta07_norm','theta08_norm','theta09_norm','theta10_norm','theta11_norm',...
    'theta12_norm','theta13_norm','theta14_norm','theta15_norm','theta16_norm','theta17_norm','theta18_norm'},...
    'Vconstraints',zeros(18,1),'Mcoeff',-1*eye(18),'Lrejection',true);

theta01 = Function('Sexpression',['<&theta01_norm&>*' num2str(OrderOfMagnitude(1))]);
theta02 = Function('Sexpression',['<&theta02_norm&>*' num2str(OrderOfMagnitude(2))]);
theta03 = Function('Sexpression',['<&theta03_norm&>*' num2str(OrderOfMagnitude(3))]);
theta04 = Function('Sexpression',['<&theta04_norm&>*' num2str(OrderOfMagnitude(4))]);
theta05 = Function('Sexpression',['<&theta05_norm&>*' num2str(OrderOfMagnitude(5))]);
theta06 = Function('Sexpression',['<&theta06_norm&>*' num2str(OrderOfMagnitude(6))]);
theta07 = Function('Sexpression',['<&theta07_norm&>*' num2str(OrderOfMagnitude(7))]);
theta08 = Function('Sexpression',['<&theta08_norm&>*' num2str(OrderOfMagnitude(8))]);
theta09 = Function('Sexpression',['<&theta09_norm&>*' num2str(OrderOfMagnitude(9))]);
theta10 = Function('Sexpression',['<&theta10_norm&>*' num2str(OrderOfMagnitude(10))]);
theta11 = Function('Sexpression',['<&theta11_norm&>*' num2str(OrderOfMagnitude(11))]);
theta12 = Function('Sexpression',['<&theta12_norm&>*' num2str(OrderOfMagnitude(12))]);
theta13 = Function('Sexpression',['<&theta13_norm&>*' num2str(OrderOfMagnitude(13))]);
theta14 = Function('Sexpression',['<&theta14_norm&>*' num2str(OrderOfMagnitude(14))]);
theta15 = Function('Sexpression',['<&theta15_norm&>*' num2str(OrderOfMagnitude(15))]);
theta16 = Function('Sexpression',['<&theta16_norm&>*' num2str(OrderOfMagnitude(16))]);
theta17 = Function('Sexpression',['<&theta17_norm&>*' num2str(OrderOfMagnitude(17))]);
theta18 = Function('Sexpression',['<&theta18_norm&>*' num2str(OrderOfMagnitude(18))]);

Xinput = Input('CXmembers',{XrvsetLargeVariance,theta01,theta02,theta03,theta04,theta05,theta06,theta07,...
    theta08,theta09,theta10,theta11,theta12,theta13,theta14,theta15,theta16,theta17,theta18},...
    'CSmembers',{'Xrvset','theta01','theta02','theta03','theta04','theta05','theta06','theta07',...
    'theta08','theta09','theta10','theta11','theta12','theta13','theta14','theta15','theta16','theta17','theta18'});

%% Test model
Xmodel = Model('Xinput',Xinput,'Xevaluator',Xeval);
% Xout_ini = Xmodel.apply(TinputDefault);
