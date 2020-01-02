Spath = fileparts(which('AIRMOD_Bayesian_part1a_LHS.m'));
OpenCossan.setWorkingPath(fullfile(Spath,'workfolder'));
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
    'Ssolverbinary','/usr/software/Nastran/20122/bin/nast20122',... '/Apps/msc/MSC_Nastran/20131/bin/nast20131',...
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
    'CSqueues',{'all.q',''},'CShostnames',{'cossan.cfd.liv.ac.uk',''},...
    'Vconcurrent',[25 0],'LremoteInjectExtract',true);
%% Prior
% Define the prior random variables. The mean values are assigned according
% to table 2.

% We set the lowerbound to 95% under the mean and upperbound to 100% over
% the mean. The magic numbers come from:
[mean_factor,v] = unifstat(0.05,2.0); % normalized bounds
CoV = sqrt(v)/mean_factor;

% support stiffness
theta01 = RandomVariable('Sdistribution','uniform','mean',mean_factor*1.8e3,'cov',CoV); TinputDefault.theta01 = 1.8e3;
theta02 = RandomVariable('Sdistribution','uniform','mean',mean_factor*7.5e3,'cov',CoV); TinputDefault.theta02 = 7.5e3;
theta03 = RandomVariable('Sdistribution','uniform','mean',mean_factor*1.3e2,'cov',CoV); TinputDefault.theta03 = 1.3e2;
theta04 = RandomVariable('Sdistribution','uniform','mean',mean_factor*7.0e1,'cov',CoV); TinputDefault.theta04 = 7.0e1;
theta05 = RandomVariable('Sdistribution','uniform','mean',mean_factor*7.0e1,'cov',CoV); TinputDefault.theta05 = 7.0e1;

% joint stiffness
theta06 = RandomVariable('Sdistribution','uniform','mean',mean_factor*1.0e7,'cov',CoV); TinputDefault.theta06 = 1.0e7;
theta07 = RandomVariable('Sdistribution','uniform','mean',mean_factor*1.0e9,'cov',CoV); TinputDefault.theta07 = 1.0e9;
theta14 = RandomVariable('Sdistribution','uniform','mean',mean_factor*2.0e7,'cov',CoV); TinputDefault.theta14 = 2.0e7;
theta15 = RandomVariable('Sdistribution','uniform','mean',mean_factor*2.0e7,'cov',CoV); TinputDefault.theta15 = 2.0e7;
theta16 = RandomVariable('Sdistribution','uniform','mean',mean_factor*7.0e6,'cov',CoV); TinputDefault.theta16 = 7.0e6;
theta17 = RandomVariable('Sdistribution','uniform','mean',mean_factor*5.0e7,'cov',CoV); TinputDefault.theta17 = 5.0e7;
theta18 = RandomVariable('Sdistribution','uniform','mean',mean_factor*5.0e7,'cov',CoV); TinputDefault.theta18 = 5.0e7;

% masses (sensor cables, screws and glue)
theta08 = RandomVariable('Sdistribution','uniform','mean',mean_factor*2.0e-1,'cov',CoV); TinputDefault.theta08 = 2.0e-1;
theta09 = RandomVariable('Sdistribution','uniform','mean',mean_factor*1.86e-1,'cov',CoV); TinputDefault.theta09 = 1.86e-1;
theta10 = RandomVariable('Sdistribution','uniform','mean',mean_factor*1.86e-1,'cov',CoV); TinputDefault.theta10 = 1.86e-1;
theta11 = RandomVariable('Sdistribution','uniform','mean',mean_factor*1.5e-2,'cov',CoV); TinputDefault.theta11 = 1.5e-2;
theta12 = RandomVariable('Sdistribution','uniform','mean',mean_factor*1.5e-2,'cov',CoV); TinputDefault.theta12 = 1.5e-2;
theta13 = RandomVariable('Sdistribution','uniform','mean',mean_factor*1.5e-2,'cov',CoV); TinputDefault.theta13 = 1.5e-2;

Xrvset = RandomVariableSet('Cmembers',{'theta01','theta02','theta03','theta04',...
    'theta05','theta06','theta07','theta08','theta09','theta10','theta11',...
    'theta12','theta13','theta14','theta15','theta16','theta17','theta18'});

Xinput = Input('CXmembers',{Xrvset},'CSmembers',{'Xrvset'});

%% Test model
Xmodel = Model('Xinput',Xinput,'Xevaluator',Xeval);
% Xout_ini = Xmodel.apply(TinputDefault);
save fullmodel.mat Xmodel
%% Monte-carlo 
Xlhs = LatinHypercubeSampling('Nsamples',5000,'Nbatches',10);
Xout_LHS = Xlhs.apply(Xmodel);
SbatchFolder = Xout_LHS.SbatchFolder;

save fullmodel.mat Xmodel Xlhs SbatchFolder