Spath = fileparts(which('AIRMOD_Bayesian_part2_model_updating.m'));
OpenCossan.setWorkingPath(fullfile(Spath,'workfolder'));

%% experimental data
FE_frequencies = [0.29 0.56 0.82 2.14 5.6500 15.1100 33.3100 33.6200 35.3900 44.6600 47.2100 52.9100...
    60.5900 67.6900 102.5900 128.6200 132.0800 145.9100 206.7300 225.7300 ....
    261.5300 262.6400 278.7100 320.1500 321.6400 324.1200 336.3100 341.1500...
    343.5500 359.5400];

test_frequencies = [0.23 0.65 0.83 2.17 5.50,14.91,31.96,32.33,34.38,43.89,46.71,51.88,58.59,...
    65.93,100.05,124.56,129.38,141.47,205.59,219.07,254.73,255.02,272.08,...
    303.96,304.32,313.68,328.55,331.18,336.21,348.68];

%% load model
load(fullfile(Spath,'fullmodel.mat'))
Xinput = Xmodel.Xinput;
Xinput = Xinput.sample('Nsamples',1);

Nmodes = 30; % nr. of extracted modes
CeigfNames = cell(1,Nmodes);
for i = 1:Nmodes
    CeigfNames{i} = ['freq' num2str(i)];
end

%% load LHS results

% SbatchFolder = fullfile(Spath, 'LHS_10000');

Xout_LHS = SimulationData.load('SfileName',fullfile(SbatchFolder,['SimulationData_batch_1_of_' num2str(Xlhs.Nbatches) '.mat']));
for i=1:Xlhs.Nbatches
    Xout_LHS = Xout_LHS.merge(SimulationData.load('SfileName',fullfile(SbatchFolder,['SimulationData_batch_' num2str(i) '_of_' num2str(Xlhs.Nbatches) '.mat'])));
end
Xout_LHS = Xout_LHS.removeData('Sname','eigenfrequencies'); % remove Dataseries for speed

%% separate data in calibration and validation
Minput = Xout_LHS.getValues('CSnames',Xmodel.Cinputnames);
Toutput = Xout_LHS.Tvalues;
NcalibrationSamples = floor(Xout_LHS.Nsamples*0.8);

% assign input and output data
XcalibrationInput = Xinput; XvalidationInput = Xinput;
XcalibrationInput.Xsamples.MsamplesPhysicalSpace = Minput(1:NcalibrationSamples,:);
XvalidationInput.Xsamples.MsamplesPhysicalSpace = Minput(NcalibrationSamples:end,:);
XcalibrationOutput = Xout_LHS; XvalidationOutput = Xout_LHS;
XcalibrationOutput.Tvalues = Toutput(1:NcalibrationSamples);
XvalidationOutput.Tvalues = Toutput(NcalibrationSamples:end);

Xnn(Nmodes) = NeuralNetwork;
Xnn_no_data(Nmodes) = NeuralNetwork;
if isempty(gcp)
    parpool(12)
end
parfor imode = 1:Nmodes
    Xnn(imode) = NeuralNetwork('Stype','HyperbolicTangent',...
        'VhiddenNodes',[14 6],...
        'XFullModel',Xmodel,...
        'Cinputnames',Xmodel.Cinputnames,...
        'Coutputnames',CeigfNames(imode),...
        'XcalibrationInput',XcalibrationInput,...
        'XcalibrationOutput',XcalibrationOutput,...
        'XvalidationInput',XvalidationInput,...
        'XvalidationOutput',XvalidationOutput);
    Xnn(imode) = Xnn(imode).calibrate();
    Xnn(imode) = Xnn(imode).validate();
    Xnn_no_data(imode) = Xnn(imode);
    Xnn_no_data(imode).XcalibrationInput=[];Xnn_no_data(imode).XcalibrationOutput=[];
    Xnn_no_data(imode).XvalidationInput=[];Xnn_no_data(imode).XvalidationOutput=[];
end

save metamodel.mat Xnn 
save metamodel_light Xnn_no_data