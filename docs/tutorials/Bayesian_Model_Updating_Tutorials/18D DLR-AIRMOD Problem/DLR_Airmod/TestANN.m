clear variable

load TestANN

XcalibrationInput.Xsamples.MsamplesPhysicalSpace = Minput(Vpermutation(1:20),:);
XvalidationInput.Xsamples.MsamplesPhysicalSpace = Minput(Vpermutation(21:30),:);
        
        XcalibrationOutput=SimulationData('Mvalues',Moutput(1:20,:),'Cnames',CeigfNames);
        XvalidationOutput=SimulationData('Mvalues',Moutput(21:30,:),'Cnames',CeigfNames);
        
                    parfor j=1:2
                Xnn(j) = NeuralNetwork('Stype','HyperbolicTangent',...
                    'VhiddenNodes',VhiddenNodes{j},...
                    'XFullModel',XfullModel,...
                    'Cinputnames',XfullModel.Xinput.CnamesFunction,...
                    'Coutputnames',CeigfNames,...
                    'XcalibrationInput',XcalibrationInput,...
                    'XcalibrationOutput',XcalibrationOutput,...
                    'XvalidationInput',XvalidationInput,...
                    'XvalidationOutput',XvalidationOutput,...
                    'Vnormminmax',[-1 1]);
                    end
                    
Xnn(1).XvalidationInput


%
for n=1:
XcorrInput.Xsamples.MsamplesPhysicalSpace = Minput;
XcorrOutput=Xnn.apply(XcorrInput);