% ANN are on the structure Tnn

%Tnn(length(VhiddenNodes)).net=[];
%Tnn(length(VhiddenNodes)).r=[];
%Tnn(length(VhiddenNodes)).cpuTime=[];

%% ANN using FANN

nOutput=size(MoutputValidation,2);
load('Xout_FE_metamodel_validationGM_noNeg')
XoutFE=Xout_FE_metamodel_validationGM_noNeg;

CmembersNoNorm= {'theta01','theta02','theta03','theta04',...
    'theta05','theta06','theta07','theta08','theta09','theta10','theta11',...
    'theta12','theta13','theta14','theta15','theta16','theta17','theta18'};

Minput = XoutFE.getValues('CSnames',CmembersNoNorm); % Get values of normalised inputs (from GMD)
Moutput=XoutFE.getValues('CSnames',CeigfNames);
% Remove NaN and Inf
indexNoFiniteInput=find(~isfinite(Minput));
indexNoFiniteOutput=find(~isfinite(Moutput));
indexNoFinite=unique([indexNoFiniteOutput;indexNoFiniteInput]);
Minput(mod(indexNoFinite,size(Minput,1)),:)=[];
Moutput(mod(indexNoFinite,size(Moutput,1)),:)=[];
        
NsamplesTraining=size(Minput,1);
NcalibrationSamples = floor(NsamplesTraining*0.9);
Vpermutation=randperm(NsamplesTraining);       

% assign input and output data
Xinput = XfullModel.Xinput.sample('Nsamples',1);
XcalibrationInput = Xinput; XvalidationInput = Xinput;
        
MinputTraining=Minput(Vpermutation(1:NcalibrationSamples),:);
MinputValidation=Minput(Vpermutation(NcalibrationSamples+1:end),:);
        
XcalibrationInput.Xsamples.MsamplesPhysicalSpace = MinputTraining;
XvalidationInput.Xsamples.MsamplesPhysicalSpace = MinputValidation;
        
MoutputTraining=Moutput(Vpermutation(1:NcalibrationSamples),:);
MoutputValidation=Moutput(Vpermutation(NcalibrationSamples+1:end),:);
XcalibrationOutput=SimulationData('Mvalues',MoutputTraining ,'Cnames',CeigfNames);
XvalidationOutput=SimulationData('Mvalues',MoutputValidation,'Cnames',CeigfNames);
       
Vnodes=[14 9];
tic;
for imode=1:nOutput
    Xnn(imode) = NeuralNetwork('Stype','HyperbolicTangent',...
                    'VhiddenNodes',Vnodes,...
                    'XFullModel',XfullModel,...
                    'Cinputnames',CmembersNoNorm,...
                    'Coutputnames',CeigfNames(imode),...
                    'XcalibrationInput',XcalibrationInput,...
                    'XcalibrationOutput',XcalibrationOutput,...
                    'XvalidationInput',XvalidationInput,...
                    'XvalidationOutput',XvalidationOutput,...
                    'Vnormminmax',[-.9 .9]);
    Xnn(imode) = Xnn(imode).calibrate();
    Xnn(imode) = Xnn(imode).validate();
    Xnn_no_data(imode) = Xnn(imode);
    Xnn_no_data(imode).XcalibrationInput=[];Xnn_no_data(imode).XcalibrationOutput=[];
    Xnn_no_data(imode).XvalidationInput=[];Xnn_no_data(imode).XvalidationOutput=[];
    fprintf('FANN %i: Performance: %f\n', ...
                        imode,Xnn(imode).VvalidationError)
end
Tnn(length(VhiddenNodes)+1).cpuTime=toc;

for imode=1:nOutput
Tnn(length(VhiddenNodes)+1).Rsq(n)=Xnn(imode).VvalidationError;  
RsqFann(imode)=Xnn(imode).VvalidationError;
end

VscoreTotFANN=sum(RsqFann);
[VscoreWorstFANN, IndexMinFANN]=min(RsqFann);
[~, IndexMaxFANN]=max(RsqFann);

% Plot validation ANN from FANN
Xnn(IndexMinFANN).plotregression;
saveas(gcf,['ValidationANN_FANN_' num2str(IndexMinFANN)],'fig')
saveas(gcf,['ValidationANN_FANN_' num2str(IndexMinFANN)],'pdf')

Xnn(IndexMaxFANN).plotregression;
saveas(gcf,['ValidationANN_FANN_' num2str(IndexMaxFANN)],'fig')
saveas(gcf,['ValidationANN_FANN_' num2str(IndexMaxFANN)],'pdf')

% Save FANN 
save FANN_25072016 Xnn Xnn_no_data
%Xnn contains the ANN from FANNimagesc(CorrMatOutputFE);




%% Plot Correlation outputs
% Moutput: full output from FE
XinputFANN = XfullModel.Xinput.sample('Nsamples',1);
XinputFANN.Xsamples.MsamplesPhysicalSpace = Minput;

outputFANN = zeros(length(TinputFANN),nOutput);
for inn = 1:nOutput
    outputFANN(:,inn)  = [Xnn(inn).McalibrationOutput; Xnn(inn).MvalidationOutput];
end

%save Mycmap Mycmap

CorrMatOutputFE=corrcoef(Moutput);
imagesc(CorrMatOutputFE);
colormap(Mycmap);colorbar;
colormap(Mycmap);colorbar;
xlabel('\omega')
ylabel('\omega')
saveas(gcf,['CorrelationOutputFE'],'fig')
saveas(gcf,['CorrelationOutputFE'],'pdf')

figure
CorrMatOutputFANN=corrcoef(outputFANN);
imagesc(CorrMatOutputFANN);
colormap(Mycmap);colorbar;
xlabel('\omega')
ylabel('\omega')
title('ANN (FANN)')
saveas(gcf,['CorrelationOutputANN'],'fig')
saveas(gcf,['CorrelationOutputANN'],'pdf')


%% Compute the R^2 of the ANN
for j=1:length(Tnn)-1
    ANNoutputs = Tnn(j).net(MinputValidation');
    for n=1:nOutput
        Y=MoutputValidation(:,n);
        Y_calc=ANNoutputs(n,:)';
        Tnn(j).Rsq(n) = 1 - sum((Y - Y_calc).^2)/sum((Y - mean(Y)).^2);
    end
    VscoreTot(j)=sum(Tnn(j).Rsq(n));
    [VscoreWorst(j), IndexMin]=min(Tnn(j).Rsq);
        
    hf=figure;
    scatter(MoutputValidation(:,IndexMin),ANNoutputs(IndexMin,:))
    hold on
    plot(MoutputValidation(:,IndexMin),MoutputValidation(:,IndexMin))
    box
    xlabel('FE output')
    ylabel('ANN output')
    title(['Validation ' CeigfNames{IndexMin}, ', ANN: ' num2str(VhiddenNodes{IndexMin})])
    grid on
    saveas(gcf,['ValidationANN_' num2str(j)],'fig')
    saveas(gcf,['ValidationANN_' num2str(j)],'pdf')
end  

VscoreTot=zeros(length(Tnn),1);
VscoreWorst=zeros(length(Tnn),1);
Vcpu=zeros(length(Tnn),1);
for n=1:length(Tnn)-1
   VscoreTot(n)=sum(Tnn(n).r)/length(Tnn(n).r);
   VscoreWorst(n)=min(Tnn(n).r);
   Vcpu(n)=Tnn(n).cpuTime; 
end
    
VscoreTot(length(VhiddenNodes)+1)=sum(RsqFann)/length(RsqFann);
VscoreWorst(length(VhiddenNodes)+1)=min(RsqFann);
Vcpu(length(VhiddenNodes)+1)=Tnn(length(VhiddenNodes)+1).cpuTime;


figure
yyaxis left
bar([VscoreTot VscoreWorst])
xlabel('ANN configuration')
ylabel('Performance')
legend({'Mean performace ($\sum_i(R_i^2)/\sum_i$)' 'Worst performace ($R_i^2$)' 'Cputime'},'Interpreter','latex')

yyaxis right
plot([1:length(Tnn)],Vcpu);

legend({'Mean performace ($\sum_i(R_i^2)/\sum_i$)' 'Worst performace ($R_i^2$)' 'Cputime'},'Interpreter','latex')

saveas(gcf,'PerformanceANN_Bar','fig')
saveas(gcf,'PerformanceANN_Bar','pdf')
    
