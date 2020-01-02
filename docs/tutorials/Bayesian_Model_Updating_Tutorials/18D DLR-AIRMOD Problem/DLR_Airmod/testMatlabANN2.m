net = feedforwardnet([20 20]);

net.trainParam.showWindow=false;
net.trainParam.showCommandLine=true;

plotfit

[net,tr] = train(net,Minput',Moutput');

houseTargets=Moutput';

% houseOutputs = net(Minput');
% trOut = houseOutputs(tr.trainInd);
% vOut = houseOutputs(tr.valInd);
% tsOut = houseOutputs(tr.testInd);
% trTarg = houseTargets(tr.trainInd);
% vTarg = houseTargets(tr.valInd);
% tsTarg = houseTargets(tr.testInd);
% plotregression(trTarg,trOut,'Train',vTarg,vOut,'Validation',...
% tsTarg,tsOut,'Testing')M1V=Minput(:,1);
M1TO=Moutput(:,1);



ANNoutputs = net(Minput');
performance = perform(net,Moutput',ANNoutputs)
[r,m,b] = regression(Moutput',ANNoutputs);

plotfit(net,Minput(:,1)',Moutput(:,1)')

net = feedforwardnet(18,20,14);


M1T=Minput(:,1);
M1TO=Moutput(:,1);
%M1V=Minput(:,1);
%M1VO=Moutput(:,1);