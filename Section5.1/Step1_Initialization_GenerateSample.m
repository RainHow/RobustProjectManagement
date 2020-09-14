clear;
clc;

%******************************************
instance_type = 1;
%******************************************
K = 3;
dis_n = 6; dis_alpha = 1; dis_beta = 3;
%******************************************


%%
% independent
outputfile = ['DATA/DC' num2str(instance_type) '/Sample_Training.mat'];
SampleSize = 10^7;

prob_ind = betarnd(dis_alpha,dis_beta,K,SampleSize);
DeltaSample = binornd(dis_n,prob_ind);

save(outputfile, 'SampleSize', 'DeltaSample');

%%
% State Generating
outputfile = ['DATA/DC' num2str(instance_type) '/StateTesting.mat'];

REALIZATION = linspace(0,dis_n,dis_n+1);

[X1,X2,X3] = ndgrid( REALIZATION, REALIZATION, REALIZATION);
State = [X1(:) X2(:) X3(:)];

StateSize = size(State,1);

prob = zeros(StateSize,1);

for i = 1:StateSize
    temp = ismember(DeltaSample',State(i,:),'rows');
    prob(i) = sum(temp>0)/SampleSize    
    i   
end

save(outputfile, 'StateSize', 'State','prob');
