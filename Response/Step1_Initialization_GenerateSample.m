clear;
clc;

%******************************************
instance_type = 2;
%******************************************
K = 7;
dis_n = 4; dis_alpha = 1; dis_beta = 3;
%******************************************


%%
% independent
outputfile = ['DATA/DC' num2str(instance_type) '/Sample_Training_3.mat'];
SampleSize = 10^7;

prob_ind = betarnd(dis_alpha,dis_beta,K,SampleSize);
DeltaSample = binornd(dis_n,prob_ind);

save(outputfile, 'SampleSize', 'DeltaSample');

%%
% Training data with 0.01 correlation 
outputfile = ['DATA/DC' num2str(instance_type) '/Sample_Training_3.mat'];
SampleSize = 10^6;


prob_com = betarnd(dis_alpha,dis_beta,1,SampleSize);
Common_Sample = binornd(dis_n,prob_com);

prob_ind = betarnd(dis_alpha,dis_beta,K,SampleSize);
Independent_Sample = binornd(dis_n,prob_ind);

Common_Chosen = binornd(1,0.1*ones(K,SampleSize));
Independent_Chosen = 1 - Common_Chosen;

DeltaSample = Common_Chosen.* (ones(K,1)*Common_Sample) + Independent_Chosen.*Independent_Sample;

save(outputfile, 'SampleSize', 'DeltaSample');
