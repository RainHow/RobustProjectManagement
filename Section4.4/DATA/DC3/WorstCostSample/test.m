clear
clc;

load Cost_aro1.mat
cost_sample_aro = cost_sample;
probability_aro = probability;

load Cost_dist1.mat
cost_sample_dist = cost_sample;
probability_dist = probability;

clear probability cost_sample

Data_aro = [ cost_sample_aro probability_aro ];
Data_dist = [ cost_sample_dist probability_dist ];

mu_aro = cost_sample_aro'*probability_aro;

mu_dist = cost_sample_dist'*probability_dist;


% tau = 234.64;
% alpha = 20;
% cost_tau_aro = cost_sample_aro - tau;
% probability_dist'*max(-alpha,cost_tau_aro)


% cost_tau_dist = cost_sample_dist - tau;
% probability_dist'*max(-alpha,cost_tau_dist)


