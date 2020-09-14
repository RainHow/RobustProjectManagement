clc;
clear;


%%
% CompareTime
instance_type = 2

Time_1 = [];
for i = 1:20
    load(['DATA/DC' num2str(instance_type) '/Instance_Solution/result_aro' num2str(i) '.mat']);
    load(['DATA/DC' num2str(instance_type) '/Instance_Solution/result_aro_rev' num2str(i) '.mat']);
    Time_1 = [Time_1; T_aro T_aro_rev];
    
end

Time_2 = [];
for i =21:39
    load(['DATA/DC' num2str(instance_type) '/Instance_Solution/result_aro' num2str(i) '.mat']);
    load(['DATA/DC' num2str(instance_type) '/Instance_Solution/result_aro_rev' num2str(i) '.mat']);
    Time_2 = [Time_2; T_aro T_aro_rev];  
end

Time_3 = [];
for i =41:45
    load(['DATA/DC' num2str(instance_type) '/Instance_Solution/result_aro' num2str(i) '.mat']);
    load(['DATA/DC' num2str(instance_type) '/Instance_Solution/result_aro_rev' num2str(i) '.mat']);
    Time_3 = [Time_3; T_aro T_aro_rev];  
end

eval(['save(''DATA/DC' num2str(instance_type) '/TimeAnalysis.mat'', ''Time_*'');']);