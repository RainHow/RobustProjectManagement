clc;
clear;

instance_type = 1;

% Index = [1:100; 361:460; 721:820; 1081:1180; 1441:1540];  %100
% Index = [1:150; 361:510; 721:870; 1081:1230; 1441:1590];  %150
% Index = [1:200; 361:560; 721:920; 1081:1280; 1441:1640];  %200
% Index = [1:250; 361:610; 721:970; 1081:1330; 1441:1690];  %250
 Index = [1:300; 361:660; 721:1020; 1081:1380; 1441:1740];  %300


Time = zeros(size(Index,2),5);
Gamma_all = zeros(size(Index,2),5);

load(['DATA/DC' num2str(instance_type) '/InstanceAnalysis.mat' ]);
NGT = [];

for i = 1:5
    %eval([ 'Index_temp = Index_' num2str(i) ';' ])
    Index_temp = Index(i,:);
    k = 1;
    for j = Index_temp
        load(['DATA/DC' num2str(instance_type) '/Instance_Solution/result_aro' num2str(j) '.mat']);
        Time(k,i) = T_aro;
        Gamma_all(k,i) = Gamma(j);
        NGT = [NGT; N, Gamma(j), T_aro];
        k = k+1;
    end
end

eval(['save(''DATA/DC' num2str(instance_type) '/NGT.mat'', ''NGT'',''Time'',''Gamma_all'');']);