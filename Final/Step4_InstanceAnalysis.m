% find the gamma coefficient of instances
clear;
clc;

%***************************************
instance_type = 1;      
%***************************************

length = 1800;

N = zeros(length,1);      % initial "record matrix"
N_pr = zeros(length,1);
LongestChain = zeros(length,1);

Gamma = zeros(length,1);     % benchmark  Gamma0=0.4

PRECEDENCE = zeros(0,2);

for i = 1 : length
    instance_information = load(['DATA/DC' num2str(instance_type) '/MAT_for_Optimization/instance' num2str(i) '.mat']);
    N(i) = instance_information.N;
    N_pr(i) = instance_information.N_pr;
    PRECEDENCE = [PRECEDENCE; double(instance_information.Precedence)];  % processed data
end


% evaluate the longest path
parfor i = 1 : length
    
    Precedence = PRECEDENCE(sum(N_pr(1:i-1))+1 : sum(N_pr(1:i)),:);  % i_th instance
    
    rome_begin;
    h = rome_model('makespan');
    T_completion = newvar(N(i),1,'nonneg');
    rome_minimize( T_completion(end) );   % final spot
    for j = 1 : N_pr(i)
        source_node = Precedence(j,1);
        destination_node = Precedence(j,2);
        rome_constraint( T_completion(destination_node) >= T_completion(source_node) + 1 );
    end
    h.solve('CPLEX');
    
    LongestChain(i) = h.objective;
    
    Gamma(i) = (LongestChain(i)-1) / (N(i)-1);    % gamma coefficient
    
end

save(['DATA/DC' num2str(instance_type) '/InstanceAnalysis.mat'], 'N', 'N_pr', 'LongestChain', 'Gamma');