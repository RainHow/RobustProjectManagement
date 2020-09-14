clear;
clc;

%****************************
instance_type = 2;
number_instance = 240;
%****************************

Instance_Index = [1:number_instance];

eval(['load(''data/DC' num2str(instance_type) '/InstanceAnalysis.mat'');']);
%eval(['load(''data/DC' num2str(instance_type) '/CostAnalysis.mat'');']);

eval(['load(''data/DC' num2str(instance_type) '/CostAnalysis_withTarget.mat'');']);


Gamma_sort = sortrows([Gamma,Instance_Index']);
Small_Gamma_ind = Gamma_sort(1:number_instance/2,2);               % fisrt 120 gamma
Large_Gamma_ind = Gamma_sort(number_instance/2+1:number_instance,2);             % last 120 gamma

type_collection = cell(1,4);
type_collection{1} = 'CCG';
type_collection{2} = 'BD';
type_collection{3} = 'EO';
type_collection{4} = 'UM';


% ********* one-way ANOVA *********
P_Value_N = zeros(11,1);
P_Value_Gamma = zeros(11,1);
Mu_N = zeros(11,2);
Mu_Gamma = zeros(11,2);
for j = 1:4
    type = type_collection{j};
    eval([ 'P_Value_N_' type '= zeros(11,1);' ])
    eval([ 'P_Value_Gamma_' type '= zeros(11,1);' ])
    eval([ 'Mu_N_' type '= zeros(11,2);' ])
    eval([ 'Mu_Gamma_' type '= zeros(11,2);'])    
    for i = 1 : 11
        
        eval([ 'Performance = P_' type '(:,i);' ])      % need change
        
        Performance_N = [Performance(1:number_instance/2) Performance(number_instance/2+1:end)];
        Performance_Gamma = [Performance(Small_Gamma_ind) Performance(Large_Gamma_ind)];
        
        eval([ 'Mu_N_' type '(i,:) = [mean(Performance_N(:,1)) mean(Performance_N(:,2))];' ])
        eval([ 'Mu_Gamma_' type '(i,:) = [mean(Performance_Gamma(:,1)) mean(Performance_Gamma(:,2))];' ])
        eval([ 'P_Value_N_' type '(i) = anova1(Performance_N,[],''off'');' ])
        eval([ 'P_Value_Gamma_' type '(i) = anova1(Performance_Gamma,[],''off'');' ])
    end
end


% ********* two-way ANOVA *********
P_Value = zeros(8,3);
for i = 1 : 11
    
    Performance = P_EO(:,i);
    
    Small_N = [ones(number_instance/2,1); zeros(number_instance/2,1)];
    Small_Gamma = zeros(number_instance,1);
    Small_Gamma(Small_Gamma_ind) = 1;
    
    group = {Small_N Small_Gamma};
    
    p = anovan(Performance, group, 'model', 2, 'display', 'off');
    
    P_Value(i,:) = p';
    
end

% eval(['save(''DATA/DC' num2str(instance_type) '/AnovaAnalysis_NoTarget.mat'', ''P_Value'', ''P_Value_*'', ''Mu_*'');']);