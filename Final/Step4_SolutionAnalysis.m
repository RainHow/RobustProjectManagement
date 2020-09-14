clear;
clc;

%*******************************
instance_type = 2;
number_instance = 120;
%*******************************

solution_type_collection = cell(1,5);
solution_type_collection{1} = 'aro';
solution_type_collection{2} = 'ccg';
solution_type_collection{3} = 'bd';
solution_type_collection{4} = 'eo';
solution_type_collection{5} = 'um';


outputfile = ['DATA/DC' num2str(instance_type) '/SolutionAnalysis.mat'];
% TimeAll = zeros(number_instance,5);
% IteAll = zeros(number_instance,5);


for solution_approach_index = [1 2 3 5]
    solution_approach = solution_type_collection{solution_approach_index};
    eval([ 'Ite_' solution_approach '= zeros(' num2str(number_instance) ',1);' ])
    eval([ 'Time_' solution_approach '= zeros(' num2str(number_instance) ',1);' ])
    
    for i = 1 : number_instance
        eval(['load(''DATA/DC' num2str(instance_type) '/Instance_Solution/result_' solution_approach num2str(i) '.mat'');']);
        eval([ 'Ite_' solution_approach '(' num2str(i) ') = Iteration_' solution_approach ';' ])
        eval([ 'Time_' solution_approach '(' num2str(i) ') = T_' solution_approach ';' ])
    end
    
end

%save(outputfile, 'Ite_*', 'Time_*', 'solution_type_collection');