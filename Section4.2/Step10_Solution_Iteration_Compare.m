clear;
clc;

instance_type = 2;
number_instance = 120;

solution_type_collection = cell(1,5);
solution_type_collection{1} = 'aro';
solution_type_collection{2} = 'ccg';
solution_type_collection{3} = 'bd';
solution_type_collection{4} = 'eo';
solution_type_collection{5} = 'um';

Ite_mean = [];

for j = 1:5
    solution_type = solution_type_collection{j};
    eval([ 'Ite_' solution_type '=zeros(120,1);' ]);
    for i = 1:number_instance
            load(['DATA/DC' num2str(instance_type) '/Instance_Solution/result_' solution_type num2str(i) '.mat']);
            eval([ 'Ite_' solution_type '(i)=Iteration_' solution_type ';' ])       
        
    end
    eval([ 'Ite_mean = [Ite_mean,mean(Ite_' solution_type ')];' ])
end



eval(['save(''DATA/DC' num2str(instance_type) '/IterationAnalysis.mat'', ''Ite_*'');']);

