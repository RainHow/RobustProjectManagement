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

Ite_mean = []; Time_mean = []; Time_master_mean = [];

for j = 1:5
    solution_type = solution_type_collection{j};
    eval([ 'Ite_' solution_type '=zeros(120,1);' ]);
    eval([ 'Time_' solution_type '=zeros(120,1);' ]);

    for i = 1:number_instance
            load(['DATA/DC' num2str(instance_type) '/Instance_Solution/result_' solution_type num2str(i) '.mat']);
            eval([ 'Ite_' solution_type '(i)=Iteration_' solution_type ';' ])  
            eval([ 'Time_' solution_type '(i)=T_' solution_type ';' ])        
    end
    
    eval([ 'Ite_mean = [Ite_mean,mean(Ite_' solution_type ')];' ])
    eval([ 'Time_mean = [Time_mean,mean(Time_' solution_type ')];' ])

end

for j = 1:3
    solution_type = solution_type_collection{j};
    eval([ 'Time_master_' solution_type '=zeros(120,1);' ]);
    for i = 1:number_instance
            load(['DATA/DC' num2str(instance_type) '/Instance_Solution/result_' solution_type num2str(i) '.mat']);
            eval([ 'Time_master_' solution_type '(i)=mean(T_master_' solution_type ');' ])               
    end    
    eval([ 'Time_master_mean = [Time_master_mean,mean(Time_master_' solution_type ')];' ])
end


eval(['save(''DATA/DC' num2str(instance_type) '/SolutionAnalysis.mat'', ''Ite_*'', ''Time_*'');']);

