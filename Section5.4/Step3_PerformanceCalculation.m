clear;
clc;
%************************************
instance_type = 2;
%************************************

solution_type_collection = cell(1,5);
solution_type_collection{1} = 'aro';
solution_type_collection{2} = 'ccg_new';
solution_type_collection{3} = 'bd';
solution_type_collection{4} = 'eo';
solution_type_collection{5} = 'um';

load(['DATA/DC' num2str(instance_type) '/StateTesting.mat']);
State = State';

for instance_index = 1:240
    instance_information = load(['DATA/DC' num2str(instance_type) '/MAT_for_Optimization/instance' num2str(instance_index) '.mat']);
    for solution_approach_index = 2
        solution_approach = solution_type_collection{solution_approach_index};
        solution_information = load(['DATA/DC' num2str(instance_type) '/Instance_Solution/result_' solution_approach num2str(instance_index) '.mat']);
        eval(['x = double(solution_information.x_' solution_approach ');'])
        eval(['y = double(solution_information.y_' solution_approach ');'])
        
        cost_sample = zeros(StateSize,1);
        parfor m = 1 : StateSize
            cost_sample(m) = instance_information.Ux'*x + instance_information.Uy'*y + second_stage_cost(x, y, State(:,m), instance_information.N, instance_information.N_pr, instance_information.Q, instance_information.bar_n, instance_information.Gf, instance_information.Gc, instance_information.d, instance_information.AF, instance_information.ACC, instance_information.AV, instance_information.AC, instance_information.AX, instance_information.AY, instance_information.AZ, instance_information.Z0, instance_information.Z, instance_information.B0);
        end
%         save(['DATA/DC' num2str(instance_type) '/CostSample/Cost_' solution_approach num2str(instance_index) '.mat'], 'cost_sample');
        save(['DATA/DC' num2str(instance_type) '/CostSample/Cost_ccg' num2str(instance_index) '.mat'], 'cost_sample');

    end
end

