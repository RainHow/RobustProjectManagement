clear;
clc;
%************************************
instance_type = 3;
%************************************

solution_type_collection = cell(1,2);
solution_type_collection{1} = 'aro';
solution_type_collection{2} = 'dist';


for instance_index = 1:120
    
    instance_information = load(['DATA/DC' num2str(instance_type) '/MAT_for_Optimization/instance' num2str(instance_index) '.mat']);
    
    dist_information = load(['DATA/DC' num2str(instance_type) '/WorstCase/worst_dist' num2str(instance_index) '.mat']);
    State = dist_information.DELTA_dist_prob';
    probability = dist_information.prob_dist;

    for solution_approach_index = 1:2
        solution_approach = solution_type_collection{solution_approach_index};
        solution_information = load(['DATA/DC' num2str(instance_type) '/WorstCase/worst_' solution_approach num2str(instance_index) '.mat']);
        eval(['x = double(solution_information.x_' solution_approach ');'])
        eval(['y = double(solution_information.y_' solution_approach ');'])
        
        E = size(State,2);
        cost_sample = zeros(E,1);
        
        parfor m = 1 : E
            cost_sample(m) = instance_information.Ux'*x + instance_information.Uy'*y + second_stage_cost(x, y, State(:,m), instance_information.N, instance_information.N_pr, instance_information.Q, instance_information.bar_n, instance_information.Gf, instance_information.Gc, instance_information.d, instance_information.AF, instance_information.ACC, instance_information.AV, instance_information.AC, instance_information.AX, instance_information.AY, instance_information.AZ, instance_information.Z0, instance_information.Z, instance_information.B0);
        end
        
        Temp = [cost_sample,probability];
        Temp = sortrows(Temp,1);
        
        cost_sample = Temp(:,1); probability = Temp(:,2);
        
        save(['DATA/DC' num2str(instance_type) '/WorstSample_DIST/Cost_' solution_approach num2str(instance_index) '.mat'], 'cost_sample', 'probability');
    end
end

