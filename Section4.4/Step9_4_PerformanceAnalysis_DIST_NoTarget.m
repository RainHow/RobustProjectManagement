clear;
clc;

instance_type = 3;
number_instance = 120;
number_type = 11;

FontSize = 36;
BigFontSize = 50;

solution_type_collection = cell(1,2);
solution_type_collection{1} = 'aro';
solution_type_collection{2} = 'dist';


P_ARO = zeros(number_instance,number_type);
P_DIST = zeros(number_instance,number_type);

for i = 1:number_instance
    load(['DATA/DC' num2str(instance_type) '/Instance_Solution/result_aro' num2str(i) '.mat']);
    
    for solution_approach_index = 1:2
        solution_approach = solution_type_collection{solution_approach_index};
        eval(['load(''DATA/DC' num2str(instance_type) '/WorstSample_DIST/Cost_' solution_approach num2str(i) '.mat'');']);
        eval(['cost_' solution_approach '=cost_sample;']);
        eval(['prob_' solution_approach '=probability;']);
        
        cost_sample_prob = [cost_sample, probability];
        cost_sample_prob = sortrows(cost_sample_prob,1);
        cost_sample_prob = [ cost_sample_prob , cumsum(cost_sample_prob(:,2)) ];
        eval(['cost_' solution_approach '_prob=cost_sample_prob;']);

        clear cost_sample cost_sample_prob
    end
    
    p_best = [cost_aro(1); cost_dist(1)];    % should be like this
    p_worst = [cost_aro(end); cost_dist(end)];
    
    p_mean = [sum(cost_aro.*prob_aro); sum(cost_dist.*prob_dist)];
    p_squd = [sum((cost_aro.^2).*prob_aro); sum((cost_dist.^2).*prob_dist) ];    
    p_variance = p_squd - p_mean.^2;
      
    for solution_approach_index = 1:2
        solution_approach = solution_type_collection{solution_approach_index};        
        eval([ 'quantile_index_' solution_approach '=max(find(cost_' solution_approach '_prob(:,3)<0.90));' ]);
        
        eval([ 'quantile_' solution_approach ' = cost_' solution_approach '_prob(quantile_index_' solution_approach '+1:end, [1 2] );' ]);
        
        eval([ 'quantile_' solution_approach '(1,2) = 0.1 - sum(quantile_' solution_approach '(2:end,2));' ])

        eval([ 'quantile_' solution_approach '(:,2) = quantile_' solution_approach '(:,2)/ sum(quantile_' solution_approach '(:,2));' ]);

    end
    
    p_var10 = [cost_aro_prob(quantile_index_aro+1,1); cost_dist_prob(quantile_index_dist+1,1) ];
    p_cvar10 = [ sum(quantile_aro(:,1).*quantile_aro(:,2)); sum(quantile_dist(:,1).*quantile_dist(:,2)) ];
    clear quantile_*
    
    for solution_approach_index = 1:2
        solution_approach = solution_type_collection{solution_approach_index};        
        eval([ 'quantile_index_' solution_approach '=max(find(cost_' solution_approach '_prob(:,3)<0.95));' ]);
        
        eval([ 'quantile_' solution_approach ' = cost_' solution_approach '_prob(quantile_index_' solution_approach '+1:end, [1 2] );' ]);
        
        eval([ 'quantile_' solution_approach '(1,2) = 0.05 - sum(quantile_' solution_approach '(2:end,2));' ])

        eval([ 'quantile_' solution_approach '(:,2) = quantile_' solution_approach '(:,2)/ sum(quantile_' solution_approach '(:,2));' ]);

    end
    
    p_var5 = [cost_aro_prob(quantile_index_aro+1,1); cost_dist_prob(quantile_index_dist+1,1)] ;
    p_cvar5 = [ sum(quantile_aro(:,1).*quantile_aro(:,2)); sum(quantile_dist(:,1).*quantile_dist(:,2)) ];
    clear quantile_*
    
    for solution_approach_index = 1:2
        solution_approach = solution_type_collection{solution_approach_index};
        eval([ 'quantile_index_' solution_approach '=min(find(cost_' solution_approach '_prob(:,1)>target));' ]); 
        eval([ 'quantile_' solution_approach ' = cost_' solution_approach '_prob(quantile_index_' solution_approach ':end, [1 2] );' ]);
        eval([ 'quantile_' solution_approach ' =[ quantile_' solution_approach ', quantile_' solution_approach '(:,2)/ sum(quantile_' solution_approach '(:,2))];' ]);
        eval([ 'quantile_' solution_approach '(:,1) = quantile_' solution_approach '(:,1) - target;' ])

    end
    
    
    p_vp = [ sum(quantile_aro(:,2));sum(quantile_dist(:,2)) ];
    p_el = [ sum(quantile_aro(:,1).*quantile_aro(:,2)); sum(quantile_dist(:,1).*quantile_dist(:,2)) ];
    p_cel = [ sum(quantile_aro(:,1).*quantile_aro(:,3)); sum(quantile_dist(:,1).*quantile_dist(:,3)); ];
    clear quantile_*
    
    
    
%     p_best = p_best/p_best(1);
%     p_worst = p_worst/p_worst(1);
%     p_mean = p_mean/p_mean(1);
%     p_variance = p_variance/p_variance(1);
%     p_var10 = p_var10/p_var10(1);
%     p_cvar10 = p_cvar10/p_cvar10(1);
%     p_var5 = p_var5/p_var5(1);
%     p_cvar5 = p_cvar5/p_cvar5(1);
%     p_vp = p_vp/p_vp(1);
%     p_el = p_el/p_el(1);
%     p_cel = p_cel/p_cel(1);
    
    performance_num = [p_best p_worst p_mean p_variance p_var10 p_cvar10 p_var5 p_cvar5 p_vp p_el p_cel];
    
    P_ARO(i,:) = performance_num(1,:);
    P_DIST(i,:) = performance_num(2,:);
    
    
end

P_DIST_1best = P_DIST(:,1);
P_DIST_2worst = P_DIST(:,2);
P_DIST_3mean = P_DIST(:,3);
P_DIST_4variance = P_DIST(:,4);
P_DIST_5var10 = P_DIST(:,5);
P_DIST_6cvar10 = P_DIST(:,6);
P_DIST_7var5 = P_DIST(:,7);
P_DIST_8cvar5 = P_DIST(:,8);
P_DIST_9vp = P_DIST(:,9);
P_DIST_10el = P_DIST(:,10);
P_DIST_11cel = P_DIST(:,11);
P_DIST_aggregated = mean(P_DIST,1);

P_ARO_1best = P_ARO(:,1);
P_ARO_2worst = P_ARO(:,2);
P_ARO_3mean = P_ARO(:,3);
P_ARO_4variance = P_ARO(:,4);
P_ARO_5var10 = P_ARO(:,5);
P_ARO_6cvar10 = P_ARO(:,6);
P_ARO_7var5 = P_ARO(:,7);
P_ARO_8cvar5 = P_ARO(:,8);
P_ARO_9vp = P_ARO(:,9);
P_ARO_10el = P_ARO(:,10);
P_ARO_11cel = P_ARO(:,11);
P_ARO_aggregated = mean(P_ARO,1);

eval(['save(''DATA/DC' num2str(instance_type) '/WorstCostAnalysis_DIST_NoTarget.mat'', ''P_*'');']);