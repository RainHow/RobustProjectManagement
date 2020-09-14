clear;
clc;

%***************************************
instance_type = 2;
number_instance = 240;
number_type = 11;
%***************************************

FontSize = 36;
BigFontSize = 50;

solution_type_collection = cell(1,5);
solution_type_collection{1} = 'aro';
solution_type_collection{2} = 'ccg';
solution_type_collection{3} = 'bd';
solution_type_collection{4} = 'eo';
solution_type_collection{5} = 'um';

P_ARO = zeros(number_instance, number_type);
P_CCG = zeros(number_instance, number_type);
P_BD = zeros(number_instance, number_type);
P_EO = zeros(number_instance, number_type);
P_UM = zeros(number_instance, number_type);

load(['DATA/DC' num2str(instance_type) '/StateTesting.mat']);

for i = 1:240
    
    load(['DATA/DC' num2str(instance_type) '/Instance_Solution/result_aro' num2str(i) '.mat']);
    
    for solution_approach_index = 1:5
        solution_approach = solution_type_collection{solution_approach_index};
        eval(['load(''DATA/DC' num2str(instance_type) '/CostSample/Cost_' solution_approach num2str(i) '.mat'');']);
        eval(['cost_' solution_approach '=cost_sample;']);
        
        cost_sample_prob = [cost_sample, prob];
        cost_sample_prob = sortrows(cost_sample_prob,1);
        cost_sample_prob = [ cost_sample_prob , cumsum(cost_sample_prob(:,2)) ];
        eval(['cost_' solution_approach '_prob=cost_sample_prob;']);

        clear cost_sample cost_sample_prob
    end
    
    p_best = [cost_aro(1); cost_ccg(1); cost_bd(1); cost_eo(1); cost_um(1)];    % should be like this
    p_worst = [cost_aro(end); cost_ccg(end); cost_bd(end); cost_eo(end); cost_um(end)];
    
    p_mean = [sum(cost_aro.*prob); sum(cost_ccg.*prob); sum(cost_bd.*prob); sum(cost_eo.*prob); sum(cost_um.*prob)];
    p_squd = [sum((cost_aro.^2).*prob); sum((cost_ccg.^2).*prob); sum((cost_bd.^2).*prob); sum((cost_eo.^2).*prob); sum((cost_um.^2).*prob)];    
    p_variance = p_squd - p_mean.^2;
      
    for solution_approach_index = 1:5
        solution_approach = solution_type_collection{solution_approach_index};        
        eval([ 'quantile_index_' solution_approach '=max(find(cost_' solution_approach '_prob(:,3)<0.90));' ]);
        eval([ 'quantile_' solution_approach ' = cost_' solution_approach '_prob(quantile_index_' solution_approach '+1:end, [1 2] );' ]);
       
        eval([ 'quantile_' solution_approach '(1,2) = 0.1 - sum(quantile_' solution_approach '(2:end,2));' ])

        eval([ 'quantile_' solution_approach '(:,2) = quantile_' solution_approach '(:,2)/ sum(quantile_' solution_approach '(:,2));' ]);

    end
    
    p_var10 = [cost_aro_prob(quantile_index_aro,1); cost_ccg_prob(quantile_index_ccg,1); cost_bd_prob(quantile_index_bd,1); cost_eo_prob(quantile_index_eo,1); cost_um_prob(quantile_index_um,1)];
    p_cvar10 = [ sum(quantile_aro(:,1).*quantile_aro(:,2)); sum(quantile_ccg(:,1).*quantile_ccg(:,2)); sum(quantile_bd(:,1).*quantile_bd(:,2)); sum(quantile_eo(:,1).*quantile_eo(:,2));sum(quantile_um(:,1).*quantile_um(:,2))];
    clear quantile_*
    
    for solution_approach_index = 1:5
        solution_approach = solution_type_collection{solution_approach_index};        
        eval([ 'quantile_index_' solution_approach '=max(find(cost_' solution_approach '_prob(:,3)<0.95));' ]);
        
        eval([ 'quantile_' solution_approach ' = cost_' solution_approach '_prob(quantile_index_' solution_approach '+1:end, [1 2] );' ]);
        eval([ 'quantile_' solution_approach '(1,2) = 0.05 - sum(quantile_' solution_approach '(2:end,2));' ])

        eval([ 'quantile_' solution_approach '(:,2) = quantile_' solution_approach '(:,2)/ sum(quantile_' solution_approach '(:,2));' ]);

    end
    
    p_var5 = [cost_aro_prob(quantile_index_aro,1); cost_ccg_prob(quantile_index_ccg,1); cost_bd_prob(quantile_index_bd,1); cost_eo_prob(quantile_index_eo,1); cost_um_prob(quantile_index_um,1)] ;
    p_cvar5 = [ sum(quantile_aro(:,1).*quantile_aro(:,2)); sum(quantile_ccg(:,1).*quantile_ccg(:,2)); sum(quantile_bd(:,1).*quantile_bd(:,2)); sum(quantile_eo(:,1).*quantile_eo(:,2));sum(quantile_um(:,1).*quantile_um(:,2))];
    clear quantile_*
    
    for solution_approach_index = 1:5
        solution_approach = solution_type_collection{solution_approach_index};
        eval([ 'quantile_index_' solution_approach '=min(find(cost_' solution_approach '_prob(:,1)>target));' ]); 
        eval([ 'quantile_' solution_approach ' = cost_' solution_approach '_prob(quantile_index_' solution_approach ':end, [1 2] );' ]);
        eval([ 'quantile_' solution_approach ' =[ quantile_' solution_approach ', quantile_' solution_approach '(:,2)/ sum(quantile_' solution_approach '(:,2))];' ]);
        eval([ 'quantile_' solution_approach '(:,1) = quantile_' solution_approach '(:,1) - target;' ])

    end
    
    p_vp = [ sum(quantile_aro(:,2));sum(quantile_ccg(:,2));sum(quantile_bd(:,2));sum(quantile_eo(:,2));sum(quantile_um(:,2)); ];
    p_el = [ sum(quantile_aro(:,1).*quantile_aro(:,2)); sum(quantile_ccg(:,1).*quantile_ccg(:,2)); sum(quantile_bd(:,1).*quantile_bd(:,2)); sum(quantile_eo(:,1).*quantile_eo(:,2));sum(quantile_um(:,1).*quantile_um(:,2))];
    p_cel = [ sum(quantile_aro(:,1).*quantile_aro(:,3)); sum(quantile_ccg(:,1).*quantile_ccg(:,3)); sum(quantile_bd(:,1).*quantile_bd(:,3)); sum(quantile_eo(:,1).*quantile_eo(:,3));sum(quantile_um(:,1).*quantile_um(:,3))];
    clear quantile_*
      
    
    p_best = p_best/p_best(1);
    p_worst = p_worst/p_worst(1);
    p_mean = p_mean/p_mean(1);
    p_variance = p_variance/p_variance(1);
    p_var10 = p_var10/p_var10(1);
    p_cvar10 = p_cvar10/p_cvar10(1);
    p_var5 = p_var5/p_var5(1);
    p_cvar5 = p_cvar5/p_cvar5(1);
    p_vp = p_vp/p_vp(1);
    p_el = p_el/p_el(1);
    p_cel = p_cel/p_cel(1);
    
    performance_num = [p_best p_worst p_mean p_variance p_var10 p_cvar10 p_var5 p_cvar5 p_vp p_el p_cel];
    
    P_ARO(i,:) = performance_num(1,:);
    P_CCG(i,:) = performance_num(2,:);
    P_BD(i,:) = performance_num(3,:);
    P_EO(i,:) = performance_num(4,:);
    P_UM(i,:) = performance_num(5,:);
 
end

P_CCG_1best = P_CCG(:,1);
P_CCG_2worst = P_CCG(:,2);
P_CCG_3mean = P_CCG(:,3);
P_CCG_4variance = P_CCG(:,4);
P_CCG_5var10 = P_CCG(:,5);
P_CCG_6cvar10 = P_CCG(:,6);
P_CCG_7var5 = P_CCG(:,7);
P_CCG_8cvar5 = P_CCG(:,8);
P_CCG_9vp = P_CCG(:,9);
P_CCG_10el = P_CCG(:,10);
P_CCG_11cel = P_CCG(:,11);
P_CCG_aggregated = mean(P_CCG,1);

P_BD_1best = P_BD(:,1);
P_BD_2worst = P_BD(:,2);
P_BD_3mean = P_BD(:,3);
P_BD_4variance = P_BD(:,4);
P_BD_5var10 = P_BD(:,5);
P_BD_6cvar10 = P_BD(:,6);
P_BD_7var5 = P_BD(:,7);
P_BD_8cvar5 = P_BD(:,8);
P_BD_9vp = P_BD(:,9);
P_BD_10el = P_BD(:,10);
P_BD_11cel = P_BD(:,11);
P_BD_aggregated = mean(P_BD,1);

P_EO_1best = P_EO(:,1);
P_EO_2worst = P_EO(:,2);
P_EO_3mean = P_EO(:,3);
P_EO_4variance = P_EO(:,4);
P_EO_5var10 = P_EO(:,5);
P_EO_6cvar10 = P_EO(:,6);
P_EO_7var5 = P_EO(:,7);
P_EO_8cvar5 = P_EO(:,8);
P_EO_9vp = P_EO(:,9);
P_EO_10el = P_EO(:,10);
P_EO_11cel = P_EO(:,11);
P_EO_aggregated = mean(P_EO,1);

P_UM_1best = P_UM(:,1);
P_UM_2worst = P_UM(:,2);
P_UM_3mean = P_UM(:,3);
P_UM_4variance = P_UM(:,4);
P_UM_5var10 = P_UM(:,5);
P_UM_6cvar10 = P_UM(:,6);
P_UM_7var5 = P_UM(:,7);
P_UM_8cvar5 = P_UM(:,8);
P_UM_9vp = P_UM(:,9);
P_UM_10el = P_UM(:,10);
P_UM_11cel = P_UM(:,11);
P_UM_aggregated = mean(P_UM,1);

eval(['save(''DATA/DC' num2str(instance_type) '/CostAnalysis_withTarget.mat'', ''P_*'');']);