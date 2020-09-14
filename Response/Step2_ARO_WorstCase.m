% Solving problem by ARO algorithm
clear;
clc;

ratio = 1.2;
instance_type = 1;
Ratio = [0.2, 1, 10, 10000];

for i =1:5
    load(['DATA/DC' num2str(instance_type) '/MAT_for_Optimization/instance1.mat']);
    load(['DATA/DC' num2str(instance_type) '/Instance_Solution/result_eo1.mat']);
    target = obj_eo * ratio;

    for j = 1:4
        
        Cor_temp = Cor * Ratio(j);
        
        outputfile = ['DATA/DC' num2str(instance_type) '/Instance_Solution/result_aro' num2str(i) '_' num2str(j) '.mat'];       
        errorfile = ['DATA/DC' num2str(instance_type) '/Error/DC1Error1.mat'];        
        
        whole_solution_t_start = tic;
        [ obj_aro, x_aro, y_aro, DELTA_aro, r_aro, r0_aro, vartheta_aro, rc_aro, Iteration_aro, T_master_aro, T_bilinear_aro ] = ARO(A_xy, b_xy, N, N_pr, M, K, Q, bar_n, Ux, Uy, Gf, Gc, d, AF, ACC, AV, AC, AX, AY, AZ, B0, Z, Z0, Mu, Sigma,Cor_temp, L, REALIZATION,Am,bm,target, errorfile );
        T_aro = toc(whole_solution_t_start);
        
        E = size(DELTA_aro,1);
        cost = zeros(E,1);
        for k =1:E
            Delta = DELTA_aro(k,:)';
            temp = Ux'*x_aro+Uy'*y_aro + second_stage_cost( x_aro, y_aro, Delta, N, N_pr, Q, bar_n, Gf, Gc, d, AF, ACC, AV, AC, AX, AY, AZ, Z0, Z, B0);
            cost(k) = max( Am(1)*temp-Am(1)*target +bm(1)*obj_aro, Am(2)*temp-Am(2)*target +bm(2)*obj_aro   );
        end
        
        ABS = abs( sum((DELTA_aro - ones(E,1)*Mu')./(ones(E,1)*sqrt(Sigma)'),2) );
        
        rome_begin;
        h = rome_model('Worstcase_Dist');
        p = newvar(E,1,'nonneg');
        rome_minimize(p'*ones(E,1));
        rome_constraint( p'*cost == r0_aro + Mu'*r_aro + Sigma'*vartheta_aro + Cor*rc_aro );     % important
        rome_constraint( sum(p) == 1);
        rome_constraint( DELTA_aro'*p == Mu );
        rome_constraint( (DELTA_aro.^2)'*p <= Sigma  );
        rome_constraint( ABS'*p <= Cor );
        h.solve('CPLEX');
        
        p_sol = h.eval(p);
        
        index = find(p_sol>0);
        DELTA_aro_prob = DELTA_aro(index,:);
        prob_aro = p_sol(index);
        
        save(outputfile, 'obj_aro', 'x_aro', 'y_aro', 'Cor_temp', 'DELTA_aro', 'p_sol', 'DELTA_aro_prob','prob_aro');
               
    end
end
