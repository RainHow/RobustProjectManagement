
clear;
clc;

%*********************************
instance_type = 3;
%*********************************

for i = 121:240
    
    outputfile = ['DATA/DC' num2str(instance_type) '/WorstCase/worst_aro' num2str(i) '.mat'];
    load(['DATA/DC' num2str(instance_type) '/MAT_for_Optimization/instance' num2str(i) '.mat']);   
    load(['DATA/DC' num2str(instance_type) '/Instance_Solution/result_aro' num2str(i) '.mat']);
    
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
    
    save(outputfile, 'obj_aro', 'x_aro', 'y_aro', 'DELTA_aro_prob', 'prob_aro');
    
end
