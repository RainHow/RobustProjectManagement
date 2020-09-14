
clear;
clc;

%*********************************
instance_type = 3;
%*********************************

for i = 121:240
     
    outputfile = ['DATA/DC' num2str(instance_type) '/WorstCase/worst_dist' num2str(i) '.mat'];
    
    load(['DATA/DC' num2str(instance_type) '/MAT_for_Optimization/instance' num2str(i) '.mat']);
    load(['DATA/DC' num2str(instance_type) '/Instance_Solution/result_dist' num2str(i) '.mat']);
    load(['DATA/DC' num2str(instance_type) '/StateTesting.mat']);
    
    [ obj_dist, r_dist, r0_dist, vartheta_dist, rc_dist, h_value, cost, E] = ARO_given_x_y(A_xy, b_xy, x_dist,y_dist, N, N_pr, M, K, Q, bar_n, Ux, Uy, Gf, Gc, d, AF, ACC, AV, AC, AX, AY, AZ, B0, Z, Z0, Mu, Sigma,Cor, L, State, Am,bm,target );
    
    ABS = abs( sum((State - ones(E,1)*Mu')./(ones(E,1)*sqrt(Sigma)'),2) );
    
    rome_begin;
    h = rome_model('Worstcase_Dist');
    p = newvar(E,1,'nonneg');
    rome_minimize(p'*ones(E,1));
    rome_constraint( p'*cost == r0_dist + Mu'*r_dist + Sigma'*vartheta_dist + Cor*rc_dist );     % important
    rome_constraint( sum(p) == 1);
    rome_constraint( State'*p == Mu );
    rome_constraint( (State.^2)'*p <= Sigma  );
    rome_constraint( ABS'*p <= Cor );
    h.solve('CPLEX');
    
    p_sol = h.eval(p);
    
    index = find(p_sol>0);
    DELTA_dist_prob = State(index,:);
    prob_dist = p_sol(index);
    
    save(outputfile, 'obj_dist', 'x_dist', 'y_dist', 'DELTA_dist_prob','prob_dist');
    
end
