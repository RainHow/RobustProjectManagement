% Searching Sampling Solution by ARO algorithm

clear;
clc;

%*********************************
instance_type = 3;
%*********************************
ratio = 1.2;

for i = 134:240
    
    outputfile = ['DATA/DC' num2str(instance_type) '/Instance_Solution/result_dist' num2str(i) '.mat'];
    
    load(['DATA/DC' num2str(instance_type) '/MAT_for_Optimization/instance' num2str(i) '.mat']);
    load(['DATA/DC' num2str(instance_type) '/Instance_Solution/result_eo' num2str(i) '.mat']);
    load(['DATA/DC' num2str(instance_type) '/StateTesting.mat'])
    
    target = obj_eo * ratio;
    
    whole_solution_t_start = tic;
    [ obj_dist, x_dist, y_dist, Time  ] = ARO_given_dist(A_xy, b_xy, N, N_pr, M, K, Q, bar_n, Ux, Uy, Gf, Gc, d, AF, ACC, AV, AC, AX, AY, AZ, B0, Z, Z0, Am,bm,target, State, prob );
    T_dist = toc(whole_solution_t_start);
    
    save(outputfile, 'x_dist', 'y_dist', 'T_dist', 'target');
    clear obj_eo target
end
