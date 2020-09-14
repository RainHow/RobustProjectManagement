
clear;
clc;

load(['DATA/DC2/MAT_for_Optimization/instance1.mat']);
load(['DATA/DC2/Instance_Solution/result_aro1.mat']);
load([ 'DATA/DC2/StateTesting.mat' ])

%%
REALIZATION = [0,6];
[X1,X2,X3,X4,X5] = ndgrid( REALIZATION, REALIZATION, REALIZATION, REALIZATION, REALIZATION);
DELTA_aro = [X1(:) X2(:) X3(:) X4(:) X5(:)];

%%
%  methed A    (using State)
E_a = size(State,1);
cost_a = zeros(E_a,1);
for i = 1:E_a
    Delta = State(i,:)';
    temp = Ux'*x_aro+Uy'*y_aro + second_stage_cost( x_aro, y_aro, Delta, N, N_pr, Q, bar_n, Gf, Gc, d, AF, ACC, AV, AC, AX, AY, AZ, Z0, Z, B0);
    cost_a(i) = max( Am(1)*temp-Am(1)*target +bm(1)*obj_aro, Am(2)*temp-Am(2)*target +bm(2)*obj_aro   );
end

ABS_a = abs( sum((State - ones(E_a,1)*Mu')./(ones(E_a,1)*sqrt(Sigma)'),2) )

rome_begin;
h = rome_model('Worstcase_Dist');

p = newvar(E_a,1,'nonneg');

rome_minimize(p'*ones(E_a,1));

rome_constraint( p'*cost_a == r0_aro + Mu'*r_aro + Sigma'*vartheta_aro + Cor*rc_aro );     % important
rome_constraint( sum(p) == 1);
rome_constraint( State'*p == Mu );
rome_constraint( (State.^2)'*p <= Sigma  );
rome_constraint( ABS_a'*p <= Cor );

relaxation_t_start = tic;
h.solve('CPLEX');
relaxation_solution_time = toc(relaxation_t_start);

p_sol_a = h.eval(p)

%%
% method b
E_b = size(DELTA_aro,1);
cost_b = zeros(E_b,1);
for i =1:E_b
    Delta = DELTA_aro(i,:)';
    temp = Ux'*x_aro+Uy'*y_aro + second_stage_cost( x_aro, y_aro, Delta, N, N_pr, Q, bar_n, Gf, Gc, d, AF, ACC, AV, AC, AX, AY, AZ, Z0, Z, B0);
    cost_b(i) = max( Am(1)*temp-Am(1)*target +bm(1)*obj_aro, Am(2)*temp-Am(2)*target +bm(2)*obj_aro   );
end

ABS_b = abs( sum((DELTA_aro - ones(E_b,1)*Mu')./(ones(E_b,1)*sqrt(Sigma)'),2) )

rome_begin;
h = rome_model('Worstcase_Dist');

p = newvar(E_b,1,'nonneg');

rome_minimize(p'*ones(E_b,1));

rome_constraint( p'*cost_b == r0_aro + Mu'*r_aro + Sigma'*vartheta_aro + Cor*rc_aro );     % important
rome_constraint( sum(p) == 1);
rome_constraint( DELTA_aro'*p == Mu );
rome_constraint( (DELTA_aro.^2)'*p <= Sigma  );
rome_constraint( ABS_b'*p <= Cor );

relaxation_t_start = tic;
h.solve('CPLEX');
relaxation_solution_time = toc(relaxation_t_start);

p_sol_b = h.eval(p)

%%
index_a = find(p_sol_a>0);
State(index_a,:)

index_b = find(p_sol_b>0);
DELTA_aro(index_b,:)
