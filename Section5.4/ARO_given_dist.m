function [ obj_sol, x_sol, y_sol, dist_solution_time  ] = ARO_dist(A_xy, b_xy, N, N_pr, M, K, Q, bar_n, Ux, Uy, Gf, Gc, d, AF, ACC, AV, AC, AX, AY, AZ, B0, Z, Z0, Am,bm,target, State, prob )
% Adding corrleation matrix

E = size(State,1);
limit_time = 600;

rome_begin;
h = rome_model('ARO_given_dist');

alpha = newvar(1,1,'nonneg');
x = newvar(M,1,'binary');
y = newvar(N,1,'binary');
aux = newvar(E,1);

Wf = newvar(N_pr*E,1,'nonneg');     % put all w_f^i into one column vector
Wc = newvar(bar_n*E,1,'nonneg');    % put all w_c^i into one column vector
V = newvar(Q*E,1,'nonneg');         % put all v^i into one column vector
C = newvar(N*E,1,'nonneg');         % put all c^i into one column vector

rome_minimize(alpha);

rome_constraint( prob'*aux <= 0 );
for i = 1 : E
    rome_constraint( aux(i) >= Am(1)*( Ux'*x+Uy'*y + Gf'*Wf((i-1)*N_pr+1:i*N_pr)+Gc'*Wc((i-1)*bar_n+1:i*bar_n) + d'*V((i-1)*Q+1:i*Q) -target) + bm(1)*alpha );
    rome_constraint( aux(i) >= Am(1)*( -target) + bm(1)*alpha );
    rome_constraint( AF*Wf((i-1)*N_pr+1:i*N_pr) + ACC*Wc((i-1)*bar_n+1:i*bar_n) + AV*V((i-1)*Q+1:i*Q) + AC*C((i-1)*N+1:i*N) >= AX*x + AY*y + AZ*(Z0+Z*(State(i,:))') + B0 );
end

rome_constraint( A_xy * [x;y] <= b_xy);

dist_t_start = tic;
h.solve('CPLEX',limit_time);
dist_solution_time = toc(dist_t_start);


obj_sol = h.objective;
x_sol = h.eval(x);
y_sol = h.eval(y);


end




