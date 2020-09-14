function [ obj_dist, r, r0, vartheta, rc, h_value, cost, E] = ARO_given_x_y(A_xy, b_xy, x_dist,y_dist, N, N_pr, M, K, Q, bar_n, Ux, Uy, Gf, Gc, d, AF, ACC, AV, AC, AX, AY, AZ, B0, Z, Z0, Mu, Sigma,Cor, L, State, Am,bm,target )
% Adding corrleation matrix

E = size(State,1);
h_value = zeros(E,1);
cost = zeros(E,1);

for i =1:E
    Delta = State(i,:)';
    h_value(i) = Ux'*x_dist+Uy'*y_dist + second_stage_cost( x_dist, y_dist, Delta, N, N_pr, Q, bar_n, Gf, Gc, d, AF, ACC, AV, AC, AX, AY, AZ, Z0, Z, B0);
end

[ obj_dist, r0, r, vartheta, rc ] = Relaxation_Solution_aro (A_xy, b_xy, x_dist, y_dist, N, N_pr, M, K, Q, bar_n, Ux, Uy, Gf, Gc, d, AF, ACC, AV, AC, AX, AY, AZ, B0, Z, Z0, Mu, Sigma, Cor, State, E, h_value, Am, bm, target);

for j =1:E
    cost(i) = max( Am(1)*h_value(j)-Am(1)*target +bm(1)*obj_dist, Am(2)*h_value(j)-Am(2)*target +bm(2)*obj_dist   );
end

end


function [ obj, r0_sol, r_sol, vartheta_sol, rc_sol] = Relaxation_Solution_aro (A_xy, b_xy,x, y, N, N_pr, M, K, Q, bar_n, Ux, Uy, Gf, Gc, d, AF, ACC, AV, AC, AX, AY, AZ, B0, Z, Z0, Mu, Sigma,Cor, DELTA, E, h_value, Am, bm, target )

CHI = DELTA.^2;
ABS = abs( sum((DELTA - ones(E,1)*Mu')./(ones(E,1)*sqrt(Sigma)'),2) );

rome_begin;
h = rome_model('ARO_subproblem');

alpha = newvar(1,1);
r0 = newvar(1,1);    % support
r = newvar(K,1);     % mean
rc = newvar(1,1,'nonneg');    % correlation
vartheta = newvar(K,1,'nonneg');    % second moment

rome_minimize(alpha);
rome_constraint(alpha >= 0 );     % important
rome_constraint(r0 + Mu'*r + Sigma'*vartheta + Cor*rc <= 0);
for i = 1 : E
    rome_constraint( Am(1)*h_value(i)-Am(1)*target+bm(1)*alpha <= r0 + DELTA(i,:)*r + CHI(i,:)*vartheta + ABS(i,:)*rc);
    rome_constraint( Am(2)*h_value(i)-Am(2)*target+bm(2)*alpha <= r0 + DELTA(i,:)*r + CHI(i,:)*vartheta + ABS(i,:)*rc);

end

relaxation_t_start = tic;
h.solve('CPLEX');
relaxation_solution_time = toc(relaxation_t_start);

r0_sol = h.eval(r0);
r_sol = h.eval(r);
vartheta_sol = h.eval(vartheta);
rc_sol = h.eval(rc);

obj = h.objective;

end


function [ obj_bilinear, Delta, binary_solution_time  ] = Bilinear_aro( x, y, r, vartheta,rc, N, N_pr, K, bar_n, Q, Gf, Gc, d, AF, ACC, AV, AC, AX, AY, AZ, B0, Z, Z0, L, Realization, Chi, Abs, limit_time, am)

M = 1.0e+7;

% construct the vector R and Vartheta
% R: column vector with size sum(L)*1, such that the first L_1 elements are all with value r_1, the next L_2 elements are all with value r_2...
% Vartheta: defined similarly to R
R = zeros(sum(L),1);
Vartheta = zeros(sum(L),1);
%Rc = zeros(sum(L),1);

for k = 1 : K
    R(sum(L(1:k-1))+1:sum(L(1:k))) = r(k)*ones(L(k),1);
    Vartheta(sum(L(1:k-1))+1:sum(L(1:k))) = vartheta(k)*ones(L(k),1);
    %    Rc(sum(L(1:k-1))+1:sum(L(1:k))) = rc*ones(L(k),1);
end


B = (AZ*Z)';
c = AX*x + AY*y + AZ*Z0 + B0;

rome_begin;
h = rome_model('Bilinear_subproblem');

p = newvar(size(AF,1),1,'nonneg');
Phi = newvar(sum(L),1,'binary');
Zeta = newvar(sum(L),1);
aux = newvar(1,1);    % represent abs()

rome_maximize( am*(c'*p + sum(Zeta)) - (R.*Realization)'*Phi - (Vartheta.*Chi)'*Phi - rc*aux);

rome_constraint( c'*p + sum(Zeta) - (R.*Realization)'*Phi - (Vartheta.*Chi)'*Phi - rc*aux <= 1000000000 ); % prevent the unbounded case
rome_constraint( aux >= Abs'*Phi );
rome_constraint( aux >= -Abs'*Phi);

rome_constraint( [AF ACC AV AC]'*p <= [Gf; Gc; d; zeros(N,1)] );
for k = 1 : K
    Index_k = sum(L(1:k-1))+1 : sum(L(1:k));
    rome_constraint( sum ( Phi(Index_k) ) == 1 );
    rome_constraint( Zeta(Index_k) <= (Realization(Index_k)*B(k,:))*p + M*(ones(L(k),1)-Phi(Index_k)) );
    rome_constraint( Zeta(Index_k) <= M*Phi(Index_k) );
end

binary_t_start = tic;
h.solve('CPLEX',limit_time);
binary_solution_time = toc(binary_t_start);
if binary_solution_time>=limit_time-1
    error('check')
end

% round it to avoid the computational error
Phi_sol = round(h.eval(Phi));
p_sol = h.eval(p);
% if error is too much, there must be some problem
if sum(abs(Phi_sol-h.eval(Phi))) > 1e-4
    error('check')
end

% use the exact Delta, solve the problem again to get the exact value
Chosen_realization = find(Phi_sol);
Delta = Realization(Chosen_realization);
chi_current = Chi(Chosen_realization);
abs_current = Abs(Chosen_realization);

sum(abs_current)

c_2 = second_stage_cost( x, y, Delta, N, N_pr, Q, bar_n, Gf, Gc, d, AF, ACC, AV, AC, AX, AY, AZ, Z0, Z, B0);

obj_bilinear = am*c_2-r'*Delta-vartheta'*chi_current-rc*abs(sum(abs_current));

end



