function [ obj, x, y, DELTA, Iteration ] = EO( A_xy, b_xy, N, N_pr, M, K, Q, bar_n, Ux, Uy, Gf, Gc, d, AF, ACC, AV, AC, AX, AY, AZ, B0, Z, Z0, Mu, Sigma, Cor, L, REALIZATION )
%ARO Summary of this function goes here
%   Solve the project planning problem using adaptive robust optimization

limit_extreme_points = 1000;
limit_time = 800;

% put all realizations into a column vector
Realization = zeros(sum(L),1);
Abs = Realization;
for k = 1 : K
    Realization(sum(L(1:k-1))+1:sum(L(1:k))) = REALIZATION(k,1:L(k))';
    Abs(sum(L(1:k-1))+1:sum(L(1:k))) = (REALIZATION(k,1:L(k))' - Mu(k)*ones(L(k),1))/sqrt(Sigma(k));
end
Chi = Realization.^2;

% Algorithm BD
stop_criterion = 0;
DELTA = REALIZATION(:,1)'; % each row is a (delta^i)', size: E*K

while stop_criterion == 0 && size(DELTA,1)<=limit_extreme_points
    
    [ obj, x, y, r0, r, vartheta, rc ] = Relaxation_Solution_EO ( A_xy, b_xy, N, N_pr, M, K, Q, bar_n, Ux, Uy, Gf, Gc, d, AF, ACC, AV, AC, AX, AY, AZ, B0, Z, Z0, Mu, Sigma, Cor, DELTA, limit_time );
       
    % **************** CHECK THE CONSTRAINT ****************
    [ obj_bilinear, Delta ] = Bilinear_EO( x, y, r, vartheta, rc, N, N_pr, K, bar_n, Q, Gf, Gc, d, AF, ACC, AV, AC, AX, AY, AZ, B0, Z, Z0, L, Realization, Chi, Abs, limit_time );
    
    obj_bilinear
    RHS = r0 - Ux'*x - Uy'*y + 1e-5
    
    c_2 = second_stage_cost( x, y, Delta, N, N_pr, Q, bar_n, Gf, Gc, d, AF, ACC, AV, AC, AX, AY, AZ, Z0, Z, B0);
    
    if obj_bilinear > r0 - Ux'*x - Uy'*y + 1e-5
        if ismember(Delta',DELTA,'rows')
            error('check')
        end
        DELTA = [DELTA; Delta']
    else
        stop_criterion = 1;
    end
    
end


Iteration = size(DELTA,1);

if size(DELTA,1) >= limit_extreme_points
    obj = NaN;
    x = NaN;
    y = NaN;
end

end


function [ obj, x_sol, y_sol, r0_sol, r_sol, vartheta_sol, rc_sol ] = Relaxation_Solution_EO ( A_xy, b_xy, N, N_pr, M, K, Q, bar_n, Ux, Uy, Gf, Gc, d, AF, ACC, AV, AC, AX, AY, AZ, B0, Z, Z0, Mu, Sigma, Cor, DELTA, limit_time )

E = size(DELTA,1);
CHI = DELTA.^2;
ABS = abs( sum((DELTA - ones(E,1)*Mu')./(ones(E,1)*sqrt(Sigma)'),2));

rome_begin;

h = rome_model('BD_subproblem');

x = newvar(M,1,'binary');
y = newvar(N,1,'binary');
r0 = newvar(1,1);
r = newvar(K,1);
rc = newvar(1,1,'nonneg');       % correlation
vartheta = newvar(K,1,'nonneg');
Wf = newvar(N_pr*E,1,'nonneg');     % put all w_f^i into one column vector
Wc = newvar(bar_n*E,1,'nonneg');    % put all w_c^i into one column vector
V = newvar(Q*E,1,'nonneg');         % put all v^i into one column vector
C = newvar(N*E,1,'nonneg');         % put all c^i into one column vector

rome_minimize( r0 + Mu'*r + Sigma'*vartheta + Cor*rc);

rome_constraint( r0 + Mu'*r + Sigma'*vartheta + Cor*rc >= -10000);    % maybe unbounded

for i = 1 : E
    rome_constraint( Ux'*x + Uy'*y + Gf'*Wf((i-1)*N_pr+1:i*N_pr) + Gc'*Wc((i-1)*bar_n+1:i*bar_n) + d'*V((i-1)*Q+1:i*Q) - DELTA(i,:)*r - CHI(i,:)*vartheta - ABS(i,:)*rc <= r0 );
    rome_constraint( AF*Wf((i-1)*N_pr+1:i*N_pr) + ACC*Wc((i-1)*bar_n+1:i*bar_n) + AV*V((i-1)*Q+1:i*Q) + AC*C((i-1)*N+1:i*N) >= AX*x + AY*y + AZ*(Z0+Z*(DELTA(i,:))') + B0 );
end

rome_constraint( A_xy * [x;y] <= b_xy);

relaxation_t_start = tic;
h.solve('CPLEX',limit_time);
relaxation_solution_time = toc(relaxation_t_start);
if relaxation_solution_time>=limit_time-1
    error('check')
end


x_sol = h.eval(x);
y_sol = h.eval(y);
r0_sol = h.eval(r0);
r_sol = h.eval(r);
vartheta_sol = h.eval(vartheta);
rc_sol = h.eval(rc);

obj = h.objective;

end

function [ obj_bilinear, Delta ] = Bilinear_EO( x, y, r, vartheta, rc, N, N_pr, K, bar_n, Q, Gf, Gc, d, AF, ACC, AV, AC, AX, AY, AZ, B0, Z, Z0, L, Realization, Chi, Abs, limit_time )

M = 1.0e+6;

% construct the vector R and Vartheta
% R: column vector with size sum(L)*1, such that the first L_1 elements are all with value r_1, the next L_2 elements are all with value r_2...
% Vartheta: defined similarly to R
R = zeros(sum(L),1);
Vartheta = zeros(sum(L),1);
for k = 1 : K
    R(sum(L(1:k-1))+1:sum(L(1:k))) = r(k)*ones(L(k),1);
    Vartheta(sum(L(1:k-1))+1:sum(L(1:k))) = vartheta(k)*ones(L(k),1);
end


B = (AZ*Z)';
c = AX*x + AY*y + AZ*Z0 + B0;

rome_begin;
h = rome_model('Bilinear_subproblem');

p = newvar(size(AF,1),1,'nonneg');
Phi = newvar(sum(L),1,'binary');
Zeta = newvar(sum(L),1);
aux = newvar(1,1);      % represent abs()

rome_maximize( c'*p + sum(Zeta) - (R.*Realization)'*Phi - (Vartheta.*Chi)'*Phi - rc*aux);

rome_constraint( c'*p + sum(Zeta) - (R.*Realization)'*Phi - (Vartheta.*Chi)'*Phi -rc*aux <= 1000000000 ); % prevent the unbounded case
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
    error('something wrong')
end

% use the exact Delta, solve the problem again to get the exact value
Chosen_realization = find(Phi_sol);
Delta = Realization(Chosen_realization);
chi_current = Chi(Chosen_realization);
abs_current = Abs(Chosen_realization);


c_2 = second_stage_cost( x, y, Delta, N, N_pr, Q, bar_n, Gf, Gc, d, AF, ACC, AV, AC, AX, AY, AZ, Z0, Z, B0);

obj_bilinear = c_2-r'*Delta-vartheta'*chi_current-rc*abs(sum(abs_current));

end