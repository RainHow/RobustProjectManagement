function [ obj, x, y, DELTA, r, r0, vartheta, rc, Iteration, master_time, bilinear_time ] = ARO(A_xy, b_xy, N, N_pr, M, K, Q, bar_n, Ux, Uy, Gf, Gc, d, AF, ACC, AV, AC, AX, AY, AZ, B0, Z, Z0, Mu, Sigma,Cor, L, REALIZATION,Am,bm,target, errorfile )
% Adding corrleation matrix

limit_extreme_points = 1000;  
limit_time = 200;

% put all realizations into a column vector
Realization = zeros(sum(L),1);
Abs = Realization;
for k = 1 : K
    Realization(sum(L(1:k-1))+1:sum(L(1:k))) = REALIZATION(k,1:L(k))';
    Abs(sum(L(1:k-1))+1:sum(L(1:k))) = (REALIZATION(k,1:L(k))' - Mu(k)*ones(L(k),1))/sqrt(Sigma(k));
end
Chi = Realization.^2;   % Realization.^2  

% Algorithm BD
stop_criterion = 0;
DELTA = [REALIZATION(:,1)'];  % each row is a (delta^i)', size: E*K; DELTA are the constriants we choose.
Alpha = []; Diff = []; x_list = []; y_list = []; master_time = []; bilinear_time = [];

while stop_criterion == 0 && size(DELTA,1)<=limit_extreme_points
    
    [ obj, x, y, r0, r, vartheta, rc, time ] = Relaxation_Solution_aro (A_xy, b_xy, N, N_pr, M, K, Q, bar_n, Ux, Uy, Gf, Gc, d, AF, ACC, AV, AC, AX, AY, AZ, B0, Z, Z0, Mu, Sigma,Cor, DELTA, Am, bm,limit_time,target);
    Alpha = [Alpha,obj]
    x_list = [x_list, x]; y_list = [y_list, y]; master_time = [master_time,time];
    
    violate = 0; i = 1;
    while i <= length(Am) && violate==0     
        [ obj_bilinear, Delta, time_bi ] = Bilinear_aro( x, y, r, vartheta,rc, N, N_pr, K, bar_n, Q, Gf, Gc, d, AF, ACC, AV, AC, AX, AY, AZ, B0, Z, Z0, L, Realization, Chi, Abs, limit_time, Am(i));
        Delta = Delta'; 
        bilinear_time = [bilinear_time, time_bi];
        if obj_bilinear > r0 + Am(i)*(target- Ux'*x - Uy'*y)-bm(i)*obj + 1e-3
            Diff = [Diff; obj_bilinear - (r0 + Am(i)*(target- Ux'*x - Uy'*y)-bm(i)*obj), obj_bilinear];
            if ismember(Delta,DELTA,'rows')
                RHS = r0 + Am(i)*(target- Ux'*x - Uy'*y)-bm(i)*obj;
                save(errorfile, 'DELTA','Delta','obj_bilinear','RHS');
                i = i+1;
            else 
                DELTA = [DELTA; Delta];
                violate = 1;
                fprintf('Violate: %d\n',i);
            end
        else
            i = i+1;
        end    
    end

%     DELTA;
    
    if violate==0
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


function [ obj, x_sol, y_sol, r0_sol, r_sol, vartheta_sol, rc_sol, relaxation_solution_time ] = Relaxation_Solution_aro (A_xy, b_xy, N, N_pr, M, K, Q, bar_n, Ux, Uy, Gf, Gc, d, AF, ACC, AV, AC, AX, AY, AZ, B0, Z, Z0, Mu, Sigma,Cor, DELTA, Am, bm,limit_time,target )

E = size(DELTA,1);
CHI = DELTA.^2;
ABS = abs( sum((DELTA - ones(E,1)*Mu')./(ones(E,1)*sqrt(Sigma)'),2) )

rome_begin;
h = rome_model('ARO_subproblem');

alpha = newvar(1,1);
x = newvar(M,1,'binary');
y = newvar(N,1,'binary');
r0 = newvar(1,1);    % support
r = newvar(K,1);     % mean
rc = newvar(1,1,'nonneg');    % correlation
vartheta = newvar(K,1,'nonneg');    % second moment

Wf = newvar(N_pr*E,1,'nonneg');     % put all w_f^i into one column vector
Wc = newvar(bar_n*E,1,'nonneg');    % put all w_c^i into one column vector
V = newvar(Q*E,1,'nonneg');         % put all v^i into one column vector
C = newvar(N*E,1,'nonneg');         % put all c^i into one column vector

rome_minimize(alpha);

rome_constraint(alpha >= 0 );     % important
rome_constraint(r0 + Mu'*r + Sigma'*vartheta + Cor*rc <= 0);
for i = 1 : E
    rome_constraint( Am(1)* (Ux'*x + Uy'*y + Gf'*Wf((i-1)*N_pr+1:i*N_pr) + Gc'*Wc((i-1)*bar_n+1:i*bar_n) + d'*V((i-1)*Q+1:i*Q) ) - DELTA(i,:)*r - CHI(i,:)*vartheta - ABS(i,:)*rc <= r0 - bm(1)*alpha + Am(1)*target);
    rome_constraint( Am(2)* (Ux'*x + Uy'*y + Gf'*Wf((i-1)*N_pr+1:i*N_pr) + Gc'*Wc((i-1)*bar_n+1:i*bar_n) + d'*V((i-1)*Q+1:i*Q) ) - DELTA(i,:)*r - CHI(i,:)*vartheta - ABS(i,:)*rc <= r0 - bm(2)*alpha + Am(2)*target);
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
abs_current = Abs(Chosen_realization)

sum(abs_current)

c_2 = second_stage_cost( x, y, Delta, N, N_pr, Q, bar_n, Gf, Gc, d, AF, ACC, AV, AC, AX, AY, AZ, Z0, Z, B0);

obj_bilinear = am*c_2-r'*Delta-vartheta'*chi_current-rc*abs(sum(abs_current));

end



