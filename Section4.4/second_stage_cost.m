function [ second_stage_cost ] = second_stage_cost( x, y, Delta, N, N_pr, Q, bar_n, Gf, Gc, d, AF, ACC, AV, AC, AX, AY, AZ, Z0, Z, B0)

rome_begin;
h = rome_model('minimal_cost');

Wf = newvar(N_pr,1,'nonneg');
Wc = newvar(bar_n,1,'nonneg');
V = newvar(Q,1,'nonneg');
C = newvar(N,1,'nonneg');

rome_minimize( Gf'*Wf + Gc'*Wc + d'*V );
rome_constraint(AF*Wf + ACC*Wc + AV*V + AC*C >= AX*x + AY*y + AZ*(Z0+Z*Delta) + B0);

h.solve('CPLEX');
second_stage_cost = h.objective;

end

