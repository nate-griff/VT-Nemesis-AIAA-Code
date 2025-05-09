CM = 2.25;
rho_e = 2.69826e-4;
rho_e = rho_e * 1000;
c_star = 4285.48;

c = 1.25 * 0.0254;
b = 2 * 0.0254;
S = c * b;

q_inf = .5 * rho_e * c_star^2;

M = CM * S * q_inf * c
M_ftlb = M * 0.737562149277