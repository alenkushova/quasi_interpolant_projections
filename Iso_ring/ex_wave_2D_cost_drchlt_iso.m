% DESCRIBE PROBLEM HERE:

clear 
close all

%% 0) LOAD OVERKILL SOLUTION:

filename     = '2D_dirichlet_u0_128x128_degree_6x6_eig4_constant_c_iso.mat'; % c=2;

myVar        = {'geometry_H','msh_H','space_H','u_0','omega','space_u1','space_u2','Du0_x1','Du0_x2'};
problem_data = load(filename,myVar{:});
load(filename,myVar{:});

%% 1) PHYSICAL DATA OF THE PROBLEM
% Physical domain
problem_data.geo_name = 'geo_ring.txt';

% Diffusion coefficient:
 problem_data.c_diff         = @(x, y) 2*ones(size(x));
 problem_data.grad_c_diff    = @(x, y) cat (1, ...        % \nabla c
                       reshape (0*x, [1, size(x)]), ...
                       reshape (0*x, [1, size(x)]));     
problem_data.dx1_c = @(x, y) 0*x; % dx1_c  
problem_data.dx2_c = @(x, y) 0*x; % dx2_c  

c  = problem_data.c_diff; dc = problem_data.grad_c_diff;

% Type of boundary conditions for each side of the domain
problem_data.nmnn_sides   = []; % Neumann 
problem_data.drchlt_sides = [1 2 3 4];  % Dirichlet
problem_data.periodic_sides = [];     % Periodic

% Source and boundary terms
problem_data.g = @(x, y) zeros (size (x));% for velocity.
problem_data.h = @(x, y, ind) zeros (size (x)); % for pressures.

%% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data
method_data.degree       = [3 3];  % Degree of the splines
method_data.regularity   = [2 2];  % Regularity of the splines
method_data.nsub         = [32 32];  % Number of subdivisions
method_data.nquad        = [4 4];  % Points for the Gaussian quadrature rule
%method_data.BC           = [2 2];  % Boundary continuity   
method_data.T            = 1;      % Final time (Time step * Numero di elementi in tempo)
method_data.k            = 1/2000; % Time step così varia solo lui.
T = method_data.T; k = method_data.k;


%% 3) CALL TO THE SOLVER
tic
    [geometry, msh, space_p, L, space_v, P] = ...
                      solve_wave_dirichlet_G_2D_iso (problem_data, method_data);
toc


%% 4) POST-PROCESSING
N = 60;
vtk_pts = {linspace(0, 1, N), linspace(0, 1, N)};
[ep_0, FP] = sp_eval (L(:,1), space_p, geometry, vtk_pts); % velocità iniziale
[ev_0, FV] = sp_eval (P(:,1), space_v, geometry, vtk_pts); % pressione iniziale
[ep_T, ~] = sp_eval (L(:,end), space_p, geometry, vtk_pts);% velocità finale
[ev_T, ~] = sp_eval (P(:,end), space_v, geometry, vtk_pts);% pressione finale
[XP, YP]  = deal (squeeze(FP(1,:,:)), squeeze(FP(2,:,:)));  %meshgrid velocità 
[XV, YV]  = deal (squeeze(FV(1,:,:)), squeeze(FV(2,:,:)));  %meshgrid pressione

%numerical solution
figure ('Units', 'pixels', 'Position', [100 200 1000 350]) ;
subplot (1,2,1)
s_iniziale = pcolor(XV,YV,ev_0); hold on;
set(s_iniziale, 'EdgeColor', 'none'); colorbar;
p1 = reshape(ep_0(1,:,:),N,N); p2 = reshape(ep_0(2,:,:),N,N);
quiver(XP,YP,p1,p2), hold off, axis ([0 2 0 2])
title ('Numerical solution at T = 0')

subplot (1,2,2)
s_finale = pcolor(XV,YV,ev_T); hold on;
set(s_finale, 'EdgeColor', 'none'); colorbar;
p1 = reshape(ep_T(1,:,:),N,N); p2 = reshape(ep_T(2,:,:),N,N);
quiver(XP,YP,p1,p2), hold off, axis ([0 2 0 2])
title ('Numerical solution at T = 1')

% exact solution
[G_U0, ~] = sp_eval (u_0, space_H, geometry_H, vtk_pts, 'gradient');
exact_p1  = @(t) c(XP,YP).*squeeze(G_U0(1,:,:))*real(exp(1i*omega*(t)));
exact_p2  = @(t) c(XP,YP).*squeeze(G_U0(2,:,:))*real(exp(1i*omega*(t)));
[U0  , ~] = sp_eval (u_0, space_H, geometry_H, vtk_pts);
exact_v = @(t) U0*real(1i*omega*exp(1i*omega*(t)));

figure ('Units', 'pixels', 'Position', [100 200 1000 350]) ;
subplot (1,2,1)
s_iniziale = pcolor(XV,YV,exact_v(0)); hold on;
set(s_iniziale, 'EdgeColor', 'none'); colorbar;
quiver(XP,YP,exact_p1(0),exact_p2(0)), hold off, axis ([0 2 0 2])
title ('Exact solution at T = 0')

subplot (1,2,2)
s_finale = pcolor(XV,YV,exact_v(T)); hold on;
set(s_finale, 'EdgeColor', 'none'); colorbar;
quiver(XP,YP,exact_p1(T),exact_p2(T)), hold off, axis ([0 2 0 2])
title ('Exact solution at T = 1')

%% error with overkill solution 
err_p = []; err_v = [];
D_space_H    = sp_vector ({space_u1, space_u2}, msh); 
sp_const_p_H = space_p.constructor(msh_H);
sp_const_v_H = space_v.constructor(msh_H);
sp_const_D_H = D_space_H.constructor(msh_H);
%errore in spazio
for n = 1 : round(T/k) : round(T/k)+1
    % errore in norma L2 per soluzione p(x,t) per ogni istante t:
    error_l2_p = my_l2_error (sp_const_p_H, msh_H, L(:,n),...
            [Du0_x1; Du0_x2], sp_const_D_H,@(x,y) c(x, y)*real(exp(1i*omega*(n*k-k))) );
    err_p = [err_p, error_l2_p];

    % errore in norma L2 per soluzione v(x,t) per ogni istante t:
    error_l2_v = my_l2_error (sp_const_v_H, msh_H, P(:,n),...
            u_0, space_H, @(x, y) real(1i*omega*exp(1i*omega*(n*k-k))));
    err_v = [err_v, error_l2_v];
end
%errore in tempo: norma L infinito
err_INF_p = max(abs(err_p)) 
err_INF_v = max(abs(err_v)) 
%errore in tempo: norma L2
err_L2_p  = (k/2*sum((err_p(1:end-1)).^2+(err_p(2:end)).^2))^(1/2)
err_L2_v  = (k/2*sum((err_v(1:end-1)).^2+(err_v(2:end)).^2))^(1/2)

%% save solutions
n = method_data.nsub(1);
filename = ['test_c_cost_' num2str(n) 'x' num2str(n) ...
                'x' num2str(1/k) 'drchlt_iso.mat'];
save(filename)

