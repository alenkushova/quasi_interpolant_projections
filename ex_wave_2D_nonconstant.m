% DESCRIBE PROBLEM HERE:

clear; close all;


%% 0) LOAD OVERKILL SOLUTION:

filename     = '2D_u0_128x128_degree_6x6_eig4_sine_c.mat'; % sin*sin;

myVar        = {'geometry_H','msh_H','space_H','u_0','omega','space_u1','space_u2','Du0_x1','Du0_x2'};
problem_data = load(filename,myVar{:});
load(filename,myVar{:});

%% 1) PHYSICAL DATA OF THE PROBLEM
% Physical domain
problem_data.geo_name = 'geo_square.txt';

% Diffusion coefficient:
problem_data.c_diff = @(x, y) sin(2*pi*x).*sin(2*pi*y)+2; % c 
problem_data.grad_c = @(x, y) cat (1, ...        % \nabla c
            reshape (2*pi*cos (2*pi*x).*sin (2*pi*y), [1, size(x)]), ...
            reshape (2*pi*sin (2*pi*x).*cos (2*pi*y), [1, size(y)]));
        
problem_data.dx1_c = @(x, y) 2*pi*cos (2*pi*x).*sin (2*pi*y); % dx1_c  
problem_data.dx2_c = @(x, y) 2*pi*sin (2*pi*x).*cos (2*pi*y); % dx2_c  

c  = problem_data.c_diff; dc = problem_data.grad_c;

%% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data
method_data.degree       = [3 3];  % Degree of the splines
method_data.regularity   = [2 2];  % Regularity of the splines
method_data.nsub         = [8 8];  % Number of subdivisions
method_data.nquad        = [4 4];  % Points for the Gaussian quadrature rule
method_data.periodic_dir = [1 2];  % Boundary continuity   
method_data.BC           = [2 2];  % Boundary continuity   
method_data.T            = 60;      % Final time (Time step * Numero di elementi in tempo)
method_data.k            = 1.e-3;   % Time step così varia solo lui.
T = method_data.T; k = method_data.k;

%% 3) CALL TO THE SOLVER: G = Galerkin, QI = quasi-interpolant
tic
    [geometry, msh, space_p, L, space_v, P] = ...
                      solve_wave_periodic_QI_2D (problem_data, method_data);
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

figure ('Units', 'pixels', 'Position', [100 200 1000 350]) ;
subplot (1,2,1)
s_iniziale = pcolor(XV,YV,ev_0); hold on;
set(s_iniziale, 'EdgeColor', 'none'); colorbar;
p1 = reshape(ep_0(1,:,:),N,N); p2 = reshape(ep_0(2,:,:),N,N);
quiver(XP,YP,p1,p2), hold off, axis ([0 1 0 1])
title ('Numerical solution at T = 0')

subplot (1,2,2)
s_finale = pcolor(XV,YV,ev_T); hold on;
set(s_finale, 'EdgeColor', 'none'); colorbar;
p1 = reshape(ep_T(1,:,:),N,N); p2 = reshape(ep_T(2,:,:),N,N);
quiver(XP,YP,p1,p2), hold off, axis ([0 1 0 1])
title ('Numerical solution at T = 1')

%% error with overkill solution???
err_p = []; err_v = [];
% constructor mesh 
D_space_H    = sp_vector ({space_u1, space_u2}, msh); 
sp_const_p_H = space_p.constructor(msh_H);
sp_const_v_H = space_v.constructor(msh_H);
sp_const_D_H = D_space_H.constructor(msh_H);

%errore in spazio
for n = 1 : round(T/k)+1 %step "1:200:round(T/k)+1" per h conv.
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
err_INF_p = max(abs(err_p));  
err_INF_v = max(abs(err_v)); 
%errore in tempo: norma L2
err_L2_p  = (k/2*sum((err_p(1:end-1)).^2+(err_p(2:end)).^2))^(1/2);
err_L2_v  = (k/2*sum((err_v(1:end-1)).^2+(err_v(2:end)).^2))^(1/2);

%% energy check
counter = 1;
for j = [1, round(T/k)+1]
[nrg_p] = sp_l2_error (space_p, msh, L(:,j), @(x,y) [0*x, 0*y]);
[nrg_v] = sp_l2_error (space_v, msh, P(:,j), @(x,y) 0*x);
nrg(counter) = (nrg_p^2+nrg_v^2)/2;
counter = counter + 1;
end
energy = abs(nrg(1)-nrg(end))/abs(nrg(1))
%% save solutions
n = method_data.nsub(1);
filename = ['test_c_sine_energy_' num2str(n) 'x' num2str(n) 'x' num2str(round(T/k)) '.mat'];
save(filename)