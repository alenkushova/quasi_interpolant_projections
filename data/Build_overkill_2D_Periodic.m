clear problem_data  
clear method_data

% 1) PHYSICAL DATA OF THE PROBLEM
% Physical domain
problem_data.geo_name = 'geo_square.txt';

% Type of boundary conditions for each side of the domain
problem_data.nmnn_sides     = [];
problem_data.drchlt_sides   = [];

% diffusion coefficient
problem_data.c_diff         = @(x, y) sin(2*pi*x).*sin(2*pi*y)+2;
problem_data.grad_c_diff    = @(x, y) cat (1, ...        % \nabla c
                      reshape (2*pi*cos (2*pi*x).*sin (2*pi*y), [1, size(x)]), ...
                      reshape (2*pi*sin (2*pi*x).*cos (2*pi*y), [1, size(x)]));

%  problem_data.c_diff         = @(x, y) 2*ones(size(x));
%  problem_data.grad_c_diff    = @(x, y) cat (1, ...        % \nabla c
%                        reshape (0*x, [1, size(x)]), ...
%                        reshape (0*x, [1, size(x)]));                   

c  = problem_data.c_diff;
dc = problem_data.grad_c_diff;

problem_data.val = 4; %scelta dell' i-esimo autovalore più piccolo!

% 2) CHOICE OF THE DISCRETIZATION PARAMETER FOR (type)-HELMHOLTZ EQ.
method_data.degree_H     = [6 6];     % Degree of the splines
method_data.regularity_H = [5 5];     % Regularity of the splines
method_data.nsub_H       = [128 128]; % Number of subdivisions 
            % --> DETERMINA QUANTO ACCURATA è LA SOLUZIONE OVERKILLED.
method_data.nquad_H      = [7 7];     % Points for the Gaussian quadrature rule
method_data.periodic_directions = [1 2];
% dof of overkilled solution:
tic
[geometry_H, msh_H, space_H, u_0, omega, space_u1, space_u2, Du0_x1, Du0_x2]...
                = Eigenf_2D_Periodic (problem_data, method_data);
toc

% save overkill solution
n1 = method_data.nsub_H(1);    n2 = method_data.nsub_H(2);
d1 = method_data.degree_H(1);  d2 = method_data.degree_H(2);
filename = ['2D_periodic_u0_', num2str(n1),'x',num2str(n2),'_degree_',...
    num2str(d1),'x',num2str(d2),'.mat'];
save(filename, 'geometry_H', 'msh_H', 'space_H', 'u_0', 'omega', 'space_u1', 'space_u2', 'Du0_x1', 'Du0_x2');

%questi parametri mi servirà usarli allora nel solver, costituiscono
%la scelta dei dati iniziali e quindi anche della soluzione esatta!
problem_data.geometry_H = geometry_H;
problem_data.msh_H      = msh_H;
problem_data.space_H    = space_H;
problem_data.u_0        = u_0;
problem_data.omega      = omega;

    
