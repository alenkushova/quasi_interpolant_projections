%% Dirichelet approximation

u   = @(x,y) sin(2*pi*x).*y.^2;
ux  = @(x,y) 2*pi*cos(2*pi*x).*y.^2;
uy  = @(x,y) 2*sin(2*pi*x).*y;
Grad_u = @(x, y) cat (1, ...
           reshape (ux(x,y), [1, size(x)]), ...
           reshape (uy(x,y), [1, size(x)])); %for sp_vec evaluation.
       
% Physical domain
geo_name = 'geo_square.txt';
% Method data
degree       = [3 3];   % Degree of the splines
regularity   = [2 2];   % Regularity of the splines
nsub         = [32 32]; % Number of subdivisions
nquad        = [4 4];   % Points for the Gaussian quadrature rule
drch_dir     = [1 2];   % Dirchlet boundaries
prdc_dir     = [];   % Periodic boundaries

geometry  = geo_load (geo_name);
[knots, zeta]  = kntrefine (geometry.nurbs.knots, nsub-1, degree, regularity);
% fondamentale chiamare quest fun. ogni volta che hai periodicità!
knots = kntunclamp(knots, degree, regularity, prdc_dir); 

% Construct msh structure
rule     = msh_gauss_nodes (nquad);
[qn, qw] = msh_set_quad_nodes (zeta, rule);
msh      = msh_cartesian (zeta, qn, qw, geometry);

% Construct SPACE structures as in de Rham complex
[knots_H1, degrees_H1]       = knt_derham(knots, degree, 'H1');
[knots_hdiv, degrees_hdiv]   = knt_derham(knots, degree, 'Hdiv');
[knots_L2, degrees_L2]       = knt_derham(knots, degree, 'L2');

space    = sp_bspline (knots_H1, degrees_H1, msh,...
            'grad-preserving', prdc_dir); % space S_{p, p}.       
scalar_spaces = cell(msh.ndim,1);
for idim = 1:msh.ndim
  scalar_spaces{idim} = sp_bspline(knots_hdiv{idim}, degrees_hdiv{idim},...
                          msh, 'grad-preserving', prdc_dir); %di nuovo!
end
space_p1 = scalar_spaces{1};                                % space S_{p, p-1}=Sp1.     
space_p2 = scalar_spaces{2};                                %   "   S_{p-1, p}=Sp2.  
space_p  = sp_vector (scalar_spaces, msh, 'div-preserving');%   "    Sp1 x Sp2.
space_v  = sp_bspline (knots_L2, degrees_L2, msh,...
            'integral-preserving',prdc_dir);                %   "   S_{p-1, p-1}.

%% Proiezione con dati al bordo di dirichlet, vediamo se funziona:
[Lam1, Lam2] = Lyche_c_2D_Drch (ux, uy, space); 
Lam  = cat(1,Lam1,Lam2);
Phi = Lyche_2D_Drch (u, space_v);
%% Plot
N = 60;
vtk_pts   = {linspace(0, 1, N), linspace(0, 1, N)};
[ep, FP]  = sp_eval (Lam(:,1), space_p, geometry, vtk_pts); % velocità 
[ev, FV]  = sp_eval (Phi(:,1), space_v, geometry, vtk_pts); % pressione 
[XP, YP]  = deal (squeeze(FP(1,:,:)), squeeze(FP(2,:,:)));  % meshgrid velocità 
[XV, YV]  = deal (squeeze(FV(1,:,:)), squeeze(FV(2,:,:)));  % meshgrid pressione

figure(1)
s = pcolor(XV,YV,ev); hold on;
set(s, 'EdgeColor', 'none'); colorbar;
p1 = reshape(ep(1,:,:),N,N); p2 = reshape(ep(2,:,:),N,N);
quiver(XP,YP,p1,p2), hold off, axis ([0 1 0 1])
title ('Quasi interpolant approximation with Drch. = 0')

%% Interpolation error:
errl2_Gradu = sp_l2_error (space_p, msh, Lam, Grad_u)
errl2_u   = sp_l2_error (space_v, msh, Phi, u)
