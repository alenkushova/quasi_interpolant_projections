% Description here to be changed 

function [geometry, msh, space_p, L, space_v, P] = ...
             solve_wave_dirichlet_G_2D_iso (problem_data, method_data)
% Extract the fields from the data structures into local variables
data_names = fieldnames (problem_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= problem_data.(data_names{iopt});']);
end
data_names = fieldnames (method_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= method_data.(data_names{iopt});']);
end

geometry  = geo_load (geo_name);
degelev  = max (degree - (geometry.nurbs.order-1), 0);
nurbs    = nrbdegelev (geometry.nurbs, degelev);
[rknots, zeta, nknots] = kntrefine (nurbs.knots, nsub-1, nurbs.order-1, regularity);

nurbs = nrbkntins (nurbs, nknots);
geometry = geo_load (nurbs); % (check continuity ask again rafa)


% check for periodic directions
if (exist('periodic_sides', 'var'))
  periodic_directions = unique(ceil(periodic_sides./2), 'legacy');
  rknots = kntunclamp(rknots, degree, regularity, periodic_directions); %fondamentale
else
  periodic_directions = [];
end

% Construct msh structure
rule     = msh_gauss_nodes (nquad);
[qn, qw] = msh_set_quad_nodes (zeta, rule);
msh      = msh_cartesian (zeta, qn, qw, geometry);

% Construct SPACE structures as in de Rham complex
[knots_H1, degrees_H1]       = knt_derham(rknots, degree, 'H1');
[knots_hdiv, degrees_hdiv]   = knt_derham(rknots, degree, 'Hdiv');
[knots_L2, degrees_L2]       = knt_derham(rknots, degree, 'L2');

space    = sp_bspline (knots_H1, degrees_H1, msh); % space S_{p, p}.       
scalar_spaces = cell(msh.ndim,1);
for idim = 1:msh.ndim
  scalar_spaces{idim} = sp_bspline(knots_hdiv{idim}, degrees_hdiv{idim},...
                          msh, 'grad-preserving', periodic_directions); %di nuovo!
end
space_p1 = scalar_spaces{1};                                              % space S_{p, p-1}.     
space_p2 = scalar_spaces{2};                                              %   "   S_{p-1, p}.  
space_p  = sp_vector (scalar_spaces, msh, 'div-preserving');              %   "    Sp1 x Sp2.
space_v  = sp_bspline (knots_L2, degrees_L2, msh, 'integral-preserving',periodic_directions); %   "   S_{p-1, p-1}.

%% INTRODUCIAMO LE MATRICI: (here it's ok)
M_DIV = op_u_v_tp (space_p, space_p, msh);% mass matrix in S_{g,g-1} x S_{g-1,g}
M_L2  = op_u_v_tp (space_v, space_v, msh);% mass matrix in S_{g-1,g-1}
A_2   = my_op_div_v_q_tp (space_p, space_v, msh, c_diff);
A_p1  = op_u_v_tp (space_p1, space_v, msh, dx1_c);
A_p2  = op_u_v_tp (space_p2, space_v, msh, dx2_c);
A_1   = [A_p1 A_p2];
A     = k/2*(A_1 + A_2);
MAT = [M_DIV, A.'; - A, M_L2];

%% PROIEZIONE DATI INIZIALI PER P e V (Controlla costruzione msh)
spazio = space_H.constructor(msh);
cambio = op_u_v_tp (spazio, space_v, msh);
rhs_v  = real(1i*omega*cambio*u_0);
spazio = space_u1.constructor(msh);
cambio = op_u_v_tp (spazio, space_p1, msh, c_diff);
rhs_p1 = cambio*Du0_x1;
spazio = space_u2.constructor(msh);
cambio = op_u_v_tp (spazio, space_p2, msh, c_diff);
rhs_p2 = cambio*Du0_x2;
rhs_p  = cat(1,rhs_p1,rhs_p2);
l    = M_DIV \ rhs_p;
p    = M_L2  \ rhs_v;
u = [l;p];

% intialize solutions
L      = [l zeros(space_p.ndof,round(T/k))];
P      = [p zeros(space_v.ndof,round(T/k))];

%% Solve the problem:
% NB: homogeneous Dirichlet boundary conditions on u becomes now some "do
% nothing" boundary conditions on v (v = u_t = 0 on boundary of Dirichlet)
% so at this point we just need to solve the problem.

for n = 1 : round(T/k)
    rhs = MAT'*u;% termine noto
    u   = MAT\rhs; % risolviamo il sistema  
    L(:,n+1) = u(1:space_p.ndof);
    P(:,n+1) = u(space_p.ndof+1:end);
end

end