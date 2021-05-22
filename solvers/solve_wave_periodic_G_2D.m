function [geometry, msh, sp, L, space_v, P] = ... 
             solve_wave_periodic_G_2D(problem_data, method_data)
                      %SEMBRA BUONO!!!
                      
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

[knots, zeta]  = kntrefine (geometry.nurbs.knots, nsub-1, degree, regularity);
knots = kntunclamp(knots, degree, regularity, periodic_dir); % fondamentale chiamare quest fun. ogni volta che hai periodicità!

% Construct msh structure
rule     = msh_gauss_nodes (nquad);
[qn, qw] = msh_set_quad_nodes (zeta, rule);
msh      = msh_cartesian (zeta, qn, qw, geometry);

% Construct SPACE structures as in de Rham complex
[knots_H1, degrees_H1]       = knt_derham(knots, degree, 'H1');
[knots_hdiv, degrees_hdiv]   = knt_derham(knots, degree, 'Hdiv');
[knots_L2, degrees_L2]       = knt_derham(knots, degree, 'L2');

space    = sp_bspline (knots_H1, degrees_H1, msh, 'grad-preserving', periodic_dir); % space S_{p, p}.       
scalar_spaces = cell(msh.ndim,1);
for idim = 1:msh.ndim
  scalar_spaces{idim} = sp_bspline(knots_hdiv{idim}, degrees_hdiv{idim},...
                          msh, 'grad-preserving', periodic_dir); %di nuovo!
end
space_p1 = scalar_spaces{1};                                              % space S_{p, p-1}.     
space_p2 = scalar_spaces{2};                                              %   "   S_{p-1, p}.  
sp       = sp_vector (scalar_spaces, msh, 'div-preserving');              %   "    Sp1 x Sp2.
space_v  = sp_bspline (knots_L2, degrees_L2, msh, 'integral-preserving',periodic_dir); %   "   S_{p-1, p-1}.

%% INTRODUCIAMO LE MATRICI:
M_DIV = op_u_v_tp (sp, sp, msh);% mass matrix in S_{g,g-1} x S_{g-1,g}
M_L2  = op_u_v_tp (space_v, space_v, msh);% mass matrix in S_{g-1,g-1}
A_2   = my_op_div_v_q_tp (sp, space_v, msh, c_diff);
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
Lam    = M_DIV \ rhs_p;
Phi    = M_L2  \ rhs_v;
L      = [Lam zeros(sp.ndof,round(T/k))];
P      = [Phi zeros(space_v.ndof,round(T/k))];

dof = cat(1,Lam,Phi);

for n = 1 : round(T/k)
% risolvo la prima equazione:   %CE UN ERRORE!!!! QUI!!!
    rhs      = MAT'*dof;% termine noto
    dof      = MAT\rhs; % risolviamo il sistema  
    L(:,n+1) = dof(1:sp.ndof);
    P(:,n+1) = dof(sp.ndof+1:end);
end
end