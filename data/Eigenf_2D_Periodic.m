function [geometry_H, msh_H, space_H, u_0, omega, space_u1, space_u2, Du0_x1, Du0_x2] = ...
              Eigenf_2D_Periodic(problem_data, method_data)
% Extract the fields from the data structures into local variables
data_names = fieldnames (problem_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= problem_data.(data_names{iopt});']);
end
data_names = fieldnames (method_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= method_data.(data_names{iopt});']);
end

% Construct geometry structure
geometry_H  = geo_load (geo_name);

[knots_H, zeta_H]  = kntrefine (geometry_H.nurbs.knots, nsub_H-1, degree_H, regularity_H);
knots_H = kntunclamp(knots_H, degree_H, regularity_H, periodic_directions);
%per costruire spazio per Grad(u0):
[knots_hcurl, degrees_hcurl] = knt_derham(knots_H, degree_H, 'Hcurl'); %ci serve lui per grad_u0

% Construct msh structure
rule_H     = msh_gauss_nodes (nquad_H);
[qn_H, qw_H] = msh_set_quad_nodes (zeta_H, rule_H);
msh_H      = msh_cartesian (zeta_H, qn_H, qw_H, geometry_H);
 
% Construct space structure
space_H  = sp_bspline(knots_H, degree_H, msh_H,[], periodic_directions); % nota argomento in piu' con direzioni periodiche
scalar_spaces = cell(msh_H.ndim,1);
for idim = 1:msh_H.ndim
  scalar_spaces{idim} = sp_bspline(knots_hcurl{idim}, degrees_hcurl{idim}, msh_H,[], periodic_directions); %di nuovo!
end
space_u1  = scalar_spaces{1};                              %   "   S_{p-1, p}     
space_u2  = scalar_spaces{2};                              %   "   S_{p, p-1}     
%D_space_H = sp_vector (scalar_spaces, msh, 'div-preserving'); % space S{p,p}.

% Compute matrices' values
K = op_gradu_gradv_tp (space_H, space_H, msh_H, @(x,y) c_diff(x,y).^2);
M = op_u_v_tp (space_H, space_H, msh_H);

% Ora dobbiamo risolvere il problema agli autovalori:
[V,D] = eigs(K,M,val,'smallestabs');
u_0 = V(:,val);

%ora fissiamo anche la norma.
a = sqrt(1/(2*u_0'*M*u_0));
u_0 = a*u_0; %adesso è una soluzione unica. 

%l'autovalore a cui è associata l'autofunzione è D(val,val)
%ma io voglio restituire omega in realtà!
omega = sqrt(D(val,val));

%proviamo a calcolare anche grad_u_0 
d1 = degree_H(1);  d2 = degree_H(2);
alpha1 = d1/(space_H.knots{1}(d1+2) - space_H.knots{1}(2));% peso der. dir. 1
alpha2 = d2/(space_H.knots{2}(d2+2) - space_H.knots{2}(2));% peso der. dir. 2

mat_u0 = reshape(u_0, space_H.ndof_dir(1),space_H.ndof_dir(2));
dx1_u0 = alpha1*[(mat_u0(2:end,:)- mat_u0(1:end-1,:)); mat_u0(1,:)-mat_u0(end,:)];
dx2_u0 = alpha2*[(mat_u0(:,2:end)- mat_u0(:,1:end-1)) (mat_u0(:,1)-mat_u0(:,end))];

Du0_x1   = reshape(dx1_u0, space_u1.ndof, 1);
Du0_x2   = reshape(dx2_u0, space_u2.ndof, 1);
end


