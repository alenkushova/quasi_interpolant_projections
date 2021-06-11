function [geometry_H, msh_H, space_H, u_0,...
    omega, space_u1, space_u2, Du0_x1, Du0_x2] =...
                        Eigenf_2D_Dirichlet_iso (problem_data, method_data)
                    
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
degelev  = max (degree_H - (geometry_H.nurbs.order-1), 0);
nurbs    = nrbdegelev (geometry_H.nurbs, degelev);
[rknots_H, zeta_H, nknots_H] = kntrefine (nurbs.knots, nsub_H-1,...
                                            nurbs.order-1, regularity_H);
nurbs = nrbkntins (nurbs, nknots_H);
geometry_H = geo_load (nurbs); % (check continuity ask again rafa)

% check for periodic directions
if (exist('periodic_sides', 'var'))
  periodic_directions = unique(ceil(periodic_sides./2), 'legacy');
  rknots_H = kntunclamp(rknots_H, degree_H, regularity_H,...
                            periodic_directions); %fondamentale
else
  periodic_directions = [];
end

[knots_hcurl, degrees_hcurl] = knt_derham(rknots_H, degree_H, 'Hcurl');

% Construct msh structure
rule_H     = msh_gauss_nodes (nquad_H);
[qn_H, qw_H] = msh_set_quad_nodes (zeta_H, rule_H);
msh_H      = msh_cartesian (zeta_H, qn_H, qw_H, geometry_H);
 
% Construct space structure
space_H  = sp_bspline(rknots_H, degree_H, msh_H,'grad-preserving',...
    periodic_directions); % nota argomento in piu' con direzioni periodiche
scalar_spaces = cell(msh_H.ndim,1);
for idim = 1:msh_H.ndim
  scalar_spaces{idim} = sp_bspline(knots_hcurl{idim},... %di nuovo!
      degrees_hcurl{idim}, msh_H, 'grad-preserving', periodic_directions);
end
space_u1  = scalar_spaces{1};  %   "   S_{p-1, p}     
space_u2  = scalar_spaces{2};  %   "   S_{p, p-1}     
%D_space_H = sp_vector (scalar_spaces, msh, 'div-preserving'); % space S{p,p}.

% Compute matrices' values
K = op_gradu_gradv_tp (space_H, space_H, msh_H, @(x,y) c_diff(x,y).^2);
M = op_u_v_tp (space_H, space_H, msh_H);

u_0 = zeros (space_H.ndof, 1);
drchlt_dofs = [];
for iside = 1:numel (drchlt_sides)
  drchlt_dofs = union (drchlt_dofs, space_H.boundary(drchlt_sides(iside)).dofs);
end
int_dofs = setdiff (1:space_H.ndof, drchlt_dofs);
% Ora dobbiamo risolvere il problema agli autovalori:
[V,D] = eigs( (K(int_dofs,int_dofs)), (M(int_dofs,int_dofs)),val,'smallestabs');
u_0(int_dofs) = V(:,val);

%ora fissiamo anche la norma.
a = sqrt(1/(2*u_0'*M*u_0));
u_0 = a*u_0; %adesso è una soluzione unica. 

%l'autovalore a cui è associata l'autofunzione è D(val,val)
%ma io voglio restituire omega in realtà!
omega = sqrt(D(val,val));

%proviamo a calcolare anche grad_u_0 
d1 = method_data.degree_H(1);  d2 = method_data.degree_H(2);
%g/(xi_{i+g+1}- \xi{i+1})
step_h1 = space_H.knots{1}(1+d1+1:end-1) - space_H.knots{1}(2:end-d1-1); 
%g/(xi_{i+g+1}- \xi{i+1})
step_h2 = space_H.knots{2}(1+d2+1:end-1) - space_H.knots{2}(2:end-d2-1); 
alpha1 = d1./step_h1;
alpha2 = d2./step_h2;
mat_u0 = reshape(u_0, space_H.ndof_dir(1),space_H.ndof_dir(2));
dx1_u0 = alpha1'.*(mat_u0(2:end,:)- mat_u0(1:end-1,:));
dx2_u0 = alpha2.*(mat_u0(:,2:end)- mat_u0(:,1:end-1));
Du0_x1   = reshape(dx1_u0, space_u1.ndof, 1);
Du0_x2   = reshape(dx2_u0, space_u2.ndof, 1);
end