function [geometry, msh_H, space_H, u_0, omega] = ...
              Eigenf_2D(problem_data, method_data)
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
geometry  = geo_load (geo_name);

domain_elev_H      = nrbdegelev(geometry.nurbs, degree_H - 1); % how many times we need to increase the order (assuming we start from linear)
[knots_H, zeta_H]  = kntrefine (domain_elev_H.knots, nsub_H-1, degree_H, regularity_H);
domain_clamped_H   = nrbkntins(domain_elev_H, {knots_H{1}(degree_H(1)+2 : end-degree_H(1)-1), knots_H{2}(degree_H(2)+2 : end-degree_H(2)-1)});
domain_unclamped_H = nrbunclamp(domain_clamped_H, BC_H);

% Construct msh structure
rule_H     = msh_gauss_nodes (nquad_H);
[qn_H, qw_H] = msh_set_quad_nodes (zeta_H, rule_H);
msh_H      = msh_cartesian (zeta_H, qn_H, qw_H, geometry);
 
% Construct space structure
knt   = domain_unclamped_H.knots;
knt_u1 = knt; knt_u2 = knt;
knt_u1(1) = num2cell (knt{1}(2:end-1),2); 
knt_u2(2) = num2cell (knt{2}(2:end-1),2);

space_H   = sp_bspline (knt, degree_H, msh_H);
space_u1  = sp_bspline (knt_u1, degree_H-[1,0], msh_H); %   "   S_{p-1, p}     
space_u2  = sp_bspline (knt_u2, degree_H-[0,1], msh_H); %   "   S_{p, p-1}     
D_space_H = sp_vector ([{space_u1} {space_u2}], msh_H);%   "   s2 \times s3.

% Compute matrices' values
[rows_K, cols_K, values_K] = op_gradu_gradv_tp (space_H, space_H, msh_H,@(x,y) c_diff(x,y).^2);
[rows_M, cols_M, values_M] = op_u_v_tp (space_H, space_H, msh_H);
Full_M = sparse(rows_M, cols_M, values_M);

% INDICIZZIAMO PER INCOLLAMENTO: 
indices = (1:space_H.ndof)';
indices = reshape(indices, space_H.ndof_dir(1), space_H.ndof_dir(2));

% incolliamo indici
indices(:,end-(BC_H):end) = indices(:,1:BC_H +1);
indices(end-(BC_H):end,:) = indices(1:BC_H +1,:);

% ripristiniamo dimensione
indices = reshape(indices, space_H.ndof, 1);

% APPLICHIAMO CONDIZIONI AL BORDO PERIODICHE (con reg. massima):
rows_M = indices(rows_M);
cols_M = indices(cols_M);
rows_K = indices(rows_K);
cols_K = indices(cols_K);

% number of reduced d.o.f. (Posso usare BC al posto di 3 e 2)
R_dof = space_H.ndof - (BC_H(1)+1)*space_H.ndof_dir(2) - (BC_H(2)+1)*(space_H.ndof_dir(1)-(BC_H(1)+1));

% assemblaggio:
M = sparse(rows_M, cols_M, values_M);
K = sparse(rows_K, cols_K, values_K);

% delete zero rows and colums from sparse matrices M, K.
M( ~any(M,2), : ) = [];  %rows of M
M( :, ~any(M,1) ) = [];  %columns of M
K( ~any(K,2), : ) = [];  %rows of K
K( :, ~any(K,1) ) = [];  %columns of K

% Ora dobbiamo risolvere il problema agli autovalori:
[V,D] = eigs(K,M,val,'smallestabs');
u_0 = V(:,val);

%:
ind_reset  = (1:R_dof);
ind_reset  = reshape (ind_reset, space_H.ndof_dir(1) - (BC_H(1)+1), ...
                            space_H.ndof_dir(2) - (BC_H(2)+1));
ind_reset  = [ind_reset  ind_reset(:,1:BC_H+1)];
ind_reset  = [ind_reset; ind_reset(1:BC_H+1,:)];
ind_reset_vec  = reshape(ind_reset, space_H.ndof, 1);

%ripristiniamo la dimensione di u_0:
u_0 = u_0(ind_reset_vec);

%ora fissiamo anche la norma.
a = sqrt(1/(2*u_0'*Full_M*u_0));
u_0 = a*u_0; %adesso è una soluzione unica. 

%l'autovalore a cui è associata l'autofunzione è D(val,val)
%ma io voglio restituire omega in realtà!
omega = sqrt(D(val,val));
end
