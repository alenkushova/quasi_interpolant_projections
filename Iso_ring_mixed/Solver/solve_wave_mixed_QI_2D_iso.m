% Description here!
function [geometry, msh, space_p, L, space_v, P] = ...
             solve_wave_mixed_QI_2D_iso (problem_data, method_data)
%% Extract the fields from the data structures into local variables
data_names = fieldnames (problem_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= problem_data.(data_names{iopt});']);
end
data_names = fieldnames (method_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= method_data.(data_names{iopt});']);
end

geometry  = geo_load (geo_name);
[knots, zeta] = kntrefine (geometry.nurbs.knots, nsub-1, degree, regularity);

%apply periodic directions if they exist.
periodic_directions = periodic_sides;
knots = kntunclamp(knots, degree, regularity, periodic_directions); %fondamentale
% Construct msh structure
rule     = msh_gauss_nodes (nquad);
[qn, qw] = msh_set_quad_nodes (zeta, rule);
msh      = msh_cartesian (zeta, qn, qw, geometry);

% Construct SPACE structures as in de Rham complex
[knots_H1, degrees_H1]       = knt_derham(knots, degree, 'H1');
[knots_hdiv, degrees_hdiv]   = knt_derham(knots, degree, 'Hdiv');
[knots_L2, degrees_L2]       = knt_derham(knots, degree, 'L2');

space    = sp_bspline (knots_H1, degrees_H1, msh, 'grad-preserving', periodic_directions); % space S_{p, p}.       
scalar_spaces = cell(msh.ndim,1);
for idim = 1:msh.ndim
  scalar_spaces{idim} = sp_bspline(knots_hdiv{idim}, degrees_hdiv{idim},...
                          msh, 'grad-preserving', periodic_directions); %di nuovo!
end                                      %   "   S_{p-1, p}.  
space_p  = sp_vector (scalar_spaces, msh, 'div-preserving');              %   "    Sp1 x Sp2.
space_v  = sp_bspline (knots_L2, degrees_L2, msh, 'integral-preserving',periodic_directions); %   "   S_{p-1, p-1}.  

%% neeed the space of p1 and p2 over the parametric domain.
geo_para = geo_load('geo_square.txt');
[~, zeta_para] = kntrefine (geo_para.nurbs.knots,...
                                        nsub-1, degree, regularity);
% Construct msh structure
rule_para  = msh_gauss_nodes (nquad);
[qnp, qwp] = msh_set_quad_nodes (zeta_para, rule_para);
msh_para   = msh_cartesian (zeta_para, qnp, qwp, geo_para);%forse solo lei basta.

scalar_spaces = cell(msh.ndim,1);
for idim = 1:msh.ndim
  scalar_spaces{idim} = sp_bspline(knots_hdiv{idim}, degrees_hdiv{idim},...
                          msh_para, 'grad-preserving', periodic_directions); %di nuovo!
end                                      %   "   S_{p-1, p}.  
spazio_p  = sp_vector (scalar_spaces, msh_para, 'div-preserving');              %   "    Sp1 x Sp2.
space_p1 = spazio_p.scalar_spaces{1};
space_p2 = spazio_p.scalar_spaces{2};
%% INTRODUCIAMO LE MATRICI: (here it's ok)
M_DIV = op_u_v_tp (space_p, space_p, msh);% mass matrix in S_{g,g-1} x S_{g-1,g}
M_L2  = op_u_v_tp (space_v, space_v, msh);% mass matrix in S_{g-1,g-1}
A = op_divu_divv_tp (space_p, space_p, msh); % (div(Bi),div(Bj))
B = op_div_v_q_tp (space_p, space_v, msh); % manca coeff c !

%% PROIEZIONE DATI INIZIALI PER P e V (Controlla costruzione msh)
spazio = space_H.constructor(msh);
cambio = op_u_v_tp (spazio, space_v, msh);
rhs_v  = real(1i*omega*cambio*u_0);

D_space_H = sp_vector ({space_u1, space_u2}, msh_H,'curl-preserving'); 
spazio    = D_space_H.constructor(msh);
cambio = op_u_v_tp (spazio, space_p, msh, c_diff);
rhs_p  = cambio*cat(1,Du0_x1,Du0_x2);

l    = M_DIV \ rhs_p;
p    = M_L2  \ rhs_v;

% intialize solutions
L      = [l zeros(space_p.ndof,round(T/k))];
P      = [p zeros(space_v.ndof,round(T/k))];

%% PROIEZIONE FUNZIONI DI BASE... (c*Bi)
% dovremo usare [X2,Y2] e [X1,Y1] risp per la prima e la seconda componente
% delle funzioni che vogliamo proiettare 

N = numel(space.knots{1}(space.degree(1)+1:end-space.degree(1)))-1;
pnt       = linspace (0, 1, 2*N +1);
pnt_m     = linspace (0, 1, 4*N +1);
[X1, Y1] = meshgrid (pnt_m, pnt); % now map this to the ring with F
[X2, Y2] = meshgrid (pnt, pnt_m); % now map this to the ring with F
iso_points1 = msh.map([reshape(X1,numel(X1),1) reshape(Y1,numel(Y1),1)]);
iso_points2 = msh.map([reshape(X2,numel(X2),1) reshape(Y2,numel(Y2),1)]);
FX1 = reshape(iso_points1(1,:,:),size(X1,1),size(X1,2));
FY1 = reshape(iso_points1(2,:,:),size(Y1,1),size(Y1,2));
FX2 = reshape(iso_points2(1,:,:),size(X2,1),size(X2,2));
FY2 = reshape(iso_points2(2,:,:),size(Y2,1),size(Y2,2));
vtk_pts_1 = {pnt_m, pnt};
vtk_pts_2 = {pnt, pnt_m};

Theta = zeros(space_p.ndof,space_p.ndof);
for i = 1 : space_p.ndof
    cfs_1 = (1:space_p1.ndof) == i;
    cfs_2 = (space_p1.ndof +1:space_p.ndof) == i;
    [basis_1, ~] = sp_eval(cfs_1, space_p1, msh_para, vtk_pts_2);
    [basis_2, ~] = sp_eval(cfs_2, space_p2, msh_para, vtk_pts_1);
    b1 = c_diff(FX2,FY2).* basis_1';
    b2 = c_diff(FX1,FY1).* basis_2';
    [theta1, theta2] = Lyche_c_2D_Mixed (b1,b2, space);
    Theta(:,i) = cat(1, theta1, theta2);                
end
filename = ['mixed_Theta_' num2str(nsub(1)) '_constant_iso.mat']
save(filename, 'Theta')

% per le derivate: (controlla BENE!!!!)
d1 = space.degree(1); d2 = space.degree(2); %gradi
step_h1 = knots{1}(d1+2:end-1) - knots{1}(2:end-d1-1);
step_h2 = abs(knots{2}(degree(1)+2)- knots{2}(2));
alpha1 = d1./step_h1;% peso derivate rispetto a direzione 1
alpha2 = d2/step_h2;% peso derivate rispetto a direzione 2
indices_d2 = (1:space_p2.ndof)';
indices_d2 = reshape(indices_d2, space_p2.ndof_dir(1), space_p2.ndof_dir(2));
indices_d2 = [indices_d2 indices_d2(:,1)];

% COEFFICIENTI DI \pi_g(c*p_h)
cph  = Theta* L(:,1);
cph1 = reshape(cph(1:space_p1.ndof), space_p1.ndof_dir(1), space_p1.ndof_dir(2));
cph2 = reshape(cph(space_p1.ndof+1:end), space_p2.ndof_dir(1), space_p2.ndof_dir(2));
cph2 = cph2(indices_d2); 

% derivate, quindi d.o.f. negli spazi S_{p,p-1} e S_{p-1,p} (calcolo div)
DLam1 = alpha1'.*(cph1(2:end,:)-cph1(1:end-1,:));
DLam2 = alpha2*(cph2(:,2:end)-cph2(:,1:end-1)); %direzione periodica
DL1   = reshape(DLam1, space_v.ndof, 1); % reshape to colum vector
DL2   = reshape(DLam2, space_v.ndof, 1); % reshape to colum vector
Div_L = [DL1+DL2, zeros(numel(DL1), round(T/k))];

%IMPOSTIAMO MATRICI DEL SISTEMA DA RISOLVERE 
A = Theta'*A*Theta;
B = k*(B*Theta)';

% assemblaggio porblema più condizioni al bordo.
mat = M_DIV-(k/2)^2*A;
% apply BC to MAT and assemble to solve then ----> (MAT*Lam = rhs);
MAT = M_DIV + (k/2)^2*A; 

for n = 1 : round(T/k)
% risolvo la prima equazione:  
    rhs = mat*L(:,n) - B*P(:,n);% termine noto
    L(:,n+1) = MAT\rhs; % risolviamo il sistema ridotto %NUOVI Lambda.  

    %__ QUI SOTTO SEMBRA OK...   
    % calcola div di questa L 
    
    cph  = Theta* L(:,n+1); %questo full_lam + L(:,n+1) già
    cph1 = reshape(cph(1:space_p1.ndof), space_p1.ndof_dir(1), space_p1.ndof_dir(2));
    cph2 = reshape(cph(space_p1.ndof+1:end), space_p2.ndof_dir(1), space_p2.ndof_dir(2));
    cph2 = cph2(indices_d2); 

    % derivate, quindi d.o.f. negli spazi S_{p,p-1} e S_{p-1,p} (calcolo div)
    DLam1 = alpha1'.*(cph1(2:end,:)-cph1(1:end-1,:));
    DLam2 = alpha2*(cph2(:,2:end)-cph2(:,1:end-1));
    DL1   = reshape(DLam1, space_v.ndof, 1); % reshape to colum vector
    DL2   = reshape(DLam2, space_v.ndof, 1); % reshape to colum vector
        
    %DIVERGENZA
    Div_L(:,n+1) = DL1 + DL2;
    
    % risolvo la seconda equazione:
    P(:,n+1) = P(:,n) + (k/2)*(Div_L(:,n+1) + Div_L(:,n)); %NUOVO Phi      
end

end