function [geometry, msh, sp, L, space_v, P] = ... 
             solve_wave_periodic_QI_2D_opt(problem_data, method_data)

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
knots = kntunclamp(knots, degree, regularity, periodic_dir); % fondamentale chiamare quest fun. ogni volta che hai periodicitÃ !

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

% INTRODUCIAMO LE MATRICI:
M = op_u_v_tp (sp, sp, msh);% mass matrix
A = op_divu_divv_tp (sp, sp, msh); % (div(Bi),div(Bj))
B = op_div_v_q_tp (sp, space_v, msh); % manca coeff c !

% PROIEZIONE DATI INIZIALI PER P    (Controlla costruzione msh)
N = numel(space.knots{1}(space.degree(1)+1:end-space.degree(1)))-1;
pnt       = linspace (0, 1, 2*N +1);
pnt_m     = linspace (0, 1, 4*N +1);
[X1, Y1] = meshgrid (pnt_m, pnt);
[X2, Y2] = meshgrid (pnt, pnt_m);
vtk_pts_1 = {pnt_m, pnt};
vtk_pts_2 = {pnt, pnt_m};
% VALUTIAMO IL GRADIENTE DI U_0 PER INT RISP DIR 1
[grad_u0_1, ~]   = sp_eval (u_0, space_H, geometry_H, vtk_pts_1, 'gradient'); % valutato su punti_m x punti
% [X1, Y1]  = deal (squeeze(f1(1,:,:)), squeeze(f1(2,:,:)));

% VALUTIAMO IL GRADIENTE DI U_0 PER INT RISP DIR 2
[grad_u0_2, ~]   = sp_eval (u_0, space_H, geometry_H, vtk_pts_2, 'gradient'); % valutato su punti x punti_m
% [X2, Y2]  = deal (squeeze(f2(1,:,:)), squeeze(f2(2,:,:)));
p1 = c_diff(X2,Y2).*(squeeze(grad_u0_2(1,:,:))'); % seconda comp di p risp x1.
p2 = c_diff(X1,Y1).*(squeeze(grad_u0_1(2,:,:))'); % prima comp di p la integro risp x2.

% PROIEZIONE DATI INIZIALI PER V
Z_dir_1   = space_v.knots{1}(space_v.degree(1)+1:end-space_v.degree(1)); % nodi interni in dir 1 
Zm_dir_1  = (Z_dir_1(1:end-1)+Z_dir_1(2:end))/2; % più punti medi 
u1        = sort([Z_dir_1 Zm_dir_1]);% ora sono vertici e punti medi ordinati 
Z_dir_2   = space_v.knots{2}(space_v.degree(2)+1:end-space_v.degree(2)); % nodi interni in dir 2 
Zm_dir_2  = (Z_dir_2(1:end-1)+Z_dir_2(2:end))/2; % più punti medi 
u2        = sort([Z_dir_2 Zm_dir_2]);% ora sono vertici e punti medi ordinati 
vtk_pts_3 = {u1, u2};
[v, ~]    = sp_eval (u_0, space_H, geometry_H, vtk_pts_3);
%[X3, Y3]  = deal (squeeze(f3(1,:,:)), squeeze(f3(2,:,:))); %forse non servono.

v_2proj = [v(:,end-2*space_v.degree(1):end) v(:,2:end-1) v(:,1:2*space_v.degree(1))];                  % in dir 1 
v_2proj = [v_2proj(end-2*space_v.degree(2):end,:); v_2proj(2:end-1,:); v_2proj(1:2*space_v.degree(2),:)]; % in dir 2


% Proiezione dati iniziali: (CONTROLLA!!) sembra buona però!
[Lam1, Lam2] = Lyche_c_2D_Periodic (p1, p2, space);
Phi          = Lyche_2D_Periodic (v_2proj*real(1i*omega*exp(1i*omega*0)), space_v);

L  = [cat(1,Lam1,Lam2) zeros(sp.ndof,round(T/k))];
P  = [Phi zeros(space_v.ndof,round(T/k))];

% PROIEZIONE FUNZIONI DI BASE... (c*Bi)
% dovremo usare [X2,Y2] e [X1,Y1] risp per la prima e la seconda componente
% delle funzioni che vogliamo proiettare 

 load('Theta_64_sine.mat');
%load('Theta_64_constant.mat');

% Theta = zeros(sp.ndof,sp.ndof);
% tic
% msh_aux2 = msh_auxiliary (geometry, msh, vtk_pts_2);
% msh_aux1 = msh_auxiliary (geometry, msh, vtk_pts_1);
% sp_aux_p1 = space_p1.constructor (msh_aux2);
% sp_aux_p2 = space_p2.constructor (msh_aux1);
% msh_aux2_struct = msh_precompute (msh_aux2);
% msh_aux1_struct = msh_precompute (msh_aux1);
% sp_aux_p1_struct = sp_precompute (sp_aux_p1, msh_aux2_struct);
% sp_aux_p2_struct = sp_precompute (sp_aux_p2, msh_aux1_struct);
% for i = 1 : sp.ndof
%     cfs_1 = (1:space_p1.ndof) == i;
%     cfs_2 = (space_p1.ndof +1:sp.ndof) == i;
%     [basis_1, ~] = sp_eval_msh(cfs_1, sp_aux_p1_struct, msh_aux2_struct);
%     [basis_2, ~] = sp_eval_msh(cfs_2, sp_aux_p2_struct, msh_aux1_struct);
%     basis_1 = reshape (basis_1, cellfun(@numel,vtk_pts_2));
%     basis_2 = reshape (basis_2, cellfun(@numel,vtk_pts_1));
%     b1 = c_diff(X2,Y2).* basis_1';
%     b2 = c_diff(X1,Y1).* basis_2';
%     [theta1, theta2] = Lyche_c_2D_Periodic(b1,b2, space);
%     Theta(:,i) = cat(1, theta1, theta2);                
% end
%toc
% per le derivate: (controlla BENE!!!!)
d1 = space.degree(1); d2 = space.degree(2); %gradi
step_h = abs(knots{1}(degree(1)+2)- knots{1}(2));
alpha1 = d1/step_h;% peso der. dir. 1
alpha2 = d2/step_h;% peso der. dir. 2
indices_d1 = (1:space_p1.ndof)';
indices_d2 = (1:space_p2.ndof)';
indices_d1 = reshape(indices_d1, space_p1.ndof_dir(1), space_p1.ndof_dir(2));
indices_d2 = reshape(indices_d2, space_p2.ndof_dir(1), space_p2.ndof_dir(2));
indices_d1 = [indices_d1;indices_d1(1,:)];
indices_d2 = [indices_d2 indices_d2(:,1)];

% COEFFICIENTI DI \pi_g(c*p_h)
cph  = Theta* L(:,1);
cph1 = reshape(cph(1:space_p1.ndof), space_p1.ndof_dir(1), space_p1.ndof_dir(2));
cph2 = reshape(cph(space_p1.ndof+1:end), space_p2.ndof_dir(1), space_p2.ndof_dir(2));
cph1 = cph1(indices_d1);
cph2 = cph2(indices_d2); 

% derivate, quindi d.o.f. negli spazi S_{p,p-1} e S_{p-1,p} (calcolo div)
DLam1 = alpha1*(cph1(2:end,:)-cph1(1:end-1,:));
DLam2 = alpha2*(cph2(:,2:end)-cph2(:,1:end-1));
DL1   = reshape(DLam1, space_p1.ndof, 1); % reshape to colum vector
DL2   = reshape(DLam2, space_p2.ndof, 1); % reshape to colum vector
Div_L = [DL1+DL2, zeros(numel(DL1), round(T/k))];

%IMPOSTIAMO MATRICI DEL SISTEMA DA RISOLVERE 
A = Theta'*A*Theta;
B = k*(B*Theta)';

% assemblaggio porblema più condizioni al bordo.
mat = M-(k/2)^2*A;
% apply BC to MAT and assemble to solve then ----> (MAT*Lam = rhs);
MAT = M + (k/2)^2*A; 

%RISOLUZIONE
for n = 1 : round(T/k)
% risolvo la prima equazione:   %CE UN ERRORE!!!! QUI!!!
    rhs = mat*L(:,n) - B*P(:,n);% termine noto
    L(:,n+1) = MAT\rhs; % risolviamo il sistema ridotto %NUOVI Lambda.  

    %__ QUI SOTTO SEMBRA OK...   
    % calcola div di questa L 
    
    cph  = Theta* L(:,n+1); %questo full_lam + L(:,n+1) già
    cph1 = reshape(cph(1:space_p1.ndof), space_p1.ndof_dir(1), space_p1.ndof_dir(2));
    cph2 = reshape(cph(space_p1.ndof+1:end), space_p2.ndof_dir(1), space_p2.ndof_dir(2));
    cph1 = cph1(indices_d1);
    cph2 = cph2(indices_d2); 

    % derivate, quindi d.o.f. negli spazi S_{p,p-1} e S_{p-1,p} (calcolo div)
    DLam1 = alpha1*(cph1(2:end,:)-cph1(1:end-1,:));
    DLam2 = alpha2*(cph2(:,2:end)-cph2(:,1:end-1));
    DL1   = reshape(DLam1, space_p1.ndof, 1); % reshape to colum vector
    DL2   = reshape(DLam2, space_p2.ndof, 1); % reshape to colum vector
        
    %DIVERGENZA
    Div_L(:,n+1) = DL1 + DL2;
    
    % risolvo la seconda equazione:
    P(:,n+1) = P(:,n) + (k/2)*(Div_L(:,n+1) + Div_L(:,n)); %NUOVO Phi      
end

end

function msh_out = msh_auxiliary (geometry, msh, pts)

ndim = numel (msh.breaks);
endpoints = zeros (2, ndim);
if (isfield (geometry, 'nurbs'))
    nurbs = geometry.nurbs;
    for idim=1:ndim
      endpoints(:,idim) = nurbs.knots{idim}([nurbs.order(idim), end-nurbs.order(idim)+1]);
    end
    clear nurbs
else
    for idim=1:ndim
      endpoints(:,idim) = [msh.breaks{idim}(1) msh.breaks{idim}(end)];
    end
end
for jj = 1:ndim
    pts{jj} = pts{jj}(:)';
    if (numel (pts{jj}) > 1)
      brk{jj} = [endpoints(1,jj), pts{jj}(1:end-1) + diff(pts{jj})/2, endpoints(2,jj)];
    else
      brk{jj} = endpoints(:,jj).';
    end
end
msh_out = msh_cartesian (brk, pts, [], geometry, 'boundary', false);
end