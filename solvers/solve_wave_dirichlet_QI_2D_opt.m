function [geometry, msh, sp, L, space_v, P] = ... 
             solve_wave_dirichlet_QI_2D_opt(problem_data, method_data)

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

% check for periodic directions
if (exist('periodic_sides', 'var'))
  periodic_directions = unique(ceil(periodic_sides./2), 'legacy');
  knots = kntunclamp(knots, degree, regularity, periodic_directions); %fondamentale
else
  periodic_directions = [];
end

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
end
space_p1 = scalar_spaces{1};                                              % space S_{p, p-1}.     
space_p2 = scalar_spaces{2};                                              %   "   S_{p-1, p}.  
sp       = sp_vector (scalar_spaces, msh, 'div-preserving');              %   "    Sp1 x Sp2.
space_v  = sp_bspline (knots_L2, degrees_L2, msh, 'integral-preserving',periodic_directions); %   "   S_{p-1, p-1}.  

% INTRODUCIAMO LE MATRICI:
M = op_u_v_tp (sp, sp, msh);% mass matrix
A = op_divu_divv_tp (sp, sp, msh); % (div(Bi),div(Bj))
B = op_div_v_q_tp (sp, space_v, msh); % manca coeff c !

% PROIEZIONE DATI INIZIALI PER P  
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
p1 = c_diff(X2,Y2).*(squeeze(grad_u0_2(1,:,:))'); % seconda comp di p la integro risp x1.
p2 = c_diff(X1,Y1).*(squeeze(grad_u0_1(2,:,:))'); % prima comp di p la integro risp x2.

% PROIEZIONE DATI INIZIALI PER V
Z_dir_1   = space_v.knots{1}(space_v.degree(1) + 1:end-space_v.degree(1)); % nodi interni in dir 1 
Zm_dir_1  = (Z_dir_1(1:end-1)+Z_dir_1(2:end))/2; % più punti medi 
u1        = sort([Z_dir_1 Zm_dir_1]);% ora sono vertici e punti medi ordinati 
Z_dir_2   = space_v.knots{2}(space_v.degree(2)+1:end-space_v.degree(2)); % nodi interni in dir 2 
Zm_dir_2  = (Z_dir_2(1:end-1)+Z_dir_2(2:end))/2; % più punti medi 
u2        = sort([Z_dir_2 Zm_dir_2]);% ora sono vertici e punti medi ordinati 
vtk_pts_3 = {u1, u2};
[v, ~]    = sp_eval (u_0, space_H, geometry_H, vtk_pts_3);

%[X3, Y3]  = deal (squeeze(f3(1,:,:)), squeeze(f3(2,:,:))); %forse non servono.
% questi si incollavano per caso periodico!
v_2proj = [v(:,end-2*space_v.degree(1):end) v(:,2:end-1) v(:,1:2*space_v.degree(1))];                  % in dir 1 
v_2proj = [v_2proj(end-2*space_v.degree(2):end,:); v_2proj(2:end-1,:); v_2proj(1:2*space_v.degree(2),:)]; % in dir 2

% Proiezione dati iniziali: (CONTROLLA!!) sembra buona però!
[Lam1, Lam2] = Lyche_c_2D_Drch_new (p1, p2, space);
Phi          = Lyche_2D_Drch (v_2proj*real(1i*omega*exp(1i*omega*0)), space_v);

L  = [cat(1,Lam1,Lam2) zeros(sp.ndof,round(T/k))];
P  = [Phi zeros(space_v.ndof,round(T/k))];

% PROIEZIONE FUNZIONI DI BASE... (c*Bi)
% dovremo usare [X2,Y2] e [X1,Y1] risp per la prima e la seconda componente
% delle funzioni che vogliamo proiettare 

load(['drchlt_Theta_' num2str(nsub(1)) '_sine.mat']);
%load('drchlt_Theta_' num2str(nsub(1)) '_constant.mat');

% per le derivate: (controlla BENE!!!!)
d1 = space.degree(1); d2 = space.degree(2); %gradi
step_h1 = knots{1}(d1+2:end-1) - knots{1}(2:end-d1-1);
step_h2 = knots{2}(d2+2:end-1) - knots{2}(2:end-d2-1);
alpha1 = d1./step_h1;% peso derivate rispetto a direzione 1
alpha2 = d2./step_h2;% peso derivate rispetto a direzione 2

% COEFFICIENTI DI \pi_g(c*p_h)
cph  = Theta* L(:,1);
cph1 = reshape(cph(1:space_p1.ndof), space_p1.ndof_dir(1), space_p1.ndof_dir(2));
cph2 = reshape(cph(space_p1.ndof+1:end), space_p2.ndof_dir(1), space_p2.ndof_dir(2));

% derivate, quindi d.o.f. negli spazi S_{p,p-1} e S_{p-1,p} (calcolo div)
DLam1 = alpha1'.*(cph1(2:end,:)-cph1(1:end-1,:));
DLam2 = alpha2.*(cph2(:,2:end)-cph2(:,1:end-1));
DL1   = reshape(DLam1, space_v.ndof, 1); % reshape to colum vector
DL2   = reshape(DLam2, space_v.ndof, 1); % reshape to colum vector
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
% risolvo la prima equazione:  
    rhs = mat*L(:,n) - B*P(:,n);% termine noto
    L(:,n+1) = MAT\rhs; % risolviamo il sistema ridotto %NUOVI Lambda.  

    %__ QUI SOTTO SEMBRA OK...   
    % calcola div di questa L 
    
    cph  = Theta* L(:,n+1); %questo full_lam + L(:,n+1) già
    cph1 = reshape(cph(1:space_p1.ndof), space_p1.ndof_dir(1), space_p1.ndof_dir(2));
    cph2 = reshape(cph(space_p1.ndof+1:end), space_p2.ndof_dir(1), space_p2.ndof_dir(2));

    % derivate, quindi d.o.f. negli spazi S_{p,p-1} e S_{p-1,p} (calcolo div)
    DLam1 = alpha1'.*(cph1(2:end,:)-cph1(1:end-1,:));
    DLam2 = alpha2.*(cph2(:,2:end)-cph2(:,1:end-1));
    DL1   = reshape(DLam1, space_v.ndof, 1); % reshape to colum vector
    DL2   = reshape(DLam2, space_v.ndof, 1); % reshape to colum vector
        
    %DIVERGENZA
    Div_L(:,n+1) = DL1 + DL2;
    
    % risolvo la seconda equazione:
    P(:,n+1) = P(:,n) + (k/2)*(Div_L(:,n+1) + Div_L(:,n)); %NUOVO Phi      
end

end

