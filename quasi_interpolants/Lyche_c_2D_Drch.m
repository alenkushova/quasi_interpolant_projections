function [Lam1, Lam2] = Lyche_c_2D_Drch(F,G,space)

% proietteremeo F = (f1,f2) (salvata come array descritto sotto)
% prima proietto f1 con Lyche_c_1D essendo f1 una matrice così salvata:
%
%
%       f1 =  |  f1(x_1,y_1) | f1(x_2,y_1) |  ... 
%             |  f1(x_1,y_2) | f1(x_2,y_2) |  ...
%             |  f1(x_1,y_3) | f1(x_2,y_3) |  ...
%             |  ...             ...
%
% dobbiamo allora usare Lyche_c_1D lungo ogni colonna, così avremo dei 
% d.o.f. colonna per colonna! poi dobbiamo usare la proiezione Lyche_1D
% sui coefficienti trovati
% N.B. usiamo cavalieri-simpson! formuale esatte per spline di grado 3.


%space.degree must have 2 components and must be as follows:
% lungo dir 1 posso usare solo p = 2, 3. 
if space.degree(1) ~= 2 && space.degree(1) ~= 3 
    msg='in dir 1, use only p = 2 oppure p = 3';
    error(msg);
end

% lungo dir 2 posso usare solo p = 2, 3.
if space.degree(2) ~= 2 && space.degree(2) ~= 3 
    msg='in dir 2, use only p = 2 oppure p = 3';
    error(msg);
end

% vorrei proiettare in S_{p,p} = S_p \otimes S_p (Tensor product of same p)
if space.degree(1) ~= space.degree(2) 
    msg='use the same degree in each direction';
    error(msg);
end

% quando f é una funzione/array:
if isa(F,'function_handle') && isa(G,'function_handle')
    % internal knots (prima direzione, tanto seconda è uniforme a prima)
    N = numel(space.knots{1}(space.degree(1)+1:end-space.degree(1)))-1;
    pnt      = linspace (0, 1, 2*N +1);
    pnt_m    = linspace (0, 1, 4*N +1);
    [X1, Y1] = meshgrid (pnt_m, pnt); % x = righe, y = colonne. qui integro in x.
    [X2, Y2] = meshgrid (pnt, pnt_m); % x = righe, y = colonne. qui integro in y.
    f1 = F(X2,Y2); 
    f2 = G(X1,Y1); % handle function
else
    f1 = F; %array over single knots and mid-points. (as in line 25-31) 
    f2 = G;
end

% PROIEZIONI:  \pi_g \otimes \pi_{g-1}^c aplicata a f1
for i = 1 : size(f1,2) % proiettiamo le colonne 
    lam1(:,i) = Lyche_c_1D(f1(:,i), space.sp_univ(2)); 
end
%per costruzione 
coeff1 = [zeros(size(lam1,1), space.sp_univ(1).degree*2)+lam1(:,1) lam1 zeros(size(lam1,1), space.sp_univ(1).degree*2)+lam1(:,end)];
% ora che abbiamo tamto coefficienti facciamo proiezione nell'altro senso.
for i = 1 : size(lam1,1) % proiettiamo le righe
    Lam1(i,:) = Lyche_1D(coeff1(i,:), space.sp_univ(1)); %chiaramente lungo y
end
Lam1 = reshape(Lam1', space.ndof-space.ndof_dir(1), 1); % CONTROLLA che il reshape sia giusto.

% ora rieptiamo per f2 ma prima si usa Lyche_1D poi la proiezione che
% commuta.

% PROIEZIONI:  \pi_{g-1}^c \otimes \pi_g applicata ad f2
%per costruzione 
coeff2 = [zeros(space.sp_univ(2).degree*2, size(f2,2))+f2(1,:); f2 ;zeros(space.sp_univ(2).degree*2, size(f2,2))+f2(end,:)];
for i = 1 : size(f2,2) % proiettiamo le colonne 
    lam2(:,i) = Lyche_1D(coeff2(:,i), space.sp_univ(2)); 
end
% ora che abbiamo tamto coefficienti facciamo proiezione nell'altro senso.
for i = 1 : size(lam2,1) % proiettiamo le righe
    Lam2(i,:) = Lyche_c_1D(lam2(i,:), space.sp_univ(1)); %chiaramente lungo y
end
Lam2 = reshape(Lam2', space.ndof-space.ndof_dir(1), 1); % CONTROLLA che il reshape sia giusto.
end