function [Lam1, Lam2] = Lyche_c_2D_Mixed(F,G,space)

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

if isa(F,'function_handle') && isa(G,'function_handle')
    % internal knots (prima direzione, tanto seconda è uniforme a prima)
    N = numel(space.knots{1}(space.degree(1)+1:end-space.degree(1)))-1;
    pnt      = linspace (0, 1, 2*N +1);
    pnt_m    = linspace (0, 1, 4*N +1);
    [X1, Y1] = meshgrid (pnt_m, pnt); % x = righe, y = colonne. qui integro in x
    [X2, Y2] = meshgrid (pnt, pnt_m); % x = righe, y = colonne quindi integro in y.
    f1 = F(X2,Y2); 
    f2 = G(X1,Y1); % handle function
else
    f1 = F; %array over single knots and mid-points. (as in line 25-31) 
    f2 = G;
end

% PROIEZIONI:  (\pi_g drchlt) \otimes (\pi_{g-1}^c prdc), aplicata a f1
% prima proietto le colonne con (\pi_{g-1}^c prdc)
for i = 1 : size(f1,2) % proiettiamo le colonne 
    lam1(:,i) = Lyche_c_1D_Periodic(f1(:,i), space.sp_univ(2)); 
end
% ora assemblo il vettore per la proiezione delle righe,
% cioè sto aggiungendo delle colonne qui:
coeff1 = [zeros(size(lam1,1), space.sp_univ(1).degree*2)+lam1(:,1) lam1 zeros(size(lam1,1), space.sp_univ(1).degree*2)+lam1(:,end)];
% ora che abbiamo i giusti coefficienti facciamo proiezione nell'altro senso.
% proiettiamo ora le righe con (\pi_g drchlt)
for i = 1 : size(lam1,1) % proiettiamo le righe
    Lam1(i,:) = Lyche_1D(coeff1(i,:), space.sp_univ(1)); %chiaramente lungo y
end
Lam1 = reshape(Lam1', space.ndof, 1); % CONTROLLA che il reshape sia giusto.

%Adesso dobbiamo proiettare l'altro campo. Quindi ora uso f2
% PROIEZIONI:  (\pi_{g-1}^c drchlt) \otimes (\pi_g prdc) applicata ad f2
% come prima parto con proiettare le colonne con (\pi_g prdc)
coeff2 = [f2(end - space.degree(1)*2 : end,:); f2(2:end,:); f2(2:space.degree(1)*2+1,:)];
for i = 1 : size(f2,2) % proiettiamo le colonne 
    lam2(:,i) = Lyche_1D_Periodic(coeff2(:,i), space.sp_univ(2)); 
end
% adesso posso subito proiettare le righe con la proiezione (\pi_{g-1}^c drchlt)
for i = 1 : size(lam2,1) % proiettiamo le righe
    Lam2(i,:) = Lyche_c_1D_new(lam2(i,:), space.sp_univ(1)); %chiaramente lungo y
end
Lam2 = reshape(Lam2', space.ndof-space.ndof_dir(2), 1); % CONTROLLA che il reshape sia giusto.
% abbiamo adesso Lam1 e Lam2 che sono i dof dei campi scalari f1 e f2. 
% complessivamente avremo Lam = [Lam1; Lam2] dof di f =(f1 f2) = (F,G).
end