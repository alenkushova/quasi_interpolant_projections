function dof = Lyche_2D_Drch(f, space)
% This code projects the scalar field f in the discrete space of H^1
% with the notation of De Rham 2D spline complex the discrete space is 
% S_{p,p} and the Projection is called \Pi^0.
% 
%    dof = Lyche_2D_Drch (f, space);
%
%   Input :  f      = function handle or vector of evaluation of f on 
%                      knots and midpoints;
%            space  = space scalar object of the discrete De Rham complex.
%
%   Output:  dof    = degrees of freedom of the projection.
%
%__________________________________________________________________________           

% controlliamo che siamo nel caso polinomi di gradi 2 e 3
% lungo la prima direzione univariata.
if space.degree(1) ~= 2 && space.degree(1) ~= 3 
    msg='in dir 1, use only p = 2 oppure p = 3';
    error(msg);
end
% e lungo la seconda direzione univariata.
if space.degree(2) ~= 2 && space.degree(2) ~= 3 
    msg='in dir 2, use only p = 2 oppure p = 3';
    error(msg);
end

% vorrei proiettare in S_{p,p} = S_p \otimes S_p (Tensor product of same p)
if space.degree(1) ~= space.degree(2) 
    msg='use the same degree in each direction';
    error(msg);
end

% quando f é una funzione/vettore:
if isa(f,'function_handle')
    %dir 1
    Z_dir_1 = space.knots{1}; % nodi interni in dir 1 
    Zm_dir_1 = (Z_dir_1(1:end-1)+Z_dir_1(2:end))/2; % più punti medi 
    u1 = sort([Z_dir_1 Zm_dir_1]);% ora sono vertici e punti medi ordinati 
    %dir 2
    Z_dir_2 = space.knots{2}; % nodi interni in dir 2 
    Zm_dir_2 = (Z_dir_2(1:end-1)+Z_dir_2(2:end))/2; % più punti medi 
    u2 = sort([Z_dir_2 Zm_dir_2]);% ora sono vertici e punti medi ordinati 

    %mesh di valutazione di f in 2D
    
    [XX, YY] = meshgrid(u1,u2); % ordinamento lessico grafico in matrice 
    % riferimento dal left-down corner come primo punto
    F = f(XX,YY); % handle function
else 
    F = f; %array over single knots and mid-points.
end

% proiettiamo prima le colonne:
for i = 1 : size(F,2) % proiettiamo le colonne 
    lam1(:,i) = Lyche_1D(F(:,i), space.sp_univ(2)); 
end
% ora proiettiamo le righe:
for i = 1 : size(lam1,1) % proiettiamo le righe
    Lam1(i,:) = Lyche_1D(lam1(i,:), space.sp_univ(1)); %chiaramente lungo y
end
dof = reshape(Lam1', space.ndof, 1); % CONTROLLA che il reshape sia giusto.
end

