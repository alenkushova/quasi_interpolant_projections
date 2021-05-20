function Lam = Lyche_2D_Periodic(f, space)
%
% This code projects the scalar field f in the discrete space of H^1
% with the notation of De Rham 2D spline complex the discrete space is 
% S_{p,p} and the Projection is called \Pi^0.
% 
%    dof = Lyche_2D_Periodic (f, space);
%
%   Input :  f      = function handle or vector of evaluation of f on 
%                      knots and midpoints;
%            space  = space scalar object of the discrete De Rham complex.
%
%   Output:  dof    = degrees of freedom of the projection.
%
%__________________________________________________________________________
% questo è un caso particolare in cui p = 2 oppure p = 3 per usare risp le
% formule esplicite con 3 e 5 nodi. Questo dovrebbe velocizzare i conti
% anzi che fare tanti backslash.
% Caso 2D per f:\Omega \subset \mathbb{R}^2 --> \mathbb{R}
%__________________________________________________________________________
% N.B. Quando devo passare un array F, che è valutato solo su
% nodi e punti medi interni, devo operare prima in questo modo
% per incollare i valori periodici di F alle 2 estremità del knot vector
%
% full_F = [F(:,end-2*space.degree:end) F(:,2:end-1) F(:,1:2*space.degree)];                  % in dir 1 
% full_F = [full_F(end-2*space.degree:end,:); full_F(2:end-1,:); full_F(1:2*space.degree,:)]; % in dir 2
%
% solo dopo chiamo     " Lam = Lyche_2D_Periodic(full_F, space) " 

% space.degree must have 2 components and must be as follows:
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

Lam = zeros(space.ndof,1);
stancil_2 = [1/4 -1 1/4;-1 4 -1;1/4 -1 1/4]; 
stancil_3 = 1/36*[1 -8 20 -8 1].*[1 -8 20 -8 1]';
h = 1;

for j = 1 : space.ndof_dir(2)
    for i = 1 : space.ndof_dir(1)
        k = space.degree*2 -1;
        Floc = F(2*j+1: 2*j+k,2*i+1: 2*i+k);% termine noto locale.
        if (space.degree(1)==2) % se abbiamo grado = 2 
            Lam(h) = sum(sum(stancil_2.*Floc));
        elseif (space.degree(1)==3) % se abbiamo grado = 3 
            Lam(h) = sum(sum(stancil_3.*Floc));
          %  Lam(i) = 1/6*(Floc(1)-8*Floc(2)+20*Floc(3)-8*Floc(4)+Floc(5));
        end
        h = h+1;
    end
end
end
