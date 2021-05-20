function Lam = Lyche_1D_Periodic(f, space)
% controllo dati inseriti correttamente
% questo è un caso particolare in cui p = 2 oppure p = 3 per usare risp le
% formule esplicite con 3 e 5 nodi. Questo dovrebbe velocizzare i conti
% anzi che fare tanti backslash.


if space.degree ~= 2 && space.degree ~= 3 
    msg='use only p = 2 oppure p = 3';
    error(msg);
end

if isa(space.knots, 'cell' )
    Z = space.knots{:}; 
else
    Z = space.knots;
end

Zm = (Z(1:end-1)+Z(2:end))/2; % più punti medi
u = sort([Z Zm]);% ora sono vertici e punti medi ordinati 

% quando f é una funzione/vettore:
if isa(f,'function_handle')
    F = f(u); % handle function
else 
    F = f; %vector over single knots and mid-points.
end

Lam = zeros(space.ndof,1);

for i = 1 : space.ndof
    k = space.degree*2 -1;
    Floc = F(2*i+1: 2*i+k);% termine noto locale.
    if (space.degree==2) % se abbiamo grado = 2 
        Lam(i) = 1/2*(-Floc(1)+4*Floc(2)-Floc(3));
    elseif (space.degree==3) % se abbiamo grado = 3 
        Lam(i) = 1/6*(Floc(1)-8*Floc(2)+20*Floc(3)-8*Floc(4)+Floc(5));
    end
end
end