function Lam = Lyche_c_1D_Periodic(f, space)
% controllo dati inseriti correttamente
% questo è un caso particolare in cui p = 2 oppure p = 3 per usare risp le
% formule esplicite con 3 e 5 nodi. Questo dovrebbe velocizzare i conti
% anzi che fare tanti backslash.

if space.degree ~= 2 && space.degree ~= 3 
    msg='use only p = 2 oppure p = 3';
    error(msg);
end

if isa(space.knots, 'cell' )
    knots = space.knots{:};  
else
    knots = space.knots;
end

Z = knots(space.degree+1:end-space.degree); % nodi interni
Zm = (Z(1:end-1)+Z(2:end))/2; % più punti medi 
v  = sort([Z Zm]);% ora sono vertici e punti medi ordinati 
vm = (v(1:end-1)+v(2:end))/2;
u  = sort([v,vm]);
% quando f é una funzione/vettore:
if isa(f,'function_handle')
    F = f(u); % handle function
else 
    F = f; %vector over single knots and mid-points.
end

%calcoliamo la sua media 
s    = numel(F);
h1   = 1/(s-1);
int1  = 0;
F1 = int1;
for i = 1 : (s-1)/2 % il primo so già che fa zero. 
    int1 = int1 + h1/3*(F(2*i-1) +4*F(2*i)+F(2*i+1)); % C-S quad.
    F1(i+1) = int1;
end
%togliamo la media alla integranda F
F = F - int1;
%calcolo integrale 
s    = numel(F);
h1   = 1/(s-1);
int  = 0;
F1 = int;
for i = 1 : (s-1)/2 % il primo so già che fa zero. 
    int = int + h1/3*(F(2*i-1) +4*F(2*i)+F(2*i+1)); % C-S quad.
    F1(i+1) = int;
end

f1 = [F1(end - space.degree*2 : end) F1(2:end) F1(2:space.degree*2+1)];
dof = Lyche_1D_Periodic (f1, space); % call to projector
d = space.degree;
step_h = knots(1+d+1) - knots(2); %g/(xi_{i+g+1}- \xi{i+1})
alpha = d/step_h;
Lam = alpha*[(dof(2:end)-dof(1:end-1)); (dof(1)-dof(end))] + int1;
end