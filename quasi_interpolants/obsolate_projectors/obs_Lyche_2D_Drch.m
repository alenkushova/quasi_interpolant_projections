function Lam = Lyche_2D_Drch(f, space)
% controllo dati inseriti correttamente
% questo è un caso particolare in cui p = 2 oppure p = 3 per usare risp le
% formule esplicite con 3 e 5 nodi. Questo dovrebbe velocizzare i conti
% anzi che fare tanti backslash.
% Caso 2D per f:\Omega \subset \mathbb{R}^2 --> \mathbb{R}

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

%incolliamo i valori di F 
Lam = zeros(space.ndof_dir(2), space.ndof_dir(1));
stancil_2 = [1/4 -1 1/4;-1 4 -1;1/4 -1 1/4]; 
stancil_3 = 1/36*[1 -8 20 -8 1].*[1 -8 20 -8 1]';

for j = space.degree(2) : space.ndof_dir(2)-(space.degree(2)-1)
    for i = space.degree(1) : space.ndof_dir(1)-(space.degree(2)-1)
        k = space.degree*2 -1;
        Floc = F(2*j+1: 2*j+k,2*i+1: 2*i+k);% termine noto locale.
        
        if (space.degree(1)==2) % se abbiamo grado = 2 
            Lam(j,i) = sum(sum(stancil_2.*Floc));
            if (j == space.degree(2))
                if(i == space.degree(1))
                    Lam(j-1,i-1) = sum(sum([1 0 0]'.*[1 0 0].*Floc));
                    Lam(j-1,i)   = sum(sum([1 0 0]'.*[-1/2 2 -1/2].*Floc));
                    Lam(j,i-1)   = sum(sum([-1/2 2 -1/2]'.*[1 0 0].*Floc));
                elseif(i == space.ndof_dir(1)-(space.degree(1)-1))
                    Lam(j-1,i)   = sum(sum([1 0 0]'.*[-1/2 2 -1/2].*Floc));
                    Lam(j-1,i+1) = sum(sum([1 0 0]'.*[0 0 1].*Floc));
                    Lam(j,i+1)   = sum(sum([-1/2 2 -1/2]'.*[0 0 1].*Floc));
                else 
                    Lam(j-1,i)   = sum(sum([1 0 0]'.*[-1/2 2 -1/2].*Floc));
                end
            elseif(j == space.ndof_dir(2)-(space.degree(2)-1))
                if(i == space.degree(1))
                    Lam(j+1,i-1) = sum(sum([0 0 1]'.*[1 0 0].*Floc));
                    Lam(j+1,i)   = sum(sum([0 0 1]'.*[-1/2 2 -1/2].*Floc));
                    Lam(j,i-1)   = sum(sum([-1/2 2 -1/2]'.*[1 0 0].*Floc));
                elseif(i == space.ndof_dir(1)-(space.degree(1)-1))
                    Lam(j+1,i)   = sum(sum([0 0 1]'.*[-1/2 2 -1/2].*Floc));
                    Lam(j+1,i+1) = sum(sum([0 0 1]'.*[0 0 1].*Floc));
                    Lam(j,i+1)   = sum(sum([-1/2 2 -1/2]'.*[0 0 1].*Floc));
                else
                    Lam(j+1,i)   = sum(sum([0 0 1]'.*[-1/2 2 -1/2].*Floc));
                end
            elseif(i == space.degree(1))
                Lam(j,i-1) = sum(sum([-1/2 2 -1/2]'.*[1 0 0].*Floc));           
            elseif(i == space.ndof_dir(1)-(space.degree(1)-1))
                Lam(j,i+1) = sum(sum([-1/2 2 -1/2]'.*[0 0 1].*Floc));           
            end
        elseif (space.degree(1)==3) % se abbiamo grado = 3
            Lam(j,i) = sum(sum(stancil_3.*Floc));
            if (j == space.degree(2))
                if(i == space.degree(1))
                    Lam(j-2,i-2) = sum(sum([1 0 0 0 0]'.*[1 0 0 0 0].*Floc));
                    Lam(j-2,i-1) = sum(sum([1 0 0 0 0]'.*1/18.*[-5 40 -24 8 -1].*Floc));
                    Lam(j-2,i)   = sum(sum([1 0 0 0 0]'.*1/6.*[1 -8 20 -8 1].*Floc));
                    Lam(j-1,i-2) = sum(sum(1/18*[-5 40 -24 8 -1]'.*[1 0 0 0 0].*Floc));
                    Lam(j-1,i-1) = sum(sum(1/18*[-5 40 -24 8 -1]'.*1/18.*[-5 40 -24 8 -1].*Floc));
                    Lam(j-1,i)   = sum(sum(1/18*[-5 40 -24 8 -1]'.*1/6.*[1 -8 20 -8 1].*Floc));
                    Lam(j,i-2)   = sum(sum(1/6*[1 -8 20 -8 1]'.*[1 0 0 0 0].*Floc));
                    Lam(j,i-1)   = sum(sum(1/6*[1 -8 20 -8 1]'.*1/18.*[-5 40 -24 8 -1].*Floc));
                elseif(i == space.ndof_dir(1)-(space.degree(1)-1))
                    Lam(j-2,i+2) = sum(sum([1 0 0 0 0]'.*[0 0 0 0 1].*Floc));
                    Lam(j-2,i+1) = sum(sum([1 0 0 0 0]'.*1/18.*[-1 8 -24 40 -5].*Floc));
                    Lam(j-2,i)   = sum(sum([1 0 0 0 0]'.*1/6.*[1 -8 20 -8 1].*Floc));
                    Lam(j-1,i+2) = sum(sum(1/18*[-5 40 -24 8 -1]'.*[0 0 0 0 1].*Floc));
                    Lam(j-1,i+1) = sum(sum(1/18*[-5 40 -24 8 -1]'.*1/18.*[-1 8 -24 40 -5].*Floc));
                    Lam(j-1,i)   = sum(sum(1/18*[-5 40 -24 8 -1]'.*1/6.*[1 -8 20 -8 1].*Floc));
                    Lam(j,i+2)   = sum(sum(1/6*[1 -8 20 -8 1]'.*[0 0 0 0 1].*Floc));
                    Lam(j,i+1)   = sum(sum(1/6*[1 -8 20 -8 1]'.*1/18.*[-1 8 -24 40 -5].*Floc));
                else 
                    Lam(j-2,i)   = sum(sum([1 0 0 0 0]'.*1/6.*[1 -8 20 -8 1].*Floc));
                    Lam(j-1,i)   = sum(sum(1/18*[-5 40 -24 8 -1]'.*1/6.*[1 -8 20 -8 1].*Floc));
                end
            elseif(j == space.ndof_dir(2)-(space.degree(2)-1))
                if(i == space.degree(1))
                    Lam(j+2,i-2) = sum(sum([0 0 0 0 1]'.*[1 0 0 0 0].*Floc));
                    Lam(j+2,i-1) = sum(sum([0 0 0 0 1]'.*1/18.*[-5 40 -24 8 -1].*Floc));
                    Lam(j+2,i)   = sum(sum([0 0 0 0 1]'.*1/6.*[1 -8 20 -8 1].*Floc));
                    Lam(j+1,i-2) = sum(sum(1/18*[-1 8 -24 40 -5]'.*[1 0 0 0 0].*Floc));
                    Lam(j+1,i-1) = sum(sum(1/18*[-1 8 -24 40 -5]'.*1/18.*[-5 40 -24 8 -1].*Floc));
                    Lam(j+1,i)   = sum(sum(1/18*[-1 8 -24 40 -5]'.*1/6.*[1 -8 20 -8 1].*Floc));
                    Lam(j,i-2)   = sum(sum(1/6*[1 -8 20 -8 1]'.*[1 0 0 0 0].*Floc));
                    Lam(j,i-1)   = sum(sum(1/6*[1 -8 20 -8 1]'.*1/18.*[-5 40 -24 8 -1].*Floc));
                elseif(i == space.ndof_dir(1)-(space.degree(1)-1))
                    Lam(j+2,i+2) = sum(sum([0 0 0 0 1]'.*[0 0 0 0 1].*Floc));
                    Lam(j+2,i+1) = sum(sum([0 0 0 0 1]'.*1/18.*[-1 8 -24 40 -5].*Floc));
                    Lam(j+2,i)   = sum(sum([0 0 0 0 1]'.*1/6.*[1 -8 20 -8 1].*Floc));
                    Lam(j+1,i+2) = sum(sum(1/18*[-1 8 -24 40 -5]'.*[0 0 0 0 1].*Floc));
                    Lam(j+1,i+1) = sum(sum(1/18*[-1 8 -24 40 -5]'.*1/18.*[-1 8 -24 40 -5].*Floc));
                    Lam(j+1,i)   = sum(sum(1/18*[-1 8 -24 40 -5]'.*1/6.*[1 -8 20 -8 1].*Floc));
                    Lam(j,i+2)   = sum(sum(1/6*[1 -8 20 -8 1]'.*[0 0 0 0 1].*Floc));
                    Lam(j,i+1)   = sum(sum(1/6*[1 -8 20 -8 1]'.*1/18.*[-1 8 -24 40 -5].*Floc));
                else 
                    Lam(j+2,i)   = sum(sum([0 0 0 0 1]'.*1/6.*[1 -8 20 -8 1].*Floc));
                    Lam(j+1,i)   = sum(sum(1/18*[-1 8 -24 40 -5]'.*1/6.*[1 -8 20 -8 1].*Floc));
                end
            elseif(i == space.degree(1))
                Lam(j,i-2) = sum(sum(1/6.*[1 -8 20 -8 1]'.*[1 0 0 0 0].*Floc));           
                Lam(j,i-1) = sum(sum(1/6.*[1 -8 20 -8 1]'.*1/18.*[-5 40 -24 8 -1].*Floc));           
            elseif(i == space.ndof_dir(1)-(space.degree(1)-1))
                Lam(j,i+2) = sum(sum(1/6.*[1 -8 20 -8 1]'.*[0 0 0 0 1].*Floc));           
                Lam(j,i+1) = sum(sum(1/6.*[1 -8 20 -8 1]'.*1/18.*[-1 8 -24 40 -5].*Floc));           
            end
        end
    end
end
Lam = reshape(Lam',1,space.ndof)';

%% N.B. possiamo velocizzare forse il tutto fecendo proprio prodotto tensore con il caso Lyche_1D in ambo le direzioni.
% in questo caso fanno comunque la stessa cosa, cambia quando faccio il
% prodotto tensore tra Lyche_1D e Lyche_c_1D.
end

