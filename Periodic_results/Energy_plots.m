function Energy_plots(nsub, timestep)
%% this script performs the graphics for the evolution of energy conservation.
T = 30; % because energy simulations were run till T = 30.
filename = ['test_energy\test_c_sine_energy_' num2str(nsub) 'x' ...
                num2str(nsub) 'x' num2str(timestep) '.mat'];
load( filename, 'space_p', 'space_v', 'msh', 'L', 'P')
% compute evolution of energy conservation:
counter = 1; energy = []; nrg = [];
for j = 1 + (timestep / T) : timestep / T : timestep +1
[nrg_p] = sp_l2_error (space_p, msh, L(:,j), @(x,y) [0*x, 0*y]);
[nrg_v] = sp_l2_error (space_v, msh, P(:,j), @(x,y) 0*x);
nrg(counter) = (nrg_p^2+nrg_v^2)/2;
energy = [energy abs(nrg(1)-nrg(counter))/abs(nrg(1))];
counter = counter +1;
end
%% plotting the graph
figure
semilogy(1:T, energy,'--ks','MarkerEdgeColor','r',...
    'MarkerFaceColor','r','LineWidth',1.5)
grid on 
title('Energy evolution');
legend('Relative error', 'Location', 'NorthWest' )
xlabel('Time'); ylabel('log_{10} of |E_0 - E_t|/ |E_0|');
