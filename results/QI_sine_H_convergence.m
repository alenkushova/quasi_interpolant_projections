%% here we perform the h- convergence analysis:
% Let's start by the constant case.
ERR_INF = []; ERR_L2 = []; h = [];
for i = 1 : 4
    N = 2^(i+2); h = [h 1/N];
    % load the workspace
    s = ['Projections\sine\test_c_sine_'...
            num2str(N) 'x' num2str(N) 'x2000.mat'];    load(s);
    % ricompute the errors if you want different evaluations here:
%     err_p = []; err_v = [];
%     D_space_H    = sp_vector ({space_u1, space_u2}, msh); 
%     sp_const_p_H = space_p.constructor(msh_H);
%     sp_const_v_H = space_v.constructor(msh_H);
%     sp_const_D_H = D_space_H.constructor(msh_H);
%     %errore in spazio
%     for n = 1:200:201 %1 : 200 : round(T/k)+1
%         % errore in norma L2 per soluzione p(x,t) per ogni istante t:
%         error_l2_p = my_l2_error (sp_const_p_H, msh_H, L(:,n),...
%                 [Du0_x1; Du0_x2], sp_const_D_H,@(x,y) c(x, y)*real(exp(1i*omega*(n*k-k))) );
%         err_p = [err_p, error_l2_p];
% 
%         % errore in norma L2 per soluzione v(x,t) per ogni istante t:
%         error_l2_v = my_l2_error (sp_const_v_H, msh_H, P(:,n),...
%                 u_0, space_H, @(x, y) real(1i*omega*exp(1i*omega*(n*k-k))));
%         err_v = [err_v, error_l2_v];
%     end
%     %errore in tempo: norma L infinito
%     err_INF_p = max(abs(err_p)); 
%     err_INF_v = max(abs(err_v)); 
%     %errore in tempo: norma L2
%     err_L2_p  = (k/2*sum((err_p(1:end-1)).^2+(err_p(2:end)).^2))^(1/2);
%     err_L2_v  = (k/2*sum((err_v(1:end-1)).^2+(err_v(2:end)).^2))^(1/2);

    
        % save errors: this saves the error i already computed
    ERR_INF = [ERR_INF [err_INF_p; err_INF_v]];
    ERR_L2  = [ERR_L2 [err_L2_p; err_L2_v]];
end
%% Analyze the errors into a plot!
figure
subplot(1,2,1)
loglog(h,ERR_INF(1,:),'-*',h,ERR_INF(2,:),'-*',h,h.^3*30,'Linewidth', 1.5)
title('h-convergence in 2D')
grid on
legend('|| p-p_h ||_{\infty}', '|| v-v_h ||_{\infty}',...
    'h^3','Location','northwest')
subplot(1,2,2)
loglog(h,ERR_L2(1,:),'-*',h,ERR_L2(2,:),'-*',h,h.^3,'Linewidth', 1.5)
title('h-convergence in 2D')
grid on
legend('|| p-p_h ||_2', '|| v-v_h ||_2',...
    'h^3','Location','northwest')