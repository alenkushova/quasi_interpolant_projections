%% here we perform the k-convergence study:
% Let's start by the constant case.
ERR_INF = []; ERR_L2 = []; h = [];
for i = 1 : 4
    N = 2^(i+2); h = [h 1/N];
    % load the workspace
    s = ['Projections\sine\test_c_sine_'...
            num2str(N) 'x' num2str(N) 'x' num2str(N) '.mat'];    load(s);    
        % save errors: this saves the errors already computed
    ERR_INF = [ERR_INF [err_INF_p; err_INF_v]];
    ERR_L2  = [ERR_L2 [err_L2_p; err_L2_v]];
end
%% Analyze the errors into a plot!
figure
subplot(1,2,1)
loglog(h,ERR_INF(1,:),'-*',h,ERR_INF(2,:),'-*',h,h.^2*900,'Linewidth', 1.5)
title('k-convergence in 2D')
grid on;
legend('|| p-p_h ||_{\infty}', '|| v-v_h ||_{\infty}',...
    'h^3','Location','northwest')

subplot(1,2,2)
loglog(h,ERR_L2(1,:),'-*',h,ERR_L2(2,:),'-*',h,h.^2*400,'Linewidth', 1.5)
title('k-convergence in 2D')
grid on; 
legend('|| p-p_h ||_2', '|| v-v_h ||_2', 'h^3',...
    'Location','northwest')
