
%% adding calculated quantities of interest
results = calc_EfUf_corr(results,k_e);

%% plotting
num_plots_fig1 = 5;

figure();
ax1 = subplot(num_plots_fig1,1,1);
plot(results.y_CL{k_e,1}.')
hold on
plot(results.y_CL{k_e,2}.')
plot(results.y_CL{k_e,3}.')
plot(results.ref);
legend('OL','CL','oracle','reference')
ylabel('y')
grid on;
xline(results.Nbar,HandleVisibility='off')

ax2 = subplot(num_plots_fig1,1,2);
plot(results.u_CL{k_e,1}.')
hold on
plot(results.u_CL{k_e,2}.')
plot(results.u_CL{k_e,3}.')
plot(results.du_CL{k_e,1}.');
legend('OL','CL','oracle','disturbance')
xline(results.Nbar,HandleVisibility='off')
ylabel('u')
grid on;

ax3 = subplot(num_plots_fig1,1,3);
plot(results.eObX{k_e,1}); hold on
plot(results.eObX{k_e,2});
xline(results.Nbar,HandleVisibility='off')
grid on
ylabel('$\max|\Gamma x_k-\hat{L}z|$','Interpreter','latex')

ax4 = subplot(num_plots_fig1,1,4);
semilogy(results.eGu{k_e,1}); hold on
semilogy(results.eGu{k_e,2});
xline(results.Nbar,HandleVisibility='off')
grid on
ylabel('$||\hat{G}_u-G_u||_\mathrm{F}^2$','Interpreter','latex')

ax5 = subplot(num_plots_fig1,1,5);
semilogy(results.EfUfnorm{k_e,1}); hold on;
semilogy(results.EfUfnorm{k_e,2});
xline(results.Nbar,HandleVisibility='off')
grid on
ylabel('$||E_fU_f^\top||_\mathrm{F}^2/f_\mathrm{ID}$','Interpreter','latex')

% ax5 = subplot(num_plots_fig1,1,5);
% semilogy(results.eLu{k_e,1}); hold on
% semilogy(results.eLu{k_e,2});
% legend('OL','CL','oracle')
% if isfield(results,'Nbar')
%     xline(results.Nbar,HandleVisibility='off')
% end
% grid on
% ylabel('$||\hat{L}_u-L_u||_\mathrm{F}^2$','Interpreter','latex')
% 
% ax6 = subplot(6,1,6);
% semilogy(results.eLy{k_e,1}); hold on
% semilogy(results.eLy{k_e,2});
% if isfield(results,'Nbar')
%     xline(results.Nbar,HandleVisibility='off')
% end
% grid on
% ylabel('$||\hat{L}_y-L_y||_\mathrm{F}^2$','Interpreter','latex')

linkaxes([ax1,ax2,ax3,ax4,ax5],'x');