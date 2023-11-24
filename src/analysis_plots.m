figure();
ax1 = subplot(6,1,1);
plot(results.y_CL{1,k_e,1}.')
hold on
plot(results.y_CL{1,k_e,2}.')
plot(ref);
ylabel('y')
grid on;

ax2 = subplot(6,1,2);
plot(results.u_CL{1,k_e,1}.')
hold on
plot(results.u_CL{1,k_e,2}.')
plot(results.du_CL{1,k_e}.');
legend('OL','CL','disturbance')
ylabel('u')
grid on;

ax3 = subplot(6,1,3);
plot(results.eObX{1,k_e,1}); hold on
plot(results.eObX{1,k_e,2});
grid on
ylabel('$\max|\Gamma x_k-\hat{L}z|$','Interpreter','latex')

ax4 = subplot(6,1,4);
semilogy(results.eGu{1,k_e,1}); hold on
semilogy(results.eGu{1,k_e,2});
grid on
ylabel('$||\hat{G}_u-G_u||_\mathrm{F}^2$','Interpreter','latex')

ax5 = subplot(6,1,5);
semilogy(results.eLu{1,k_e,1}); hold on
semilogy(results.eLu{1,k_e,2});
legend('OL','CL')
grid on
ylabel('$||\hat{L}_u-L_u||_\mathrm{F}^2$','Interpreter','latex')

ax6 = subplot(6,1,6);
semilogy(results.eLy{1,k_e,1}); hold on
semilogy(results.eLy{1,k_e,2});
grid on
ylabel('$||\hat{L}_y-L_y||_\mathrm{F}^2$','Interpreter','latex')

linkaxes([ax1,ax2,ax3,ax4,ax5,ax6],'x');