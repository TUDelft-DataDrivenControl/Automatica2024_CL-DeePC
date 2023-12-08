figure();
ax1 = subplot(6,1,1);
plot(results.y_CL{k_e,1}.')
hold on
plot(results.y_CL{k_e,2}.')
plot(results.ref);
ylabel('y')
grid on;
if isfield(results,'Nbar')
    xline(results.Nbar,HandleVisibility='off')
end

ax2 = subplot(6,1,2);
plot(results.u_CL{k_e,1}.')
hold on
plot(results.u_CL{k_e,2}.')
plot(results.du_CL{k_e,1}.');
legend('OL','CL','disturbance')
if isfield(results,'Nbar')
    xline(results.Nbar,HandleVisibility='off')
end
ylabel('u')
grid on;

ax3 = subplot(6,1,3);
plot(results.eObX{k_e,1}); hold on
plot(results.eObX{k_e,2});
if isfield(results,'Nbar')
    xline(results.Nbar,HandleVisibility='off')
end
grid on
ylabel('$\max|\Gamma x_k-\hat{L}z|$','Interpreter','latex')

ax4 = subplot(6,1,4);
semilogy(results.eGu{k_e,1}); hold on
semilogy(results.eGu{k_e,2});
if isfield(results,'Nbar')
    xline(results.Nbar,HandleVisibility='off')
end
grid on
ylabel('$||\hat{G}_u-G_u||_\mathrm{F}^2$','Interpreter','latex')

ax5 = subplot(6,1,5);
semilogy(results.eLu{k_e,1}); hold on
semilogy(results.eLu{k_e,2});
legend('OL','CL')
if isfield(results,'Nbar')
    xline(results.Nbar,HandleVisibility='off')
end
grid on
ylabel('$||\hat{L}_u-L_u||_\mathrm{F}^2$','Interpreter','latex')

ax6 = subplot(6,1,6);
semilogy(results.eLy{k_e,1}); hold on
semilogy(results.eLy{k_e,2});
if isfield(results,'Nbar')
    xline(results.Nbar,HandleVisibility='off')
end
grid on
ylabel('$||\hat{L}_y-L_y||_\mathrm{F}^2$','Interpreter','latex')

linkaxes([ax1,ax2,ax3,ax4,ax5,ax6],'x');