function plot_run(results,ref,Nbar_all,fignum,k_N,k_e)
    
    figure(fignum);
    clf;
    
    Nbar = Nbar_all(k_N);
    num_c = size(results.u_CL,3);

    u_ol = results.u_OL{k_N,k_e};
    y_ol = results.y_OL{k_N,k_e};
    Label = results.CzLabel;
    Color = {'#DC3220','#005AB5'};

    % subplot with outputs & reference
    ax1 = subplot(2,1,1);
    xline(ax1,2*Nbar+0.5,'k--','HandleVisibility','off'); hold on;
    plot(ax1,Nbar+1:Nbar+length(ref),ref,'--',...
        'DisplayName','reference',...
        'color', [.5 .5 .5], ... grey
        'linewidth', 1.5);%      thicker line
    ylabel(ax1,'$y_k$','interpreter','latex');
    grid on
    
    % subplot with inputs
    ax2 = subplot(2,1,2);
    xline(ax2,2*Nbar+0.5,'k--'); hold on;
    yline(ax2,-15,'r--');
    yline(ax2, 15,'r--');
%     ylim(ax2,[-16,16]);
    ylabel(ax2,'$u_k$','interpreter','latex')
    grid on
    
    lineref1 = cell(num_c,1);
    lineref2 = cell(num_c,1);
    Styles = struct('LineStyle',cell(num_c,1),'Color',cell(num_c,1),'LineWidth',cell(num_c,1));
    for k_c = num_c:-1:1
        
        u_cl = results.u_CL{k_N,k_e,k_c};
        y_cl = results.y_CL{k_N,k_e,k_c};
    
        u = [u_ol,u_cl];
        y = [y_ol,y_cl];
                
        lineref1{k_c} = plot(ax1,1:length(y),y,'DisplayName', Label{k_c});
        lineref1{k_c}.LineWidth = 1;
        lineref1{k_c}.Color = Color{k_c};
        Styles(k_c).LineStyle = lineref1{k_c}.LineStyle;
        Styles(k_c).Color     = lineref1{k_c}.Color;
        Styles(k_c).LineWidth = lineref1{k_c}.LineWidth;

        styles2use = namedargs2cell(Styles(k_c));
        lineref2{k_c} = plot(ax2,1:length(u),u,styles2use{:});
    end
    legend(ax1);
    
    linkaxes([ax1 ax2],'x')
    xlim([Nbar,length(u)])
end