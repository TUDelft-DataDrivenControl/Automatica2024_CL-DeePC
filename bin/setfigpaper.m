function setfigpaper(varargin)
%
% SETFIGPAPER applies paper style to your figure.
% Call this function after making all the details of the figure: plots, legends, labels, etc.
%
% OPTIONS:
%
%   - 'Width'       Simple: one integer. With aspect ratio: array with two
%                   elements, the width and the aspect ratio.
%                   Examples: ('Width',15) or ('Width',[15 0.7])
%                   Default: [19.7556 0.75]
%   - 'FontSize'    The font size in all texts in the figure.
%                   Examples: ('FontSize',10)
%                   Default: 10
%   - 'FontName'    Font family.
%                   Examples: ('FontName','Arial')
%                   Default: 'Helvetica'
%   - 'Interpreter' Text interpreter syntax, 'Latex', 'tex' or 'none'.
%                   If you don't want to modify the interpreter the value must be [] (empty value).
%                   Examples: ('Interpreter','Latex')
%                   Default: []
%   - 'Grayscale'   Change all colours to grayscale. [boolean]
%                   Examples: ('Grayscale',1)
%                   Default: 0
%   - 'LineWidth'   Set line width only for axes, not for data.
%                   Examples: ('LineWidth',2)
%                   Default: 0.5
%   - 'Figure'      To apply the function to a specific figure.
%                   Examples: ('Figure',gcf)
%                   Default: Last figure created
%
% EXAMPLES WITH NAME-VALUE PAIR ARGUMENTS:
%
%   SETFIGPAPER
%   SETFIGPAPER('Width',20,'Interpreter','Latex')
%   SETFIGPAPER('Width',[20,0.65],'Interpreter','Latex')
%   SETFIGPAPER('Width',[20,0.65],'FontSize',13)
%   SETFIGPAPER('Grayscale',true)
%   ...
%
% EXAMPLES WITH VALUES ONLY (same examples as above):
%
%   SETFIGPAPER
%   SETFIGPAPER(20,10,'Helvetica','Latex')
%   SETFIGPAPER([20,0.65],10,'Helvetica','Latex')
%   SETFIGPAPER([20,0.65],13)
%   SETFIGPAPER([15 0.85],10,'Helvetica','Latex',true)
%   ...
%
% TESTED WITH:
%   - Line Plots:   plot, plot3, stairs, errorbar, area, stackedplot, loglog,
%                   semilogy, semilogx, fplot, fimplicit, triplot.
%   - Polar Plots:  polarplot, polarhistogram, polarscatter, compass, ezpolar.
%   - Surf/Mesh:    surf, surfc, surfl, ribbon, pcolor, fsurf, fimplicit3, mesh,
%                   meshc, meshz, waterfall, fmesh, trimesh, trisurf.
%   - Distribution: histogram, histogram2, pie (view troubleshooting), pie3 (view troubleshooting),
%                   scatter, scatter3, scatterhistogram, spy, plotmatrix, heatmap,
%                   wordcloud, parallelplot.
%   - Discrete:     bar, bar3, barh, bar3h, pareto (view troubleshooting), stem,
%                   stem3, scatter, scatter3, stairs.
%   - Geographic:   geobubble, geoplot, geoscatter.
%   - Contour:      contour, contourf, contour3, contourslice, fcontour.
%   - Vector field: quiver, quiver3, feather.
%   - Volume:       streamline, streamslice, streamparticles, streamribbon,
%                   streamtube, coneplot.
%   - Animation:    animatedline (add 'setfigpaper' after 'drawnow'),
%                   comet (only when the animation ends),
%                   comet3 (only when the animation ends).
%   - Images        image, imagesc.
%
% TROUBLESHOOTING:
%   - pie and pie3: If you use latex interpreter with default labels (percentage) you
%                   need to fix the plot before call setfigpaper. Solution:
%                   p = pie([1 3 0.5 2.5 2]);
%                   txt = findobj(p,'Type','Text');
%                   arrayfun(@(x) set(x,'String',strrep(x.String,'%','\%')),txt,...
%                   'UniformOutput',false)
%                   setfigpaper
%   - pareto:       If you use latex interpreter with default ticklabels (percentage) you
%                   need to fix the plot before call setfigpaper. Solution:
%                   [H,ax]  = pareto([200 120 555 608], {'Fred','Ginger','Norman','Max'});
%                   ax(2).YTickLabel = arrayfun(@(x) strrep(x,'%','\%'),ax(2).YTickLabel);
%                   setfigpaper
%
% Author: Jose M. Requena Plens, 2020. (info@jmrplens.com   | joreple@upv.es)
% Original concept: Noe Jimenez, 2014. (noe.jimenez@csic.es | nojigon@i3m.upv.es)
%

% ========================================================================================
% Inputs
p = inputParser;
addOptional(p,'Width',[19.7556 0.75],...        % Default: [19.7556 0.75]
    @(x) length(x)<=2);
addOptional(p,'FontSize',10,...                 % Default: 10
    @(x) isnumeric(x) && isscalar(x));
addOptional(p,'FontName','Helvetica',@ischar);	% Default: Helvetica
addOptional(p,'Interpreter',[],...
    @(x) ischar(x) || isempty(x));              % Default: []
addOptional(p,'Grayscale',0,...                 % Default: false or 0
    @(x) x==0 || x==1);
addOptional(p,'LineWidth',0.5,...               % Default: 0.5
    @(x) isnumeric(x) && isscalar(x) && x>0);
addOptional(p,'Figure',gcf,...                  % Default: Actual gcf
    @(x) isobject(x));
parse(p,varargin{:});
% Unpacking
fwidth    = p.Results.Width;
fsize     = p.Results.FontSize;
font      = p.Results.FontName;
type      = p.Results.Interpreter;
lwidth    = p.Results.LineWidth;
grayscale = p.Results.Grayscale;
fig       = p.Results.Figure;
if length(fwidth) == 1
    aspect = 0.75;
else
    aspect = fwidth(2);
end
if isempty(type)
    notype = 1; % Unmodified interpreter
else
    validatestring(type,{'tex','latex','none'});
    notype = 0; % Modified interpreter
end
% Common Colors
Black = [0 0 0];
White = [1 1 1];

% ========================================================================================
% General
% Set initial properties
set(fig,'WindowStyle','normal');
set(fig,'paperunits','centimeters');
set(fig,'units','centimeters');
% Set new size and new position
set(fig,'PaperSize',[fwidth(1) fwidth(1)*aspect]);
pos=get(fig,'position'); % Get actual position
set(fig,'Position',[pos(1) pos(2) fwidth(1) fwidth(1)*aspect]);
% General color properties
set(fig,'color',White); % General Background
set(findobj(fig,'color',[0.15 0.15 0.15]),'color',Black); % Gray to black
% General text
txt = findobj(fig, 'Type', 'Text');
arrayfun(@(x) set(x,'Fontname',font,'FontSize', fsize), txt) ;
if notype == 0
    arrayfun(@(x) set(x,'Interpreter', type), txt) ;
end

set(0,'ShowHiddenHandles','on'); % For hidden objects (e.g. annotations)

% ========================================================================================
% Common axes
haxis=findobj(fig,'type','axes');
for iaxis=1:length(haxis)
    % Axes linewidth
    set(haxis(iaxis),'linewidth',lwidth);
    % Color Properties
    set(haxis(iaxis),'gridcolor',Black,'minorgridcolor',Black,...
        'xcolor',Black,'ycolor',Black,'zcolor',Black);
    if ~strcmp(haxis(iaxis).Color,'none')
        set(haxis(iaxis),'color',White); % Background
    end
    % Details
    set(haxis(iaxis),'box','on');
    % Text properties
    ht=findobj(haxis(iaxis),'type','text'); % Get all text
    set(ht,'fontname',font,'fontsize',fsize,'color',Black);
    set(haxis(iaxis),'fontname',font,'fontsize',fsize)
    set(haxis(iaxis),'XMinorTick','on','YMinorTick','on','ZMinorTick','on')
    if notype == 0
        set(ht,'interpreter',type);
        set(haxis(iaxis),'TickLabelInterpreter',type)
    end
    % Axis labels
    set([haxis(iaxis).XLabel haxis(iaxis).YLabel haxis(iaxis).ZLabel],...
        'fontname',font,'fontsize',fsize,'color',Black);
    if notype == 0
        set([haxis(iaxis).XLabel haxis(iaxis).YLabel haxis(iaxis).ZLabel],...
            'interpreter',type);
    end
    % Title
    set(haxis(iaxis).Title,'fontname',font,'fontsize',fsize,'color',Black);
    if notype == 0
        set(haxis(iaxis).Title,'interpreter',type);
    end
    % Text for nodes for graphplots
    if ~verLessThan('matlab','9.6') % Only for Matlab 2018 and later
        grap = findobj(haxis(iaxis),'Type','graphplot');
        set(grap,'nodefontname',font,'nodefontsize',fsize,...
            'EdgeFontSize',fsize-2,'EdgeFontName',font)
        if notype == 0
            set(grap,'interpreter',type)
        end
    end
    % Grayscale
    if grayscale == 1
        set(fig,'Colormap', rgb2gray(get(fig, 'Colormap')));
        haxis(iaxis).ColorOrder = rgb2gray(haxis(iaxis).ColorOrder);
        
        for jj = 1:length(haxis(iaxis).Children)
            for ii = length(haxis(iaxis).YAxis) % 1 or 2 Yaxis
                if length(haxis(iaxis).YAxis)==2
                    if ii==1
                        yyaxis right
                    elseif ii==2
                        yyaxis left
                    end
                end
                try %#ok<TRYNC>
                    haxis(iaxis).Children(jj).Color = ...
                        rgb2gray(haxis(iaxis).Children(jj).Color);
                end
                try %#ok<TRYNC>
                    haxis(iaxis).Children(jj).EdgeColor = ...
                        rgb2gray(haxis(iaxis).Children(jj).EdgeColor);
                end
                try %#ok<TRYNC>
                    haxis(iaxis).Children(jj).FaceColor = ...
                        rgb2gray(haxis(iaxis).Children(jj).FaceColor);
                end
                try %#ok<TRYNC>
                    haxis(iaxis).Children(jj).MarkerFaceColor = ...
                        rgb2gray(haxis(iaxis).Children(jj).MarkerFaceColor);
                end
                try %#ok<TRYNC>
                    haxis(iaxis).Children(jj).MarkerEdgeColor = ...
                        rgb2gray(haxis(iaxis).Children(jj).MarkerEdgeColor);
                end
                try %#ok<TRYNC>
                    haxis(iaxis).Children(jj).CData = ...
                        rgb2gray(haxis(iaxis).Children(jj).CData);
                end
                try %#ok<TRYNC>
                    haxis(iaxis).Children(jj).NodeColor = ...
                        rgb2gray(haxis(iaxis).Children(jj).NodeColor);
                end
            end
        end
        % Two axes
        if length(haxis(iaxis).YAxis)==2
            haxis(iaxis).YAxis(1).Color = rgb2gray(haxis(iaxis).YAxis(1).Color);
            haxis(iaxis).YAxis(2).Color = rgb2gray(haxis(iaxis).YAxis(2).Color);
        end
    end
    % Other texts (e.g. values in countour plot)
    drawnow
    warning('off','MATLAB:structOnObject');warning('off','MATLAB:hg:EraseModeIgnored')
    warning('off','MATLAB:hg:willberemoved')
    ax = struct(haxis(iaxis).Children);
    warning('on','MATLAB:structOnObject');warning('on','MATLAB:hg:EraseModeIgnored')
    warning('off','MATLAB:hg:willberemoved')
    if ~isempty(ax)
        for jj = 1:length(ax.NodeChildren)
            if isa(ax.NodeChildren(jj),'matlab.graphics.primitive.world.Text')
                ax.NodeChildren(jj).Font.Name = font;
                ax.NodeChildren(jj).Font.Size  = fsize;
                if notype == 0
                    ax.NodeChildren(jj).Interpreter= type;
                end
            end
        end
    end
end

% ========================================================================================
% Polar axes
polaxis=findobj(fig,'type','polaraxes');
for iaxis=1:length(polaxis)
    % Polar Text
    set(polaxis(iaxis).ThetaAxis,'fontname',font,'fontsize',fsize,'color',Black);
    set(polaxis(iaxis).ThetaAxis.Label,'fontname',font,'fontsize',fsize,'color',Black);
    if notype == 0
        set(polaxis(iaxis).ThetaAxis.Label,'interpreter',type);
    end
    % Radial Text
    set(polaxis(iaxis).RAxis,'fontname',font,'fontsize',fsize,'color',Black);
    set(polaxis(iaxis).RAxis.Label,'fontname',font,'fontsize',fsize,'color',Black);
    if notype == 0
        set(polaxis(iaxis).RAxis.Label,'interpreter',type);
    end
    % Title
    set(polaxis(iaxis).Title,'fontname',font,'fontsize',fsize,'color',Black);
    if notype == 0
        set(polaxis(iaxis).Title,'interpreter',type);
    end
    % General axis
    set(polaxis(iaxis),...
        'LineWidth',lwidth,'Color',White,'GridColor',Black,'MinorGridColor',Black,...
        'RColor',Black,'ThetaColor',Black)
    if notype == 0
        set(polaxis(iaxis),'TickLabelInterpreter',type)
    end
    % Grayscale
    if grayscale == 1
        set(fig,'Colormap', rgb2gray(get(fig, 'Colormap')));
        polaxis(iaxis).ColorOrder = rgb2gray(polaxis(iaxis).ColorOrder);
        for jj = 1:length(polaxis(iaxis).Children)
            try %#ok<TRYNC>
                polaxis(iaxis).Children(jj).Color = ...
                    rgb2gray(polaxis(iaxis).Children(jj).Color);
            end
            try %#ok<TRYNC>
                polaxis(iaxis).Children(jj).EdgeColor = ...
                    rgb2gray(polaxis(iaxis).Children(jj).EdgeColor);
            end
            try %#ok<TRYNC>
                polaxis(iaxis).Children(jj).FaceColor = ...
                    rgb2gray(polaxis(iaxis).Children(jj).FaceColor);
            end
            try %#ok<TRYNC>
                polaxis(iaxis).Children(jj).MarkerFaceColor = ...
                    rgb2gray(polaxis(iaxis).Children(jj).MarkerFaceColor);
            end
            try %#ok<TRYNC>
                polaxis(iaxis).Children(jj).MarkerEdgeColor = ...
                    rgb2gray(polaxis(iaxis).Children(jj).MarkerEdgeColor);
            end
            try %#ok<TRYNC>
                polaxis(iaxis).Children(jj).CData = ...
                    rgb2gray(polaxis(iaxis).Children(jj).CData);
            end
            try %#ok<TRYNC>
                polaxis(iaxis).Children(jj).NodeColor = ...
                    rgb2gray(polaxis(iaxis).Children(jj).NodeColor);
            end
        end
    end
end

% ========================================================================================
% Stacked axes
staxis=findobj(fig,'type','Stacked');
for iaxis=1:length(staxis)
    % Axis
    for jj = 1:length(staxis(iaxis).AxesProperties)
        warning('off','MATLAB:structOnObject')
        ax = struct(staxis(iaxis).AxesProperties(jj));
        warning('on','MATLAB:structOnObject')
        % Axis lines
        set(ax.Axes,'linewidth',lwidth,'XColor',Black,'YColor',Black,'ZColor',Black,...
            'MinorGridColor',Black,'GridColor',Black);
        % Text
        ht = findobj(ax.Axes,'Type','Text');
        set (ht,'fontname',font,'fontsize',fsize,'color',Black)
        if notype == 0
            set (ht,'interpreter',type)
        end
        % Axis texts (tickslabels)
        set(ax.Axes,'fontname',font,'fontsize',fsize)
        if notype == 0
            set(ax.Axes,'TickLabelInterpreter',type)
        end
    end
    % Legend
    warning('off','MATLAB:structOnObject')
    bx = struct(staxis(iaxis));
    warning('on','MATLAB:structOnObject')
    set(bx.LegendHandle,'Color',White,'EdgeColor',Black,'FontName',font,'FontSize',fsize,...
        'LineWidth',lwidth)
    if notype == 0
        set(bx.LegendHandle,'Interpreter',type)
    end
    % Grayscale
    if grayscale == 1
        staxis(iaxis).Color = rgb2gray(staxis(iaxis).Color);
        try %#ok<TRYNC>
            staxis(iaxis).MarkerFaceColor = ...
                rgb2gray(staxis(iaxis).MarkerFaceColor);
        end
        try %#ok<TRYNC>
            staxis(iaxis).MarkerEdgeColor = ...
                rgb2gray(staxis(iaxis).MarkerEdgeColor);
        end
        for jj = 1:length(staxis.LineProperties)
            try %#ok<TRYNC>
                staxis(iaxis).LineProperties(jj).MarkerFaceColor = ...
                    rgb2gray(staxis(iaxis).LineProperties(jj).MarkerFaceColor);
            end
            try %#ok<TRYNC>
                staxis(iaxis).LineProperties(jj).MarkerEdgeColor = ...
                    rgb2gray(staxis(iaxis).LineProperties(jj).MarkerEdgeColor);
            end
        end
    end
end

% ========================================================================================
% ScatterHistogram
scataxis=findobj(fig,'type','scatterhistogram');
for iaxis=1:length(scataxis)
    % Grayscale
    if grayscale == 1
        scataxis(iaxis).Color = rgb2gray(scataxis(iaxis).Color);
    end
    warning('off','MATLAB:structOnObject')
    ax = struct(scataxis(iaxis));
    warning('on','MATLAB:structOnObject')
    % Axis linewidth
    set(ax.Axes,'linewidth',lwidth);
    % Title
    set(ax.Axes.Title,'fontname',font,'fontsize',fsize,'color',Black);
    if notype == 0
        set(ax.Axes.Title,'interpreter',type);
    end
    % Axis labels
    set([ax.Axes.XLabel ax.Axes.YLabel ax.Axes.ZLabel],...
        'fontname',font,'fontsize',fsize,'color',Black);
    if notype == 0
        set([ax.Axes.XLabel ax.Axes.YLabel ax.Axes.ZLabel],'interpreter',type);
    end
    % Axis texts (tickslabels)
    set(ax.Axes,'fontname',font,'fontsize',fsize)
    if notype == 0
        set(ax.Axes,'TickLabelInterpreter',type)
    end
    % Legend
    set(ax.NodeChildren(1),...
        'fontname',font,'fontsize',fsize);
    if notype == 0
        set(ax.NodeChildren(1),'interpreter',type);
    end
    set(ax.NodeChildren(1),'EdgeColor',Black,'LineWidth',lwidth);
    for jj = 2:length(ax.NodeChildren)
        axnode = ax.NodeChildren(jj);
        % Axis linewidth
        set(axnode,'linewidth',lwidth);
        % Axis color
        set(axnode,'XColor',Black,'YColor',Black,'ZColor',Black)
        % Axis labels
        set([axnode.XLabel axnode.YLabel axnode.ZLabel],...
            'fontname',font,'fontsize',fsize,'color',Black);
        if notype == 0
            set([axnode.XLabel axnode.YLabel axnode.ZLabel],'interpreter',type);
        end
        % Title
        set(axnode.Title,'fontname',font,'fontsize',fsize,'color',Black);
        if notype == 0
            set(axnode.Title,'interpreter',type);
        end
    end
    
end

% ========================================================================================
% Heatmap
heataxis=findobj(fig,'type','Heatmap');
for iaxis=1:length(heataxis)
    % Main text
    set(heataxis(iaxis),...
        'fontname',font,'fontsize',fsize,'FontColor',Black);
    % Axis
    warning('off','MATLAB:structOnObject')
    ax = struct(heataxis(iaxis));
    warning('on','MATLAB:structOnObject')
    % Axis line
    set(ax.Axes,'linewidth',lwidth,'GridColor',Black,'MinorGridColor',Black);
    % Title
    set(ax.Axes.Title,'fontname',font,'fontsize',fsize,'color',Black);
    if notype == 0
        set(ax.Axes.Title,'interpreter',type);
    end
    % Axis labels
    set([ax.Axes.XLabel ax.Axes.YLabel ax.Axes.ZLabel],...
        'fontname',font,'fontsize',fsize,'color',Black);
    if notype == 0
        set([ax.Axes.XLabel ax.Axes.YLabel ax.Axes.ZLabel],...
            'interpreter',type);
    end
    % Axis texts (tickslabels)
    set(ax.Axes,'fontname',font,'fontsize',fsize)
    if notype == 0
        set(ax.Axes,'TickLabelInterpreter',type)
    end
    % Axis color
    set(ax.Axes,'XColor',Black,'YColor',Black,'ZColor',Black)
    % Inner Data
    set(ax.Heatmap,'fontname',font,'fontsize',fsize);
    if notype == 0
        if ~verLessThan('matlab','9.6') % Only for Matlab 2018 and later
            set(ax.Heatmap,'Interpreter',type);
        end
    end
    set(ax.Heatmap,'LineColor',Black,'LineWidth',lwidth);
    % Colorbar
    colbar = findobj(ax.NodeChildren,'Type','Colorbar');
    for ii = 1:length(colbar)
        set(colbar(ii),'fontname',font,'fontsize',fsize)
        set(colbar(ii).Label,'fontname',font,'fontsize',fsize)
        if notype == 0
            set(colbar(ii),'TickLabelInterpreter',type)
            set(colbar(ii).Label,'Interpreter',type)
        end
    end
    % Grayscale
    if grayscale == 1
        set(fig,'Colormap', rgb2gray(get(fig, 'Colormap')));
        set(heataxis(iaxis),'Colormap', rgb2gray(get(heataxis(iaxis), 'Colormap')));
    end
end

% ========================================================================================
% WordCloud
wordaxis=findobj(fig,'type','wordcloud');
for iaxis=1:length(wordaxis)
    % Grayscale
    if grayscale == 1
        wordaxis(iaxis).HighlightColor = rgb2gray(wordaxis(iaxis).HighlightColor);
        wordaxis(iaxis).Color = rgb2gray(wordaxis(iaxis).Color);
    end
    warning('off','MATLAB:structOnObject')
    ax=struct(wordaxis(iaxis));
    warning('on','MATLAB:structOnObject')
    % Title
    set(ax.Axes.Title,'fontname',font,'fontsize',fsize,'color',Black);
    if notype == 0
        set(ax.Axes.Title,'interpreter',type);
    end
    % Words style
    ht=findobj(ax.Axes,'Type','Text');
    set (ht,'fontname',font)
    if notype == 0
        set (ht,'interpreter',type)
    end
end

% ========================================================================================
% Parallelplot
paraxis=findobj(fig,'type','Parallel');
for iaxis=1:length(paraxis)
    % Grayscale
    if grayscale == 1
        paraxis(iaxis).Color = rgb2gray(paraxis(iaxis).Color);
    end
    % Main text
    set(paraxis(iaxis),...
        'fontname',font,'fontsize',fsize);
    warning('off','MATLAB:structOnObject')
    ax = struct(paraxis(iaxis));
    warning('on','MATLAB:structOnObject')
    % Axis linewidth
    set(ax.Axes,'linewidth',lwidth);
    % Title
    set(ax.Axes.Title,'fontname',font,'fontsize',fsize,'color',Black);
    if notype == 0
        set(ax.Axes.Title,'interpreter',type);
    end
    % Axis labels
    set([ax.Axes.XLabel ax.Axes.YLabel ax.Axes.ZLabel],...
        'fontname',font,'fontsize',fsize,'color',Black);
    if notype == 0
        set([ax.Axes.XLabel ax.Axes.YLabel ax.Axes.ZLabel],...
            'interpreter',type);
    end
    % Axis texts (tickslabels)
    set(ax.Axes,'fontname',font,'fontsize',fsize)
    set(ax.YRulers,'fontname',font,'fontsize',fsize)
    if notype == 0
        set(ax.Axes,'TickLabelInterpreter',type)
        set(ax.YRulers,'TickLabelInterpreter',type)
    end
    % Legend
    set(ax.LegendHandle.Legend,'fontname',font,'fontsize',fsize); % Text
    if notype == 0
        set(ax.LegendHandle.Legend,'interpreter',type);
    end
    set(ax.LegendHandle.Legend,'Color',White,'EdgeColor',Black,'LineWidth',lwidth) % Box
    if strcmpi(type,'latex') % Replace for <undefined> legend
        arrayfun(@(x) set(x,'String',strrep(x.String,'<','$<$')),ax.LegendHandle.Legend)
        arrayfun(@(x) set(x,'String',strrep(x.String,'>','$>$')),ax.LegendHandle.Legend)
    end
end

% ========================================================================================
% GeoBubble
geobuaxis=findobj(fig,'type','geobubble');
for iaxis=1:length(geobuaxis)
    warning('off','MATLAB:structOnObject')
    ax = struct(geobuaxis(iaxis));
    warning('on','MATLAB:structOnObject')
    % Axis
    set(ax.Axes,'linewidth',lwidth,'GridColor',Black,'MinorGridColor',Black,...
        'AxisColor',Black);
    % Title
    set(ax.Axes.Title,'fontname',font,'fontsize',fsize,'color',Black);
    if notype == 0
        set(ax.Axes.Title,'interpreter',type);
    end
    % Axis labels
    set([ax.Axes.LongitudeLabel ax.Axes.LatitudeLabel],...
        'fontname',font,'fontsize',fsize,'color',Black);
    if notype == 0
        set([ax.Axes.LongitudeLabel ax.Axes.LatitudeLabel],...
            'interpreter',type);
    end
    if strcmpi(type,'latex') % Replace degree symbol
        ax.Axes.LongitudeAxis.TickLabels = ...
            arrayfun(@(x) strrep(x,char(176),'$^{\circ}$'),ax.Axes.LongitudeAxis.TickLabels);
        ax.Axes.LatitudeAxis.TickLabels = ...
            arrayfun(@(x) strrep(x,char(176),'$^{\circ}$'),ax.Axes.LatitudeAxis.TickLabels);
    end
    set([ax.Axes.LongitudeAxis ax.Axes.LatitudeAxis],...
        'fontname',font,'fontsize',fsize);
    if notype == 0
        set([ax.Axes.LongitudeAxis ax.Axes.LatitudeAxis],...
            'TickLabelInterpreter',type);
    end
    % Scalebar
    set(ax.Axes.Scalebar,'fontname',font,'fontsize',fsize-2)
    % Legend
    leg = findobj(ax.NodeChildren,'Type','legend');
    set(leg,'fontname',font,'fontsize',fsize,'Linewidth',lwidth,'EdgeColor',Black)
    if notype == 0
        set(leg,'interpreter',type)
    end
    if strcmpi(type,'latex') % Replace for <undefined> legend
        arrayfun(@(x) set(x,'String',strrep(x.String,'<','$<$')),leg)
        arrayfun(@(x) set(x,'String',strrep(x.String,'>','$>$')),leg)
    end
    % Size Legend
    sizleg = findobj(ax.NodeChildren,'Type','Axes');
    set(sizleg,'fontname',font,'fontsize',fsize,'Linewidth',lwidth,...
        'GridColor',Black,'MinorGridColor',Black,...
        'XColor',Black,'YColor',Black,'ZColor',Black)
    if notype == 0
        set(sizleg,'TickLabelInterpreter',type)
    end
    set(sizleg.Children,'Linewidth',lwidth)
    ht = findobj(sizleg,'Type','Text');
    set (ht,'fontname',font,'fontsize',fsize)
    if notype == 0
        set (ht,'interpreter',type)
    end
    % Grayscale
    if grayscale == 1
        geosca = findobj(ax.Axes.Children,'Type','Scatter');
        for ii=1:length(geosca)
            geosca(ii).CData = rgb2gray(geosca(ii).CData);
        end
        sizleg.Colormap = rgb2gray(sizleg.Colormap);
        geosca2 = findobj(sizleg.Children,'Type','Scatter');
        for ii=1:length(geosca2)
            geosca2(ii).CData = rgb2gray(geosca2(ii).CData);
            geosca2(ii).MarkerEdgeColor = rgb2gray(geosca2(ii).MarkerEdgeColor);
            geosca2(ii).MarkerFaceColor = rgb2gray(geosca2(ii).MarkerFaceColor);
        end
        % Maps
        if (strcmp(geobuaxis(iaxis).Basemap,'colorterrain') ||...
                strcmp(geobuaxis(iaxis).Basemap,'landcover') )
            geobuaxis(iaxis).Basemap = 'grayterrain';
        end
        if strcmp(geobuaxis(iaxis).Basemap,'bluegreen')
            geobuaxis(iaxis).Basemap = 'grayland';
        end
    end
end

% ========================================================================================
% GeoPlot
geoaxis=findobj(fig,'type','geo');
if ~isempty(geobuaxis); geoaxis=[]; end
for iaxis=1:length(geoaxis)
    % Main axis
    set(geoaxis(iaxis),'LineWidth',lwidth,'GridColor',Black,'MinorGridColor',Black,...
        'AxisColor',Black);
    warning('off','MATLAB:structOnObject')
    ax = struct(geoaxis(iaxis));
    warning('on','MATLAB:structOnObject')
    % Title
    set(ax.Title,'fontname',font,'fontsize',fsize,'color',Black);
    if notype == 0
        set(ax.Title,'interpreter',type);
    end
    % Axis labels
    set([ax.LongitudeLabel ax.LatitudeLabel],...
        'fontname',font,'fontsize',fsize,'color',Black);
    if notype == 0
        set([ax.LongitudeLabel ax.LatitudeLabel],'interpreter',type);
    end
    if strcmpi(type,'latex') % Replace degree symbol
        ax.LongitudeAxis.TickLabels = ...
            arrayfun(@(x) strrep(x,char(176),'$^{\circ}$'),ax.LongitudeAxis.TickLabels);
        ax.LatitudeAxis.TickLabels = ...
            arrayfun(@(x) strrep(x,char(176),'$^{\circ}$'),ax.LatitudeAxis.TickLabels);
    end
    set([ax.LongitudeAxis ax.LatitudeAxis],'fontname',font,'fontsize',fsize);
    if notype == 0
        set([ax.LongitudeAxis ax.LatitudeAxis],'TickLabelInterpreter',type);
    end
    % Scalebar
    set(ax.Scalebar,'fontname',font,'fontsize',fsize-2)
    % Grayscale
    if grayscale == 1
        % Lines
        geolines = findobj(ax.Children,'Type','line');
        for ii=1:length(geolines)
            geolines(ii).Color = rgb2gray(geolines(ii).Color);
        end
        % Scatters
        geoscat = findobj(ax.Children,'Type','scatter');
        for ii=1:length(geoscat)
            if strcmp(geoscat.MarkerEdgeColor,'flat')
                geoscat.MarkerEdgeColor = 'black';
            end
        end
        % Maps
        if (strcmp(geoaxis(iaxis).Basemap,'colorterrain') ||...
                strcmp(geoaxis(iaxis).Basemap,'landcover') )
            geoaxis(iaxis).Basemap = 'grayterrain';
        end
        if strcmp(geoaxis(iaxis).Basemap,'bluegreen')
            geoaxis(iaxis).Basemap = 'grayland';
        end
    end
end

drawnow % Refresh

% ========================================================================================
% Legend General
leg=findobj(fig,'Type','Legend');
set(leg,'fontname',font,'fontsize',fsize); % Text
if notype == 0
    set(leg,'interpreter',type);
end
set(leg,'Color',White,'EdgeColor',Black,'LineWidth',lwidth) % Box

% ========================================================================================
% Colorbar General
colbar = findobj(fig,'Type','Colorbar');
for ii = 1:length(colbar)
    set(colbar(ii),'fontname',font,'fontsize',fsize)
    set(colbar(ii).Label,'fontname',font,'fontsize',fsize)
    set(colbar(ii).Title,'fontname',font,'fontsize',fsize)
    set(colbar(ii),'LineWidth',lwidth,'Color',Black)
    if notype == 0
        set(colbar(ii),'TickLabelInterpreter',type)
        set(colbar(ii).Label,'Interpreter',type)
         set(colbar(ii).Title,'Interpreter',type)
    end
end

% ========================================================================================
% Annotation
annot=findobj(fig,'Type','Annotation');
for ii = 1:length(annot)
    % TextArrow
    tearr = findobj(annot(ii),'Type','TextArrow');
    set(tearr,'fontname',font,'fontsize',fsize,'linewidth',lwidth)
    if notype == 0
        set(tearr,'interpreter',type)
    end
    % TextBox
    tebox = findobj(annot(ii),'Type','TextBox');
    set(tebox,'fontname',font,'fontsize',fsize,'linewidth',lwidth)
    if notype == 0
        set(tebox,'interpreter',type)
    end
    % Arrow
    arr = findobj(annot(ii),'Type','Arrow');
    set(arr,'LineWidth',lwidth)
    % Ellipse
    ell = findobj(annot(ii),'Type','Ellipse');
    set(ell,'LineWidth',lwidth)
    % Rectangle
    recta = findobj(annot(ii),'Type','rect');
    set(recta,'LineWidth',lwidth)
    % Line
    linee = findobj(annot(ii),'Type','lin');
    set(linee,'LineWidth',lwidth)
    % Double Arrow
    dobarr = findobj(annot(ii),'Type','DoubleEndArrow');
    set(dobarr,'LineWidth',lwidth)
    % Grayscale
    if grayscale == 1
        % TextArrow
        for k = 1:length(tearr)
            tearr(k).Color = rgb2gray(tearr(k).Color);
            if isnumeric(tearr(k).TextBackgroundColor)
                tearr(k).TextBackgroundColor = rgb2gray(tearr(k).TextBackgroundColor);
            end
            tearr(k).TextColor = rgb2gray(tearr(k).TextColor);
            if isnumeric(tearr(k).TextEdgeColor)
                tearr(k).TextEdgeColor = rgb2gray(tearr(k).TextEdgeColor);
            end
        end
        % TextBox
        for k = 1:length(tebox)
            if isnumeric(tebox(k).BackgroundColor)
                tebox(k).BackgroundColor = rgb2gray(tebox(k).BackgroundColor);
            end
            tebox(k).Color = rgb2gray(tebox(k).Color);
            tebox(k).EdgeColor = rgb2gray(tebox(k).EdgeColor);
        end
        % Arrow
        for k = 1:length(arr)
            arr(k).Color = rgb2gray(arr(k).Color);
        end
        % Ellipse
        for k = 1:length(ell)
            ell(k).Color = rgb2gray(ell(k).Color);
            if isnumeric(ell(k).FaceColor)
                ell(k).FaceColor = rgb2gray(ell(k).FaceColor);
            end
        end
        % Rectangle
        for k = 1:length(recta)
            recta(k).Color = rgb2gray(recta(k).Color);
            if isnumeric(recta(k).FaceColor)
                recta(k).FaceColor = rgb2gray(recta(k).FaceColor);
            end
        end
        % Line
        for k = 1:length(linee)
            linee(k).Color = rgb2gray(linee(k).Color);
        end
        % Double Arrow
        for k = 1:length(dobarr)
            dobarr(k).Color = rgb2gray(dobarr(k).Color);
        end
    end
end

% ========================================================================================
% Grayscale for old 'text' command
if grayscale == 1
    for ii = 1:length(txt) % 'txt' is defined at the begining
        try %#ok<TRYNC>
            txt(ii).Color = rgb2gray(txt(ii).Color);
        end
        try %#ok<TRYNC>
            txt(ii).EdgeColor = rgb2gray(txt(ii).EdgeColor);
        end
        try %#ok<TRYNC>
            txt(ii).BackgroundColor = rgb2gray(txt(ii).BackgroundColor);
        end
    end
end

% ========================================================================================
% End process
drawnow
set(0,'ShowHiddenHandles','off');
set(fig, 'PaperPosition', [0 0 fwidth(1) fwidth(1)*aspect]);
set(fig, 'PaperSize', [fwidth(1) fwidth(1)*aspect]);
movegui(fig,"center")