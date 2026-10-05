function fig = plot_chance_prl_figure3(report, opts)
%PLOT_CHANCE_PRL_FIGURE3 Overlay computed classifications on approximate reference.
%   Dashed lines are visually digitized from the supplied PRL scan, not fits
%   to the numerical output. Reference-reading uncertainty is about 0.5 beta
%   percentage points and 0.003 in indentation; no exact original data exist.
    arguments
        report (1,1) struct
        opts.ExportPath (1,1) string = ""
        opts.Visible (1,1) string = "off"
    end
    fig=figure('Color','w','Visible',char(opts.Visible)); ax=axes(fig); hold(ax,'on');
    plot(ax,[.095 .11 .15 .185 .22 .25 .28 .30 .31], ...
        [1.6 1.6 1.6 1.75 2 2.4 3 3.7 4.6],'--','Color',[.55 .55 .55], ...
        'LineWidth',1.3,'DisplayName','PRL boundary (approximate reading)');
    plot(ax,[.10 .11 .14 .18 .20 .24 .26 .29 .303 .31], ...
        [24 22 17 13.4 12.4 11.1 10.7 9 7.2 5.8],'--', ...
        'Color',[.55 .55 .55],'LineWidth',1.3,'HandleVisibility','off');
    missing=0;
    for j=1:numel(report.points)
        p=report.points{j}; if isempty(p), continue; end
        if ~isfinite(p.betaPercent), missing=missing+1; continue; end
        switch string(p.status)
            case "stable", marker='o'; face='none'; color=[.05 .35 .75];
            case "unstable", marker='o'; face=[.8 .15 .1]; color=face;
            case "marginal", marker='d'; face='none'; color=[.5 .2 .5];
            otherwise, marker='x'; face='none'; color=[.3 .3 .3];
        end
        plot(ax,p.indent,p.betaPercent,marker,'Color',color,'MarkerFaceColor',face, ...
            'MarkerSize',7,'LineWidth',1.2,'HandleVisibility','off');
    end
    plot(ax,NaN,NaN,'o','Color',[.05 .35 .75],'MarkerSize',7,'DisplayName','Computed stable (sampled)');
    plot(ax,NaN,NaN,'o','Color',[.8 .15 .1],'MarkerFaceColor',[.8 .15 .1], ...
        'MarkerSize',7,'DisplayName','Computed unstable');
    plot(ax,NaN,NaN,'x','Color',[.3 .3 .3],'MarkerSize',7,'DisplayName','Unresolved');
    xlabel(ax,'Inboard indentation d/(2a)'); ylabel(ax,'Total-field volume beta [%]');
    xlim(ax,[.08 .35]); ylim(ax,[0 26]); grid(ax,'on');
    title(ax,sprintf('Chance PRL Fig. 3 test; R/a=4, b/a=1.386; %d without beta',missing));
    legend(ax,'Location','northeast','FontSize',9);
    if strlength(opts.ExportPath)>0, exportgraphics(fig,opts.ExportPath,'Resolution',180); end
end
