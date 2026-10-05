function [fig, info] = plot_indented_cshape_beta(scan, trace, opts)
%PLOT_INDENTED_CSHAPE_BETA Sampled stable/unstable points and found boundaries.
%   Unresolved roots are shown as brackets; failures with no beta cannot
%   be placed on this plane and are counted explicitly in the title.
    arguments
        scan (1,1) struct
        trace (1,1) struct = struct('roots',struct([]))
        opts.XAxis (1,1) string {mustBeMember(opts.XAxis,["iOut","indent"])} = "iOut"
        opts.Visible (1,1) string {mustBeMember(opts.Visible,["on","off"])} = "on"
        opts.ExportPath (1,1) string = ""
    end
    fig=figure('Color','w','Visible',char(opts.Visible));
    ax=axes(fig); hold(ax,'on');
    colors=lines(numel(scan.aspectRatio));
    hidden=0;
    for ia=1:numel(scan.aspectRatio)
        color=colors(ia,:);
        for ii=1:numel(scan.indentation)
            for ip=1:numel(scan.pressureAmplitude)
                p=scan.points{ip,ii,ia};
                if isempty(p), continue; end
                if ~isfinite(p.betaPercent), hidden=hidden+1; continue; end
                x=p.(char(opts.XAxis));
                switch string(p.status)
                    case "stable", marker='o'; face='none';
                    case "unstable", marker='o'; face=color;
                    case "marginal", marker='d'; face='none';
                    otherwise, marker='x'; face='none';
                end
                plot(ax,x,p.betaPercent,marker,'Color',color, ...
                    'MarkerFaceColor',face,'MarkerSize',5,'HandleVisibility','off');
            end
        end
        % Match branches by transition direction and pressure rank. NaNs
        % at missing indentations break lines rather than bridging gaps.
        roots=trace.roots;
        if ~isempty(roots)
            roots=roots(abs([roots.aspectRatio]-scan.aspectRatio(ia))<1.e-10);
            kinds=unique(string({roots.kind}));
            for kind=kinds(:).'
                rows=cell(numel(scan.indentation),1); maximum=0;
                for ii=1:numel(rows)
                    selected=roots(abs([roots.indent]-scan.indentation(ii))<1.e-10 & ...
                        string({roots.kind})==kind);
                    selected=selected([selected.converged]);
                    [~,order]=sort([selected.pressureAmplitude]);
                    rows{ii}=selected(order); maximum=max(maximum,numel(selected));
                end
                for branch=1:maximum
                    xs=nan(numel(rows),1); ys=xs;
                    for ii=1:numel(rows)
                        if numel(rows{ii})>=branch
                            r=rows{ii}(branch);
                            xs(ii)=r.(char(opts.XAxis)); ys(ii)=r.betaPercent;
                        end
                    end
                    style='-';
                    if abs(scan.aspectRatio(ia)-4)<1.e-10, style='--'; end
                    plot(ax,xs,ys,style,'Color',color,'LineWidth',1.4,'HandleVisibility','off');
                    plot(ax,xs,ys,'.','Color',color,'MarkerSize',12,'HandleVisibility','off');
                end
            end
            for r=roots(~[roots.converged])
                x=r.(char(opts.XAxis));
                plot(ax,[x,x],r.betaBracket,':','Color',color,'LineWidth',2,'HandleVisibility','off');
            end
        end
        plot(ax,NaN,NaN,'-','Color',color,'LineWidth',1.4, ...
            'DisplayName',sprintf('R_g/a_g = %.4g',scan.aspectRatio(ia)));
    end
    plot(ax,NaN,NaN,'ko','MarkerFaceColor','none','DisplayName','Stable (sampled)');
    plot(ax,NaN,NaN,'ko','MarkerFaceColor','k','DisplayName','Unstable');
    plot(ax,NaN,NaN,'kd','DisplayName','Marginal');
    plot(ax,NaN,NaN,'kx','DisplayName','Unresolved');
    if opts.XAxis=="iOut", xlabel(ax,'Outboard indentation d_{out}/(2a_g)');
    else, xlabel(ax,'Construction indentation I'); end
    ylabel(ax,'Volume beta [%]'); grid(ax,'on');
    title(ax,sprintf('Ideal infinite-n ballooning; %d unresolved points without beta',hidden));
    legend(ax,'Location','best');
    info=struct('numberPointsWithoutBeta',hidden, ...
        'boundaryMatching','transition direction and pressure rank; missing rows break curves');
    if strlength(opts.ExportPath)>0
        folder=fileparts(opts.ExportPath);
        if strlength(folder)>0 && ~isfolder(folder), mkdir(folder); end
        exportgraphics(fig,opts.ExportPath,'Resolution',180);
    end
end
