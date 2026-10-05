function fig = plot_dense_indented_cshape_boundary(scan,trace,opts)
%PLOT_DENSE_INDENTED_CSHAPE_BOUNDARY Actual refined first-instability brackets.
    arguments
        scan (1,1) struct
        trace (1,1) struct
        opts.XAxis (1,1) string {mustBeMember(opts.XAxis,["iOut","indent"])} = "iOut"
    end
    x=nan(numel(scan.indentation),1); y=x; lower=x; upper=x;
    for ii=1:numel(x)
        roots=trace.roots;
        if isempty(roots), continue; end
        r=roots(abs([roots.indent]-scan.indentation(ii))<1.e-10 & ...
            string({roots.kind})=="stable_to_unstable" & [roots.converged]);
        if isempty(r), continue; end
        [~,order]=sort([r.pressureAmplitude]); r=r(order(1));
        x(ii)=r.(char(opts.XAxis)); lower(ii)=min(r.betaBracket); upper(ii)=max(r.betaBracket);
        y(ii)=(lower(ii)+upper(ii))/2;
    end
    fig=figure('Visible','off','Color','w','Position',[100,100,1000,650]);
    ax=axes(fig,'Position',[.12,.20,.82,.67]);
    errorbar(ax,x,y,y-lower,upper-y,'o-','LineWidth',1.5,'MarkerSize',6, ...
        'Color',[.1,.25,.45],'MarkerFaceColor',[.1,.25,.45],'CapSize',8);
    if opts.XAxis=="iOut", xlabel(ax,'Outboard indentation, d_{out}/(2a_g)');
    else, xlabel(ax,'Construction indentation, I'); end
    ylabel(ax,'First critical \beta_{av} [%]');
    title(ax,'Indented C-shape: first ballooning stability boundary');
    subtitle(ax,sprintf('%d/%d converged transition brackets',sum(isfinite(y)),numel(y)));
    grid(ax,'on'); ax.GridAlpha=.15; ax.FontSize=12; box(ax,'on');
    finiteX=x(isfinite(x));
    if ~isempty(finiteX)
        padding=.04*max(.01,max(finiteX)-min(finiteX));
        xlim(ax,[min(finiteX)-padding,max(finiteX)+padding]);
    end
    annotation(fig,'textbox',[.12,.035,.82,.075],'String', ...
        ['Markers: midpoint of the actual beta bracket; bars: bracket endpoints. ', ...
         'Pressure brackets refined to 0.5%. Lines join the computed brackets. ', ...
         '11 surfaces x 8 phases, 96 cells/turn.'], ...
        'EdgeColor','none','FontSize',10,'Color',[.3,.3,.3]);
end
