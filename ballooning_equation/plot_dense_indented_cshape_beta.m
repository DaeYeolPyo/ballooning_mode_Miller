function fig = plot_dense_indented_cshape_beta(scan,trace,opts)
%PLOT_DENSE_INDENTED_CSHAPE_BETA Actual beta averages and sampled stability.
    arguments
        scan (1,1) struct
        trace (1,1) struct = struct('roots',struct([]))
        opts.XAxis (1,1) string {mustBeMember(opts.XAxis,["iOut","indent"])} = "iOut"
    end
    fig=figure('Visible','off','Color','w','Position',[100,100,1100,700]);
    ax=axes(fig,'Position',[.09,.19,.70,.69]); hold(ax,'on');
    statuses=["stable","unstable","marginal","unknown"];
    colors=[.10,.48,.30;.78,.19,.18;.85,.58,.05;.30,.30,.30];
    names=["Stable (sampled)","Unstable","Marginal band","Unresolved"];
    markers={'o','o','d','x'};
    plotted=0;
    for k=1:numel(statuses)
        x=[]; y=[];
        for j=1:numel(scan.points)
            p=scan.points{j};
            if isempty(p) || string(p.status)~=statuses(k) || ~isfinite(p.betaPercent), continue; end
            x(end+1)=p.(char(opts.XAxis)); y(end+1)=p.betaPercent; %#ok<AGROW>
        end
        plotted=plotted+numel(x);
        if k==4
            scatter(ax,x,y,25,colors(k,:),markers{k},'LineWidth',1.2,'DisplayName',names(k));
        else
            scatter(ax,x,y,22,colors(k,:),markers{k},'filled','DisplayName',names(k));
        end
    end
    if ~isempty(trace.roots)
        roots=trace.roots;
        for kind=unique(string({roots.kind}))
            ys=nan(numel(scan.indentation),1); xs=ys;
            for ii=1:numel(scan.indentation)
                r=roots(abs([roots.indent]-scan.indentation(ii))<1.e-10 & ...
                    string({roots.kind})==kind);
                r=r([r.converged]);
                if isempty(r), continue; end
                [~,order]=sort([r.pressureAmplitude]); r=r(order(1));
                xs(ii)=r.(char(opts.XAxis)); ys(ii)=r.betaPercent;
            end
            plot(ax,xs,ys,'k-','LineWidth',1.6,'DisplayName',strrep(char(kind),'_',' '));
        end
    end
    if opts.XAxis=="iOut", xlabel(ax,'Outboard indentation, d_{out}/(2a_g)');
    else, xlabel(ax,'Construction indentation, I'); end
    ylabel(ax,'\beta_{av} = 2\mu_0\langle p\rangle_V / \langle B^2\rangle_V [%]');
    allX=zeros(numel(scan.indentation),1);
    for ii=1:numel(allX)
        shape=normalize_indented_cshape_boundary(scan.config.AspectRatio, ...
            scan.config.Kappa,scan.config.MinorRadius,scan.config.Delta, ...
            scan.indentation(ii),NBoundary=scan.config.Grid.nBoundary);
        allX(ii)=shape.(char(opts.XAxis));
    end
    margin=.04*max(.01,max(allX)-min(allX));
    xlim(ax,[min(allX)-margin,max(allX)+margin]);
    title(ax,'Indented C-shape: computed ideal infinite-n ballooning map');
    subtitle(ax,sprintf('R_g/a_g = %.4f, \\kappa_g = %.4f, \\delta = %.2f, B_{ref} = %.1f T', ...
        scan.config.AspectRatio,scan.config.Kappa,scan.config.Delta,scan.config.BReference));
    ax.FontSize=12; ax.LineWidth=1; grid(ax,'on'); ax.GridAlpha=.12; box(ax,'on');
    ylim(ax,[0,max(1,1.04*max(scan.betaPercent,[],'all','omitnan'))]);
    legend(ax,'Position',[.81,.62,.18,.23],'FontSize',10);
    annotation(fig,'textbox',[.09,.02,.88,.07],'String',sprintf( ...
        '%d/%d attempted; %d beta values; %d without beta. %d surfaces x %d phases; %d cells/turn.\nMarkers are actual solves. The curve joins refined zero-eigenvalue brackets.', ...
        sum(scan.done,'all'),numel(scan.done),plotted,sum(scan.done,'all')-plotted, ...
        numel(scan.config.Ballooning.RadialLabels), ...
        numel(scan.config.Ballooning.Theta0),scan.config.Ballooning.CellsPerTurn), ...
        'EdgeColor','none','FontSize',10,'Color',[.3,.3,.3]);
end
