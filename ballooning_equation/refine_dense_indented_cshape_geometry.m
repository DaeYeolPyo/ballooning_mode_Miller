function scan = refine_dense_indented_cshape_geometry(columnFile,opts)
%REFINE_DENSE_INDENTED_CSHAPE_GEOMETRY Re-evaluate stored native equilibria.
%   The physical equilibrium, beta, and acceptance tolerances are retained.
%   Actual geometry options are stored with every updated point. Run only
%   after the process writing this column has finished.
    arguments
        columnFile (1,1) string
        opts.RefineBoundary (1,1) logical = true
        opts.AdapterNR (1,1) double = 513
        opts.NSurfacePoints (1,1) double = 4096
        opts.FourierModes (1,1) double = 64
        opts.ProjectionIterations (1,1) double = 3
    end
    setup_indented_cshape_scan;
    saved=load(columnFile,'scan'); scan=saved.scan;
    controls=scan.config.Ballooning;
    for name={'AdapterNR','NSurfacePoints','FourierModes','ProjectionIterations'}
        controls.(name{1})=opts.(name{1});
    end
    args=namedargs2cell(controls);
    for ip=1:numel(scan.pressureAmplitude)
        eq=scan.equilibria{ip}; p=scan.points{ip};
        if isempty(eq) || isempty(p), continue; end
        if p.complete && isequaln(p.ballooning.options,controls), continue; end
        timer=tic; bal=evaluate_equilibrium_ballooning(eq,args{:});
        p.ballooning=bal; p.lambdaMax=bal.lambdaMax;
        p.lambdaMaxObserved=bal.lambdaMaxObserved;
        p.status=bal.status; p.complete=bal.complete;
        if ~isempty(p.failure) && strcmp(p.failure.stage,'ballooning') && bal.complete
            p.failure=[];
        end
        p.wallTime=p.wallTime+toc(timer);
        p.geometryRefinement=struct('reason','higher resolution at unchanged quality limits', ...
            'options',controls);
        scan.points{ip}=p; scan.status(ip)=p.status; scan.lambdaMax(ip)=p.lambdaMax;
        fprintf('GEOMETRY I %.3f p0 %.6g beta %.6f: complete %d %s lambda %.8g\n', ...
            p.indent,p.pressureAmplitude,p.betaPercent,p.complete,p.status,p.lambdaMax);
        write_checkpoint(scan,columnFile);
    end
    scan.finished=all(scan.done,'all');
    scan.numberStable=sum(scan.status=="stable",'all');
    scan.numberUnstable=sum(scan.status=="unstable",'all');
    scan.numberMarginal=sum(scan.status=="marginal",'all');
    scan.numberUnknown=sum(scan.status=="unknown",'all');
    scan.geometryRefinement=controls;
    write_checkpoint(scan,columnFile);
    if opts.RefineBoundary
        refined=scan; refined.config.Ballooning=controls;
        trace=trace_indented_cshape_marginal(refined,PressureRelativeTolerance=.005,MaxIterations=12);
        trace.signature=scan.signature; trace.geometryRefinement=controls;
        [directory,name]=fileparts(columnFile);
        traceFile=fullfile(directory,strrep(name,'column_','trace_')+".mat");
        temporary=[tempname(char(directory)),'.mat'];
        save(temporary,'trace','-v7.3'); movefile(temporary,traceFile,'f');
    end
end

function write_checkpoint(scan,path)
    temporary=[tempname(char(fileparts(path))),'.mat'];
    save(temporary,'scan','-v7.3'); movefile(temporary,path,'f');
end
