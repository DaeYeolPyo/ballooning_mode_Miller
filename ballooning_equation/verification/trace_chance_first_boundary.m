function result=trace_chance_first_boundary(checkpoint,outputFile)
%TRACE_CHANCE_FIRST_BOUNDARY Actual eigenvalue crossing, independent of PRL line.
    setup_indented_cshape_scan;
    saved=load(checkpoint,'report'); source=saved.report;
    assert(numel(source.shapes)==1 && numel(source.equilibria)>=2);
    shape=source.shapes{1};
    config=make_chance_prl_figure3_config(Indentation=shape.iIn);
    config.Ballooning.RadialLabels=.05:.05:.95;
    config.Ballooning.Theta0=2*pi*((0:7)-4)/8;
    config.Ballooning.CellsPerTurn=96;
    config.Ballooning.NSurfacePoints=2048;
    config.Ballooning.FourierModes=48;
    scan=struct('aspectRatio',4,'indentation',shape.iIn,'pressureAmplitude',[], ...
        'points',{cell(2,1)},'equilibria',{cell(2,1)},'shapes',{{shape}}, ...
        'predictors',{{[]}},'config',config);
    input=build_indented_cshape_scan_input(config,shape,source.equilibria{1}.scanPressureAmplitude);
    prepared=prepare_GS_problem(input);
    prepared.qEvaluation.contourCache=build_q_contour_cache(prepared.mesh,config.QEvaluation.contourNR);
    history=struct([]); warm=[];
    for k=1:2
        warm=source.equilibria{k};
        p0=warm.scanPressureAmplitude;
        point=evaluateAtPressure(1,1,p0);
        scan.points{k}=point; scan.equilibria{k}=warm;
        scan.pressureAmplitude(k)=p0;
    end
    trace=trace_indented_cshape_marginal(scan,PressureRelativeTolerance=.01, ...
        MaxIterations=8,Evaluator=@evaluateAtPressure);
    trace.options.Evaluator=[];
    result=struct('scan',scan,'trace',trace,'history',history);
    save(outputFile,'result','-v7.3');
    if ~isempty(trace.roots), disp(trace.roots); end

    function point=evaluateAtPressure(~,~,p0)
        [point,equilibrium]=evaluate_indented_cshape_scan_point(config,shape,p0,prepared,warm);
        if ~isempty(equilibrium), warm=equilibrium; end
        if isempty(history), history=point; else, history(end+1)=point; end
        fprintf('boundary probe p0 %.7g: beta %.7g%% %s complete %d lambda %+.6e\n', ...
            p0,point.betaPercent,point.status,point.complete,point.lambdaMax);
        if ~isempty(point.failure), fprintf('%s\n',point.failure.message); end
        save(outputFile,'scan','history','-v7.3');
    end
end
