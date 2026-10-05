function report = test_indented_cshape_pipeline(opts)
%TEST_INDENTED_CSHAPE_PIPELINE Analytic geometry, beta and crossing checks.
%   RunIntegration=true additionally runs a real prescribed-q C-shape
%   equilibrium, phase/surface sweep, prepared reuse and checkpoint resume.
    arguments
        opts.RunIntegration (1,1) logical = false
    end
    root=fileparts(fileparts(mfilename('fullpath')));
    addpath(root); setup_indented_cshape_scan;
    for A=[3,4,7,10]
        for I=[0,.6,.8,1.15]
            shape=normalize_indented_cshape_boundary(A,3,.8,-.4,I,NBoundary=4096);
            b=shape.boundary;
            measuredA=(max(b.R)+min(b.R))/(max(b.R)-min(b.R));
            measuredK=(max(b.Z)-min(b.Z))/(max(b.R)-min(b.R));
            assert(abs(measuredA-A)<2.e-5 && abs(measuredK-3)<1.e-5);
            assert(numel(b.R)==4096 && shape.polygonArea>0);
        end
    end
    try
        normalize_indented_cshape_boundary(4,3,.8,-.4,2);
        error('test:ExpectedFailure','Self-contact was not rejected.');
    catch e
        assert(strcmp(e.identifier,'cshape:Boundary:Topology'));
    end
    % Independent analytic integral on R in [2,3], Z in [0,3-R].
    % psi=R-2, F=2R gives Bp^2=1/R^2 and Bphi^2=4.
    eq=manufactured_equilibrium;
    % Cached interpolation must reproduce a general cubic and constants,
    % including exact P3-node values and NaN semantics outside the mesh.
    Xq=[eq.mesh.points;2.17,.21;2.36,.1;1.9,.2;2.9,.9];
    cache=build_P3_value_cache(eq.mesh,Xq);
    cubic=@(R,Z)R.^3+2*R.^2.*Z-3*R.*Z.^2+.7*Z.^3+2*R-Z+1;
    nodal=cubic(eq.mesh.points(:,1),eq.mesh.points(:,2));
    cached=cache.operator*nodal;
    expected=cubic(Xq(:,1),Xq(:,2));
    assert(max(abs(cached(cache.inside)-expected(cache.inside)))<1.e-12);
    assert(norm(cache.operator(1:10,:)-speye(10),inf)<1.e-12);
    assert(max(abs(cache.operator(cache.inside,:)*ones(10,1)-1))<1.e-12);
    assert(all(~cache.inside(end-1:end)));
    beta=compute_volume_beta(eq,QuadratureDegree=14);
    exactVolume=7*pi/3;
    exactBp=2*pi*(3*log(3/2)-1);
    exactBeta=2*(4*pi*1.e-7)*(3.e4*exactVolume)/(exactBp+4*exactVolume);
    assert(abs(beta.volume-exactVolume)<1.e-11);
    assert(abs(beta.Bp2VolumeIntegral-exactBp)/exactBp<1.e-8);
    assert(abs(beta.beta-exactBeta)/exactBeta<1.e-9);
    eq.psi=-eq.psi; eq.global.sibry=-1;
    reversed=compute_volume_beta(eq,QuadratureDegree=14);
    assert(abs(reversed.beta-beta.beta)<1.e-13);
    polynomial=@(p)(p-1)*(p-2)*(p-3);
    seed=synthetic_scan([0,.7,1.5,2.5,3.5,4],polynomial);
    trace=trace_indented_cshape_marginal(seed, ...
        Evaluator=@(ia,ii,p)synthetic_point(p,polynomial), ...
        PressureAbsoluteTolerance=1.e-9,PressureRelativeTolerance=1.e-8, ...
        MaxIterations=40,Verbose=false);
    assert(trace.numberRoots==3 && trace.numberConverged==3);
    assert(max(abs([trace.roots.pressureAmplitude]-[1,2,3]))<1.e-7);
    monotonic=@(p)p-1;
    one=synthetic_scan([0,.7,1.5,4],monotonic);
    first=trace_indented_cshape_marginal(one,Evaluator=@(ia,ii,p)synthetic_point(p,monotonic),Verbose=false);
    assert(first.numberRoots==1 && strcmp(first.roots.kind,'stable_to_unstable'));
    stable=synthetic_scan([0,1,2],@(p)-1);
    empty=trace_indented_cshape_marginal(stable,Evaluator=@(ia,ii,p)synthetic_point(p,@(v)-1),Verbose=false);
    assert(empty.numberRoots==0);
    gap=seed; gap.points{3}.complete=false; gap.points{3}.lambdaMax=NaN;
    broken=trace_indented_cshape_marginal(gap, ...
        Evaluator=@(ia,ii,p)synthetic_point(p,polynomial),Verbose=false);
    assert(broken.numberRoots==1 && broken.numberSkippedIntervals==2);
    sampled=synthetic_scan([0,1,2,3,4],polynomial);
    zeros=trace_indented_cshape_marginal(sampled, ...
        Evaluator=@(ia,ii,p)synthetic_point(p,polynomial),Verbose=false);
    assert(zeros.numberRoots==3 && zeros.numberNewEvaluations==0);
    % A failed refinement midpoint must keep an unresolved bracket.
    blocked=trace_indented_cshape_marginal(one, ...
        Evaluator=@(ia,ii,p)unresolved_point(p), ...
        PressureAbsoluteTolerance=0,Verbose=false);
    assert(blocked.numberRoots==1 && blocked.numberConverged==0 && ...
        strcmp(blocked.roots.terminationReason,'unresolved_midpoint'));
    % A bracket can first meet tolerance on the final allowed iteration.
    narrow=synthetic_scan([0,2],@(p)p-.9);
    last=trace_indented_cshape_marginal(narrow, ...
        Evaluator=@(ia,ii,p)synthetic_point(p,@(v)v-.9), ...
        MaxIterations=1,PressureAbsoluteTolerance=1,Verbose=false);
    assert(last.numberConverged==1 && last.roots.iterations==1 && ...
        strcmp(last.roots.terminationReason,'pressure_bracket'));
    [figureHandle,~]=plot_indented_cshape_beta(seed,trace,Visible="off");
    close(figureHandle);
    report=struct('geometryCases',16,'betaRelativeError',abs(beta.beta-exactBeta)/exactBeta, ...
        'threeCrossings',trace.numberRoots,'singleCrossing',first.numberRoots, ...
        'integration',[]);
    fprintf('Analytic geometry/beta/crossing checks passed.\n');
    if opts.RunIntegration
        config=make_indented_cshape_scan_config(Indentation=.8,CentralBeta=.002, ...
            GridSize=[24,32],NBoundary=80);
        config.Ballooning.RadialLabels=[.35,.7];
        config.Ballooning.Theta0=[-pi,0];
        config.Ballooning.NSurfacePoints=256;
        config.Ballooning.FourierModes=24;
        config.Ballooning.AdapterNR=161;
        config.Ballooning.CellsPerTurn=16;
        folder=fullfile(root,'tmp','cshape_tests');
        if ~isfolder(folder), mkdir(folder); end
        config.CheckpointFile=[tempname(folder),'.mat'];
        cleanup=onCleanup(@()delete(config.CheckpointFile));
        scan=scan_indented_cshape_beta(config);
        p=scan.points{1};
        assert(p.equilibriumConverged && p.complete, ...
            'Real C-shape point failed: %s',describe_failure(p));
        assert(isfinite(p.betaPercent) && all(isfinite(p.ballooning.lambda),'all'));
        config.Resume=true;
        % Stability cached before the open-contour correction is obsolete.
        current=scan; scan.version=2;
        save(config.CheckpointFile,'scan','-v7.3');
        rejected=false;
        try
            scan_indented_cshape_beta(config);
        catch e
            rejected=strcmp(e.identifier,'cshape:Scan:ResumeMismatch');
        end
        assert(rejected,'An obsolete stability checkpoint must be rejected.');
        scan=current; save(config.CheckpointFile,'scan','-v7.3');
        resumed=scan_indented_cshape_beta(config);
        assert(isequaln(resumed.points,scan.points));
        assert(scan.version==3 && ~isempty(scan.predictors{1}) && ...
            scan.predictors{1}.solverResult.converged && ...
            scan.predictors{1}.constraintResult.converged);
        assert(all(cellfun(@isempty,scan.equilibria),'all'));
        % RetryUnknown must recompute an unfinished classification while
        % retaining the converged predictor even with StoreEquilibria=false.
        scan.status(1)="unknown"; scan.points{1}.status="unknown";
        scan.points{1}.complete=false;
        save(config.CheckpointFile,'scan','-v7.3');
        config.RetryUnknown=true;
        retried=scan_indented_cshape_beta(config);
        assert(retried.points{1}.equilibriumConverged && retried.points{1}.complete && ...
            strcmp(retried.points{1}.status,'stable'));
        config.Resume=false; config.CheckpointFile='';
        shape=scan.shapes{1};
        input=build_indented_cshape_scan_input(config,shape,config.PressureAmplitude);
        prepared=prepare_GS_problem(input);
        input.dim.R0=input.dim.R0+.1;
        try
            solve_GS(input,prepared);
            error('test:ExpectedFailure','Mismatched prepared geometry was accepted.');
        catch e
            assert(strcmp(e.identifier,'GS:Solve:PreparedMismatch'));
        end
        report.integration=p;
        fprintf('Real C-shape equilibrium, 2 surfaces x 2 phases, checkpoint resume and unknown retry passed.\n');
    end
end

function eq=manufactured_equilibrium
    L=[1,0,0;0,1,0;0,0,1;2/3,1/3,0;1/3,2/3,0; ...
        0,2/3,1/3;0,1/3,2/3;1/3,0,2/3;2/3,0,1/3;1/3,1/3,1/3];
    V=[2,0;3,0;2,1];
    P=L*V;
    eq=struct('mesh',struct('P1points',V,'P1elements',[1,2,3], ...
        'points',P,'elements',1:10),'psi',P(:,1)-2, ...
        'global',struct('simag',0,'sibry',1), ...
        'profile',struct('psiN',[0;.5;1],'pressure',3.e4*ones(3,1),'Fpol',[4;5;6]));
end

function scan=synthetic_scan(pressure,fun)
    shape=struct('aspectRatio',4,'indent',.8,'iOut',.1644);
    points=cell(numel(pressure),1);
    for j=1:numel(pressure), points{j}=synthetic_point(pressure(j),fun); end
    scan=struct('aspectRatio',4,'indentation',.8,'pressureAmplitude',pressure, ...
        'shapes',{{shape}},'points',{points});
end

function p=synthetic_point(pressure,fun)
    value=fun(pressure);
    status="stable";
    if value>0, status="unstable"; elseif value==0, status="marginal"; end
    p=struct('pressureAmplitude',pressure,'lambdaMax',value, ...
        'betaPercent',10*pressure/(1+pressure),'complete',true,'status',status, ...
        'aspectRatio',4,'indent',.8,'iOut',.1644);
end

function p=unresolved_point(pressure)
    p=synthetic_point(pressure,@(v)NaN);
    p.complete=false; p.status="unknown";
end

function message=describe_failure(p)
    message='Inspect the saved geometry/phase diagnostics.';
    if ~isempty(p.failure), message=p.failure.message; end
    if ~isempty(p.ballooning)
        items=p.ballooning.failures;
        first=find(~cellfun(@isempty,items),1);
        if ~isempty(first), message=items{first}.message; end
    end
end
