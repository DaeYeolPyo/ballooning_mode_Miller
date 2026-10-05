function scan = complete_dense_indented_cshape_column(columnFile,opts)
%COMPLETE_DENSE_INDENTED_CSHAPE_COLUMN Recover actual missing pressure points.
%   Neighbor equilibria are initial guesses only. Every accepted point is
%   solved at its own boundary/pressure with the original GS/q tolerances.
%   RefineMesh increases numerical resolution without changing the physical
%   problem. Per-point solverConfiguration records the controls actually used.
%   Do not run concurrently with another writer of the same column.
    arguments
        columnFile (1,1) string
        opts.NeighborFile (1,1) string = ""
        opts.ProbeFile (1,1) string = ""
        opts.RefineMesh (1,1) logical = false
        opts.RefineBoundary (1,1) logical = true
    end
    setup_indented_cshape_scan;
    saved=load(columnFile,'scan'); scan=saved.scan;
    config=scan.config;
    config.Ballooning.NSurfacePoints=4096;
    config.Ballooning.FourierModes=64; config.Ballooning.ProjectionIterations=3;
    if opts.RefineMesh
        config.Grid.NR=64; config.Grid.NZ=96;
        config.QEvaluation.contourNR=400;
    end
    config.Continuation.Enabled=false;
    config.Continuation.MaximumQIterationRestarts=1;
    config.PQConstraint.outer.maxIterations=24;
    controls=struct('Grid',config.Grid,'Picard',config.Picard, ...
        'PQConstraint',config.PQConstraint,'QEvaluation',config.QEvaluation, ...
        'Ballooning',config.Ballooning,'Continuation',config.Continuation);
    neighbor=[];
    if strlength(opts.NeighborFile)>0
        saved=load(opts.NeighborFile,'scan'); neighbor=saved.scan;
        if ~isequal(neighbor.pressureAmplitude,scan.pressureAmplitude)
            error('cshape:Dense:NeighborPressures','Neighbor pressure controls must match.');
        end
    end
    if strlength(opts.ProbeFile)>0
        probe=load(opts.ProbeFile);
        if probe.point.equilibriumConverged
            ip=probe.ip;
            if probe.point.pressureAmplitude~=scan.pressureAmplitude(ip) || ...
                    probe.point.indent~=scan.indentation
                error('cshape:Dense:ProbeMismatch','Recovery probe has a different physical point.');
            end
            probe.point.solverConfiguration=struct('Grid',probe.config.Grid, ...
                'Picard',probe.config.Picard,'PQConstraint',probe.config.PQConstraint, ...
                'QEvaluation',probe.config.QEvaluation,'Ballooning',probe.config.Ballooning, ...
                'Continuation',probe.config.Continuation);
            scan=store_point(scan,ip,probe.point,probe.equilibrium);
            scan.predictors{1}=probe.equilibrium;
            write_checkpoint(scan,columnFile);
        end
    end
    prepared=[]; warm=scan.predictors{1};
    for ip=1:numel(scan.pressureAmplitude)
        if scan.done(ip) && scan.points{ip}.equilibriumConverged, continue; end
        p0=scan.pressureAmplitude(ip);
        if isempty(prepared)
            input=build_indented_cshape_scan_input(config,scan.shapes{1},p0);
            prepared=prepare_GS_problem(input);
            prepared.qEvaluation.contourCache=build_q_contour_cache( ...
                prepared.mesh,prepared.qEvaluation.contourNR);
        end
        seeds={warm,[]}; labels=["previous accepted pressure","neighbor at same pressure"];
        nSeeds=1;
        if ~isempty(neighbor) && ~isempty(neighbor.equilibria{ip})
            seeds{2}=neighbor.equilibria{ip}; nSeeds=2;
        end
        history=cell(nSeeds,1);
        for j=1:nSeeds
            [point,equilibrium,predictor]=evaluate_indented_cshape_scan_point( ...
                config,scan.shapes{1},p0,prepared,seeds{j});
            history{j}=struct('seed',labels(j),'equilibriumAttempts',point.equilibriumAttempts);
            if point.equilibriumConverged, break; end
        end
        point.solverConfiguration=controls;
        point.recoverySeedAttempts=history(1:j);
        scan=store_point(scan,ip,point,equilibrium);
        if point.equilibriumConverged
            warm=predictor; scan.predictors{1}=warm;
        end
        fprintf('RECOVERY I %.3f p0 %.6g beta %.6f: complete %d %s lambda %.8g\n', ...
            point.indent,p0,point.betaPercent,point.complete,point.status,point.lambdaMax);
        if ~isempty(point.failure), fprintf('  %s\n',point.failure.message); end
        write_checkpoint(scan,columnFile);
    end
    scan.numericalRecoveryConfiguration=controls;
    scan.finished=all(scan.done,'all');
    scan.numberStable=sum(scan.status=="stable",'all');
    scan.numberUnstable=sum(scan.status=="unstable",'all');
    scan.numberMarginal=sum(scan.status=="marginal",'all');
    scan.numberUnknown=sum(scan.status=="unknown",'all');
    write_checkpoint(scan,columnFile);
    scan=refine_dense_indented_cshape_geometry(columnFile,RefineBoundary=false);
    if opts.RefineBoundary
        refined=scan; refined.config=config;
        trace=trace_indented_cshape_marginal(refined,PressureRelativeTolerance=.005,MaxIterations=12);
        trace.signature=scan.signature; trace.numericalRecoveryConfiguration=controls;
        [directory,name]=fileparts(columnFile);
        traceFile=fullfile(directory,strrep(name,'column_','trace_')+".mat");
        temporary=[tempname(char(directory)),'.mat'];
        save(temporary,'trace','-v7.3'); movefile(temporary,traceFile,'f');
    end
    fprintf('RECOVERY FINISHED I %.3f: %d/%d equilibria, %d unknown\n', ...
        scan.indentation,sum(isfinite(scan.betaPercent)),numel(scan.done),scan.numberUnknown);
end

function scan=store_point(scan,ip,p,eq)
    scan.points{ip}=p; scan.equilibria{ip}=eq; scan.done(ip)=true;
    scan.betaPercent(ip)=p.betaPercent; scan.lambdaMax(ip)=p.lambdaMax; scan.status(ip)=p.status;
end

function write_checkpoint(scan,path)
    temporary=[tempname(char(fileparts(path))),'.mat'];
    save(temporary,'scan','-v7.3'); movefile(temporary,path,'f');
end
