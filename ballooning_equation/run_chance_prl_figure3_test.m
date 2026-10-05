function report = run_chance_prl_figure3_test(config, outputDirectory, initialEquilibria)
%RUN_CHANCE_PRL_FIGURE3_TEST Actual GS/q and ballooning reproduction experiment.
%   Uses the paper's Eq. (2), not an indented-Miller approximation. Pressure
%   is varied until measured total-field volume beta reaches each target.
%   Only GS/q-converged equilibria enter the stability calculation.
    root=fileparts(mfilename('fullpath')); setup_indented_cshape_scan;
    if nargin<1 || isempty(config), config=make_chance_prl_figure3_config; end
    if nargin<2, outputDirectory=fullfile(root,'output','chance_prl_figure3'); end
    if nargin<3 || isempty(initialEquilibria)
        initialEquilibria=cell(numel(config.Indentation),1);
    end
    if ~iscell(initialEquilibria) || numel(initialEquilibria)~=numel(config.Indentation)
        error('chance:Test:Seeds','Use one seed cell per indentation.');
    end
    if ~isfolder(outputDirectory), mkdir(outputDirectory); end
    file=fullfile(outputDirectory,'reproduction.mat');
    if isfile(file)
        saved=load(file,'report'); report=saved.report;
        if report.version~=2
            error('chance:Test:OldCheckpoint','This checkpoint predates the contour/eigenvalue corrections; use a separate folder.');
        end
        if ~isequaln(report.config,config)
            error('chance:Test:CheckpointMismatch','Choose a separate folder for another configuration.');
        end
    else
        report=struct('config',config,'indentation',config.Indentation, ...
            'targetBetaPercent',config.TargetBetaPercent, ...
            'shapes',{cell(numel(config.Indentation),1)}, ...
            'points',{cell(numel(config.TargetBetaPercent),numel(config.Indentation))}, ...
            'equilibria',{cell(numel(config.TargetBetaPercent),numel(config.Indentation))}, ...
            'predictors',{initialEquilibria(:)}, ...
            'betaHistory',{cell(numel(config.TargetBetaPercent),numel(config.Indentation))}, ...
            'done',false(numel(config.TargetBetaPercent),numel(config.Indentation)), ...
            'source',config.PaperReference,'version',2);
    end
    magPressure=config.BReference^2/(2*(4*pi*1.e-7));
    for ii=1:numel(config.Indentation)
        if all(report.done(:,ii)), continue; end
        shape=chance_prl_bean_boundary(config.Indentation(ii),A=config.BeanA, ...
            B=config.BeanB,AspectRatio=config.AspectRatio,Kappa=config.Kappa, ...
            NBoundary=config.Grid.nBoundary);
        report.shapes{ii}=shape; warm=report.predictors{ii};
        if ~isempty(warm)
            assert(abs(warm.reference.R0-shape.Rgeo)<1.e-8*shape.Rgeo && ...
                abs(warm.reference.aMinor-shape.aMinor)<1.e-8*shape.aMinor, ...
                'chance:Test:SeedGeometry','The seed must use the same bean boundary.');
            validate_paper_profiles(warm,config,warm.scanPressureAmplitude);
        end
        input=build_indented_cshape_scan_input(config,shape,.02*magPressure);
        prepared=prepare_GS_problem(input);
        prepared.qEvaluation.contourCache=build_q_contour_cache(prepared.mesh, ...
            prepared.qEvaluation.contourNR);
        for ip=1:numel(config.TargetBetaPercent)
            if report.done(ip,ii), continue; end
            target=config.TargetBetaPercent(ip); timer=tic;
            fprintf('\nChance Eq.(2): d/(2a) %.5f, target volume beta %.4g%%\n',shape.iIn,target);
            if isempty(warm)
                pressure=3*target/100*magPressure;
            else
                beta=compute_volume_beta(warm,QuadratureDegree=config.BetaQuadratureDegree);
                pressure=warm.scanPressureAmplitude*(target/beta.betaPercent)^.85;
            end
            history=struct([]); equilibrium=[]; beta=[]; lastFailure=[];
            for iteration=1:config.MaximumBetaIterations
                [equilibrium,attempts,predictor]=solve_indented_cshape_equilibrium( ...
                    config,shape,pressure,prepared,warm);
                if ~isempty(predictor), warm=predictor; end
                item=struct('pressureAmplitude',pressure,'betaPercent',NaN, ...
                    'relativeBetaError',NaN,'equilibriumAttempts',attempts);
                if isempty(equilibrium)
                    lastFailure=attempts.lastFailure;
                    history=append(history,item); break
                end
                validate_paper_profiles(equilibrium,config,pressure);
                beta=compute_volume_beta(equilibrium,QuadratureDegree=config.BetaQuadratureDegree);
                item.betaPercent=beta.betaPercent;
                item.relativeBetaError=(beta.betaPercent-target)/target;
                history=append(history,item);
                fprintf('  beta step %d: p0 %.7g, measured %.6g%%, q max %.3e\n', ...
                    iteration,pressure,beta.betaPercent,equilibrium.constraintResult.finalMaxAbsRelativeQ);
                if abs(item.relativeBetaError)<=config.BetaRelativeTolerance, break; end
                % Secant inside a known beta bracket; otherwise use a damped
                % pressure-ratio prediction to avoid aggressive overshooting.
                low=find([history.betaPercent]<target,1,'last');
                high=find([history.betaPercent]>target,1,'last');
                if ~isempty(low) && ~isempty(high)
                    a=history(low); b=history(high);
                    fraction=(target-a.betaPercent)/(b.betaPercent-a.betaPercent);
                    fraction=min(max(fraction,.15),.85);
                    pressure=a.pressureAmplitude+fraction*(b.pressureAmplitude-a.pressureAmplitude);
                else
                    ratio=min(max(target/beta.betaPercent,.5),2);
                    pressure=pressure*ratio^.8;
                end
                report.predictors{ii}=warm;
                report.betaHistory{ip,ii}=history;
                save(file,'report','-v7.3');
            end
            point=struct('aspectRatio',shape.aspectRatio,'indent',shape.iIn,'iIn',shape.iIn, ...
                'pressureAmplitude',NaN,'betaPercent',NaN,'lambdaMax',NaN, ...
                'lambdaMaxObserved',NaN,'complete',false,'status',"unknown", ...
                'targetBetaPercent',target,'targetReached',false, ...
                'solverDiagnostics',[],'ballooning',[],'failure',lastFailure,'wallTime',NaN);
            if ~isempty(equilibrium)
                point.pressureAmplitude=equilibrium.scanPressureAmplitude;
                point.betaPercent=beta.betaPercent;
                point.targetReached=abs(beta.betaPercent-target)<=config.BetaRelativeTolerance*target;
                point.solverDiagnostics=struct('GSResidual',equilibrium.solverResult.finalRelativeResidual, ...
                    'qMaximumRelativeError',equilibrium.constraintResult.finalMaxAbsRelativeQ, ...
                    'qRmsLogError',equilibrium.constraintResult.finalRmsLogQ);
                args=namedargs2cell(config.Ballooning);
                bal=evaluate_equilibrium_ballooning(equilibrium,args{:});
                point.ballooning=bal; point.lambdaMax=bal.lambdaMax;
                point.lambdaMaxObserved=bal.lambdaMaxObserved;
                point.status=bal.status; point.complete=bal.complete;
            end
            point.wallTime=toc(timer);
            report.points{ip,ii}=point; report.equilibria{ip,ii}=equilibrium;
            report.predictors{ii}=warm; report.betaHistory{ip,ii}=history;
            report.done(ip,ii)=true;
            fprintf('  RESULT d/(2a) %.5f: beta %.6g%%, lambda %+.6e, %s, complete %d (%.1fs)\n', ...
                shape.iIn,point.betaPercent,point.lambdaMax,point.status,point.complete,point.wallTime);
            if ~isempty(point.failure), fprintf('  %s\n',point.failure.message); end
            if ~isempty(point.ballooning) && ~point.complete
                first=find(~cellfun(@isempty,point.ballooning.failures),1);
                if ~isempty(first), fprintf('  %s\n',point.ballooning.failures{first}.message); end
            end
            save(file,'report','-v7.3');
        end
    end
    report.finished=all(report.done,'all');
    save(file,'report','-v7.3');
    plot_chance_prl_figure3(report,ExportPath=string(fullfile(outputDirectory,'figure3_comparison.png')));
end

function validate_paper_profiles(eq,config,p0)
    assert(eq.solverResult.converged && eq.constraintResult.converged);
    u=eq.profile.psiN(:); pressure=p0*(1-u.^2).^2;
    assert(norm(eq.profile.pressure-pressure,inf)<5.e-5*p0);
    assert(norm(eq.profile.dp_dpsiN+4*p0*u.*(1-u.^2),inf)<1.e-12*max(1,p0));
    target=polyval(fliplr(config.QCoefficients),eq.qResult.psiN(:));
    active=eq.qResult.psiN<=config.PQConstraint.FFprime.fitRange(2);
    assert(max(abs(eq.qResult.q(active)./target(active)-1)) ...
        <=config.PQConstraint.outer.qRelativeTolerance+2.e-6);
end

function value=append(value,item)
    if isempty(value), value=item; else, value(end+1)=item; end
end
