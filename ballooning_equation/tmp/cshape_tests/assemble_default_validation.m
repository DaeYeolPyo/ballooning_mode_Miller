root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(root); addpath(fullfile(root,'verification')); setup_indented_cshape_scan;
folder=fullfile(root,'tmp','cshape_tests');
config=make_indented_cshape_scan_config;
data06=load(fullfile(folder,'full_shape_06_evaluation.mat'));
data08=load(fullfile(folder,'full_shape_08.mat'),'scan');
data10=load(fullfile(folder,'full_shape_10.mat'),'scan');
valid08=load(fullfile(folder,'partial_default_validation.mat'),'report');
valid10=load(fullfile(folder,'full_shape_10_validation.mat'),'report');
points=cell(3,3); equilibria=cell(3,3); shapes=cell(3,1);
for ii=1:3
    shapes{ii}=normalize_indented_cshape_boundary(config.AspectRatio, ...
        config.Kappa,config.MinorRadius,config.Delta,config.Indentation(ii), ...
        NBoundary=config.Grid.nBoundary);
end
for ip=1:3
    value=data06.results{ip}; bal=value.ballooning; shape=shapes{1};
    points{ip,1}=struct('aspectRatio',shape.aspectRatio,'indent',shape.indent, ...
        'iOut',shape.iOut,'pressureAmplitude',config.PressureAmplitude(ip), ...
        'betaPercent',value.beta.betaPercent,'lambdaMax',bal.lambdaMax, ...
        'complete',bal.complete,'status',bal.status,'ballooning',bal);
    equilibria{ip,1}=value.equilibrium;
    points{ip,3}=valid10.report.points{ip};
    equilibria{ip,3}=data10.scan.equilibria{ip};
end
for ip=1:2
    points{ip,2}=valid08.report.points{ip,1};
    equilibria{ip,2}=data08.scan.equilibria{ip};
end
candidate=fullfile(folder,'high_shape_prediction_08.mat');
if isfile(candidate)
    high=load(candidate,'equilibrium');
    if ~isempty(high.equilibrium)
        one=config; one.Indentation=.8; one.PressureAmplitude=config.PressureAmplitude(3);
        verified=verify_indented_cshape_equilibria(one,{high.equilibrium});
        points{3,2}=verified.points{1}; equilibria{3,2}=high.equilibrium;
    end
end
if isempty(points{3,2})
    shape=shapes{2};
    points{3,2}=struct('aspectRatio',shape.aspectRatio,'indent',shape.indent, ...
        'iOut',shape.iOut,'pressureAmplitude',config.PressureAmplitude(3), ...
        'betaPercent',NaN,'lambdaMax',NaN,'complete',false,'status',"unknown");
end
maximumResidual=0; maximumQError=0; maximumQRms=0; count=0;
for ii=1:3
    for ip=1:3
        eq=equilibria{ip,ii}; if isempty(eq), continue; end
        input=build_indented_cshape_scan_input(config,shapes{ii},config.PressureAmplitude(ip));
        assert(eq.solverResult.converged && eq.constraintResult.converged);
        assert(norm(eq.profile.dp_dpsiN-input.profile.pprime,inf) ...
            <1.e-12*max(1,norm(input.profile.pprime,inf)));
        assert(abs(eq.reference.R0-shapes{ii}.Rgeo)<1.e-12 && ...
            abs(eq.reference.aMinor-shapes{ii}.aMinor)<1.e-12);
        assert(eq.solverResult.finalRelativeResidual<config.Picard.residualTolerance);
        assert(eq.constraintResult.finalMaxAbsRelativeQ<=config.PQConstraint.outer.qRelativeTolerance);
        assert(eq.constraintResult.finalRmsLogQ<=config.PQConstraint.outer.qRmsLogTolerance);
        u=eq.qResult.psiN; target=interp1(input.profile.q.psiN,input.profile.q.q,u,'pchip');
        active=u<=config.PQConstraint.FFprime.fitRange(2);
        assert(max(abs(eq.qResult.q(active)./target(active)-1)) ...
            <=config.PQConstraint.outer.qRelativeTolerance);
        bal=points{ip,ii}.ballooning;
        assert(bal.complete && bal.options.NSurfacePoints==1024 && bal.options.FourierModes==48);
        assert(isequal(bal.psiN(:),config.Ballooning.RadialLabels(:)) && ...
            isequal(bal.theta0(:),config.Ballooning.Theta0(:)));
        maximumResidual=max(maximumResidual,eq.solverResult.finalRelativeResidual);
        maximumQError=max(maximumQError,eq.constraintResult.finalMaxAbsRelativeQ);
        maximumQRms=max(maximumQRms,eq.constraintResult.finalRmsLogQ);
        count=count+1;
    end
end
report=struct('config',config,'points',{points},'equilibria',{equilibria}, ...
    'numberRequested',9,'numberValidated',count,'maximumGSResidual',maximumResidual, ...
    'maximumQRelativeError',maximumQError,'maximumQRmsLogError',maximumQRms, ...
    'phaseSolvesPerEquilibrium',88);
save(fullfile(folder,'default_validation.mat'),'report','-v7.3');
scan=struct('aspectRatio',config.AspectRatio,'indentation',config.Indentation, ...
    'pressureAmplitude',config.PressureAmplitude,'points',{points},'shapes',{shapes});
[fig,~]=plot_indented_cshape_beta(scan,Visible='off', ...
    ExportPath=string(fullfile(folder,'default_seed_validation.png')));
close(fig);
fprintf('Validated %d/9 native GS/q + 88-phase points. Max GS %.6e, q max %.6e, q RMS %.6e.\n', ...
    count,maximumResidual,maximumQError,maximumQRms);
