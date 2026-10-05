function report = verify_indented_cshape_equilibria(config, equilibria)
%VERIFY_INDENTED_CSHAPE_EQUILIBRIA Check saved GS results and resample ballooning.
%   A cell array in pressure x indentation x aspect-ratio order can be
%   passed from scan.equilibria. Reuses GS solves when only surface/phase
%   resolution changes. Refuses unconverged or mismatched pressure inputs.
    setup_indented_cshape_scan;
    nP=numel(config.PressureAmplitude); nI=numel(config.Indentation);
    nA=numel(config.AspectRatio);
    assert(numel(equilibria)==nP*nI*nA,'Missing equilibria.');
    equilibria=reshape(equilibria,nP,nI,nA);
    points=cell(nP,nI,nA);
    maximumResidual=0; maximumQError=0; maximumQRms=0;
    for ia=1:nA
        for ii=1:nI
            shape=normalize_indented_cshape_boundary(config.AspectRatio(ia), ...
                config.Kappa,config.MinorRadius,config.Delta,config.Indentation(ii), ...
                NBoundary=config.Grid.nBoundary);
            for ip=1:nP
                eq=equilibria{ip,ii,ia};
                assert(~isempty(eq) && eq.solverResult.converged && eq.constraintResult.converged, ...
                    'A requested GS/q equilibrium is missing or unconverged.');
                input=build_indented_cshape_scan_input(config,shape,config.PressureAmplitude(ip));
                assert(norm(eq.profile.dp_dpsiN-input.profile.pprime,inf) ...
                    <=1.e-12*max(1,norm(input.profile.pprime,inf)),'Pressure input mismatch.');
                assert(abs(eq.reference.R0-shape.Rgeo)<1.e-12 && ...
                    abs(eq.reference.aMinor-shape.aMinor)<1.e-12 && ...
                    abs(eq.reference.B0-config.BReference)<1.e-12,'Geometry/reference mismatch.');
                assert(isequal(eq.qGeometry.interiorPsiN(:),config.QEvaluation.interiorPsiN(:)), ...
                    'The requested q-evaluation surfaces differ from the saved solve.');
                u=eq.qResult.psiN(:);
                target=interp1(input.profile.q.psiN,input.profile.q.q,u,'pchip');
                fit=config.PQConstraint.FFprime.fitRange;
                active=u>=fit(1) & u<=fit(2);
                assert(max(abs(eq.qResult.q(active)./target(active)-1)) ...
                    <=config.PQConstraint.outer.qRelativeTolerance,'Prescribed-q target mismatch.');
                assert(eq.solverResult.finalRelativeResidual<config.Picard.residualTolerance);
                assert(eq.constraintResult.finalMaxAbsRelativeQ<=config.PQConstraint.outer.qRelativeTolerance);
                assert(eq.constraintResult.finalRmsLogQ<=config.PQConstraint.outer.qRmsLogTolerance);
                beta=compute_volume_beta(eq,QuadratureDegree=config.BetaQuadratureDegree);
                args=namedargs2cell(config.Ballooning);
                bal=evaluate_equilibrium_ballooning(eq,args{:});
                assert(bal.complete,'A requested ballooning surface/phase solve is incomplete.');
                points{ip,ii,ia}=struct('aspectRatio',shape.aspectRatio,'indent',shape.indent, ...
                    'iOut',shape.iOut,'pressureAmplitude',config.PressureAmplitude(ip), ...
                    'betaPercent',beta.betaPercent,'lambdaMax',bal.lambdaMax, ...
                    'complete',bal.complete,'status',bal.status,'ballooning',bal);
                maximumResidual=max(maximumResidual,eq.solverResult.finalRelativeResidual);
                maximumQError=max(maximumQError,eq.constraintResult.finalMaxAbsRelativeQ);
                maximumQRms=max(maximumQRms,eq.constraintResult.finalRmsLogQ);
                fprintf('Validated I %.3f p0 %.6g: beta %.5g%%, lambda %+.6e, %s\n', ...
                    shape.indent,config.PressureAmplitude(ip),beta.betaPercent,bal.lambdaMax,bal.status);
            end
        end
    end
    report=struct('config',config,'points',{points},'numberEquilibria',numel(points), ...
        'maximumGSResidual',maximumResidual,'maximumQRelativeError',maximumQError, ...
        'maximumQRmsLogError',maximumQRms, ...
        'phaseSolvesPerEquilibrium',numel(config.Ballooning.RadialLabels)*numel(config.Ballooning.Theta0));
end
