root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(root); setup_indented_cshape_scan;
[scan,trace,summary]=collect_dense_indented_cshape_beta_scan;
assert(summary.numberComputed==374 && summary.numberWithBeta==374);
assert(summary.numberComplete==374 && summary.numberUnknown==0);
residuals=[]; qErrors=[]; closure=[];
for ii=1:numel(scan.indentation)
    for ip=1:numel(scan.pressureAmplitude)
        p=scan.points{ip,ii}; eq=scan.equilibria{ip,ii};
        assert(p.indent==scan.indentation(ii));
        assert(p.pressureAmplitude==scan.pressureAmplitude(ip));
        assert(eq.scanPressureAmplitude==p.pressureAmplitude);
        assert(eq.solverResult.converged && eq.constraintResult.converged);
        assert(p.equilibriumConverged && p.complete && isfinite(p.betaPercent));
        assert(p.solverDiagnostics.innerResidual<=scan.config.Picard.residualTolerance);
        assert(p.solverDiagnostics.qConstraintMaximumRelativeError<= ...
            scan.config.PQConstraint.outer.qRelativeTolerance);
        assert(p.solverDiagnostics.qConstraintRmsLogError<= ...
            scan.config.PQConstraint.outer.qRmsLogTolerance);
        bal=p.ballooning;
        assert(isequal(size(bal.lambda),[11,8]) && all(isfinite(bal.lambda),'all'));
        assert(p.lambdaMax==max(bal.lambda,[],'all') && bal.numberFailed==0);
        assert(scan.betaPercent(ip,ii)==p.betaPercent && scan.lambdaMax(ip,ii)==p.lambdaMax);
        assert(scan.status(ip,ii)==p.status);
        for k=1:numel(bal.surfaceDiagnostics)
            d=bal.surfaceDiagnostics{k}; b=bal.options;
            assert(d.geometry.psiNRms<=b.MaximumFluxRms);
            assert(d.geometry.psiNMaximum<=b.MaximumFluxError);
            assert(d.geometry.minimumNormalAlignment>=b.MinimumNormalAlignment);
            assert(abs(d.geometry.qCorrection-1)<=b.MaximumQCorrection);
            assert(d.thetaPsi.inputClosureDefect<=b.MaximumInputClosureDefect);
            closure(end+1)=d.thetaPsi.inputClosureDefect; %#ok<SAGROW>
        end
        residuals(end+1)=p.solverDiagnostics.innerResidual; %#ok<SAGROW>
        qErrors(end+1)=p.solverDiagnostics.qConstraintMaximumRelativeError; %#ok<SAGROW>
    end
end
assert(trace.numberConverged==17 && trace.numberRoots==17);
for r=trace.roots
    assert(r.lambdaBracket(1)<0 && r.lambdaBracket(2)>0);
    assert(diff(r.pressureBracket)<=.005*max(abs(r.pressureBracket)));
end
quality=struct('summary',summary,'maxGSResidual',max(residuals), ...
    'maxConstraintQRelativeError',max(qErrors),'maxInputClosureDefect',max(closure), ...
    'numberConvergedBrackets',trace.numberConverged);
save(fullfile(root,'output','indented_cshape_dense','quality_verification.mat'),'quality');
disp(quality);
files={'make_dense_indented_cshape_scan_config.m','scan_indented_cshape_beta.m', ...
    'examples/run_dense_indented_cshape_beta_scan.m','complete_dense_indented_cshape_column.m', ...
    'refine_dense_indented_cshape_geometry.m','collect_dense_indented_cshape_beta_scan.m', ...
    'plot_dense_indented_cshape_beta.m','plot_dense_indented_cshape_boundary.m', ...
    'solve_indented_cshape_equilibrium.m','verification/test_rhs_profile_reuse.m', ...
    '../GS_fixed/GS/fem/assemble_RHS.m','../GS_fixed/GS/fem/construct_RHS.m'};
for j=1:numel(files)
    issues=checkcode(fullfile(root,files{j}),'-id');
    fprintf('%s: %d analyzer messages\n',files{j},numel(issues));
    assert(isempty(issues));
end
fprintf('FINAL DATA QUALITY CHECK PASSED: all requested points and 17 brackets.\n');
