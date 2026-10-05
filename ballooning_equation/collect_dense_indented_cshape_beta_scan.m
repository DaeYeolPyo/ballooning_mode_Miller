function [scan,trace,summary] = collect_dense_indented_cshape_beta_scan(outputDirectory,opts)
%COLLECT_DENSE_INDENTED_CSHAPE_BETA_SCAN Assemble actual independent solves.
%   Pending/failed points are preserved. No synthetic/interpolated beta or
%   stability values enter the combined checkpoint, table, or scatter plot.
    arguments
        outputDirectory (1,1) string = ""
        opts.Export (1,1) logical = true
    end
    root=fileparts(mfilename('fullpath'));
    if strlength(outputDirectory)==0
        outputDirectory=string(fullfile(root,'output','indented_cshape_dense'));
    end
    config=make_dense_indented_cshape_scan_config;
    nP=numel(config.PressureAmplitude); nI=numel(config.Indentation);
    config.CheckpointFile=char(fullfile(outputDirectory,'beta_scan.mat'));
    signature=rmfield(config,{'CheckpointFile','Resume','Verbose','RetryUnknown'});
    scan=struct('version',3,'config',config,'signature',signature, ...
        'pressureAmplitude',config.PressureAmplitude,'indentation',config.Indentation, ...
        'aspectRatio',config.AspectRatio,'points',{cell(nP,nI)}, ...
        'shapes',{cell(nI,1)},'done',false(nP,nI), ...
        'betaPercent',nan(nP,nI),'lambdaMax',nan(nP,nI), ...
        'status',repmat("pending",nP,nI),'equilibria',{cell(nP,nI)}, ...
        'predictors',{cell(nI,1)},'columnNumericalRecovery',{cell(nI,1)});
    roots=struct([]); evaluations=0; skipped=0; numberTraced=0;
    for ii=1:nI
        file=fullfile(outputDirectory,sprintf('column_%02d.mat',ii));
        if ~isfile(file), continue; end
        data=load(file,'scan'); one=data.scan;
        expected=signature; expected.Indentation=config.Indentation(ii);
        if one.version~=3 || ~isequaln(one.signature,expected)
            error('cshape:Dense:ColumnMismatch','Column %d has a different configuration.',ii);
        end
        for name={'points','done','betaPercent','lambdaMax','status','equilibria'}
            scan.(name{1})(:,ii)=one.(name{1})(:,1);
        end
        scan.shapes{ii}=one.shapes{1}; scan.predictors{ii}=one.predictors{1};
        if isfield(one,'numericalRecoveryConfiguration')
            scan.columnNumericalRecovery{ii}=one.numericalRecoveryConfiguration;
        end
        traceFile=fullfile(outputDirectory,sprintf('trace_%02d.mat',ii));
        if isfile(traceFile)
            saved=load(traceFile,'trace'); t=saved.trace;
            if ~isfield(t,'signature') || ~isequaln(t.signature,expected)
                error('cshape:Dense:TraceMismatch','Trace %d has a different configuration.',ii);
            end
            roots=[roots,t.roots]; %#ok<AGROW>
            evaluations=evaluations+t.numberNewEvaluations;
            skipped=skipped+t.numberSkippedIntervals;
            numberTraced=numberTraced+1;
        end
    end
    scan.finished=all(scan.done,'all');
    scan.numberStable=sum(scan.status=="stable",'all');
    scan.numberUnstable=sum(scan.status=="unstable",'all');
    scan.numberMarginal=sum(scan.status=="marginal",'all');
    scan.numberUnknown=sum(scan.status=="unknown",'all');
    trace=struct('roots',roots,'numberRoots',numel(roots), ...
        'numberConverged',sum(arrayfun(@(r)r.converged,roots)), ...
        'numberSkippedIntervals',skipped,'numberNewEvaluations',evaluations, ...
        'numberColumnsTraced',numberTraced);
    rows=cell(0,22);
    for ii=1:nI
        for ip=1:nP
            p=scan.points{ip,ii}; if isempty(p), continue; end
            residual=NaN; qError=NaN; complete=p.complete; failure="";
            if ~isempty(p.solverDiagnostics)
                residual=p.solverDiagnostics.innerResidual;
                qError=p.solverDiagnostics.qConstraintMaximumRelativeError;
            end
            if ~isempty(p.failure), failure=string(p.failure.message); end
            grid=config.Grid; qe=config.QEvaluation; pq=config.PQConstraint;
            if isfield(p,'solverConfiguration')
                grid=p.solverConfiguration.Grid; qe=p.solverConfiguration.QEvaluation;
                pq=p.solverConfiguration.PQConstraint;
            end
            failedPhases=NaN; worstU=NaN; worstPhase=NaN;
            surfaceSamples=NaN; modes=NaN; projections=NaN; closure=NaN;
            if ~isempty(p.ballooning)
                failedPhases=p.ballooning.numberFailed;
                worstU=p.ballooning.worst.psiN; worstPhase=p.ballooning.worst.theta0;
                b=p.ballooning.options;
                surfaceSamples=b.NSurfacePoints; modes=b.FourierModes;
                projections=b.ProjectionIterations;
                diagnostics=p.ballooning.surfaceDiagnostics;
                measured=cellfun(@(d)d.thetaPsi.inputClosureDefect, ...
                    diagnostics(~cellfun(@isempty,diagnostics)));
                if ~isempty(measured), closure=max(measured); end
                if strlength(failure)==0 && p.ballooning.numberFailed>0
                    first=find(~cellfun(@isempty,p.ballooning.failures),1);
                    failure=string(p.ballooning.failures{first}.message);
                end
            end
            rows(end+1,:)={p.indent,p.iOut,p.pressureAmplitude,p.betaPercent, ...
                string(p.status),complete,p.lambdaMax,p.lambdaMaxObserved, ...
                residual,qError,failedPhases,worstU,worstPhase,failure, ...
                grid.NR,grid.NZ,qe.contourNR,pq.FFprime.nBasis, ...
                surfaceSamples,modes,projections,closure}; %#ok<AGROW>
        end
    end
    samples=cell2table(rows,'VariableNames',{'I','iOut','p0_Pa','beta_av_percent', ...
        'status','complete','lambdaMax','lambdaMaxObserved','GS_relative_residual', ...
        'q_max_relative_error','failed_phase_solves','worst_psiN','worst_theta0','failure', ...
        'GS_NR','GS_NZ','q_contour_NR','FFprime_basis','surface_samples', ...
        'geometry_Fourier_modes','geometry_projection_iterations','max_input_closure_defect'});
    summary=struct('numberRequested',nP*nI,'numberComputed',sum(scan.done,'all'), ...
        'numberWithBeta',sum(isfinite(scan.betaPercent),'all'), ...
        'numberComplete',sum(cellfun(@(p)~isempty(p) && p.complete,scan.points),'all'), ...
        'numberTracedColumns',numberTraced, ...
        'numberStable',scan.numberStable,'numberUnstable',scan.numberUnstable, ...
        'numberMarginal',scan.numberMarginal,'numberUnknown',scan.numberUnknown);
    if opts.Export
        if ~isfolder(outputDirectory), mkdir(outputDirectory); end
        temporary=[tempname(char(outputDirectory)),'.mat'];
        save(temporary,'scan','trace','summary','-v7.3');
        movefile(temporary,config.CheckpointFile,'f');
        writetable(samples,fullfile(outputDirectory,'samples.csv'));
        save(fullfile(outputDirectory,'marginal_trace.mat'),'trace','-v7.3');
        if ~isempty(roots)
            writetable(struct2table(roots),fullfile(outputDirectory,'marginal_brackets.csv'));
        end
        fig=plot_dense_indented_cshape_beta(scan,trace);
        exportgraphics(fig,fullfile(outputDirectory,'beta_av_indentation.png'),'Resolution',200);
        savefig(fig,fullfile(outputDirectory,'beta_av_indentation.fig'));
        close(fig);
        fig=plot_dense_indented_cshape_beta(scan,trace,XAxis="indent");
        exportgraphics(fig,fullfile(outputDirectory,'beta_av_construction_I.png'),'Resolution',200);
        close(fig);
        if ~isempty(roots)
            fig=plot_dense_indented_cshape_boundary(scan,trace);
            exportgraphics(fig,fullfile(outputDirectory,'critical_beta_indentation.png'),'Resolution',200);
            savefig(fig,fullfile(outputDirectory,'critical_beta_indentation.fig')); close(fig);
            fig=plot_dense_indented_cshape_boundary(scan,trace,XAxis="indent");
            exportgraphics(fig,fullfile(outputDirectory,'critical_beta_construction_I.png'),'Resolution',200);
            close(fig);
        end
        write_report(fullfile(outputDirectory,'scan_report.md'),scan,trace,summary,samples);
    end
    fprintf('Dense map: %d/%d attempted, %d beta values, %d complete; stable %d unstable %d marginal %d unknown %d.\n', ...
        summary.numberComputed,summary.numberRequested,summary.numberWithBeta, ...
        summary.numberComplete,summary.numberStable,summary.numberUnstable, ...
        summary.numberMarginal,summary.numberUnknown);
end

function write_report(path,scan,trace,summary,samples)
    fid=fopen(path,'w');
    if fid<0, error('cshape:Dense:Report','Could not open the scan report.'); end
    cleanup=onCleanup(@()fclose(fid));
    fprintf(fid,'# Computed original C-shape dense map\n\n');
    fprintf(fid,'Updated %s.\n\n',char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')));
    fprintf(fid,'Attempted **%d/%d** requested points; %d actual beta values; %d complete stability scans.\n\n', ...
        summary.numberComputed,summary.numberRequested,summary.numberWithBeta,summary.numberComplete);
    fprintf(fid,'Stable %d; unstable %d; marginal %d; unknown %d.\n\n', ...
        summary.numberStable,summary.numberUnstable,summary.numberMarginal,summary.numberUnknown);
    known=scan.betaPercent(isfinite(scan.betaPercent));
    if ~isempty(known)
        fprintf(fid,'Actual volume beta range: %.6f%% to %.6f%%.\n\n',min(known),max(known));
        fprintf(fid,'Maximum accepted GS relative residual: %.6g; maximum accepted q relative error: %.6g.\n\n', ...
            max(samples.GS_relative_residual,[],'omitnan'),max(samples.q_max_relative_error,[],'omitnan'));
    end
    fprintf(fid,'The map uses actual beta = 2 mu0 integral(p dV)/integral(B^2 dV), including both poloidal and toroidal field.\n\n');
    fprintf(fid,'The stability classification samples %d surfaces and %d phases per equilibrium.\n\n', ...
        numel(scan.config.Ballooning.RadialLabels),numel(scan.config.Ballooning.Theta0));
    fprintf(fid,'Numerical resolution varies only where recovery is needed. Each CSV row records the actual GS/q and surface/Fourier controls; the MAT point records full solverConfiguration and ballooning.options. Acceptance tolerances and physical profiles are unchanged.\n\n');
    fprintf(fid,'Refined %d/%d transition brackets in %d columns, using %d extra equilibrium/stability evaluations.\n\n', ...
        trace.numberConverged,trace.numberRoots,trace.numberColumnsTraced,trace.numberNewEvaluations);
    if ~isempty(trace.roots)
        fprintf(fid,'| I | d_out/(2 a_g) | transition | beta_av bracket [%%] | converged |\n');
        fprintf(fid,'|---:|---:|:---|:---|:---|\n');
        for r=trace.roots
            fprintf(fid,'| %.3f | %.6f | %s | %.6f to %.6f | %d |\n', ...
                r.indent,r.iOut,strrep(r.kind,'_',' '),r.betaBracket,r.converged);
        end
        fprintf(fid,'\n');
    end
    fprintf(fid,'These are the current solver results, not an independent validation of the high-beta ballooning formulation. See ../../indented_cshape_dense_README.md for equations, resolution, and resume commands.\n');
end
