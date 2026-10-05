function trace = trace_indented_cshape_marginal(scan, opts)
%TRACE_INDENTED_CSHAPE_MARGINAL Refine every observed eigenvalue crossing.
%   Does not require re-entry into stability. Adjacent incomplete/failed
%   points are not bridged. Evaluator(ia,ii,p0) is an optional independent
%   callback returning a scan-point structure, useful for analytic checks.
    arguments
        scan (1,1) struct
        opts.PressureRelativeTolerance (1,1) double {mustBeFinite,mustBePositive} = 1.e-3
        opts.PressureAbsoluteTolerance (1,1) double {mustBeFinite,mustBeNonnegative} = 10
        opts.MaxIterations (1,1) double {mustBeInteger,mustBePositive} = 24
        opts.Evaluator = []
        opts.Verbose (1,1) logical = true
    end
    setup_indented_cshape_scan;
    if ~isempty(opts.Evaluator) && ~isa(opts.Evaluator,'function_handle')
        error('cshape:Trace:Evaluator','Evaluator must be a function handle.');
    end
    roots=struct([]); skipped=0; evaluations=0;
    for ia=1:numel(scan.aspectRatio)
        for ii=1:numel(scan.indentation)
            shape=scan.shapes{ii,ia};
            prepared=[]; warm=[];
            if isfield(scan,'predictors'), warm=scan.predictors{ii,ia}; end
            if isempty(shape)
                skipped=skipped+max(0,numel(scan.pressureAmplitude)-1);
                continue
            end
            points=scan.points(:,ii,ia);
            % Exactly sampled zeroes need no duplicate adjacent brackets.
            for ip=1:numel(points)
                p=points{ip};
                if usable(p) && p.lambdaMax==0
                    r=root_record(shape,p,p,p,'sampled_zero',true,'sampled_zero',0);
                    roots=append(roots,r);
                end
            end
            for ip=1:numel(points)-1
                lo=points{ip}; hi=points{ip+1};
                if ~usable(lo) || ~usable(hi)
                    skipped=skipped+1; continue
                end
                if lo.lambdaMax==0 || hi.lambdaMax==0 || sign(lo.lambdaMax)==sign(hi.lambdaMax)
                    continue
                end
                if lo.lambdaMax<0, kind='stable_to_unstable'; else, kind='unstable_to_stable'; end
                if isfield(scan,'equilibria')
                    candidates=scan.equilibria(ip:ip+1,ii,ia);
                    candidates=candidates(~cellfun(@isempty,candidates));
                    if ~isempty(candidates)
                        warm=candidates{1};
                    end
                end
                best=lo;
                if abs(hi.lambdaMax)<abs(lo.lambdaMax), best=hi; end
                reason='maximum_iterations'; converged=false; count=0;
                for k=1:opts.MaxIterations
                    width=hi.pressureAmplitude-lo.pressureAmplitude;
                    tolerance=max(opts.PressureAbsoluteTolerance, ...
                        opts.PressureRelativeTolerance*max(abs([lo.pressureAmplitude,hi.pressureAmplitude])));
                    if width<=tolerance
                        converged=true; reason='pressure_bracket'; break
                    end
                    p0=(lo.pressureAmplitude+hi.pressureAmplitude)/2;
                    if isempty(opts.Evaluator)
                        if isempty(prepared)
                            input=build_indented_cshape_scan_input(scan.config,shape,p0);
                            prepared=prepare_GS_problem(input);
                            prepared.qEvaluation.contourCache=build_q_contour_cache( ...
                                prepared.mesh,prepared.qEvaluation.contourNR);
                        end
                        [mid,eq]=evaluate_indented_cshape_scan_point(scan.config,shape,p0,prepared,warm);
                        if ~isempty(eq), warm=eq; end
                    else
                        mid=opts.Evaluator(ia,ii,p0);
                    end
                    evaluations=evaluations+1; count=count+1;
                    if ~usable(mid)
                        reason='unresolved_midpoint'; break
                    end
                    if abs(mid.lambdaMax)<abs(best.lambdaMax), best=mid; end
                    if mid.lambdaMax==0
                        lo=mid; hi=mid; best=mid;
                        converged=true; reason='sampled_zero'; break
                    elseif sign(mid.lambdaMax)==sign(lo.lambdaMax)
                        lo=mid;
                    else
                        hi=mid;
                    end
                end
                if strcmp(reason,'maximum_iterations')
                    tolerance=max(opts.PressureAbsoluteTolerance, ...
                        opts.PressureRelativeTolerance*max(abs([lo.pressureAmplitude,hi.pressureAmplitude])));
                    if hi.pressureAmplitude-lo.pressureAmplitude<=tolerance
                        converged=true; reason='pressure_bracket';
                    end
                end
                % The best sampled point must lie inside the final bracket.
                if best.pressureAmplitude<lo.pressureAmplitude || best.pressureAmplitude>hi.pressureAmplitude
                    best=lo;
                    if abs(hi.lambdaMax)<abs(lo.lambdaMax), best=hi; end
                end
                r=root_record(shape,lo,hi,best,kind,converged,reason,count);
                roots=append(roots,r);
                if opts.Verbose
                    fprintf('A %.4f I %.4f %s: p0 [%.6g,%.6g], %s\n', ...
                        shape.aspectRatio,shape.indent,kind,lo.pressureAmplitude,hi.pressureAmplitude,reason);
                end
            end
        end
    end
    trace=struct('roots',roots,'numberRoots',numel(roots), ...
        'numberConverged',sum(arrayfun(@(r)r.converged,roots)), ...
        'numberSkippedIntervals',skipped,'numberNewEvaluations',evaluations,'options',opts);
end

function yes=usable(p)
    yes=isstruct(p) && isfield(p,'complete') && p.complete && ...
        isfield(p,'lambdaMax') && isfinite(p.lambdaMax) && ...
        isfinite(p.pressureAmplitude) && isfinite(p.betaPercent);
end

function r=root_record(shape,lo,hi,best,kind,converged,reason,count)
    r=struct('aspectRatio',shape.aspectRatio,'indent',shape.indent,'iOut',shape.iOut, ...
        'kind',kind,'pressureAmplitude',best.pressureAmplitude, ...
        'betaPercent',best.betaPercent,'lambda',best.lambdaMax, ...
        'pressureBracket',[lo.pressureAmplitude,hi.pressureAmplitude], ...
        'betaBracket',[lo.betaPercent,hi.betaPercent], ...
        'lambdaBracket',[lo.lambdaMax,hi.lambdaMax], ...
        'converged',converged,'terminationReason',reason,'iterations',count);
end

function out=append(out,r)
    if isempty(out), out=r; else, out(end+1)=r; end
end
