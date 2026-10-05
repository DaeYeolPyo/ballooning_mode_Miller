function trace = trace_salpha_marginal(state, shat, opts)
%TRACE_SALPHA_MARGINAL Trace the first lambda_max=0 boundary in alpha.
%
%   trace = TRACE_SALPHA_MARGINAL(state, shat)
%
% Each requested shat row is independently bracketed by AlphaBracket and
% refined by bisection.  The stable side must have lambda_max <= 0 and the
% unstable side lambda_max >= 0.  Rows without such a bracket are returned
% as NaN with status "unbracketed".

    arguments
        state (1,1) struct
        shat double
        opts.AlphaBracket (1,2) double = [0, 3]
        opts.AlphaTolerance (1,1) double = 5.e-3
        opts.MaxIterations (1,1) double = 20
        opts.Theta0 double = []
        opts.ThetaBoundary (1,1) double = 5*pi
        opts.NBalloonPoints (1,1) double = 151
        opts.EnforcePeriodicity (1,1) logical = true
        opts.Verbose (1,1) logical = true
    end

    if ~isvector(shat) || isempty(shat) || ...
            ~isreal(shat) || any(~isfinite(shat))
        error('salpha:Marginal:InvalidShat', ...
            'shat must be a nonempty finite real vector.');
    end
    shat = shat(:).';

    if any(~isfinite(opts.AlphaBracket)) || ...
            opts.AlphaBracket(2) <= opts.AlphaBracket(1)
        error('salpha:Marginal:InvalidBracket', ...
            'AlphaBracket must be a finite increasing pair.');
    end

    validateattributes(opts.AlphaTolerance, {'numeric'}, ...
        {'real', 'finite', 'scalar', 'positive'}, ...
        mfilename, 'AlphaTolerance');
    validateattributes(opts.MaxIterations, {'numeric'}, ...
        {'real', 'finite', 'scalar', 'integer', 'positive'}, ...
        mfilename, 'MaxIterations');

    nShat = numel(shat);
    marginalAlpha = nan(1, nShat);
    lambdaAtMarginal = nan(1, nShat);
    bestTheta0 = nan(1, nShat);
    bracketWidth = nan(1, nShat);
    iterationCount = zeros(1, nShat);
    evaluationCount = zeros(1, nShat);
    status = strings(1, nShat);
    message = strings(1, nShat);

    pointOptions = { ...
        'Theta0', opts.Theta0, ...
        'ThetaBoundary', opts.ThetaBoundary, ...
        'NBalloonPoints', opts.NBalloonPoints, ...
        'EnforcePeriodicity', opts.EnforcePeriodicity};

    for shatIndex = 1:nShat
        alphaLeft = opts.AlphaBracket(1);
        alphaRight = opts.AlphaBracket(2);

        try
            leftPoint = evaluate_salpha_point( ...
                state, shat(shatIndex), alphaLeft, pointOptions{:});
            rightPoint = evaluate_salpha_point( ...
                state, shat(shatIndex), alphaRight, pointOptions{:});
            evaluationCount(shatIndex) = 2;

            if leftPoint.lambdaMax > 0 || rightPoint.lambdaMax < 0
                status(shatIndex) = "unbracketed";
                message(shatIndex) = sprintf( ...
                    'lambda(left)=%+.6e, lambda(right)=%+.6e', ...
                    leftPoint.lambdaMax, rightPoint.lambdaMax);
                continue;
            end

            for iteration = 1:opts.MaxIterations
                alphaMiddle = 0.5*(alphaLeft+alphaRight);
                middlePoint = evaluate_salpha_point( ...
                    state, shat(shatIndex), alphaMiddle, ...
                    pointOptions{:});
                evaluationCount(shatIndex) = ...
                    evaluationCount(shatIndex)+1;

                if middlePoint.lambdaMax <= 0
                    alphaLeft = alphaMiddle;
                    leftPoint = middlePoint;
                else
                    alphaRight = alphaMiddle;
                    rightPoint = middlePoint;
                end

                iterationCount(shatIndex) = iteration;
                if alphaRight-alphaLeft <= opts.AlphaTolerance
                    break;
                end
            end

            % Linear interpolation within the final sign-changing bracket
            % is more accurate than selecting either bisection endpoint.
            denominator = rightPoint.lambdaMax-leftPoint.lambdaMax;
            if denominator == 0
                alphaRoot = 0.5*(alphaLeft+alphaRight);
            else
                alphaRoot = alphaLeft-leftPoint.lambdaMax ...
                    *(alphaRight-alphaLeft)/denominator;
            end

            if abs(leftPoint.lambdaMax) <= abs(rightPoint.lambdaMax)
                representativePoint = leftPoint;
            else
                representativePoint = rightPoint;
            end

            marginalAlpha(shatIndex) = alphaRoot;
            lambdaAtMarginal(shatIndex) = ...
                representativePoint.lambdaMax;
            bestTheta0(shatIndex) = ...
                representativePoint.bestTheta0;
            bracketWidth(shatIndex) = alphaRight-alphaLeft;
            status(shatIndex) = "converged";
        catch exception
            status(shatIndex) = "failed";
            message(shatIndex) = string(exception.message);
        end

        if opts.Verbose
            fprintf(['  marginal row %d/%d: shat=%+.5g, ', ...
                'alpha=%+.6f, status=%s\n'], ...
                shatIndex, nShat, shat(shatIndex), ...
                marginalAlpha(shatIndex), status(shatIndex));
        end
    end

    trace = struct();
    trace.shat = shat;
    trace.marginalAlpha = marginalAlpha;
    trace.lambdaAtMarginal = lambdaAtMarginal;
    trace.bestTheta0 = bestTheta0;
    trace.bracketWidth = bracketWidth;
    trace.iterationCount = iterationCount;
    trace.evaluationCount = evaluationCount;
    trace.status = status;
    trace.message = message;
    trace.alphaBracket = opts.AlphaBracket;
    trace.alphaTolerance = opts.AlphaTolerance;
    trace.numberConverged = nnz(status == "converged");
    trace.numberUnbracketed = nnz(status == "unbracketed");
    trace.numberFailed = nnz(status == "failed");
end
