function scan = scan_local_salpha_diagram(localEq,targetPsiN,opts)
%SCAN_LOCAL_SALPHA_DIAGRAM Scan the frozen-geometry local stability plane.
%
%   Positive lambdaMax is unstable. lambdaMax is maximized over theta0.

narginchk(2,3);
if nargin<3 || isempty(opts)
    opts = struct();
end
opts = validate_options(opts);
model = prepare_local_salpha_model(localEq,targetPsiN,struct( ...
    'nPeriods',opts.nPeriods, ...
    'poloidalModeCutoff',opts.poloidalModeCutoff));

sValues = opts.sHatValues(:);
alphaValues = opts.alphaValues(:);
theta0Values = linspace(-pi,pi,opts.nTheta0+1).';
theta0Values(end) = [];
nS = numel(sValues);
nA = numel(alphaValues);
nT = numel(theta0Values);

lambdaTheta0 = nan(nS,nA,nT);
residualTheta0 = nan(nS,nA,nT);
failMessage = strings(nS,nA,nT);

scanStart = tic;
for phaseIndex = 1:nT
    basis = assemble_local_salpha_phase_basis( ...
        model,theta0Values(phaseIndex),opts.elementsPerPeriod);
    fields = basis.coefficients;
    matrices = basis.matrices;
    free = (2:numel(basis.thetaDof)-1).';
    initialVector = ones(numel(free),1)/sqrt(numel(free));

    for sIndex = 1:nS
        sHat = sValues(sIndex);
        g = fields.g0+sHat*fields.gS+sHat^2*fields.gSS;
        f = fields.f0+sHat*fields.fS+sHat^2*fields.fSS;
        G = matrices.G0+sHat*matrices.GS+sHat^2*matrices.GSS;
        F = matrices.F0+sHat*matrices.FS+sHat^2*matrices.FSS;

        if any(g<=0) || any(f<=0)
            failMessage(sIndex,:,phaseIndex) = ...
                "Non-positive reconstructed g or f.";
            continue
        end

        for alphaIndex = 1:nA
            alpha = alphaValues(alphaIndex);
            c = alpha*fields.cA+alpha^2*fields.cAA ...
                +sHat*alpha*fields.cSA;
            C = alpha*matrices.CA+alpha^2*matrices.CAA ...
                +sHat*alpha*matrices.CSA;
            try
                [lambda,residual,vector] = leading_eigenvalue( ...
                    G,C,F,c,f,free,initialVector,opts);
                lambdaTheta0(sIndex,alphaIndex,phaseIndex) = lambda;
                residualTheta0(sIndex,alphaIndex,phaseIndex) = residual;
                initialVector = vector;
            catch ME
                failMessage(sIndex,alphaIndex,phaseIndex) = ...
                    string(ME.message);
                initialVector = ones(numel(free),1)/sqrt(numel(free));
            end
        end
    end
    if opts.verbose
        fprintf(['s-alpha phase %d/%d complete: theta0/pi=%+.3f, ' ...
                 'elapsed %.1f s\n'],phaseIndex,nT, ...
            theta0Values(phaseIndex)/pi,toc(scanStart))
    end
end

[lambdaMax,bestIndex] = max(lambdaTheta0,[],3,'omitnan');
bestTheta0 = reshape(theta0Values(bestIndex),size(lambdaMax));
allFailed = all(isnan(lambdaTheta0),3);
lambdaMax(allFailed) = NaN;
bestTheta0(allFailed) = NaN;

scan = struct();
scan.model = model;
scan.sHat = sValues;
scan.alpha = alphaValues;
scan.theta0 = theta0Values;
scan.lambdaMax = lambdaMax;
scan.lambdaTheta0 = lambdaTheta0;
scan.residualTheta0 = residualTheta0;
scan.bestTheta0 = bestTheta0;
scan.unstable = lambdaMax>0;
scan.failMessage = failMessage;
scan.options = opts;
scan.elapsedSeconds = toc(scanStart);
scan.diagnostics = struct( ...
    'nFailed',nnz(strlength(failMessage)>0), ...
    'maximumResidual',max(residualTheta0,[],'all','omitnan'), ...
    'minimumLambda',min(lambdaMax,[],'all','omitnan'), ...
    'maximumLambda',max(lambdaMax,[],'all','omitnan'), ...
    'unstableFraction',mean(scan.unstable,'all'));
end


function [lambda,residual,vector] = leading_eigenvalue( ...
        G,C,F,c,f,free,initialVector,opts)
A = 0.5*((C-G)+(C-G).');
F = 0.5*(F+F.');
Afree = A(free,free);
Ffree = F(free,free);
[~,flag] = chol(Ffree);
if flag~=0
    error('GS:ballooning:NonPositiveScanMass', ...
        'The scan mass matrix is not positive definite.');
end
upperBound = max(c./f);
shift = upperBound+max(1,abs(upperBound));
eigsOptions = struct('tol',opts.eigenTolerance, ...
    'maxit',opts.maxIterations,'disp',0,'v0',initialVector);
[vector,D,eigsFlag] = eigs(Afree,Ffree,1,shift,eigsOptions);
if eigsFlag~=0 || ~isfinite(D)
    error('GS:ballooning:ScanEigsFailure', ...
        'Sparse leading-eigenvalue solve did not converge.');
end
lambda = real(D);
vector = real(vector);
vector = vector/sqrt(vector.'*Ffree*vector);
Av = Afree*vector;
Fv = Ffree*vector;
residual = norm(Av-lambda*Fv) ...
    /max(norm(Av)+abs(lambda)*norm(Fv),eps);
end


function opts = validate_options(opts)
if ~isstruct(opts) || ~isscalar(opts)
    error('GS:ballooning:InvalidSAlphaScanOptions', ...
        'opts must be a scalar struct.');
end
defaults = struct( ...
    'sHatValues',linspace(0,4,17), ...
    'alphaValues',linspace(0,4,17), ...
    'nTheta0',8,'nPeriods',6,'elementsPerPeriod',64, ...
    'poloidalModeCutoff',24,'eigenTolerance',1e-9, ...
    'maxIterations',1500,'verbose',true);
unknown = setdiff(fieldnames(opts),fieldnames(defaults));
if ~isempty(unknown)
    error('GS:ballooning:UnknownSAlphaScanOption', ...
        'Unknown option: %s',strjoin(unknown,', '));
end
names = fieldnames(defaults);
for k = 1:numel(names)
    if ~isfield(opts,names{k})
        opts.(names{k}) = defaults.(names{k});
    end
end
validateattributes(opts.sHatValues,{'numeric'}, ...
    {'real','finite','vector','nonnegative'},mfilename,'sHatValues');
validateattributes(opts.alphaValues,{'numeric'}, ...
    {'real','finite','vector','nonnegative'},mfilename,'alphaValues');
if any(diff(opts.sHatValues)<=0) || any(diff(opts.alphaValues)<=0)
    error('GS:ballooning:NonMonotoneSAlphaGrid', ...
        'sHatValues and alphaValues must be strictly increasing.');
end
integerFields = {'nTheta0','nPeriods','elementsPerPeriod', ...
    'poloidalModeCutoff','maxIterations'};
for k = 1:numel(integerFields)
    name = integerFields{k};
    validateattributes(opts.(name),{'numeric'}, ...
        {'real','finite','scalar','integer','positive'},mfilename,name);
end
if opts.nTheta0<2
    error('GS:ballooning:TooFewTheta0', ...
        'At least two theta0 phases are required.');
end
if opts.elementsPerPeriod<2*opts.poloidalModeCutoff+2
    error('GS:ballooning:AngularGridTooCoarse', ...
        'elementsPerPeriod does not resolve the retained poloidal modes.');
end
validateattributes(opts.eigenTolerance,{'numeric'}, ...
    {'real','finite','scalar','positive','<',1}, ...
    mfilename,'eigenTolerance');
validateattributes(opts.verbose,{'logical'}, ...
    {'scalar'},mfilename,'verbose');
end
