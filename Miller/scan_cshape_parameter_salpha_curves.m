function out = scan_cshape_parameter_salpha_curves(eq, psin0, parameterName, varargin)
%SCAN_CSHAPE_PARAMETER_SALPHA_CURVES Scan marginal curves vs C-shape parameter.
%
%   out = scan_cshape_parameter_salpha_curves(eq, psin0, 'B')
%   out = scan_cshape_parameter_salpha_curves(filename, psin0, 'C', ...)
%
% The selected coefficient in
%
%   R = A + B sin(theta) + C cos(2 theta)
%   Z = G cos(theta) - H sin(2 theta)
%
% is varied while the remaining fitted coefficients and radial derivatives
% are held fixed.  For each value, scan_cshape_salpha_diagram is called and
% the marginal curve is lambda_max = 0.

    ip = inputParser;
    addParameter(ip, 'ParameterValues', [], @(x)isnumeric(x)&&isvector(x) || isempty(x));
    addParameter(ip, 'ParameterScale', [0.7, 0.9, 1.1, 1.3], @(x)isnumeric(x)&&isvector(x));
    addParameter(ip, 'SGrid', linspace(0.0, 5.0, 13), @(x)isnumeric(x)&&isvector(x));
    addParameter(ip, 'AlphaGrid', linspace(0.0, 4.0, 13), @(x)isnumeric(x)&&isvector(x));
    addParameter(ip, 'Theta0Grid', linspace(0.0, pi, 5), @(x)isnumeric(x)&&isvector(x));
    addParameter(ip, 'NThetaCoeff', 161, @(x)isnumeric(x)&&isscalar(x)&&x>=32);
    addParameter(ip, 'NGeom', 601, @(x)isnumeric(x)&&isscalar(x)&&x>=128);
    addParameter(ip, 'Dpsin', 1e-4, @(x)isnumeric(x)&&isscalar(x)&&x>0);
    addParameter(ip, 'UseBpFit', true, @(x)islogical(x)&&isscalar(x));
    addParameter(ip, 'ShearDefinition', 'flux', @(x)ischar(x)||isstring(x));
    addParameter(ip, 'ThetaB', 5*pi, @(x)isnumeric(x)&&isscalar(x)&&x>0);
    addParameter(ip, 'NBalloon', 201, @(x)isnumeric(x)&&isscalar(x)&&x>=51);
    addParameter(ip, 'NEigs', 6, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
    addParameter(ip, 'BN', [], @(x) isempty(x) || (isnumeric(x)&&isscalar(x)&&x>0));
    addParameter(ip, 'aN', [], @(x) isempty(x) || (isnumeric(x)&&isscalar(x)&&x>0));
    addParameter(ip, 'Q0', [], @(x) isempty(x) || (isnumeric(x)&&isscalar(x)));
    addParameter(ip, 'AlphaFactor', [], @(x) isempty(x) || (isnumeric(x)&&isscalar(x)&&x~=0));
    addParameter(ip, 'IncludeFPrime', true, @(x)islogical(x)&&isscalar(x));
    addParameter(ip, 'UseParallel', false, @(x)islogical(x)&&isscalar(x));
    addParameter(ip, 'AutoStartPool', true, @(x)islogical(x)&&isscalar(x));
    addParameter(ip, 'PoolType', 'threads', @(x)ischar(x)||isstring(x));
    addParameter(ip, 'NumWorkers', [], @(x) isempty(x) || (isnumeric(x)&&isscalar(x)&&x>=1));
    addParameter(ip, 'Verbose', true, @(x)islogical(x)&&isscalar(x));
    parse(ip, varargin{:});
    opt = ip.Results;

    pname = upper(string(parameterName));
    names = ["A", "B", "C", "G", "H"];
    idx = find(names == pname, 1);
    if isempty(idx)
        error('scan_cshape_parameter_salpha_curves:BadParameterName', ...
            'parameterName must be one of A, B, C, G, or H.');
    end

    if ischar(eq) || isstring(eq)
        eqInput = char(eq);
        eq = read_geqdsk(eqInput);
    else
        eqInput = '';
    end

    base = cshape_local_equilibrium(eq, psin0, ...
        'NTheta', opt.NThetaCoeff, ...
        'NGeom', opt.NGeom, ...
        'Dpsin', opt.Dpsin, ...
        'UseBpFit', opt.UseBpFit, ...
        'BN', opt.BN, ...
        'aN', opt.aN, ...
        'Q0', opt.Q0, ...
        'AlphaFactor', opt.AlphaFactor);

    baseValue = base.coeffs(idx, 1);
    if isempty(opt.ParameterValues)
        if abs(baseValue) > eps
            parameterValues = baseValue .* opt.ParameterScale(:).';
        else
            parameterValues = opt.ParameterScale(:).';
        end
    else
        parameterValues = opt.ParameterValues(:).';
    end
    parameterValues = unique(parameterValues, 'stable');

    nP = numel(parameterValues);
    scans = cell(1, nP);

    useParallel = opt.UseParallel && nP > 1;
    if useParallel
        [useParallel, poolMsg] = ensure_parallel_pool(opt, nP);
        if opt.Verbose && strlength(poolMsg) > 0
            fprintf('%s\n', poolMsg);
        end
    end

    if opt.Verbose
        fprintf('C-shape %s scan: %d values, base %s = %.8g\n', ...
            pname, nP, pname, baseValue);
    end

    if useParallel
        optPar = opt;
        optPar.Verbose = false;
        parfor iparam = 1:nP
            scans{iparam} = run_one_parameter_scan(eq, psin0, base, idx, ...
                parameterValues(iparam), optPar);
        end

        if opt.Verbose
            for iparam = 1:nP
                fprintf('  %s scan %d/%d done (%s = %.8g)\n', ...
                    pname, iparam, nP, pname, parameterValues(iparam));
            end
        end
    else
        for iparam = 1:nP
            if opt.Verbose
                fprintf('\n%s scan %d/%d: %s = %.8g\n', ...
                    pname, iparam, nP, pname, parameterValues(iparam));
            end

            scans{iparam} = run_one_parameter_scan(eq, psin0, base, idx, ...
                parameterValues(iparam), opt);
        end
    end

    out = struct();
    out.model = 'cshape-local-parameter-scan';
    out.parameterName = char(pname);
    out.parameterIndex = idx;
    out.parameterValues = parameterValues;
    out.baseValue = baseValue;
    out.base = base;
    out.scans = scans;
    out.eqInput = eqInput;
    out.psin = psin0;
    out.options = opt;
end

%==========================================================================
function scan = run_one_parameter_scan(eq, psin0, base, idx, value, opt)
    coeffOverride = base.coeffs;
    coeffOverride(idx, 1) = value;

    scan = scan_cshape_salpha_diagram(eq, psin0, ...
        'SGrid', opt.SGrid, ...
        'AlphaGrid', opt.AlphaGrid, ...
        'Theta0Grid', opt.Theta0Grid, ...
        'NThetaCoeff', opt.NThetaCoeff, ...
        'NGeom', opt.NGeom, ...
        'Dpsin', opt.Dpsin, ...
        'UseBpFit', opt.UseBpFit, ...
        'ShearDefinition', opt.ShearDefinition, ...
        'ThetaB', opt.ThetaB, ...
        'NBalloon', opt.NBalloon, ...
        'NEigs', opt.NEigs, ...
        'BN', opt.BN, ...
        'aN', opt.aN, ...
        'Q0', opt.Q0, ...
        'AlphaFactor', opt.AlphaFactor, ...
        'CoeffOverride', coeffOverride, ...
        'IncludeFPrime', opt.IncludeFPrime, ...
        'Verbose', opt.Verbose);
end

%==========================================================================
function [ok, msg] = ensure_parallel_pool(opt, nP)
    ok = false;

    if exist('parpool', 'file') ~= 2 || exist('gcp', 'file') ~= 2
        msg = "Parallel Computing Toolbox functions were not found. Falling back to serial parameter scan.";
        return;
    end

    try
        if ~license('test', 'Distrib_Computing_Toolbox')
            msg = "Parallel Computing Toolbox license is not available. Falling back to serial parameter scan.";
            return;
        end
    catch
        msg = "Could not verify Parallel Computing Toolbox license. Falling back to serial parameter scan.";
        return;
    end

    pool = gcp('nocreate');
    if isempty(pool)
        if ~opt.AutoStartPool
            msg = "No parallel pool is open and AutoStartPool=false. Falling back to serial parameter scan.";
            return;
        end

        poolType = char(opt.PoolType);
        nWorkers = opt.NumWorkers;
        if isempty(nWorkers)
            nWorkers = nP;
        else
            nWorkers = min(round(nWorkers), nP);
        end

        try
            try
                pool = parpool(poolType, nWorkers);
            catch
                pool = parpool(nWorkers);
            end
            msg = sprintf('Started parallel pool with %d workers.', pool.NumWorkers);
        catch ME
            msg = "Could not start a parallel pool (" + string(ME.message) + "). Falling back to serial parameter scan.";
            return;
        end
    else
        msg = sprintf('Using existing parallel pool with %d workers.', pool.NumWorkers);
    end

    ok = true;
end
