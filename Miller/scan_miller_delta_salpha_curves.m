function out = scan_miller_delta_salpha_curves(p, varargin)
%SCAN_MILLER_DELTA_SALPHA_CURVES Scan marginal s-alpha curves vs Miller shape.
%
%   out = scan_miller_delta_salpha_curves(p)
%   out = scan_miller_delta_salpha_curves(p,'Name',value,...)
%
% By default this is a Miller Fig. 5 style scan: delta is varied, while the
% remaining Miller shape/local-equilibrium parameters are held fixed. Use
% 'ScanParameter','kappa' and 'KappaGrid',... to scan elongation instead.
% For each shape value, scan_miller_salpha_diagram is called and the
% marginal curve is lambda=0.

    ip = inputParser;
    addParameter(ip, 'ScanParameter', 'delta', @(x)ischar(x)||isstring(x));
    addParameter(ip, 'DeltaGrid', default_delta_grid(p), @(x)isnumeric(x)&&isvector(x));
    addParameter(ip, 'KappaGrid', default_kappa_grid(p), @(x)isnumeric(x)&&isvector(x));
    addParameter(ip, 'DeltaSign', 1.0, @(x)isnumeric(x)&&isscalar(x)&&any(x == [-1, 1]));
    addParameter(ip, 'SGrid', linspace(0.0, 5.0, 13), @(x)isnumeric(x)&&isvector(x));
    addParameter(ip, 'AlphaGrid', linspace(0.0, 4.0, 13), @(x)isnumeric(x)&&isvector(x));
    addParameter(ip, 'Theta0Grid', linspace(0.0, pi, 5), @(x)isnumeric(x)&&isvector(x));
    addParameter(ip, 'NThetaCoeff', 161, @(x)isnumeric(x)&&isscalar(x)&&x>=32);
    addParameter(ip, 'NGeom', 601, @(x)isnumeric(x)&&isscalar(x)&&x>=128);
    addParameter(ip, 'FFprimeScale', 1.0, @(x)isnumeric(x)&&isscalar(x)&&isfinite(x));
    addParameter(ip, 'AlphaMapping', 'volume', @(x)ischar(x)||isstring(x));
    addParameter(ip, 'ThetaB', 5*pi, @(x)isnumeric(x)&&isscalar(x)&&x>0);
    addParameter(ip, 'NBalloon', 201, @(x)isnumeric(x)&&isscalar(x)&&x>=51);
    addParameter(ip, 'NEigs', 6, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
    addParameter(ip, 'UseParallel', false, @(x)islogical(x)&&isscalar(x));
    addParameter(ip, 'AutoStartPool', true, @(x)islogical(x)&&isscalar(x));
    addParameter(ip, 'PoolType', 'threads', @(x)ischar(x)||isstring(x));
    addParameter(ip, 'NumWorkers', [], @(x) isempty(x) || (isnumeric(x)&&isscalar(x)&&x>=1));
    addParameter(ip, 'Verbose', true, @(x)islogical(x)&&isscalar(x));
    parse(ip, varargin{:});
    opt = ip.Results;

    scanParameter = lower(string(opt.ScanParameter));
    switch scanParameter
        case "delta"
            parameterName = 'delta';
            parameterLabel = '\delta';
            parameterGrid = unique(opt.DeltaGrid(:).', 'stable');
            internalParameterGrid = opt.DeltaSign.*parameterGrid;
        case "kappa"
            parameterName = 'kappa';
            parameterLabel = '\kappa';
            parameterGrid = unique(opt.KappaGrid(:).', 'stable');
            internalParameterGrid = parameterGrid;
        otherwise
            error('scan_miller_delta_salpha_curves:BadScanParameter', ...
                'ScanParameter must be ''delta'' or ''kappa''.');
    end
    scanName = [upper(parameterName(1)), parameterName(2:end)];

    nP = numel(parameterGrid);
    scans = cell(1, nP);

    useParallel = opt.UseParallel && nP > 1;
    if useParallel
        [useParallel, poolMsg] = ensure_parallel_pool(opt, nP);
        if opt.Verbose && strlength(poolMsg) > 0
            fprintf('%s\n', poolMsg);
        end
    end

    if useParallel
        if opt.Verbose
            fprintf('Running %s scans in parallel across %d values.\n', ...
                parameterName, nP);
        end

        verboseInside = opt.Verbose;
        optPar = opt;
        optPar.Verbose = false;
        parfor id = 1:nP
            scans{id} = run_one_parameter_scan(p, parameterName, ...
                internalParameterGrid(id), optPar);
        end

        if verboseInside
            for id = 1:nP
                fprintf('  %s scan %d/%d done (%s = %.4g, internal = %.4g)\n', ...
                    parameterName, id, nP, parameterName, ...
                    parameterGrid(id), internalParameterGrid(id));
            end
        end
    else
        for id = 1:nP
            if opt.Verbose
                fprintf('\n%s scan %d/%d: %s = %.4g, internal = %.4g\n', ...
                    scanName, id, nP, parameterName, parameterGrid(id), ...
                    internalParameterGrid(id));
            end
            scans{id} = run_one_parameter_scan(p, parameterName, ...
                internalParameterGrid(id), opt);
        end
    end

    out = struct();
    out.scan_parameter = parameterName;
    out.parameter_label = parameterLabel;
    out.parameter_values = parameterGrid;
    out.internal_parameter_values = internalParameterGrid;
    if scanParameter == "delta"
        out.delta = parameterGrid;
        out.internal_delta = internalParameterGrid;
    else
        out.kappa = parameterGrid;
        out.internal_kappa = internalParameterGrid;
    end
    out.scans = scans;
    out.params = p;
    out.options = opt;
end

function deltaGrid = default_delta_grid(p)
    if isfield(p, 'delta')
        d0 = p.delta;
    else
        d0 = 0.4;
    end
    deltaGrid = [0.0, 0.2, d0, 0.6];
    deltaGrid = deltaGrid(abs(deltaGrid) < 0.98);
end

function kappaGrid = default_kappa_grid(p)
    if isfield(p, 'kappa')
        k0 = p.kappa;
    else
        k0 = 1.6;
    end
    kappaGrid = [1.2, 1.4, k0, 1.8];
end

function scan = run_one_parameter_scan(p, parameterName, parameterValue, opt)
    pScan = p;
    pScan.(parameterName) = parameterValue;

    scan = scan_miller_salpha_diagram(pScan, ...
        'SGrid', opt.SGrid, ...
        'AlphaGrid', opt.AlphaGrid, ...
        'Theta0Grid', opt.Theta0Grid, ...
        'NThetaCoeff', opt.NThetaCoeff, ...
        'NGeom', opt.NGeom, ...
        'FFprimeScale', opt.FFprimeScale, ...
        'AlphaMapping', opt.AlphaMapping, ...
        'ThetaB', opt.ThetaB, ...
        'NBalloon', opt.NBalloon, ...
        'NEigs', opt.NEigs, ...
        'Verbose', opt.Verbose);
end

function [ok, msg] = ensure_parallel_pool(opt, nD)
    ok = false;
    msg = "";

    if exist('parpool', 'file') ~= 2 || exist('gcp', 'file') ~= 2
        msg = "Parallel Computing Toolbox functions were not found. Falling back to serial delta scan.";
        return;
    end

    try
        if ~license('test', 'Distrib_Computing_Toolbox')
            msg = "Parallel Computing Toolbox license is not available. Falling back to serial delta scan.";
            return;
        end
    catch
        msg = "Could not verify Parallel Computing Toolbox license. Falling back to serial delta scan.";
        return;
    end

    pool = gcp('nocreate');
    if isempty(pool)
        if ~opt.AutoStartPool
            msg = "No parallel pool is open and AutoStartPool=false. Falling back to serial delta scan.";
            return;
        end

        poolType = char(opt.PoolType);
        nWorkers = opt.NumWorkers;
        if isempty(nWorkers)
            nWorkers = nD;
        else
            nWorkers = min(round(nWorkers), nD);
        end

        try
            try
                pool = parpool(poolType, nWorkers);
            catch
                pool = parpool(nWorkers);
            end
            msg = sprintf('Started parallel pool with %d workers.', pool.NumWorkers);
        catch ME
            msg = "Could not start a parallel pool (" + string(ME.message) + "). Falling back to serial delta scan.";
            return;
        end
    else
        msg = sprintf('Using existing parallel pool with %d workers.', pool.NumWorkers);
    end

    ok = true;
end
