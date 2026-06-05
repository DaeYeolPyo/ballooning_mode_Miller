function out = scan_miller_salpha_diagram(p, varargin)
%SCAN_MILLER_SALPHA_DIAGRAM Scan ideal-ballooning stability in Miller s-alpha space.
%
%   out = scan_miller_salpha_diagram(p)
%   out = scan_miller_salpha_diagram(p,'Name',value,...)
%
% The scan holds the non-s-alpha Miller shape parameters fixed and varies
% p.s_hat and p.alpha. At each point, it scans theta0 and keeps the largest
% ideal-ballooning eigenvalue. The marginal curve is lambda_max = 0.

    if exist('solve_ballooning_eigenvalue', 'file') ~= 2
        error('scan_miller_salpha_diagram:MissingSolver', ...
            ['solve_ballooning_eigenvalue.m is not on the MATLAB path. ', ...
             'Run addpath(fullfile(pwd,''..'',''Gaur'')) before starting the scan.']);
    end

    ip = inputParser;
    addParameter(ip, 'SGrid', linspace(0.0, 5.0, 13), @(x)isnumeric(x)&&isvector(x));
    addParameter(ip, 'AlphaGrid', linspace(0.0, 4.0, 13), @(x)isnumeric(x)&&isvector(x));
    addParameter(ip, 'Theta0Grid', linspace(0.0, pi, 5), @(x)isnumeric(x)&&isvector(x));
    addParameter(ip, 'NThetaCoeff', 161, @(x)isnumeric(x)&&isscalar(x)&&x>=32);
    addParameter(ip, 'NGeom', 601, @(x)isnumeric(x)&&isscalar(x)&&x>=128);
    addParameter(ip, 'ThetaB', 5*pi, @(x)isnumeric(x)&&isscalar(x)&&x>0);
    addParameter(ip, 'NBalloon', 201, @(x)isnumeric(x)&&isscalar(x)&&x>=51);
    addParameter(ip, 'NEigs', 6, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
    addParameter(ip, 'Verbose', true, @(x)islogical(x)&&isscalar(x));
    parse(ip, varargin{:});
    opt = ip.Results;

    sGrid = opt.SGrid(:).';
    alphaGrid = opt.AlphaGrid(:);
    theta0Grid = opt.Theta0Grid(:).';

    nS = numel(sGrid);
    nA = numel(alphaGrid);
    nT = numel(theta0Grid);

    lambdaMax = nan(nA, nS);
    lambdaMatrix = nan(nA, nS, nT);
    bestTheta0 = nan(nA, nS);
    usedFullEig = false(nA, nS, nT);
    failMessage = strings(nA, nS, nT);

    if opt.Verbose
        fprintf('Miller s-alpha scan: %d alpha x %d s x %d theta0 = %d solves\n', ...
            nA, nS, nT, nA*nS*nT);
    end

    for ia = 1:nA
        for is = 1:nS
            pScan = p;
            pScan.s_hat = sGrid(is);
            pScan.alpha = alphaGrid(ia);

            vals = nan(1, nT);
            for it = 1:nT
                try
                    bal = miller_ballooning_coefficients(pScan, ...
                        'Theta0', theta0Grid(it), ...
                        'NTheta', opt.NThetaCoeff, ...
                        'NGeom', opt.NGeom);

                    sol = solve_ballooning_eigenvalue(bal, ...
                        'Theta0', theta0Grid(it), ...
                        'ThetaB', opt.ThetaB, ...
                        'N', opt.NBalloon, ...
                        'NEigs', opt.NEigs, ...
                        'Plot', false);

                    if isfinite(sol.lambda_refined)
                        vals(it) = real(sol.lambda_refined);
                    else
                        vals(it) = real(sol.lambda);
                    end
                    usedFullEig(ia,is,it) = sol.used_full_eig;
                catch ME
                    vals(it) = NaN;
                    failMessage(ia,is,it) = string(ME.message);
                end
            end

            lambdaMatrix(ia,is,:) = vals;
            [lambdaMax(ia,is), imax] = max(vals, [], 'omitnan');
            if isfinite(lambdaMax(ia,is))
                bestTheta0(ia,is) = theta0Grid(imax);
            end
        end

        if opt.Verbose
            fprintf('  alpha row %d/%d done (alpha = %.4g)\n', ia, nA, alphaGrid(ia));
        end
    end

    out = struct();
    out.s = sGrid;
    out.alpha = alphaGrid;
    out.theta0 = theta0Grid;
    out.lambda_max = lambdaMax;
    out.lambda_theta0 = lambdaMatrix;
    out.best_theta0 = bestTheta0;
    out.used_full_eig = usedFullEig;
    out.fail_message = failMessage;
    out.params = p;
    out.options = opt;
end
