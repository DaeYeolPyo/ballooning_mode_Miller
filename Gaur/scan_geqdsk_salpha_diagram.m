function out = scan_geqdsk_salpha_diagram(eqInput, psin0, varargin)
%SCAN_GEQDSK_SALPHA_DIAGRAM Scan local ideal-ballooning stability from GEQDSK.
%
%   out = scan_geqdsk_salpha_diagram(eq, psin0)
%   out = scan_geqdsk_salpha_diagram(filename, psin0, 'Name', value, ...)
%
% The flux-surface geometry, q, F, and B are read from the GEQDSK
% equilibrium at psin0.  The local magnetic shear and pressure gradient are
% then overridden on an (s_hat, alpha) grid.  At each grid point the
% ballooning phase theta0 is scanned and the largest eigenvalue is retained.
%
% By default the scan assumes psin ~= rho^2, so
%   dq/dpsin = q*s_hat/(2*psin).
%
% The default alpha mapping is
%   alpha = -2*q^2*R0*<|grad psin|>/B_N^2 * (mu0*dp/dpsin).

    if exist('solve_ballooning_eigenvalue', 'file') ~= 2
        error('scan_geqdsk_salpha_diagram:MissingSolver', ...
            'solve_ballooning_eigenvalue.m is not on the MATLAB path.');
    end

    ip = inputParser;
    addParameter(ip, 'SGrid', linspace(0.0, 7.0, 15), @(x)isnumeric(x)&&isvector(x));
    addParameter(ip, 'AlphaGrid', linspace(0.0, 10.0, 19), @(x)isnumeric(x)&&isvector(x));
    addParameter(ip, 'Theta0Grid', linspace(0.0, pi, 5), @(x)isnumeric(x)&&isvector(x));
    addParameter(ip, 'JacMode', 4, @(x)isnumeric(x)&&isscalar(x));
    addParameter(ip, 'NThetaCoeff', 128, @(x)isnumeric(x)&&isscalar(x)&&x>=16);
    addParameter(ip, 'Dpsin', 1e-3, @(x)isnumeric(x)&&isscalar(x)&&x>0);
    addParameter(ip, 'ShearDefinition', 'flux', @(x)ischar(x)||isstring(x));
    addParameter(ip, 'AlphaFactor', [], @(x) isempty(x) || (isnumeric(x)&&isscalar(x)&&x~=0));
    addParameter(ip, 'ThetaB', 5*pi, @(x)isnumeric(x)&&isscalar(x)&&x>0);
    addParameter(ip, 'NBalloon', 151, @(x)isnumeric(x)&&isscalar(x)&&x>=51);
    addParameter(ip, 'NEigs', 4, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
    addParameter(ip, 'Verbose', true, @(x)islogical(x)&&isscalar(x));
    parse(ip, varargin{:});
    opt = ip.Results;

    if ischar(eqInput) || isstring(eqInput)
        eq = read_geqdsk(char(eqInput));
        sourceFile = char(eqInput);
    else
        eq = eqInput;
        if isfield(eq, 'filename')
            sourceFile = eq.filename;
        else
            sourceFile = '';
        end
    end

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
    alphaFactorUsed = NaN;

    if opt.Verbose
        fprintf('GEQDSK s-alpha scan at psin=%.6g: %d alpha x %d s x %d theta0 = %d solves\n', ...
            psin0, nA, nS, nT, nA*nS*nT);
        if strlength(string(sourceFile)) > 0
            fprintf('  file: %s\n', sourceFile);
        end
    end

    for ia = 1:nA
        for is = 1:nS
            vals = nan(1, nT);
            for it = 1:nT
                try
                    bal = eq_ballooning_coefficients(eq, psin0, ...
                        'JacMode', opt.JacMode, ...
                        'NTheta', opt.NThetaCoeff, ...
                        'Dpsin', opt.Dpsin, ...
                        'Theta0', theta0Grid(it), ...
                        'S_hat', sGrid(is), ...
                        'Alpha', alphaGrid(ia), ...
                        'ShearDefinition', opt.ShearDefinition, ...
                        'AlphaFactor', opt.AlphaFactor);

                    if ~isfinite(alphaFactorUsed) && isfield(bal, 'alphaFactor')
                        alphaFactorUsed = bal.alphaFactor;
                    end

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
            fprintf('  alpha row %d/%d done (alpha = %.4g)\n', ...
                ia, nA, alphaGrid(ia));
        end
    end

    out = struct();
    out.model = 'geqdsk-local';
    out.source_file = sourceFile;
    out.psin = psin0;
    out.s = sGrid;
    out.alpha = alphaGrid;
    out.theta0 = theta0Grid;
    out.lambda_max = lambdaMax;
    out.lambda_theta0 = lambdaMatrix;
    out.best_theta0 = bestTheta0;
    out.used_full_eig = usedFullEig;
    out.fail_message = failMessage;
    out.alphaFactor = alphaFactorUsed;
    out.eq = eq;
    out.options = opt;
    out.plot_title = sprintf('GEQDSK local $s$-$\\alpha$ diagram, $\\psi_N=%.3g$', psin0);
end
