function out = scan_cshape_salpha_diagram(eq, psin0, varargin)
%SCAN_CSHAPE_SALPHA_DIAGRAM Scan ideal-ballooning stability for C-shape.
%
%   out = scan_cshape_salpha_diagram(eq, psin0)
%   out = scan_cshape_salpha_diagram(filename, psin0, 'Name', value, ...)
%
% This mirrors scan_miller_salpha_diagram.m, but the local geometry is
% fitted from the target GEQDSK flux surface using
%
%   R = A + B sin(theta) + C cos(2 theta),
%   Z = G cos(theta) - H sin(2 theta).
%
% For each (s_hat, alpha), theta0 is scanned and the largest ballooning
% eigenvalue is retained.  The marginal curve is lambda_max = 0.

    if exist('solve_ballooning_eigenvalue', 'file') ~= 2
        error('scan_cshape_salpha_diagram:MissingSolver', ...
            ['solve_ballooning_eigenvalue.m is not on the MATLAB path. ', ...
             'Run addpath(fullfile(pwd,''..'',''Gaur'')) before starting the scan.']);
    end

    ip = inputParser;
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
    addParameter(ip, 'CoeffOverride', [], @(x) isempty(x) || isnumeric(x));
    addParameter(ip, 'IncludeFPrime', true, @(x)islogical(x)&&isscalar(x));
    addParameter(ip, 'Verbose', true, @(x)islogical(x)&&isscalar(x));
    parse(ip, varargin{:});
    opt = ip.Results;

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
        'AlphaFactor', opt.AlphaFactor, ...
        'CoeffOverride', opt.CoeffOverride);

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
        fprintf('C-shape s-alpha scan: %d alpha x %d s x %d theta0 = %d solves\n', ...
            nA, nS, nT, nA*nS*nT);
        fprintf('  psin=%.6g, q_profile=%.8g, q_check=%.8g\n', ...
            psin0, base.q_profile, base.q_check);
        fprintf('  coeffs [A B C G H]=[% .6g % .6g % .6g % .6g % .6g]\n', ...
            base.coeffs(:,1));
        fprintf('  alphaFactor=%.12e, shearDefinition=%s\n', ...
            base.alphaFactor, char(lower(string(opt.ShearDefinition))));
    end

    for ia = 1:nA
        for is = 1:nS
            vals = nan(1, nT);

            for it = 1:nT
                try
                    bal = cshape_ballooning_coefficients(base, ...
                        'S_hat', sGrid(is), ...
                        'Alpha', alphaGrid(ia), ...
                        'Theta0', theta0Grid(it), ...
                        'ShearDefinition', opt.ShearDefinition, ...
                        'IncludeFPrime', opt.IncludeFPrime);

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
    out.model = 'cshape-local';
    out.eqInput = eqInput;
    out.psin = psin0;
    out.s = sGrid;
    out.alpha = alphaGrid;
    out.theta0 = theta0Grid;
    out.lambda_max = lambdaMax;
    out.lambda_theta0 = lambdaMatrix;
    out.best_theta0 = bestTheta0;
    out.used_full_eig = usedFullEig;
    out.fail_message = failMessage;
    out.base = base;
    out.options = opt;
end
