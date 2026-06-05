function out = run_cshape_parameter_salpha_example(parameterName)
%RUN_CSHAPE_PARAMETER_SALPHA_EXAMPLE Run one C-shape parameter s-alpha scan.

    thisDir = fileparts(mfilename('fullpath'));
    addpath(thisDir);
    addpath(fullfile(thisDir, '..', 'Gaur'));

    geqdskFile = fullfile(thisDir, '..', 'Bishop', 'curavture_analysis', ...
        'geqdsk', 'scan_B2.5_C18_G08_H11.geqdsk');
    psin0 = 0.90;

    % Coarse first-pass scan. Increase these once the curve locations are
    % identified.
    ns = 16;
    nalpha = 16;
    ntheta0 = 12;
    sGrid = linspace(0.0, 3.0, ns);
    alphaGrid = linspace(0.0, 7.0, nalpha);
    theta0Grid = linspace(0.0, pi, ntheta0);

    eq = read_geqdsk(geqdskFile);
    base = cshape_local_equilibrium(eq, psin0, ...
        'NTheta', 161, ...
        'NGeom', 801, ...
        'UseBpFit', true);

    names = ["A", "B", "C", "G", "H"];
    pname = upper(string(parameterName));
    idx = find(names == pname, 1);
    if isempty(idx)
        error('run_cshape_parameter_salpha_example:BadParameterName', ...
            'parameterName must be one of A, B, C, G, or H.');
    end

    baseValue = base.coeffs(idx, 1);
    parameterValues = baseValue .* [0.7, 0.9, 1.1, 1.3];

    fprintf('C-shape %s parameter scan setup\n', pname);
    fprintf('  file        : %s\n', geqdskFile);
    fprintf('  psin        : %.6f\n', psin0);
    fprintf('  q profile   : %.8f\n', base.q_profile);
    fprintf('  q check     : %.8f\n', base.q_check);
    fprintf('  base coeffs : [A B C G H] = [% .8f % .8f % .8f % .8f % .8f]\n', ...
        base.coeffs(:,1));
    fprintf('  %s values    :', pname);
    fprintf(' %.8g', parameterValues);
    fprintf('\n');

    out = scan_cshape_parameter_salpha_curves(eq, psin0, char(pname), ...
        'ParameterValues', parameterValues, ...
        'SGrid', sGrid, ...
        'AlphaGrid', alphaGrid, ...
        'Theta0Grid', theta0Grid, ...
        'NThetaCoeff', 161, ...
        'NGeom', 801, ...
        'ThetaB', 5*pi, ...
        'NBalloon', 201, ...
        'NEigs', 6, ...
        'UseBpFit', true, ...
        'ShearDefinition', 'flux', ...
        'UseParallel', false, ...
        'Verbose', true);

    plot_cshape_parameter_marginal_curves(out);

    fprintf('C-shape %s s-alpha marginal-curve scan done.\n', pname);
    for iparam = 1:numel(out.parameterValues)
        lam = out.scans{iparam}.lambda_max;
        fprintf('  %s = %.8g: lambda range %.6g to %.6g\n', ...
            pname, out.parameterValues(iparam), ...
            min(lam(:), [], 'omitnan'), ...
            max(lam(:), [], 'omitnan'));
    end
end
