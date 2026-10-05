function merluc = FluxSurface_Mercier_Luc(surf, eqfunc, opts)
%FLUXSURFACE_MERCIER_LUC Build local geometry directly from a flux surface.
%
%   merluc = FLUXSURFACE_MERCIER_LUC(surf, eqfunc)
%
% Unlike Miller_Mercier_Luc, this routine does not assume an analytic
% parameterization of the surface family.  It constructs the geometry from
% one extracted flux surface and evaluates |grad(psi)| from eqfunc.  The
% resulting structure satisfies the interface consumed by
% evaluate_theta_psi and evaluate_gcf.
%
% A truncated periodic Fourier representation supplies smooth first and
% second tangential derivatives.  The surface is projected back toward its
% requested psi level between smoothing passes.  Bp is smoothed in log
% space, then given one constant correction so its contour-integral q is
% exactly the supplied surface q.
%
% Name-value options
%   FourierModes       : retained positive Fourier modes (default 32)
%   ProjectionIterations : smooth/project passes (default 2)
%   EnforceQ           : apply constant Bp correction (default true)

    arguments
        surf (1,1) struct
        eqfunc (1,1) struct
        opts.FourierModes (1,1) double = 32
        opts.ProjectionIterations (1,1) double = 2
        opts.EnforceQ (1,1) logical = true
    end

    requiredSurfaceFields = { ...
        'R', 'Z', 'psi', 'psiN', 'F', 'q'};
    assert_fields(surf, requiredSurfaceFields, 'surf');

    requiredFunctionFields = {'eval'};
    assert_fields(eqfunc, requiredFunctionFields, 'eqfunc');

    Rdata = surf.R(:);
    Zdata = surf.Z(:);

    if numel(Rdata) ~= numel(Zdata) || numel(Rdata) < 32 || ...
            any(~isfinite(Rdata)) || any(~isfinite(Zdata))
        error('localEquilibrium:FluxSurface:InvalidCurve', ...
            'surf.R and surf.Z must define a finite closed curve.');
    end

    scale = max([range(Rdata), range(Zdata), 1]);
    if hypot(Rdata(end)-Rdata(1), Zdata(end)-Zdata(1)) ...
            <= 1.e-12*scale
        Rdata(end) = [];
        Zdata(end) = [];
    end

    nPoint = numel(Rdata);
    maxAvailableMode = floor((nPoint-1)/2);

    validateattributes(opts.FourierModes, {'numeric'}, ...
        {'real', 'finite', 'scalar', 'integer', '>=', 2, ...
         '<=', maxAvailableMode}, ...
        mfilename, 'FourierModes');

    validateattributes(opts.ProjectionIterations, {'numeric'}, ...
        {'real', 'finite', 'scalar', 'integer', '>=', 0, '<=', 10}, ...
        mfilename, 'ProjectionIterations');

    signedArea = 0.5*sum( ...
        Rdata.*circshift(Zdata, -1) ...
        - circshift(Rdata, -1).*Zdata);

    if ~isfinite(signedArea) || signedArea <= 0
        error('localEquilibrium:FluxSurface:Orientation', ...
            ['The input surface must be counter-clockwise. ', ...
             'Use extract_flux_surface with Orientation="outboard_up".']);
    end

    % Alternate periodic smoothing with normal projection onto the target
    % flux level.  This removes contour-grid corners without allowing the
    % filtered curve to drift appreciably across flux surfaces.
    for iteration = 1:opts.ProjectionIterations
        Rfiltered = periodic_fourier_filter( ...
            Rdata, opts.FourierModes);
        Zfiltered = periodic_fourier_filter( ...
            Zdata, opts.FourierModes);

        value = eqfunc.eval(Rfiltered, Zfiltered);
        gradPsiSquared = value.psiR.^2 + value.psiZ.^2;

        if any(~isfinite(gradPsiSquared)) || ...
                any(gradPsiSquared <= eps)
            error('localEquilibrium:FluxSurface:InvalidGradient', ...
                'The smoothed surface contains an invalid grad(psi).');
        end

        correction = (surf.psi-value.psi)./gradPsiSquared;

        Rdata = Rfiltered + correction.*value.psiR;
        Zdata = Zfiltered + correction.*value.psiZ;
    end

    [Rcore, RtCore, RttCore] = periodic_fourier_geometry( ...
        Rdata, opts.FourierModes);
    [Zcore, ZtCore, ZttCore] = periodic_fourier_geometry( ...
        Zdata, opts.FourierModes);

    surfaceValue = eqfunc.eval(Rcore, Zcore);
    gradPsiCore = surfaceValue.gradPsi(:);

    if any(~isfinite(gradPsiCore)) || any(gradPsiCore <= 0)
        error('localEquilibrium:FluxSurface:InvalidGradient', ...
            'The final surface contains an invalid |grad(psi)|.');
    end

    dlDtCore = hypot(RtCore, ZtCore);
    if any(~isfinite(dlDtCore)) || any(dlDtCore <= 0)
        error('localEquilibrium:FluxSurface:DegenerateCurve', ...
            'The Fourier surface has a zero or invalid tangent length.');
    end

    cosuRawCore = RtCore./dlDtCore;
    sinuRawCore = -ZtCore./dlDtCore;

    rawProjection = ...
        sinuRawCore.*surfaceValue.psiR(:) ...
        + cosuRawCore.*surfaceValue.psiZ(:);

    normalSign = sign(median(rawProjection));
    if normalSign == 0
        error('localEquilibrium:FluxSurface:UndefinedNormal', ...
            'Could not align the surface normal with grad(psi).');
    end

    sinuCore = normalSign*sinuRawCore;
    cosuCore = normalSign*cosuRawCore;

    alignedProjection = ...
        sinuCore.*surfaceValue.psiR(:) ...
        + cosuCore.*surfaceValue.psiZ(:);

    if any(alignedProjection <= 0)
        error('localEquilibrium:FluxSurface:NormalMisalignment', ...
            ['The smooth geometric normal is not aligned with ', ...
             'grad(psi) everywhere. Increase the retained Fourier modes.']);
    end

    invRcRawCore = ...
        (RtCore.*ZttCore-ZtCore.*RttCore)./dlDtCore.^3;
    invRcCore = normalSign*invRcRawCore;

    % The raw elementwise P3 gradient is C0-discontinuous across element
    % edges.  Smooth log(Bp) rather than Bp so positivity is preserved.
    BpRawCore = gradPsiCore./Rcore;
    logBpCore = periodic_fourier_filter( ...
        log(BpRawCore), opts.FourierModes);
    BpCore = exp(logBpCore);

    thetaCore = (0:nPoint-1).'*(2*pi/nPoint);

    R = close_periodic(Rcore);
    Z = close_periodic(Zcore);
    Rt = close_periodic(RtCore);
    Zt = close_periodic(ZtCore);
    Rtt = close_periodic(RttCore);
    Ztt = close_periodic(ZttCore);
    dlDt = close_periodic(dlDtCore);

    cosuRaw = close_periodic(cosuRawCore);
    sinuRaw = close_periodic(sinuRawCore);
    cosu = close_periodic(cosuCore);
    sinu = close_periodic(sinuCore);

    invRcRaw = close_periodic(invRcRawCore);
    invRc = close_periodic(invRcCore);
    Bp = close_periodic(BpCore);

    thetaGeo = [thetaCore; 2*pi];

    if ~isfinite(surf.F) || ~isfinite(surf.q) || surf.q == 0
        error('localEquilibrium:FluxSurface:InvalidSurfaceConstant', ...
            'surf.F and surf.q must be finite, and q must be nonzero.');
    end

    qBeforeCorrection = surf.F/(2*pi)*trapz( ...
        thetaGeo, dlDt./(R.^2.*Bp));

    qCorrection = qBeforeCorrection/surf.q;
    if ~isfinite(qCorrection) || qCorrection <= 0
        error('localEquilibrium:FluxSurface:InvalidQCorrection', ...
            ['The Bp contour integral has a sign inconsistent with ', ...
             'the supplied F and q.']);
    end

    if opts.EnforceQ
        Bp = qCorrection*Bp;
    else
        qCorrection = 1;
    end

    Bphi = surf.F./R;
    B2 = Bp.^2 + Bphi.^2;
    B = sqrt(B2);

    H = surf.F./(surf.q*(R.^2).*Bp);
    thetaPEST = cumtrapz(thetaGeo, H.*dlDt);

    qCheck = surf.F/(2*pi)*trapz( ...
        thetaGeo, dlDt./(R.^2.*Bp));

    % Use physical psi itself as the local radial coordinate.  These
    % derivatives satisfy grad(psi).e_psi=1 at the sampled points.
    psiRCore = surfaceValue.psiR(:);
    psiZCore = surfaceValue.psiZ(:);
    gradPsiSquared = psiRCore.^2 + psiZCore.^2;
    Rpsi = close_periodic(psiRCore./gradPsiSquared);
    Zpsi = close_periodic(psiZCore./gradPsiSquared);
    jacobian = Rpsi.*Zt-Rt.*Zpsi;

    merluc = struct();

    merluc.model = 'direct-flux-surface';

    merluc.theta_geo = thetaGeo;
    merluc.theta_PEST = thetaPEST;

    merluc.R = R;
    merluc.Z = Z;
    merluc.Rt = Rt;
    merluc.Zt = Zt;
    merluc.Rtt = Rtt;
    merluc.Ztt = Ztt;
    merluc.Rr = Rpsi;
    merluc.Zr = Zpsi;
    merluc.Rpsi = Rpsi;
    merluc.Zpsi = Zpsi;

    merluc.jac = jacobian;
    merluc.dl_dt = dlDt;

    merluc.cosu = cosu;
    merluc.sinu = sinu;
    merluc.cosu_raw = cosuRaw;
    merluc.sinu_raw = sinuRaw;
    merluc.normal_sign = normalSign;
    merluc.normal_projection_raw = close_periodic(rawProjection);

    merluc.invRc = invRc;
    merluc.invRc_raw = invRcRaw;

    merluc.dpsi_dr = 1;
    merluc.Bp = Bp;
    merluc.Bphi = Bphi;
    merluc.B2 = B2;
    merluc.B = B;
    merluc.H = H;

    merluc.F = surf.F;
    merluc.q = surf.q;
    merluc.q_check = qCheck;

    merluc.diagnostics = struct();
    merluc.diagnostics.FourierModes = opts.FourierModes;
    merluc.diagnostics.ProjectionIterations = ...
        opts.ProjectionIterations;
    merluc.diagnostics.qBeforeCorrection = qBeforeCorrection;
    merluc.diagnostics.qCorrection = qCorrection;
    merluc.diagnostics.qRelativeError = ...
        (qCheck-surf.q)/max(abs(surf.q), eps);
    merluc.diagnostics.psiNRms = sqrt(mean( ...
        (surfaceValue.psiN(:)-surf.psiN).^2));
    merluc.diagnostics.psiNMaximum = max(abs( ...
        surfaceValue.psiN(:)-surf.psiN));
    merluc.diagnostics.minimumNormalAlignment = ...
        min(alignedProjection./gradPsiCore);
    merluc.diagnostics.PESTClosure = thetaPEST(end)-2*pi;
    merluc.diagnostics.BpRaw = close_periodic(BpRawCore);
end


function filtered = periodic_fourier_filter(value, maximumMode)
    [filtered, ~, ~] = periodic_fourier_geometry(value, maximumMode);
end


function [filtered, firstDerivative, secondDerivative] = ...
    periodic_fourier_geometry(value, maximumMode)

    value = value(:);
    nPoint = numel(value);

    if mod(nPoint, 2) == 0
        modeNumber = [0:nPoint/2-1, 0, -nPoint/2+1:-1].';
    else
        modeNumber = [ ...
            0:(nPoint-1)/2, -(nPoint-1)/2:-1].';
    end

    spectrum = fft(value);
    spectrum(abs(modeNumber) > maximumMode) = 0;

    filtered = real(ifft(spectrum));
    firstDerivative = real(ifft(1i*modeNumber.*spectrum));
    secondDerivative = real(ifft(-(modeNumber.^2).*spectrum));
end


function value = close_periodic(coreValue)
    value = [coreValue(:); coreValue(1)];
end


function assert_fields(value, names, variableName)
    if ~isstruct(value) || ~isscalar(value)
        error('localEquilibrium:FluxSurface:InvalidStructure', ...
            '%s must be a scalar structure.', variableName);
    end

    for k = 1:numel(names)
        name = names{k};
        if ~isfield(value, name)
            error('localEquilibrium:FluxSurface:MissingField', ...
                '%s.%s is required.', variableName, name);
        end
    end
end
