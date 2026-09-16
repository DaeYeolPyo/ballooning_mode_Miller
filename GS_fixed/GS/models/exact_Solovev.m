function model = exact_Solovev(params)
%EXACT_SOLOVEV General analytic Solov'ev equilibrium
%
%   psi(R,Z) =
%       A1/8 * R^4
%       + A2/2 * Z^2
%       + c1
%       + c2 * R^2
%       + c3 * (R^4 - 4*R^2*Z^2)
%       + c4 * (R^2*log(R) - Z^2)
%
% Project sign convention:
%
%   -DeltaStar(psi)
%       = mu0*R^2*dp/dpsi + F*dF/dpsi
%
% Therefore:
%
%   dp/dpsi   = -A1/mu0
%   F*dF/dpsi = -A2
%
% Required params fields:
%   A1, A2, c1, c2, c3, c4
%   axisBracket = [Rleft, Rright]
%   rhoMax      = maximum boundary-search distance
%
% Optional:
%   mu0
%   psiBoundary
%   boundaryScanPoints

    narginchk(1,1);

    if ~isstruct(params) || ~isscalar(params)
        error('[EXACT_SOLOVEV] params must be a scalar structure.');
    end

    requiredFields = { ...
        'A1','A2', ...
        'c1','c2','c3','c4', ...
        'axisBracket','rhoMax'};

    for k = 1:numel(requiredFields)
        name = requiredFields{k};

        if ~isfield(params,name)
            error('[EXACT_SOLOVEV] Missing parameter: %s',name);
        end
    end

    p = params;

    if ~isfield(p,'mu0')
        p.mu0 = 4*pi*1.e-7;
    end

    if ~isfield(p,'psiBoundary')
        p.psiBoundary = 0.0;
    end

    if ~isfield(p,'boundaryScanPoints')
        p.boundaryScanPoints = 500;
    end

    scalarFields = { ...
        'A1','A2','c1','c2','c3','c4', ...
        'mu0','psiBoundary','rhoMax'};

    for k = 1:numel(scalarFields)
        name = scalarFields{k};

        validateattributes(p.(name),{'numeric'}, ...
            {'real','finite','scalar'}, ...
            mfilename,['params.',name]);
    end

    validateattributes(p.mu0,{'numeric'}, ...
        {'real','finite','scalar','positive'}, ...
        mfilename,'params.mu0');

    validateattributes(p.rhoMax,{'numeric'}, ...
        {'real','finite','scalar','positive'}, ...
        mfilename,'params.rhoMax');

    validateattributes(p.boundaryScanPoints,{'numeric'}, ...
        {'real','finite','scalar','integer','>=',20}, ...
        mfilename,'params.boundaryScanPoints');

    validateattributes(p.axisBracket,{'numeric'}, ...
        {'real','finite','vector','numel',2,'positive'}, ...
        mfilename,'params.axisBracket');

    axisBracket = p.axisBracket(:).';

    if axisBracket(1) >= axisBracket(2)
        error(['[EXACT_SOLOVEV] axisBracket must satisfy ', ...
               'axisBracket(1) < axisBracket(2).']);
    end

    % =========================================================
    % Exact flux
    % =========================================================

    psi = @(R,Z) ...
          (p.A1/8).*R.^4 ...
        + (p.A2/2).*Z.^2 ...
        + p.c1 ...
        + p.c2.*R.^2 ...
        + p.c3.*(R.^4-4.*R.^2.*Z.^2) ...
        + p.c4.*(R.^2.*log(R)-Z.^2);

    % =========================================================
    % Exact first derivatives
    % =========================================================

    dpsi_dR = @(R,Z) ...
          (p.A1/2).*R.^3 ...
        + 2.*p.c2.*R ...
        + p.c3.*(4.*R.^3-8.*R.*Z.^2) ...
        + p.c4.*(2.*R.*log(R)+R) ...
        + 0.*Z;

    dpsi_dZ = @(R,Z) ...
          p.A2.*Z ...
        - 8.*p.c3.*R.^2.*Z ...
        - 2.*p.c4.*Z ...
        + 0.*R;

    % =========================================================
    % Magnetic axis
    %
    % This solution is up-down symmetric, so an ordinary
    % magnetic axis lies on Z=0.
    % =========================================================

    Raxis = fzero( ...
        @(R) dpsi_dR(R,0.0), ...
        axisBracket);

    Zaxis = 0.0;

    psiAxis = psi(Raxis,Zaxis);
    psiBoundary = p.psiBoundary;

    dpsiExact = psiBoundary-psiAxis;

    fluxScale = max([ ...
        1,abs(psiAxis),abs(psiBoundary)]);

    if abs(dpsiExact) <= 100*eps(fluxScale)
        error(['[EXACT_SOLOVEV] Axis and boundary fluxes ', ...
               'are indistinguishable.']);
    end

    % Normalized flux:
    %
    %   psiN = 0 at the magnetic axis
    %   psiN = 1 at the boundary

    psiN = @(R,Z) ...
        (psi(R,Z)-psiAxis)./dpsiExact;

    % =========================================================
    % Grad-Shafranov operator
    %
    % DeltaStar(psi) = A1*R^2 + A2
    % =========================================================

    deltaStarPsi = @(R,Z) ...
        p.A1.*R.^2+p.A2+0.*Z;

    % =========================================================
    % Physical profile derivatives
    % =========================================================

    pprimePhysical = -p.A1/p.mu0;
    FFprimePhysical = -p.A2;

    % Derivatives with respect to normalized flux
    %
    % dp/dpsiN = dp/dpsi * dpsi/dpsiN
    %           = dp/dpsi * dpsiExact

    dp_dpsiN = ...
        pprimePhysical*dpsiExact;

    FdF_dpsiN = ...
        FFprimePhysical*dpsiExact;

    % Strong-form right-hand side:
    %
    %   -DeltaStar(psi)
    %       = -A1*R^2-A2

    strongSource = @(R,Z) ...
        -p.A1.*R.^2-p.A2+0.*Z;

    % Weak-form source used by the present FEM code:
    %
    %   mu0*R*pprime + FFprime/R
    %       = -A1*R-A2/R

    weakSource = @(R,Z) ...
        -p.A1.*R-p.A2./R+0.*Z;

    % =========================================================
    % Package result
    % =========================================================

    model = struct();

    model.parameters = p;

    model.A1 = p.A1;
    model.A2 = p.A2;

    model.c1 = p.c1;
    model.c2 = p.c2;
    model.c3 = p.c3;
    model.c4 = p.c4;

    model.axisPoint   = [Raxis,Zaxis];
    model.psiAxis     = psiAxis;
    model.psiBoundary = psiBoundary;
    model.dpsi        = dpsiExact;

    if psiAxis > psiBoundary
        model.axisMode = 'max';
    else
        model.axisMode = 'min';
    end

    model.pprimePhysical  = pprimePhysical;
    model.FFprimePhysical = FFprimePhysical;

    model.pprime  = dp_dpsiN;
    model.FFprime = FdF_dpsiN;

    model.psi       = psi;
    model.psiN      = psiN;
    model.dpsi_dR   = dpsi_dR;
    model.dpsi_dZ   = dpsi_dZ;

    model.deltaStarPsi = deltaStarPsi;
    model.strongSource = strongSource;
    model.weakSource   = weakSource;

    % The LCFS is now an implicit level set. Find its first
    % intersection along each ray from the magnetic axis.

    model.boundary = @(theta) ...
        trace_Solovev_boundary( ...
            theta,psi,[Raxis,Zaxis], ...
            psiBoundary,p.rhoMax, ...
            p.boundaryScanPoints);

    model.makeProfiles = @(nw) ...
        make_Solovev_profiles( ...
            nw,dp_dpsiN,FdF_dpsiN);
end


function boundary = trace_Solovev_boundary( ...
    theta,psi,axisPoint,psiBoundary,rhoMax,nScan)

    theta = theta(:);

    nTheta = numel(theta);
    boundary = zeros(nTheta,2);

    Raxis = axisPoint(1);
    Zaxis = axisPoint(2);

    for k = 1:nTheta
        ct = cos(theta(k));
        st = sin(theta(k));

        rhoLimit = rhoMax;

        % Prevent evaluation of log(R) at R <= 0.
        if ct < 0
            positiveRLimit = ...
                (1-1.e-12)*Raxis/(-ct);

            rhoLimit = min(rhoLimit,positiveRLimit);
        end

        rhoScan = linspace(0,rhoLimit,nScan).';

        Rscan = Raxis+rhoScan.*ct;
        Zscan = Zaxis+rhoScan.*st;

        levelValue = ...
            psi(Rscan,Zscan)-psiBoundary;

        crossing = find( ...
            levelValue(1:end-1).*levelValue(2:end) <= 0, ...
            1,'first');

        if isempty(crossing)
            error([ ...
                '[EXACT_SOLOVEV] Could not locate boundary ', ...
                'at theta = %.8f rad. Increase rhoMax or ', ...
                'check whether the selected contour is closed.'], ...
                theta(k));
        end

        rhoLeft  = rhoScan(crossing);
        rhoRight = rhoScan(crossing+1);

        boundaryEquation = @(rho) ...
            psi( ...
                Raxis+rho.*ct, ...
                Zaxis+rho.*st) ...
            - psiBoundary;

        rhoBoundary = fzero( ...
            boundaryEquation, ...
            [rhoLeft,rhoRight]);

        boundary(k,1) = ...
            Raxis+rhoBoundary.*ct;

        boundary(k,2) = ...
            Zaxis+rhoBoundary.*st;
    end
end


function profiles = make_Solovev_profiles( ...
    nw,dp_dpsiN,FdF_dpsiN)

    validateattributes(nw,{'numeric'}, ...
        {'real','finite','scalar','integer','>=',2}, ...
        mfilename,'nw');

    profiles = struct();

    profiles.psiN = linspace(0,1,nw).';

    profiles.dp_dpsiN = ...
        dp_dpsiN*ones(nw,1);

    profiles.FdF_dpsiN = ...
        FdF_dpsiN*ones(nw,1);
end