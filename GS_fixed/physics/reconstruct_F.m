function F = reconstruct_F(psi,FBoundary,FFspec,fluxSpec)
%RECONSTRUCT_F Recover F=R*Bphi from F*dF/dpsi and its LCFS value.
%
%   For a tabulated GEQDSK profile, optional FFspec.FValues are used
%   directly. Otherwise F^2 is obtained from
%
%       F(psi)^2 = F_boundary^2 - 2*integral_psi^psi_b FF'(u) du.

narginchk(4,4);
validateattributes(FBoundary,{'numeric'}, ...
    {'real','finite','scalar','positive'},mfilename,'FBoundary');
if ~isstruct(FFspec) || ~isscalar(FFspec) || ~isfield(FFspec,'type')
    error('GS:physics:InvalidFFprimeSpec', ...
        'FFspec must be a scalar profile struct with type.');
end
requiredFlux = {'psiAxis','psiBoundary'};
if ~isstruct(fluxSpec) || ~isscalar(fluxSpec) || ...
        ~all(isfield(fluxSpec,requiredFlux))
    error('GS:physics:InvalidFluxSpec', ...
        'fluxSpec must contain psiAxis and psiBoundary.');
end

normalizationOptions = struct('rangePolicy','clip','tolerance',1e-8);
if isfield(fluxSpec,'rangePolicy')
    normalizationOptions.rangePolicy = fluxSpec.rangePolicy;
end
if isfield(fluxSpec,'tolerance')
    normalizationOptions.tolerance = fluxSpec.tolerance;
end
psiN = normalize_flux(psi,fluxSpec.psiAxis, ...
    fluxSpec.psiBoundary,normalizationOptions);
span = fluxSpec.psiBoundary-fluxSpec.psiAxis;
type = lower(strrep(strtrim(char(FFspec.type)),'_','-'));

switch type
    case 'constant'
        integralToBoundary = FFspec.value*span.*(1-psiN);
        F2 = FBoundary^2-2*integralToBoundary;

    case 'normalized-power'
        integralToBoundary = FFspec.axisValue*span ...
            .*(1-psiN).^(FFspec.exponent+1)/(FFspec.exponent+1);
        F2 = FBoundary^2-2*integralToBoundary;

    case 'f-power'
        F = FFspec.boundaryF+(FFspec.axisF-FFspec.boundaryF) ...
            .*(1-psiN).^FFspec.exponent;
        check_F(F);
        return

    case 'tabulated'
        if isfield(FFspec,'FValues') && ~isempty(FFspec.FValues)
            Fspec = struct('psiN',FFspec.psiN, ...
                'values',FFspec.FValues,'method',FFspec.method);
            F = evaluate_tabulated_profile(psiN,Fspec,'F');
            boundaryMismatch = abs(Fspec.values(end)-FBoundary) ...
                /max(abs(FBoundary),eps);
            if boundaryMismatch>1e-8
                error('GS:physics:InconsistentFBoundary', ...
                    ['The tabulated F edge value and FBoundary differ by ' ...
                     'more than the allowed relative tolerance.']);
            end
            check_F(F);
            return
        end
        grid = unique([linspace(0,1,4097).';FFspec.psiN(:)]);
        profileValues = evaluate_tabulated_profile(grid,FFspec,"FFprime");
        primitive = cumtrapz(grid,profileValues);
        integralGrid = span*(primitive(end)-primitive);
        integralToBoundary = interp1(grid,integralGrid,psiN,'pchip');
        F2 = FBoundary^2-2*integralToBoundary;

    otherwise
        error('GS:physics:UnsupportedFFprimeProfile', ...
            'Cannot reconstruct F for FFprime profile type %s.',type);
end

if any(F2(:)<=0) || any(~isfinite(F2(:)))
    error('GS:physics:NonPositiveFSquared', ...
        'The selected FFprime profile and FBoundary produce F^2<=0.');
end
F = sqrt(F2);
end


function check_F(F)
if ~isreal(F) || any(~isfinite(F(:))) || any(F(:)<=0)
    error('GS:physics:InvalidReconstructedF', ...
        'Reconstructed F must contain positive finite real values.');
end
end
