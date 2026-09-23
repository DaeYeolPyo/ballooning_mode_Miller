function profile = evaluate_FFprime_model(model, psiN)
%EVALUATE_FFPRIME_MODEL Evaluate a direct reduced-order FF' model.
%
%   profile = evaluate_FFprime_model(model,psiN)
%
% FF' is evaluated directly from low-order coefficients.  F and G=F^2
% are then reconstructed by integration from the prescribed boundary F.

    narginchk(2, 2);

    requiredFields = { ...
        'coefficients', 'nBasis', 'enforceZeroEdge', 'Fboundary'};

    assert_required_fields(model, requiredFields, 'model');

    psiN = psiN(:);
    coefficients = model.coefficients(:);

    if numel(coefficients) ~= model.nBasis
        error('GS:FFprimeModel:CoefficientSize', ...
            'model.coefficients must contain model.nBasis entries.');
    end

    if numel(psiN) < 2 || any(~isfinite(psiN)) || ...
            any(diff(psiN) <= 0)
        error('GS:FFprimeModel:InvalidPsiN', ...
            'psiN must be finite and strictly increasing.');
    end

    endpointTolerance = 100*eps(max(1,max(abs(psiN))));

    if abs(psiN(1)) > endpointTolerance || ...
            abs(psiN(end)-1) > endpointTolerance
        error('GS:FFprimeModel:ProfileEndpoints', ...
            'psiN must span psi_N=0 to psi_N=1.');
    end

    psiN(1) = 0;
    psiN(end) = 1;

    basis = construct_FFprime_basis( ...
        psiN, model.nBasis, model.enforceZeroEdge);

    FFprime = basis*coefficients;

    profile = reconstruct_G_from_FFprime( ...
        psiN, FFprime, model.Fboundary);

    profile.basis = basis;
    profile.coefficients = coefficients;
    profile.representation = 'direct-FFprime-reduced-basis';
end


function assert_required_fields(value, requiredFields, valueName)
%ASSERT_REQUIRED_FIELDS Validate a scalar structure interface.

    if ~isstruct(value) || ~isscalar(value)
        error('GS:FFprimeModel:InvalidStructure', ...
            '%s must be a scalar structure.', valueName);
    end

    for k = 1:numel(requiredFields)
        if ~isfield(value, requiredFields{k})
            error('GS:FFprimeModel:MissingField', ...
                '%s.%s is required.', valueName, requiredFields{k});
        end
    end
end
