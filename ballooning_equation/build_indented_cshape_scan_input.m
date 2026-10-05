function input = build_indented_cshape_scan_input(config, shape, p0, warmEquilibrium)
%BUILD_INDENTED_CSHAPE_SCAN_INPUT Use physical p0 and dp/dpsiN inputs.
    if nargin<4, warmEquilibrium=[]; end
    validateattributes(p0,{'double'},{'scalar','real','finite','nonnegative'});
    n=config.Grid.nProfile;
    u=linspace(0,1,n).';
    q=polyval(fliplr(config.QCoefficients(:).'),u);
    if any(q<=0) || any(~isfinite(q))
        error('cshape:Input:Q','Target q must be finite and positive.');
    end
    input=struct();
    input.grid=config.Grid;
    input.dim=struct('R0',shape.Rgeo,'B0',config.BReference,'aMinor',shape.aMinor);
    input.boundary=struct('R',shape.boundary.R,'Z',shape.boundary.Z,'psiBoundary',0);
    input.Picard=config.Picard;
    input.profile=struct('mode',3,'pprime',-4*p0*u.*(1-u.^2), ...
        'q',struct('psiN',u,'q',q),'pedge',0, ...
        'normalization',struct('mode','fixed_Fb','Fb',shape.Rgeo*config.BReference), ...
        'FSign',1,'initial',struct('type','cylindrical'));
    input.pqConstraint=config.PQConstraint;
    if ~isempty(warmEquilibrium)
        initialFF=warmEquilibrium.profile.FdF_dpsiN;
        if isfield(config,'PressureCompensatedWarmStart') && config.PressureCompensatedWarmStart
            uw=warmEquilibrium.profile.psiN;
            deltaDp=-4*p0*uw.*(1-uw.^2)-warmEquilibrium.profile.dp_dpsiN;
            % Initial cylindrical source balance across a pressure step.
            % The full GS/q solve still determines the accepted FFprime.
            initialFF=initialFF-(4*pi*1.e-7)*shape.Rgeo^2*deltaDp;
        end
        input.profile.initial=struct('type','FFprime', ...
            'psiN',warmEquilibrium.profile.psiN, ...
            'FFprime',initialFF);
        % A converged predictor is the relevant branch. If the new pressure
        % fails, reduce the pressure step rather than switching FF' branches.
        input.pqConstraint.FFprime.initialScales=1;
    end
    input.qEvaluation=config.QEvaluation;
    input.output=struct('verbose',false,'checkTime',false,'showPlot',false);
end
