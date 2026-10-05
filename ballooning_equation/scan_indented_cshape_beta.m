function scan = scan_indented_cshape_beta(config,opts)
%SCAN_INDENTED_CSHAPE_BETA Prescribed-q equilibria over shape and pressure.
%   Dimensions: pressure x indentation x aspect ratio. The pressure vector
%   is sorted; every requested point is attempted. No second stable region
%   is assumed. Checkpoints can be resumed using config.Resume=true.
%   InitialPredictors is an optional nIndentation-by-nAspect cell array of
%   converged equilibria used only as initial guesses for a new scan.
    arguments
        config (1,1) struct
        opts.InitialPredictors cell = {}
    end
    setup_indented_cshape_scan;
    config=validate_config(config);
    signature=rmfield(config,{'CheckpointFile','Resume','Verbose','RetryUnknown'});
    nP=numel(config.PressureAmplitude);
    nI=numel(config.Indentation); nA=numel(config.AspectRatio);
    checkpoint=char(config.CheckpointFile);
    if config.Resume
        if isempty(checkpoint) || ~isfile(checkpoint)
            error('cshape:Scan:Resume','Resume requires an existing checkpoint.');
        end
        saved=load(checkpoint,'scan');
        scan=saved.scan;
        if scan.version~=3 || ~isequaln(signature,scan.signature)
            error('cshape:Scan:ResumeMismatch','Checkpoint schema or configuration differs from this scan.');
        end
        scan.config=config;
        if config.RetryUnknown
            scan.done(scan.status=="unknown")=false;
        end
    else
        if ~isempty(checkpoint) && isfile(checkpoint)
            error('cshape:Scan:ExistingCheckpoint','Use Resume=true or choose a new checkpoint path.');
        end
        predictors=cell(nI,nA);
        if ~isempty(opts.InitialPredictors)
            if ~isequal(size(opts.InitialPredictors),[nI,nA])
                error('cshape:Scan:InitialPredictors','InitialPredictors must be nIndentation-by-nAspect.');
            end
            predictors=opts.InitialPredictors;
        end
        scan=struct('version',3,'config',config,'signature',signature, ...
            'pressureAmplitude',config.PressureAmplitude,'indentation',config.Indentation, ...
            'aspectRatio',config.AspectRatio,'points',{cell(nP,nI,nA)}, ...
            'shapes',{cell(nI,nA)},'done',false(nP,nI,nA), ...
            'betaPercent',nan(nP,nI,nA),'lambdaMax',nan(nP,nI,nA), ...
            'status',repmat("pending",nP,nI,nA),'equilibria',{cell(nP,nI,nA)}, ...
            'predictors',{predictors});
    end
    for ia=1:nA
        for ii=1:nI
            if all(scan.done(:,ii,ia)), continue; end
            shape=normalize_indented_cshape_boundary(config.AspectRatio(ia), ...
                config.Kappa,config.MinorRadius,config.Delta,config.Indentation(ii), ...
                NBoundary=config.Grid.nBoundary);
            scan.shapes{ii,ia}=shape;
            first=find(~scan.done(:,ii,ia),1);
            input=build_indented_cshape_scan_input(config,shape,config.PressureAmplitude(first));
            prepared=prepare_GS_problem(input);
            prepared.qEvaluation.contourCache=build_q_contour_cache( ...
                prepared.mesh,prepared.qEvaluation.contourNR);
            warm=scan.predictors{ii,ia};
            if isempty(warm) && first>1 && config.StoreEquilibria
                for j=first-1:-1:1
                    if ~isempty(scan.equilibria{j,ii,ia})
                        warm=scan.equilibria{j,ii,ia}; break
                    end
                end
            end
            for ip=1:nP
                if scan.done(ip,ii,ia), continue; end
                p0=config.PressureAmplitude(ip);
                [point,equilibrium,predictor]=evaluate_indented_cshape_scan_point(config,shape,p0,prepared,warm);
                % A failed solve must never overwrite the last converged predictor.
                if ~isempty(predictor), warm=predictor; end
                scan.predictors{ii,ia}=warm;
                scan.points{ip,ii,ia}=point;
                scan.betaPercent(ip,ii,ia)=point.betaPercent;
                scan.lambdaMax(ip,ii,ia)=point.lambdaMax;
                scan.status(ip,ii,ia)=point.status;
                scan.done(ip,ii,ia)=true;
                if config.StoreEquilibria, scan.equilibria{ip,ii,ia}=equilibrium; end
                if config.Verbose
                    fprintf('A %.4f I %.4f p0 %.5g: beta %.5g%%, lambda %+.5g, %s (%.1fs)\n', ...
                        shape.aspectRatio,shape.indent,p0,point.betaPercent, ...
                        point.lambdaMax,point.status,point.wallTime);
                    if ~isempty(point.failure), fprintf('  %s\n',point.failure.message); end
                    if isempty(point.failure) && ~isempty(point.ballooning) && ~point.complete
                        firstFailure=find(~cellfun(@isempty,point.ballooning.failures),1);
                        if ~isempty(firstFailure)
                            fprintf('  ballooning: %s\n',point.ballooning.failures{firstFailure}.message);
                        end
                    end
                end
                save_checkpoint(scan,checkpoint);
            end
        end
    end
    scan.finished=all(scan.done,'all');
    scan.numberStable=sum(scan.status=="stable",'all');
    scan.numberUnstable=sum(scan.status=="unstable",'all');
    scan.numberMarginal=sum(scan.status=="marginal",'all');
    scan.numberUnknown=sum(scan.status=="unknown",'all');
    save_checkpoint(scan,checkpoint);
end

function config=validate_config(config)
    if ~isfield(config,'RetryUnknown'), config.RetryUnknown=false; end
    validateattributes(config.RetryUnknown,{'logical'},{'scalar'});
    required={'AspectRatio','Kappa','MinorRadius','Delta','Indentation', ...
        'BReference','PressureAmplitude','QCoefficients','Grid','Picard', ...
        'PQConstraint','QEvaluation','Ballooning','BetaQuadratureDegree', ...
        'StoreEquilibria','CheckpointFile','Resume','Verbose'};
    if ~isstruct(config) || ~isscalar(config) || ~all(isfield(config,required))
        error('cshape:Scan:Config','Start with make_indented_cshape_scan_config.');
    end
    validateattributes(config.AspectRatio,{'double'},{'real','vector','nonempty','finite','>',1});
    validateattributes(config.Indentation,{'double'},{'real','vector','nonempty','finite','nonnegative'});
    validateattributes(config.PressureAmplitude,{'double'},{'real','vector','nonempty','finite','nonnegative'});
    for name={'Kappa','MinorRadius','BReference'}
        validateattributes(config.(name{1}),{'double'},{'real','scalar','finite','positive'});
    end
    validateattributes(config.Delta,{'double'},{'real','scalar','finite','>',-1,'<',1});
    validateattributes(config.QCoefficients,{'double'},{'real','vector','nonempty','finite'});
    for name={'AspectRatio','Indentation','PressureAmplitude'}
        key=name{1}; value=config.(key)(:).';
        if numel(unique(value))~=numel(value)
            error('cshape:Scan:Grid','Duplicate values in %s.',key);
        end
        config.(key)=sort(value);
    end
    if ~ischar(config.CheckpointFile) && ~(isstring(config.CheckpointFile) && isscalar(config.CheckpointFile))
        error('cshape:Scan:Checkpoint','CheckpointFile must be a text scalar.');
    end
    labels=config.Ballooning.RadialLabels;
    qlabels=config.QEvaluation.interiorPsiN;
    fit=config.PQConstraint.FFprime.fitRange;
    if min(labels)<min(qlabels)-1.e-12 || max(labels)>max(qlabels)+1.e-12 || ...
            min(labels)<fit(1) || max(labels)>fit(2)
        error('cshape:Scan:RadialCoverage', ...
            'Ballooning labels must lie inside the q-evaluation and constraint fit ranges.');
    end
end

function save_checkpoint(scan,path)
    if isempty(path), return; end
    folder=fileparts(path);
    if isempty(folder), folder=pwd; end
    if ~isfolder(folder), mkdir(folder); end
    temporary=[tempname(folder),'.mat'];
    cleanup=onCleanup(@()delete_if_present(temporary));
    save(temporary,'scan','-v7.3');
    [ok,message]=movefile(temporary,path,'f');
    if ~ok, error('cshape:Scan:CheckpointWrite','%s',message); end
end

function delete_if_present(path)
    if isfile(path), delete(path); end
end
