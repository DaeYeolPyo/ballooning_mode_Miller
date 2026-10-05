function run_dense_indented_cshape_beta_scan(outputDirectory,opts)
%RUN_DENSE_INDENTED_CSHAPE_BETA_SCAN Resume independent indentation columns.
%   Different MATLAB processes can run disjoint IndentationIndices. Each
%   pressure point is checkpointed by scan_indented_cshape_beta. Results are
%   assembled with collect_dense_indented_cshape_beta_scan(outputDirectory).
    arguments
        outputDirectory (1,1) string = ""
        opts.IndentationIndices double = 1:17
        opts.RefineBoundary (1,1) logical = true
        opts.SeedDirectory (1,1) string = ""
    end
    root=fileparts(fileparts(mfilename('fullpath')));
    addpath(root); setup_indented_cshape_scan;
    if strlength(outputDirectory)==0
        outputDirectory=string(fullfile(root,'output','indented_cshape_dense'));
    end
    config=make_dense_indented_cshape_scan_config;
    validateattributes(opts.IndentationIndices,{'double'}, ...
        {'vector','integer','positive','<=',numel(config.Indentation)});
    if numel(unique(opts.IndentationIndices))~=numel(opts.IndentationIndices)
        error('cshape:Dense:Indices','Each indentation index must occur once.');
    end
    if ~isfolder(outputDirectory), mkdir(outputDirectory); end
    for ii=opts.IndentationIndices(:).'
        one=config; one.Indentation=config.Indentation(ii);
        one.CheckpointFile=char(fullfile(outputDirectory,sprintf('column_%02d.mat',ii)));
        one.Resume=isfile(one.CheckpointFile);
        seeds={[]};
        if ~one.Resume && strlength(opts.SeedDirectory)>0
            seedFile=fullfile(opts.SeedDirectory,sprintf('column_%02d.mat',ii));
            if isfile(seedFile)
                data=load(seedFile,'scan');
                valid=find(~cellfun(@isempty,data.scan.equilibria(:,1)),1);
                if ~isempty(valid), seeds{1}=data.scan.equilibria{valid,1}; end
            end
        end
        scan=scan_indented_cshape_beta(one,InitialPredictors=seeds);
        traceFile=fullfile(outputDirectory,sprintf('trace_%02d.mat',ii));
        if opts.RefineBoundary
            if isfile(traceFile)
                saved=load(traceFile,'trace');
                if ~isfield(saved.trace,'signature') || ~isequaln(saved.trace.signature,scan.signature)
                    error('cshape:Dense:TraceMismatch','Trace configuration differs at column %d.',ii);
                end
            else
                trace=trace_indented_cshape_marginal(scan,PressureRelativeTolerance=.005, ...
                    MaxIterations=12);
                trace.signature=scan.signature;
                temporary=[tempname(char(outputDirectory)),'.mat'];
                save(temporary,'trace','-v7.3');
                movefile(temporary,traceFile,'f');
            end
        end
        fprintf('FINISHED COLUMN %d/%d I %.3f: %d stable, %d unstable, %d marginal, %d unknown\n', ...
            ii,numel(config.Indentation),one.Indentation,scan.numberStable, ...
            scan.numberUnstable,scan.numberMarginal,scan.numberUnknown);
    end
end
