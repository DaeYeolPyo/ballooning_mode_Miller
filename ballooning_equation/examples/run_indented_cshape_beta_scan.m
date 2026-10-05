function [scan, trace] = run_indented_cshape_beta_scan(outputDirectory)
%RUN_INDENTED_CSHAPE_BETA_SCAN Small starting scan; refine before publication.
%   Call from the examples directory or add it to the MATLAB path.
    root=fileparts(fileparts(mfilename('fullpath')));
    addpath(root);
    setup_indented_cshape_scan;
    if nargin<1, outputDirectory=fullfile(root,'output','indented_cshape'); end
    if ~isfolder(outputDirectory), mkdir(outputDirectory); end
    config=make_indented_cshape_scan_config;
    config.CheckpointFile=fullfile(outputDirectory,'beta_scan.mat');
    config.RetryUnknown=true;
    if isfile(config.CheckpointFile)
        saved=load(config.CheckpointFile,'scan');
        signature=rmfield(config,{'CheckpointFile','Resume','Verbose','RetryUnknown'});
        if saved.scan.version~=3 || ~isequaln(saved.scan.signature,signature)
            backup=fullfile(outputDirectory, ...
                ['beta_scan_previous_',char(datetime('now','Format','yyyyMMdd_HHmmssSSS')),'.mat']);
            [ok,message]=movefile(config.CheckpointFile,backup);
            if ~ok, error('cshape:Example:CheckpointBackup','%s',message); end
            fprintf('Preserved incompatible checkpoint: %s\n',backup);
        end
    end
    config.Resume=isfile(config.CheckpointFile);
    scan=scan_indented_cshape_beta(config);
    trace=trace_indented_cshape_marginal(scan);
    save(fullfile(outputDirectory,'marginal_trace.mat'),'trace','-v7.3');
    plot_indented_cshape_beta(scan,trace, ...
        ExportPath=string(fullfile(outputDirectory,'beta_stability.png')));
end
