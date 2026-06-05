function eq = read_geqdsk(filename)
%READ_GEQDSK Read an EFIT GEQDSK file into a MATLAB struct.
%
%   eq = READ_GEQDSK(filename)
%
% Output struct fields typically include:
%   eq.header
%   eq.idum
%   eq.nw, eq.nh
%   eq.rdim, eq.zdim, eq.rcentr, eq.rleft, eq.zmid
%   eq.rmaxis, eq.zmaxis, eq.simag, eq.sibry, eq.bcentr
%   eq.current
%   eq.fpol, eq.pres, eq.ffprim, eq.pprime, eq.qpsi
%   eq.psirz
%   eq.nbbbs, eq.limitr
%   eq.rbbbs, eq.zbbbs, eq.rlim, eq.zlim
%   eq.rgrid, eq.zgrid
%   eq.psin
%
% Notes:
% - This reader is designed for standard EFIT GEQDSK formatting.
% - Numeric fields are read in free-format style after converting D exponents
%   to E, so it is robust to common Fortran-style output.

    fid = fopen(filename, 'r');
    if fid == -1
        error('read_geqdsk:FileOpenError', ...
              'Could not open file: %s', filename);
    end

    cleanupObj = onCleanup(@() fclose(fid));

    %-----------------------------
    % 1) Read first line (header)
    %-----------------------------
    firstLine = fgetl(fid);
    if ~ischar(firstLine)
        error('read_geqdsk:EmptyFile', 'File appears to be empty.');
    end

    [header, idum, nw, nh] = parse_first_line(firstLine);

    eq = struct();
    eq.filename = filename;
    eq.header   = header;
    eq.idum     = idum;
    eq.nw       = nw;
    eq.nh       = nh;

    %-----------------------------------------
    % 2) Read the rest of file as raw strings
    %-----------------------------------------
    raw = textscan(fid, '%s', 'Delimiter', '\n', 'Whitespace', '');
    lines = raw{1};

    if isempty(lines)
        error('read_geqdsk:UnexpectedEOF', ...
              'No numeric content found after header.');
    end

    % Join all lines, convert Fortran D exponents to E
    txt = strjoin(lines, ' ');
    txt = regexprep(txt, '[dD]', 'E');

    nums = sscanf(txt, '%f');

    if isempty(nums)
        error('read_geqdsk:ParseError', ...
              'Failed to parse numeric data from file.');
    end

    %-----------------------------------------
    % 3) Parse numeric content
    %-----------------------------------------
    k = 1;

    % First 20 scalars (5 lines x 4 values in standard GEQDSK)
    vals20 = take(nums, k, 20); k = k + 20;

    eq.rdim   = vals20(1);
    eq.zdim   = vals20(2);
    eq.rcentr = vals20(3);
    eq.rleft  = vals20(4);
    eq.zmid   = vals20(5);

    eq.rmaxis = vals20(6);
    eq.zmaxis = vals20(7);
    eq.simag  = vals20(8);
    eq.sibry  = vals20(9);
    eq.bcentr = vals20(10);

    eq.current = vals20(11);

    % Keep the remaining values too, since some GEQDSK variants reuse/store
    % duplicated placeholders in these slots.
    eq.scalar20_raw = vals20;

    % 1D profiles
    eq.fpol   = take(nums, k, nw); k = k + nw;
    eq.pres   = take(nums, k, nw); k = k + nw;
    eq.ffprim = take(nums, k, nw); k = k + nw;
    eq.pprime = take(nums, k, nw); k = k + nw;

    % 2D psi(R,Z), usually stored as (nw x nh) in file order
    psivec = take(nums, k, nw * nh); k = k + nw * nh;
    eq.psirz = reshape(psivec, [nw, nh]);

    % q profile
    eq.qpsi = take(nums, k, nw); k = k + nw;

    % Boundary counts
    counts = take(nums, k, 2); k = k + 2;
    eq.nbbbs = round(counts(1));
    eq.limitr = round(counts(2));

    % Plasma boundary points
    if eq.nbbbs > 0
        bbb = take(nums, k, 2 * eq.nbbbs); k = k + 2 * eq.nbbbs;
        eq.rbbbs = bbb(1:2:end);
        eq.zbbbs = bbb(2:2:end);
    else
        eq.rbbbs = [];
        eq.zbbbs = [];
    end

    % Limiter points
    if eq.limitr > 0
        lim = take(nums, k, 2 * eq.limitr); k = k + 2 * eq.limitr;
        eq.rlim = lim(1:2:end);
        eq.zlim = lim(2:2:end);
    else
        eq.rlim = [];
        eq.zlim = [];
    end

    %-----------------------------------------
    % 4) Construct useful grids
    %-----------------------------------------
    eq.rgrid = linspace(eq.rleft, eq.rleft + eq.rdim, eq.nw);
    eq.zgrid = linspace(eq.zmid - eq.zdim/2, eq.zmid + eq.zdim/2, eq.nh);
    eq.psin  = (eq.qpsi * 0);  % initialize same size

    dpsi = eq.sibry - eq.simag;
    if abs(dpsi) > 0
        eq.psin = linspace(0, 1, eq.nw);
    else
        eq.psin = nan(1, eq.nw);
    end

    % Optional sanity info
    eq.psi_axis = eq.simag;
    eq.psi_bdry = eq.sibry;
end

%======================================================================
function [header, idum, nw, nh] = parse_first_line(line)
% Parse GEQDSK first line:
% [header text .........] idum nw nh
%
% Standard files usually have the last three whitespace-separated tokens
% as integers.

    tokens = regexp(strtrim(line), '\s+', 'split');
    if numel(tokens) < 3
        error('read_geqdsk:HeaderParseError', ...
              'Could not parse first line: %s', line);
    end

    nh   = str2double(tokens{end});
    nw   = str2double(tokens{end-1});
    idum = str2double(tokens{end-2});

    if any(isnan([idum, nw, nh]))
        error('read_geqdsk:HeaderParseError', ...
              'Failed to parse idum/nw/nh from first line: %s', line);
    end

    if numel(tokens) > 3
        header = strjoin(tokens(1:end-3), ' ');
    else
        header = '';
    end
end

%======================================================================
function x = take(arr, startIdx, n)
% Safely take n numbers from arr starting at startIdx.
    endIdx = startIdx + n - 1;
    if endIdx > numel(arr)
        error('read_geqdsk:UnexpectedEOF', ...
              'File ended unexpectedly while parsing numeric data.');
    end
    x = arr(startIdx:endIdx);
end