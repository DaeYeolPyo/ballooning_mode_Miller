function eq = eeq(filename)
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

    [header, idum, nw_in, nh_in] = parse_first_line(firstLine);

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
    eq = struct();
    k = 1;

    % First 20 scalars (5 lines x 4 values in standard GEQDSK)
    vals20 = take(nums, k, 20); k = k + 20;

    rwid      = vals20(1);
    zhei      = vals20(2);
    rcentr    = vals20(3);
    rleft     = vals20(4);
    zmid      = vals20(5);

    eq.R_mag  = vals20(6);
    eq.Z_mag  = vals20(7);
    eq.psi_0  = vals20(8);
    eq.psi_a  = vals20(9);
    eq.bcentr = vals20(10);

    % Keep the remaining values too, since some GEQDSK variants reuse/store
    % duplicated placeholders in these slots.
    scalar20_raw = vals20;

    % 1D profiles
    f         = take(nums, k, nw_in); k = k + nw_in;
    p         = take(nums, k, nw_in); k = k + nw_in;
    ffprim    = take(nums, k, nw_in); k = k + nw_in;
    pprime    = take(nums, k, nw_in); k = k + nw_in;

    % 2D psi(R,Z), usually stored as (nw x nh) in file order
    psivec    = take(nums, k, nw_in*nh_in); k = k + nw_in*nh_in;
    sefit_psi = reshape(psivec, [nw_in, nh_in]);

    % q profile
    q = take(nums, k, nw_in); k = k + nw_in;

    % Boundary counts
    counts = take(nums, k, 2); k = k + 2;
    nbbbs = round(counts(1));
    limitr = round(counts(2));

    % Plasma boundary points
    if nbbbs > 0
        bbb = take(nums, k, 2*nbbbs); k = k + 2*nbbbs;
        rbbbs = bbb(1:2:end);
        zbbbs = bbb(2:2:end);
    else
        rbbbs = [];
        zbbbs = [];
    end
    
    spsi_bar = linspace(0, 1, nw_in);
    sefit_R = rleft + linspace(0, rwid, nw_in);
    sefit_Z = zhei*(linspace(0, 1, nh_in) - 0.5);

    eq = setup(eq, nw_in, nh_in, spsi_bar, sefit_R, sefit_Z, ...
        f, p, q, sefit_psi, nbbbs, rbbbs, zbbbs, eq.bcentr, ...
        rleft, rwid, -0.5*zhei, zhei, 8, false);
end

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

function x = take(arr, startIdx, n)
% Safely take n numbers from arr starting at startIdx.
    endIdx = startIdx + n - 1;
    if endIdx > numel(arr)
        error('read_geqdsk:UnexpectedEOF', ...
              'File ended unexpectedly while parsing numeric data.');
    end
    x = arr(startIdx:endIdx);
end

function eq = setup(eq, nw_in, nh_in, spsi_bar, sefit_R, sefit_Z, ...
    f, p, q, sefit_psi, nbbbs, rbbbs, zbbbs, bcentr, rleft, rwid, zoffset, ...
    zhei, big, calc_aminor)

    eq.nw = nw_in*big;
    eq.nr = eq.nw;
    eq.nh = nh_in*big;
    eq.nt = eq.nh;

    eq.psi_bar = linspace(0, 1, eq.nw);
    eq.fp        = spline(spsi_bar, f, eq.psi_bar);
    eq.pressure  = spline(spsi_bar, p, eq.psi_bar);
    eq.qsf       = spline(spsi_bar, q, eq.psi_bar);

    eq_R = rleft + linspace(0, rwid, eq.nw);
    eq_Z = zoffset + linspace(0, zhei, eq.nh);

    [RR, ZZ] = meshgrid(eq_R, eq_Z);
    [sefit_RR, sefit_ZZ] = meshgrid(sefit_R, sefit_Z);
    eq.eqpsi_2d = interp2(sefit_RR, sefit_ZZ, sefit_psi,...
        RR, ZZ, "spline");

    thetab = atan2(zbbbs - eq.Z_mag, rbbbs - eq.R_mag);
    r_bound = hypot(rbbbs - eq.R_mag, zbbbs - eq.Z_mag);
    sort_abcd(thetab, r_bound, zbbbs, rbbbs);

    if thetab(1) == thetab(2)
        thetab(1) = thetab(1) + 2*pi;
        sort_abcd(thetab, r_bound, zbbbs, rbbbs);
    end

    if thetab(end-1) == thetab(end)
        thetab(end) = thetab(end) - 2*pi;
        sort_abcd(thetab, r_bound, zbbbs, rbbbs);
    end

    for i = 1:nbbbs-1
        if thetab(i+1) == thetab(i)
            thetab(i+1) = thetab(i+1) + 1.e-8;
        end
    end
    eq.thetab = thetab;
    eq.r_bound = r_bound;

    if calc_aminor
        eq.aminor = a_minor(rbbbs, zbbbs, eq.Z_mag);
    else
        eq.aminor = 1.0;
    end

    eq.B_psi = zeros(nr, nt);
    eq.diam = zeros(eq.nr, 1); eq.rc = zeros(eq.nr, 1);
    eq.r_bound = eq.r_bound/eq.aminor;
    eq.R_mag = eq.R_mag/eq.aminor;
    eq.Z_mag = eq.Z_mag/eq.aminor;
    eq_R = eq_R/eq.aminor;
    eq_Z = eq_Z/eq.aminor;
    eq.eq_R = eq_R; eq.eq_Z = eq_Z;
    [eq.R_psi, eq.Z_psi] = ndgrid(eq_R, eq_Z);

    eq.B_T = abs(bcentr);
    psi_N = eq.B_T*(eq.aminor^2);
    eq.psi_a = eq.psi_a/psi_N;
    eq.psi_0 = eq.psi_0/psi_N;
    eq.eqpsi_2d = eq.eqpsi_2d./psi_N;
    f_N = eq.B_T*eq.aminor;
    eq.fp = eq.fp/f_N;
    eq.beta = 8*pi*eq.pressure*1.e-7/(eq.B_T^2);
    eq.beta_0 = eq.beta(1);
    eq.pressure = eq.pressure/eq.pressure(1);

    eq.eqth = atan2(eq.Z_psi - eq.Z_mag, eq.R_psi - eq.R_mag);

    mask = (eq.Z_psi == eq.Z_mag) & (eq.R_psi == eq.R_mag);
    eq.eqth(mask) = 0;

    eq.eqpsi = eq.psi_bar*(eq.psi_a - eq.psi_0) + eq.psi_0;
    eq.has_full_theta_range = true;
end

function [a, b, c, d] = sort_abcd(a, b, c, d)
%SORT_ABCD Sort a in ascending order and reorder b, c, d identically.

    n = numel(a);

    if numel(b) ~= n || numel(c) ~= n || numel(d) ~= n
        error('sort_abcd:SizeMismatch', ...
              'a, b, c, and d must have the same number of elements.');
    end

    shapeA = size(a);
    shapeB = size(b);
    shapeC = size(c);
    shapeD = size(d);

    avec = a(:);
    bvec = b(:);
    cvec = c(:);
    dvec = d(:);

    [avec, idx] = sort(avec, 'ascend');

    a = reshape(avec, shapeA);
    b = reshape(bvec(idx), shapeB);
    c = reshape(cvec(idx), shapeC);
    d = reshape(dvec(idx), shapeD);
end

function a = a_minor(r, z, Z_mag)
    r = r(:);
    z = z(:);

    nz = 5;
    half_nz = floor(nz/2);
    invalidPoint = -nz;
    debug = false;

    n = numel(r);

    if numel(z) ~= n
        error('a_minor:SizeMismatch', ...
              'r and z must have the same number of elements.');
    end

    if n < nz
        error('a_minor:TooFewPoints', ...
              'Number of boundary points is less than nz = %d.', nz);
    end

    idx1 = (half_nz + 1):-1:1;

    ztmp1 = z(idx1);
    rtmp1 = r(idx1);

    [ztmp1, rtmp1] = remove_duplicate_z_for_interp(ztmp1, rtmp1);

    r1 = interp1(ztmp1, rtmp1, Z_mag, 'spline', 'extrap');

    i1 = invalidPoint;

    for i = nz:n
        if z(i) - Z_mag > 0
            i1 = i - 1;
            break;
        end
    end

    if i1 == invalidPoint
        error('a_minor:CrossingNotFound', ...
              ['Could not find point near magnetic-axis elevation ', ...
               'on low-field side.']);
    end

    k = 0;
    rtmp = zeros(nz, 1);
    ztmp = zeros(nz, 1);

    index = i1 - half_nz + k;

    if index < 1 || index > n
        error('a_minor:IndexOutOfRange', ...
              'Initial spline index is out of range.');
    end

    rtmp(1) = r(index);
    ztmp(1) = z(index);

    for ii = 2:nz
        index = i1 - half_nz + ii - 1 + k;

        if index < 1 || index > n
            error('a_minor:IndexOutOfRange', ...
                  'Spline index is out of range.');
        end

        rtmp(ii) = r(index);
        ztmp(ii) = z(index);

        dist2 = (rtmp(ii) - rtmp(ii-1))^2 + ...
                (ztmp(ii) - ztmp(ii-1))^2;

        if dist2 < 1.0e-7
            k = k + 1;

            if index + 1 > n
                error('a_minor:IndexOutOfRange', ...
                      'Duplicate-point correction index exceeds array size.');
            end

            rtmp(ii) = r(index + 1);
            ztmp(ii) = z(index + 1);
        end
    end

    [ztmp, rtmp] = remove_duplicate_z_for_interp(ztmp, rtmp);

    r2 = interp1(ztmp, rtmp, Z_mag, 'spline', 'extrap');

    a = (r2 - r1) / 2;
end

% ======================================================================
function [zout, rout] = remove_duplicate_z_for_interp(zin, rin)
%REMOVE_DUPLICATE_Z_FOR_INTERP
% MATLAB interp1 requires unique x-values. This helper removes duplicate
% z-values while preserving order.

    zin = zin(:);
    rin = rin(:);

    [zout, ia] = unique(zin, 'stable');
    rout = rin(ia);

    if numel(zout) < 2
        error('a_minor:TooFewUniqueZ', ...
              'Not enough unique z points for interpolation.');
    end

    % interp1 with spline prefers monotonic x-values.
    % The Fortran spline likely allowed arbitrary ordering through its own
    % implementation, but MATLAB interp1 wants sorted sample points.
    [zout, order] = sort(zout);
    rout = rout(order);
end