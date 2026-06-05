function geom = initialise_geometry_and_get_bmag(eq, rhoc, thetaShift, theta, bmag, bpol)
    geom = struct();

    geom.input = struct();
    geom.input.surf = struct();
    geom.input.thShift = thetaShift;
    geom.input.theta = theta;

    geom.eq = eq;

    [geom.eq, rpmin, rpmax, geom.input.surf.rmaj, B_T0, avgrmid] = initialize(geom.eq);

    rp = rpofrho(rhoc, rpmin, rpmax);
    drhodrp = drho_drp(rp, geom.input.surf.dr);
end

function soln = rpofrho(rho, rpmin, rpmax)
    xerrbi = 1.e-4; xerrsec = 1.e-8;
    a = rpmin; b = rpmax;

    % Selected irho = 3 case only.
    % It can be considered to use irho = 1, but not implemented
    % Please refer to geometry.f90 #1303-1311
    f = @(x) psifun(x, rpmin, rpmax);
    fval = rho*func*(b);

    [soln, ier] = root(f, fval, a, b, xerrbi, xerrsec);
    if ier > 0
        error('error in rpofrho');
    end
end

function [soln, ier] = root(f, fval, a, b, xerrbi, xerrsec)
    ier = 0;

    a1 = a; b1 = b;

    f1 = f(a1) - fval; f2 = f(b1) - fval;

    skip_bisection = (xerrbi <= 0);

    if ~skip_bisection && f1 * f2 > 0
        fprintf('f1 and f2 have same sign in bisection routine root\n');
        fprintf('a1,f1,b1,f2 = %.16g %.16g %.16g %.16g\n', a1, f1, b1, f2);
        fprintf('fval = %.16g\n', fval);

        ier = 1;
        skip_bisection = true;
    end

    if ~skip_bisection
        niter = 1 + floor(log(abs(b1 - a1) / xerrbi) / log(2));

        for i = 1:niter
            trial = 0.5 * (a1 + b1);
            f3 = f(trial) - fval;

            if f3 * f1 > 0
                a1 = trial;
                f1 = f3;
            else
                b1 = trial;
                f2 = f3;
            end
        end
    end

    if abs(f1) > abs(f2)
        f3 = f1;
        f1 = f2;
        f2 = f3;

        aold = a1;
        a1 = b1;
        b1 = aold;
    end

    for i = 1:10
        aold = a1;
        f3 = f2 - f1;

        if abs(f3) < 1.0e-11
            f3 = 1.0e-11;
        end

        a1 = a1 - f1 * (b1 - a1) / f3;

        f2 = f1;
        b1 = aold;
        f1 = f(a1) - fval;

        if abs(a1 - b1) < xerrsec
            break;
        end
    end

    soln = a1;

    if abs(f1) > xerrsec
        ier = ier + 2;
    end
end

function psi = psifun(rp, rpmin, rpmax)
    psi = min([1.0; max([0.0; (rp - rpmin)/(rpmax - rpmin)])]);
end

function DD = drho_drp(rp, dr, rpmin, rpmax)
    rp1 = rp*(1-dr); rp2 = rp*(1+dr);

    % Selected irho = 3 case only.
    % It can be considered to use irho = 1, but not implemented
    % Please refer to geometry.f90 #1303-1311
    rho1 = psifun(rp1, rpmin, rpmax);
    rho2 = psifun(rp2, rpmin, rpmax);

    DD = (rho2 - rho1)/(rp2 - rp1);
end

function [geom, rgrid] = rtgrid(geom, rp, theta)
    
end