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