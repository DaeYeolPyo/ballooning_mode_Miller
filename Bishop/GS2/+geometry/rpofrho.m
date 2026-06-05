function rp = rpofrho(obj, rho)
    xerrbi = 1.e-4;
    xerrsec = 1.e-8;

    a = obj.rpmin; b = obj.rpmax;

    switch obj.irho
%         case 1
%             func = @(x) obj.phi(x);
%             fval = rho*rho*func(b);
%         case 2
%             func = @(x) obj.diameter(x);
%             fval = rho*func(b);
        case 3
            func = @(x) obj.psifun(x);
            fval = rho*func(b);
%         case 4
%             func = @(x) obj.rhofun(x);
%             fval = rho*func(b);
    end

    [rp, ier] = obj.root(func, fval, a, b, xerrbi, xerrsec);
    if ier > 0
        error('error in rpofrho');
    end
end