function psi = psifun(obj, rp)
    psi = min([1.0; ...
        max([0.0; (rp - obj.rpmin)/(obj.rpmax - obj.rpmin)])]);
end