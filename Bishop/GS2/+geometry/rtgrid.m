function rgrid = rtgrid(obj, rp, theta)
    rgrid = zeros(2*obj.ntgrid + 1, 1);

    for i = 1:length(rgrid)
        rgrid(i) = obj.geom.rfun(rp, theta(i));
    end
end