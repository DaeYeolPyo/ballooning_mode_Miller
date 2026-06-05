function get_grids(obj)
end

function eik_get_grids(obj)
    
end

function salpha_get_grids(obj)
    obj.theta = linspace(-ntgrid, ntgrid, 2*ntgrid+1).*(2*pi/obj.ntheta);
end