function [theta, R, Z, lambda] = get_raw_surface_data(surface, name)
    surf = surface.surf;
    ang = surface.ang;

    R = surf.R;
    Z = surf.Z;

    switch name
        case "PEST"
            theta = ang.theta_PEST;
            lambda = zeros(size(theta));

        case "Hamada"
            theta = ang.theta_Hamada;
            lambda = ang.lambda_Hamada;

        case "Boozer"
            theta = ang.theta_Boozer;
            lambda = ang.lambda_Boozer;

        otherwise
            error('validate_sfl_remap:BadName', ...
                  'Unknown coordinate name: %s', name);
    end
end