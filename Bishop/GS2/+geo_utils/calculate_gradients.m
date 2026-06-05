function calculate_gradients(obj)
    drm = obj.derm(obj.R_psi, 'E');
    dzm = obj.derm(obj.Z_psi, 'O');

    dtmp_m = obj.derm(obj.eqpsi_2d, 'E');
    obj.dpcart = obj.eqdcart(dtmp_m, drm, dzm);
    obj.dpbish = obj.eqdbish(obj.dpcart, obj.dpcart);

    dtmp_m = obj.derm(obj.B_psi, 'E');
    obj.dbcart = obj.eqdcart(dtmp_m, drm, dzm);
    obj.dbbish = obj.eqdbish(obj.dbcart, obj.dpcart);

    dtmp_m = obj.derm(obj.eqth, 'T');
    obj.dtcart = obj.eqdcart(dtmp_m, drm, dzm);
    obj.dtbish = obj.eqdbish(obj.dtcart, obj.dpcart);
end