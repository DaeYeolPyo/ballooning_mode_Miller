function qprime = s_to_qprime(s, ints)
    Vprime = ints.Vprime;
    V = ints.V;

    qprime = s*Vprime/(2*V);
end