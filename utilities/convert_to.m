function data = convert_to(dataInG, unit)
%[data] = convert_to(dataInG, unit)
    switch unit
        case 'T'
            conv = 0.0001;
        case 'mT'
            conv = 0.1;
        case 'microT'
            conv = 100;
        case 'muT'
            conv = 100;
        case 'nT'
            conv = 100000;
        case 'G'
            conv = 1;
    end
    msg = sprintf('converting 1 G -> %i%s: ', conv, unit);
    logMsg('debug',msg,1,0);
    data = dataInG * conv;
end