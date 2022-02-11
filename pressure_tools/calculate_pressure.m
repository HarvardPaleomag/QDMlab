function fit = calculate_pressure(filePath, reference)

arguments
    filePath
    reference = 'Steele2017'
end

pShiftMap = containers.Map({'Doherty2014','Steele2017','Ivady2014'}, ...
    {14.58, 11.73, 10.3});    %MHz/GPa
pShiftMap.keys()
pShift = pShiftMap(reference);

filePath = automatic_input_ui__(filePath, 'type','file');
fit = load(filePath);

% add newstyle B111 if it does not exist
if ~isfield(fit, 'B111')
    fit.B111 = B111(fit);
end

msg = sprintf('Calculating pressure using %s calibration (%.2f MHz/GPa)', reference, pShift);
logMsg('info',msg,1,0);
fit.P = fit.B111.centerShift * 1000 /pShift;
msg = sprintf('Saving Data into old structure >> %s', filePath);
logMsg('debug',msg,1,0);
save(filePath, '-struct', 'fit');
end