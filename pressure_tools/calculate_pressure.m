function fit = calculate_pressure(filePath, reference, kwargs)

arguments
    filePath
    reference = 'Steele2017'
    kwargs.save = 1;
    kwargs.threshold = 0;
end

pShiftMap = containers.Map({'Doherty2014','Steele2017','Ivady2014'}, ...
    {14.58, 11.73, 10.3});    %MHz/GPa

pShift = pShiftMap(reference);

filePath = automatic_input_ui__(filePath, 'type','file');
fit = load(filePath);

% add newstyle B111 if it does not exist
if ~isfield(fit, 'B111')
    fit.B111 = B111(fit);
end

msg = sprintf('Calculating pressure using %s calibration (%.2f MHz/GPa)', reference, pShift);
logMsg('info',msg,1,0);
P = fit.B111.centerShift * 1000 /pShift;

chiSquares = fit.leftNeg.chiSquares + fit.leftPos.chiSquares + fit.rightNeg.chiSquares + fit.rightPos.chiSquares;

if kwargs.threshold > 0
    filter = ~isoutlier(chiSquares, 'quartiles', 'threshold', kwargs.threshold);
else
    filter = ones(size(chiSquares));
end

nFilter = numel(chiSquares)-numel(nonzeros(filter));

if nFilter > 0
    msg = sprintf('Removing %i / %i pixels because the chiSquared value were detected as outliers',...
        nFilter, numel(chiSquares));
    logMsg('waring',msg,1,0);
end

P(~filter) = nan;
P(P<0) = nan;

fit.P = P;
fit.medianP = median(fit.P(:), 'omitnan');
fit.stdP = std(fit.P(:), 'omitnan');

msg = sprintf('Saving Data into old structure >> %s', filePath);
logMsg('debug',msg,1,0);

if kwargs.save
    save(filePath, '-struct', 'fit');
end

end