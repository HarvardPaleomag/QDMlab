function binning = detect_binning(inData)

% check for differences in LED naming
if isfield(inData, 'ledImg')
    led = inData.ledImg;
    data = inData.B111ferro;
elseif isfield(data, 'newLED')
    led = data.newLED;
    data = inData.Bz;
end
    
if size(led) ~= size(data)
    binning = (size(led,1) / size(data,1));
else
    binning = 1;
end

if binning > 1
    disp(['<>   binning detected: ' num2str(binning)])
end