function show_references()
%show_references()
    global messageDisplayed

    try
        caller = dbstack(1);
        callerName =  caller.name; % get functionName
    catch
        callerName = 'QDMdataprocessing';
    end
    
    if isempty(messageDisplayed)
        messageDisplayed.(callerName) = false;
    end

    messages.QDMdataprocessing = ['QDMdataprocessing.m -- MATLAB script for performing cropping, background compensation,\n',...
                                  'interpolation, calibration, and upward continuation of QDM data.\n',...
                                  'Processed data are saved to a .mat file and images are saved to .png files.\n\n',...
                                  'Reference: R.R. Fu, E.A. Lima, M.W.R. Volk, R. Trubko (2020) "High-sensitivity\n',...
                                  'moment magnetometry with the quantum diamond microscope,"  Geochemistry,\n',...
                                  'Geophysics, Geosystems.\n\n',...
                                  '===============================================================================================\n',...
                                  '  Matlab code by Eduardo A. Lima and Roger R. Fu\n',...
                                  '  Copyright (C) 2017-2020 MIT Paleomagnetism Lab, Harvard Paleomagnetics Lab\n',...
                                  '===============================================================================================\n'];

    messages.fit_resonance = ['===================================================================================================\n',...
                              '  GPUfit: Copyright (c) 2017 Mark Bates, Adrian Przybylski, Björn Thiel, and Jan Keller-Findeisen\n',...
                              '===================================================================================================\n'];
    
    if ~isfield(messages, callerName)
        return
    end
    
    if ~isfield(messageDisplayed, callerName)
        messageDisplayed.(callerName) = false;
    end

    if ~messageDisplayed.(callerName)
        fprintf(messages.(callerName))
        messageDisplayed.(callerName) = true;
    end
end
