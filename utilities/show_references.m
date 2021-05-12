function show_references()
    global messageDisplayed

    try
        caller = dbstack(1);
        callerName =  caller.name; % get functionName
    catch
        callerName = 'QDMdataprocessing';
    end
    
    messages.QDMdataprocessing = 'QDMdataprocessing.m -- MATLAB script for performing cropping, background compensation,\ninterpolation, calibration, and upward continuation of QDM data.\nProcessed data are saved to a .mat file and images are saved to .png\nfiles.\n\nReference: R.R. Fu, E.A. Lima, M.W.R. Volk, R. Trubko (2020) "High-sensitivity\nmoment magnetometry with the quantum diamond microscope,"  Geochemistry,\nGeophysics, Geosystems.\n\n-----------------------------------------------------------------------------------\n                Matlab code by Eduardo A. Lima and Roger R. Fu\n     Copyright (C) 2017-2020 MIT Paleomagnetism Lab, Harvard Paleomagnetics Lab\n-----------------------------------------------------------------------------------\n';
    
    if ~exists_struct(messages, callerName)
        return
    end
    
    if isempty(messageDisplayed)
        messageDisplayed.(callerName) = false;
    end
    
    if ~messageDisplayed.(callerName)
        fprintf(messages.(callerName))
        messageDisplayed.(callerName) = true;
    end
end