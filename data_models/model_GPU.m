function model = model_GPU(p, freq, kwargs)
%[model] = model_GPU(p, freq; 'diamond', 'data', 'checkPlot')
% calculates a model from the fitting parameters from GPU_fit
% depending on the diamond type
% 
% N14:
% ----
%     p[1]: freq. of center peak
%     p[2]: width of peaks
%     p[3,4,5]: conrast of each peak
%     p[6]: baseline
% N15:
% ----
%     p[1]: freq. of center peak
%     p[2]: width of peaks
%     p[3,4]: conrast of each peak
%     p[5]: baseline
% 
% Returns
% -------
%     y values for each frequency
    
arguments
    p double
    freq double
    kwargs.diamond {mustBeMember(kwargs.diamond, ...
        ['N15', 'N14', 'DAC', 'singlet', 'doublet', 'triplet', 'gaussian'])} = 'N15';
    kwargs.data {mustBeNumericOrLogical} = false;
    kwargs.checkPlot (1, 1) {mustBeBoolean(kwargs.checkPlot)} = false;
end

x = freq;
switch kwargs.diamond
    case {'N14', 'triplet'}
    msg = sprintf('calculating lorentzian triplet');
    logMsg('debug',msg,1,0);

    Ahyp = 0.002158; % longitudinal hyperfine for 14N

    model = p(6)...
        -p(3)*p(2).^2./((x-p(1)+Ahyp).^2+p(2).^2)...
        -p(4)*p(2).^2./((x-p(1)).^2+p(2).^2)...
        -p(5)*p(2).^2./((x-p(1)-Ahyp).^2+p(2).^2);
    case {'N15', 'doublet'}
        msg = sprintf('calculating lorentzian doublet');
        logMsg('debug',msg,1,0);

        Ahyp = 0.0015; % 1/2 longitudinal hyperfine for 15N

        model = p(5)...
            -p(3)*p(2).^2./((x-p(1)+Ahyp).^2+p(2).^2)...
            -p(4)*p(2).^2./((x-p(1)-Ahyp).^2+p(2).^2);
        
    case {'singlet', 'DAC'}
        msg = sprintf('calculating lorentzian signlet');
        logMsg('debug',msg,1,0);
        model = p(4)...
            -p(3)*p(2).^2./((x-p(1)).^2+p(2).^2);
    case {'gaussian'}
        msg = sprintf('calculating gaussian for DAC diamond'); 
        logMsg('debug',msg,1,0); 
        model = p(1) * exp(-(x - p(2)).^2 / (2 * p(3)^2)) + p(4)-1; 
end

if kwargs.checkPlot
    figure
    hold on
    plot(freq, squeeze(kwargs.data))
    plot(freq, model)
end