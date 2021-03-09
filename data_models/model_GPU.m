function model = model_GPU(p, freq, kwargs)
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
    kwargs.diamond {mustBeMember(kwargs.diamond, ['N15', 'N14'])} = 'N14';
    kwargs.data {mustBeNumericOrLogical} = false;
    kwargs.checkPlot (1,1) {mustBeMember(kwargs.checkPlot, [1, 0])} = 0;
end

x = freq;
if strcmp(kwargs.diamond, 'N14')
    Ahyp = 0.002158; % longitudinal hyperfine for 14N

    model = p(6)...
        -p(3)*p(2).^2./((x-p(1)+Ahyp).^2+p(2).^2)...
        -p(4)*p(2).^2./((x-p(1)).^2+p(2).^2)...
        -p(5)*p(2).^2./((x-p(1)-Ahyp).^2+p(2).^2);
end
if strcmp(kwargs.diamond, 'N15')
    Ahyp = 0.0015; % 1/2 longitudinal hyperfine for 15N

    model = p(5)...
        -p(3)*p(2).^2./((x-p(1)+Ahyp).^2+p(2).^2)...
        -p(4)*p(2).^2./((x-p(1)-Ahyp).^2+p(2).^2);
end

if kwargs.checkPlot
    figure
    hold on
    plot(freq, squeeze(kwargs.data))
    plot(freq, model)
end