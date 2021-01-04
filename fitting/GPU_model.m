function model = GPU_model(p, freq, kwargs)
arguments
    p double
    freq double
    kwargs.diamond (1,1) {mustBeMember(kwargs.diamond,[14,15])} = 14;
    kwargs.data {mustBeNumericOrLogical} = false;
    kwargs.checkPlot (1,1) {mustBeMember(kwargs.checkPlot, [1, 0])} = 0;
end

x = freq;
if kwargs.diamond == 14
    Ahyp = 0.002158; % longitudinal hyperfine for 14N

    model = p(6)...
        -p(3)*p(2).^2./((x-p(1)+Ahyp).^2+p(2).^2)...
        -p(4)*p(2).^2./((x-p(1)).^2+p(2).^2)...
        -p(5)*p(2).^2./((x-p(1)-Ahyp).^2+p(2).^2);
end
kwargs
if kwargs.checkPlot
    figure
    hold on
    plot(freq, squeeze(kwargs.data))
    plot(freq, model)
end