%% N15
expData15 = load('D:\data\N15data\Jul21_2020_FOV1\run_00001.mat');
%%
[fit, initialGuess, badPixels] = fit_resonance(expData15, 4, 1, 'diamond', 'N15', 'checkPlot', 1);
%%
fits = GPU_fit('D:\data\N15data\Jul21_2020_FOV1', 4, 'diamond', 'N15');
%%
fits = QDM_lorentzian_fit('D:\data\N15data\Jul21_2020_FOV1', 4, 'diamond', 'N15', 'checkPlot', 0);
%%
data = load('D:\data\N15data\Jul21_2020_FOV1\4x4Binned\B111dataToPlot.mat');
%%
QDM_figure(data.B111para)


%% N14
expData14 = load('D:\data\mike\NRM\run_00000.mat');
%%
[binDataNorm, freq] = prepare_raw_data(expData14, 4, 1);
%%
[fit, initialGuess, badPixels] = fit_resonance(expData14, 4, 1, 'checkPlot', 1);
%%
fits = GPU_fit('D:\data\mike\NRM', 4, 'diamond', 'N14');
%%
fits = QDM_lorentzian_fit('D:\data\mike\NRM', 4, 'diamond', 'N14', 'checkPlot', 0);
%%
data = load('D:\data\mike\NRM\4x4Binned\B111dataToPlot.mat');
%%
QDM_figure(data.B111para)