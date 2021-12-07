b111 = load('/Users/mike/Dropbox/science/_projects/QDMlab_paper/data/QDM_data/NRM_MIL/4x4Binned/final_fits_(4x4).mat')
%%
pos = load('/Users/mike/Dropbox/science/_projects/QDMlab_paper/data/QDM_data/NRM_MIL/run_00000.mat')
neg = load('/Users/mike/Dropbox/science/_projects/QDMlab_paper/data/QDM_data/NRM_MIL/run_00001.mat')
%%
check_fits(b111, pos, neg)
%% 
close all

idxBad = [275,208]*4
idxGood = [291 200]*4

[lposG, rposG, lnegG, rnegG, freqList] = get_pixel_spectra(pos, neg, idxGood)
[lposB, rposB, lnegB, rnegB, freqList] = get_pixel_spectra(pos, neg, idxBad)
subplot(1,2,1)
hold on
plot(freqList(:,1), lnegB)
plot(freqList(:,1), lnegG)
xlabel('f [Hz]')
ylabel('contrast')
ylabel('contrast')
legend('bad pixel', 'good pixel')

subplot(1,2,2)
hold on
plot(freqList(:,2), rnegB)
plot(freqList(:,2), rnegG)
xlabel('f [Hz]')
ylabel('contrast')
legend('bad pixel', 'good pixel')