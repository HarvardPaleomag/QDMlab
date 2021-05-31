d = subtract_blank('nFiles', '/Users/mike/Dropbox/science/_projects/QDMlab_paper/data/QDM_data/NRM_MIL/4x4Binned/B111dataToPlot.mat', ...
                   'blankFile', '/Users/mike/Dropbox/science/_projects/QDMlab_paper/data/QDM_data/blanks/Mar6_2020_2/4x4Binned/B111dataToPlot.mat')
dkeys = d.keys();
k = dkeys{1};
d(k)
%% more testing
a = 'D:\data\mike\NRM';
b = 'D:\data\mike\NRM';
%%
ima = imread(fullfile(a, 'LED.jpg'));
imb = imread(fullfile(b, 'LED.jpg'));
%%
[transForm, refFrame] = get_image_tform(ima, imb,...
    'checkPlot', true, 'title', 'alignment');
%% bin image B by 4
imb_ = imresize(imb, 1/4, 'method', 'box');
%%
tformdata = tform_data(imb_, transForm, refFrame, 4);
%%
figure
subplot(2,2,1)
imagesc(ima)
axis 'equal'; axis xy; axis tight;
subplot(2,2,2)
imagesc(imresize(ima, 1/4, 'method', 'box'))
axis 'equal'; axis xy; axis tight;
subplot(2,2,3)
imagesc(imb)
axis 'equal'; axis xy; axis tight;
subplot(2,2,4)
imagesc(tformdata)
axis 'equal'; axis xy; axis tight;
%%
[transForm_, refFrame_] = tform_bin_down(transForm, refFrame, 4);
tformdata_ = tform_data(imb_, transForm_, refFrame_);
%%
close all
figure
subplot(2,2,1)
imagesc(ima)
axis 'equal'; axis xy; axis tight;
subplot(2,2,2)
imagesc(imresize(ima, 1/4, 'method', 'box'))
axis 'equal'; axis xy; axis tight;
subplot(2,2,3)
imagesc(imb)
axis 'equal'; axis xy; axis tight;
subplot(2,2,4)
imagesc(tformdata_)
axis 'equal'; axis xy; axis tight;
