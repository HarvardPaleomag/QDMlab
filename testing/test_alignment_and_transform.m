d = subtract_blank('D:\data\N15data\Jul21_2020_FOV1\4x4Binned', 'D:\data\N15data\Jul21_2020_FOV1\4x4Binned')
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
[transForm, refFrame] = get_image_tform2(ima, imb,...
    'checkPlot', true, 'title', 'alignment');
%% bin image B by 4
imb_ = imresize(imb, 1/4, 'method', 'box');
%%
tformdata = tform_data(imb_, transForm, refFrame);
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
