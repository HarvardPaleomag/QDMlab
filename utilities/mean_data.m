function meanData = mean_data(data)
m1 = nanmean(data,1);
m2 = nanmean(squeeze(m1),1);
meanData = squeeze(m2);