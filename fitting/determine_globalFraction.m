function [globalFraction, GF, idx] = determine_globalFraction(expData, binSize, nRes, GFs)
    arguments
        expData
        binSize
        nRes
        GFs = 0:0.05:0.8
    end
    fits = fit_resonance(expData, binSize, nRes, 'globalFraction', GFs);
    nElem = 1000;
    GF = zeros(1,nElem);

    slope = 0.0106; %determined from fitting Chi2 values of blank map
    zfs = 2.870;
    GF0 = reshape(fits.resonance(:,:,1), 1,[]);
    chi0 = reshape(fits.chiSquares(:,:,1),1,[]);
    [res, idx] = sort(abs(GF0-zfs),'descend');
    idx = idx(1:nElem);

    for i = 2:size(GFs,2)
        chi = reshape(fits.chiSquares(:,:,i), 1,[]);
        chi = chi(idx);
        chiSmaller = chi - (slope * GFs(i)) < chi0(idx);

        GF(chiSmaller) = GFs(i);
    end
    globalFraction = mean(GF(GF > 0));
%     GF = zeros(size(fits.chiSquares,[1,2]));
%     
%     slope = 0.0106; %determined from fitting Chi2 values of blank map
%     
%     for i = 2:size(GFs,2)
%         chi = fits.chiSquares(:,:,i);
%         chiSmaller = chi - (slope * GFs(i)) < fits.chiSquares(:,:,1);
%         
%         GF(chiSmaller) = GFs(i);
%     end
%     globalFraction = median(GF(GF > 0));
end