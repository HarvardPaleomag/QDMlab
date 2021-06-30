function binnedImage = BinImage(Image,binSize)

if binSize == 1
    BinnedImage = Image;
    return
end

% crop to integer number of bins
Image = Image(1:end-mod(end,binSize),1:end-mod(end,binSize));
binnedImage = imresize(Image, 1/binSize, 'method', 'box');

end