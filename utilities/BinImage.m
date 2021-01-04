function BinnedImage = BinImage(Image,BinSize)

% crop to integer number of bins
Image = Image(1:end-mod(end,BinSize),1:end-mod(end,BinSize));

sizeX = size(Image,2);
sizeY = size(Image,1);

BinnedYImage = squeeze(...
    mean(reshape(Image,BinSize,sizeY/BinSize,sizeX),1));

BinnedImage = squeeze(...
    mean(reshape(BinnedYImage',BinSize,sizeX/BinSize,sizeY/BinSize),1))';

end