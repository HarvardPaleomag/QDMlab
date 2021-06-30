function pixErr = detect_saturated(imgStack)

stdD = nanstd(imgStack);
pixErr = stdD < 1e-4; % find idx of pixels with std close to zero -> all freqs have same value

end