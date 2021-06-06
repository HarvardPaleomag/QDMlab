classdef qdmmap
    
    properties 
       dataFolder
       binSize
       calculatedBins
       b111ferro
       b111para
    end
    
    methods
        function obj = qdmmap(dataFolder, binSizes)
            maps = ODMR_to_B111('nFolders', dataFolder,'binSizes', binSizes, 'save', 0);
            obj.calculatedBins = binSizes;
            close all
            for n = 1:size(binSizes,2)
                obj.b111ferro.(sprintf('bin%i')) = maps{1,n}.('B111ferro');
                obj.b111para.(sprintf('bin%i')) = maps{1,n}.('B111para');
            end
        end
    end
end

function test()
clc
close all
% d = ODMR_to_B111('nFolders', 'D:\data\mike\NRM', 'binSizes',[16, 32], 'save', false);
% close all
map = qdmmap('D:\data\mike\NRM', [16,32,64]);
end