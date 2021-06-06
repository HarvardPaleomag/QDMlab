classdef qdmmap
    
    properties 
       dataFolder
       binSize
       map
    end
    
    methods
        function obj = qdmmap(dataFolder, binSize)
            map = ODMR_to_B111(dataFolder, binSizes)
        end
    end
end