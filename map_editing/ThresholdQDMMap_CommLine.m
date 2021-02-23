function [] = ThresholdQDMMap_CommLine(INFILE,high,low)


[filepath,name,ext]=fileparts(INFILE);
disp(filepath)
load(INFILE);

%the upper and lower bounds
ub=high;
lb=low;



[a b] = size(B111ferro);
Nreplacements=0;
for i = 1:a
    for j = 1:b
        if (B111ferro(i,j) < lb) | (B111ferro(i,j) > ub)
            B111ferro(i,j) = 0;
            Nreplacements = Nreplacements + 1;
        end
    end
end

disp([num2str(Nreplacements) ' entries in B111ferro replaced by 0'])

if exist('ledImg','var')
        save([filepath '/' 'B111dataToPlot_truncated' '.mat'],'B111ferro','B111para','ledImg','negDiff','posDiff');
    else
        save([filepath '/' 'B111dataToPlot_truncated' '.mat'],'B111ferro','B111para');
end
