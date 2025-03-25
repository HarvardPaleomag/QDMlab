function results = make_nROI_struct(infile)

% takes the output of a fit_sources_series() run and makes a new nROI
% struct that can be used as input to another run 

xMins = infile.xMin(:,1);
xMaxs = infile.xMax(:,1);
yMins = infile.yMin(:,1);
yMaxs = infile.yMax(:,1);

cellArray = cell(1, length(xMins));

for i = 1:length(xMins)
    outmatrix = [xMins(i),yMins(i),xMaxs(i)-xMins(i),yMaxs(i)-yMins(i)];
    cellArray{1,i} = outmatrix;
end

results=cellArray;

