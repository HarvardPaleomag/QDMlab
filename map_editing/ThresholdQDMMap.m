[longfilename, pathname] = uigetfile('*.mat', 'Pick a magnetic field map file');
fullfilename=[pathname longfilename];
parsedpath=strsplit(pathname,filesep);
savefilename=parsedpath(end);
savefilename=char(savefilename);

load(fullfilename); %load the B111 variables

minfield = min(min(B111ferro));
maxfield = max(max(B111ferro));

disp(['Min and Max values of the map are:  ' num2str(minfield) ' , ' num2str(maxfield)])

%asks user for the upper and lower bounds
ub=input('Upper truncation value [0.3 G]: ');
if isempty(ub)
    ub=0.3;
end

lb=input('Lower truncation value [-0.3 G]: ');
if isempty(lb)
    lb=-0.3;
end

paramean=mean(mean(B111para));
[a b] = size(B111ferro);
Nreplacements=0;
for i = 1:a
    for j = 1:b
        if (B111ferro(i,j) < lb) | (B111ferro(i,j) > ub)
            B111ferro(i,j) = 0;
            Nreplacements = Nreplacements + 1;
        end
        if ((B111para(i,j)-paramean) < lb) | ((B111para(i,j)-paramean) > ub)
            B111para(i,j) = 0;
            Nreplacements = Nreplacements + 1;
        end
    end
end

disp([num2str(Nreplacements) ' entries in B111ferro and B111para replaced by 0'])

save([pathname longfilename(1:end-4) '_truncated.mat'],'B111ferro', 'B111para', 'FitCvgBoth', 'ledImg', 'negDiff', 'posDiff', '-mat');
disp(sprintf(['\nSaving ' pathname longfilename(1:end-4) '_truncated.mat...\n']))
