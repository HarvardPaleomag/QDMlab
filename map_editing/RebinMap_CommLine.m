function [] = RebinMap_CommLine(INFILE,rebin)
% This script takes an input B111 OR Bz map and bins with by a specified
% amount

clear B111ferro;
clear Bz;
clear ledImg;
clear newLED;

%Load the correct Bz file
%Load the correct B111 or Bz file
[filepath,name,ext]=fileparts(INFILE);
load(INFILE);

%Crop a single region for source fitting
if exist('B111ferro')
    if mod(size(B111ferro,1),rebin) ~= 0
        error('The dimensions of the input maps are not divisible by the rebinning size.');
    end
    if mod(size(B111ferro,2),rebin) ~= 0
        error('The dimensions of the input maps are not divisible by the rebinning size.');
    end
    
    %B111ferro
    temp=zeros(size(B111ferro,1)/rebin, size(B111ferro,2)/rebin);
    %i are rows; j are columns
    for i=1:size(temp,1)
        rowrange=[rebin*(i-1)+1 rebin*(i-1)+rebin];
        for j=1:size(temp,2)
            colrange=[rebin*(j-1)+1 rebin*(j-1)+rebin];
            temp(i,j)=mean(mean(B111ferro(rowrange(1):rowrange(1),colrange(1):colrange(1))));
        end
    end
    B111ferro=temp;
    %B111para
    temp=zeros(size(B111para,1)/rebin, size(B111para,2)/rebin);
    %i are rows; j are columns
    for i=1:size(temp,1)
        rowrange=[rebin*(i-1)+1 rebin*(i-1)+rebin];
        for j=1:size(temp,2)
            colrange=[rebin*(j-1)+1 rebin*(j-1)+rebin];
            temp(i,j)=mean(mean(B111para(rowrange(1):rowrange(1),colrange(1):colrange(1))));
        end
    end
    B111para=temp;
    
    %{
    %FitCvg
    temp=zeros(size(FitCvgBoth,1)/rebin, size(FitCvgBoth,2)/rebin);
    %i are rows; j are columns
    for i=1:size(temp,1)
        rowrange=[rebin*(i-1)+1 rebin*(i-1)+rebin];
        for j=1:size(temp,2)
            colrange=[rebin*(j-1)+1 rebin*(j-1)+rebin];
            temp(i,j)=ceil(mean(mean(FitCvgBoth(rowrange(1):rowrange(1),colrange(1):colrange(1)))));
        end
    end
    FitCvgBoth=temp;
    %negDiff
    temp=zeros(size(negDiff,1)/rebin, size(negDiff,2)/rebin);
    %i are rows; j are columns
    for i=1:size(temp,1)
        rowrange=[rebin*(i-1)+1 rebin*(i-1)+rebin];
        for j=1:size(temp,2)
            colrange=[rebin*(j-1)+1 rebin*(j-1)+rebin];
            temp(i,j)=mean(mean(negDiff(rowrange(1):rowrange(1),colrange(1):colrange(1))));
        end
    end
    negDiff=temp;
    %posDiff
    temp=zeros(size(posDiff,1)/rebin, size(posDiff,2)/rebin);
    %i are rows; j are columns
    for i=1:size(temp,1)
        rowrange=[rebin*(i-1)+1 rebin*(i-1)+rebin];
        for j=1:size(temp,2)
            colrange=[rebin*(j-1)+1 rebin*(j-1)+rebin];
            temp(i,j)=mean(mean(posDiff(rowrange(1):rowrange(1),colrange(1):colrange(1))));
        end
    end
    %}

    posDiff=temp;
    if exist('ledImg')
        if exist('corners')
            save([filepath '/B111_' num2str(rebin) 'x' num2str(rebin) 'Binned.mat'],'B111ferro','B111para','corners','ledImg');
        else
            save([filepath '/B111_' num2str(rebin) 'x' num2str(rebin) 'Binned.mat'],'B111ferro','B111para','ledImg');
        end
    else
        if exist('corners')
            save([filepath '/B111_' num2str(rebin) 'x' num2str(rebin) 'Binned.mat'],'B111ferro','B111para','corners');
        else
            save([filepath '/B111_' num2str(rebin) 'x' num2str(rebin) 'Binned.mat'],'B111ferro','B111para');
        end
    end
    
else
    if exist('Bz')
        if mod(size(Bz,1),rebin) ~= 0
            error('The dimensions of the input maps are not divisible by the rebinning size.');
        end
        if mod(size(Bz,2),rebin) ~= 0
            error('The dimensions of the input maps are not divisible by the rebinning size.');
        end
        
        %Bz
        temp=zeros(size(Bz,1)/rebin, size(Bz,2)/rebin);
        %i are rows; j are columns
        for i=1:size(temp,1)
            rowrange=[rebin*(i-1)+1 rebin*(i-1)+rebin];
            for j=1:size(temp,2)
                colrange=[rebin*(j-1)+1 rebin*(j-1)+rebin];
                temp(i,j)=mean(mean(Bz(rowrange(1):rowrange(1),colrange(1):colrange(1))));
            end
        end
        Bz=temp;
        %Bt
        temp=zeros(size(Bt,1)/rebin, size(Bt,2)/rebin);
        %i are rows; j are columns
        for i=1:size(temp,1)
            rowrange=[rebin*(i-1)+1 rebin*(i-1)+rebin];
            for j=1:size(temp,2)
                colrange=[rebin*(j-1)+1 rebin*(j-1)+rebin];
                temp(i,j)=mean(mean(Bt(rowrange(1):rowrange(1),colrange(1):colrange(1))));
            end
        end
        Bt=temp;
        
        if exist('newLED')
            if exist('corners')
                save([filepath '/Bz_' num2str(rebin) 'x' num2str(rebin) 'Binned.mat'],'Bz','Bt','h','step','corners','newLED');
            else
                save([filepath '/Bz_' num2str(rebin) 'x' num2str(rebin) 'Binned.mat'],'Bz','Bt','h','step','newLED');
            end
        else
            if exist('corners')
                save([filepath '/Bz_' num2str(rebin) 'x' num2str(rebin) 'Binned.mat'],'Bz','Bt','h','step','corners');
            else
                save([filepath '/Bz_' num2str(rebin) 'x' num2str(rebin) 'Binned.mat'],'Bz','Bt','h','step');
            end
        end
    else
        error('No file found');
    end
end

