%% DEMO
% This demo will go through the whole process of fitting a source at
% each step of an experiment and then preparing that data to be added to
% MagIC. In this case, the data is a subset of an AF demagnetization
% of a 2G ARM applied to a slice of ALH 84001, and we fit the rightmost 
% source in the field of view.

%% SETUP
stem = ['data_txt' filesep];
paths = {[stem,'AF0G' filesep],[stem, 'AF50G_rep1' filesep],[stem,'AF50G_rep2' filesep],...
    [stem,'AF75G' filesep],[stem,'AF100G' filesep]};
Bz_file = 'Bz_uc0_uc50.mat';

%% FIT SAME SOURCE AT SEVERAL AF DEMAG STEPS
% Note that the first argument here is the source ID (what will be a single
% specimen in MagIC). This value should be the same across the entire
% experiment for a given source.

% Run this and pick the approximate location of the source in each image.

for i = 1:length(paths)
    fit_source('T1','filePath',[paths{i} Bz_file],'fitOrder',1,'cropFactor',60);
    
    close all
end

%% (OPTIONAL!) MAKE RESULTS STRUCTURE FROM DIPOLEINVERSIONS.TXT FILES
% Note: running this as-is will give the prompt 'UC = 50? [Y]/N'. This is
% the correct UC value in this case, and was estimated from the 'uc_'
% suffixes in the file name. When in doubt it is better to supply a known
% upward continuation value in microns using the 'UC' keyword.
res_tmp = txt2results(paths,'T1','fileName',Bz_file);


%% MAGIC CONVERSION

% Now, we can take the fit data we just generated and get it ready to put
% into MagIC. 

% First, we make a list of the Bz files of interest (including the Bz file
% name). Maybe here we choose to only use some of the steps in the results
% file we made (in this case, numbering will completely ignore any skipped
% steps):
data_files = paths([1,2,4]);

% We also need to supply a contribution ID (typically acquired when the
% contribution is initiated in MagIC, but this one is fake)
contrib_ID = 55555;

% Next, we set the name of the experiment to use and the name of the 
% specimen. In this case, each line describing steps of this experiment in 
% the MagIC contribution file will name the experiment name
% ‘QDM_5a_ARM-Demag-AF_##’ where ## is the step number within the
% experiment. 
experiment = 'QDM-ARM-Demag-AF'; % name of the experiment
specimen = 'T1'; % specimen name
UC = 50; % the upward continuation amount in microns (the code will attempt to guess this if it isn't supplied) 

% We also make a structure that includes the details of the experiment
% we’re adding. This includes relevant method codes and lists describing
% the details of any treatments applied.
magic.method_codes = 'LP-QDM:LP-ARM-AFD'; % method codes for a QDM measurement and AF demag of ARM
magic.treat_ac_field = [3000,25,50,75,100]/10000; % AF applied during step
magic.treat_dc_field = [2,0,0,0,0]/10000; % bias field
magic.treat_dc_field_phi = [0,0,0,0,0]; % bias field direction
magic.treat_dc_field_theta = [0,0,0,0,0]; % bias field direction

% We also add the DOI of the person who performed the experiments and
% fitting for this data:
magic.analysts = '0000-0002-7149-2693'; % DOI of the person who performed the analysis
magic.instrument = 'QDM-Antares'; % name of the QDM used to make the measurements

% It's easiest to convert this structure to a cell to supply keyword
% arguments :)
arg_array = namedargs2cell(magic);

% Now we can run results2MagIC to get everything ready for a MagIC upload!
% Here we set the keyword arguments to do everything: write the fit data to
% the contribution file (writeToContrib), make the condensed MagIC files 
% for the Bz and B111 images (makeImgs), make and save the cropped Bz map
% (saveBzFitted), and make and save the context images (saveContext). 

txt2MagIC(data_files,[stem,'demo_MagIC'],contrib_ID,...
    experiment,specimen,'UC',UC,'writeMeas',true,'writeImg',true,'fileName',Bz_file,...
    'makeImgs',true,'saveBzFitted',true,'saveContext',true,'doUC',false,arg_array{:})

% A few notes about arguments here: in this demo, we are fitting data from
% a map that has already been upward-continued rather than cropping the
% source out and upward-continuing later. I therefore specified the UC
% value here (although the code would make a guess and check with me if I hadn't
% done so) and set 'doUC' to false (this prevents the code from duplicating 
% upward-continuation while plotting the BzFitted images).


