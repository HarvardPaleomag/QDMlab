
%% DEMO
% This demo will go through the process of converting QDM image data into
% MagIC format for archiving.

%% SINGLE FILE SETUP

% Here we'll use the first step of an AF sequence as an example:
stem = ['data_ims' filesep];
path = [stem 'AF0G' filesep];
Bz_file = 'Bz_uc0.mat';
Bz_uc_file = 'Bz_uc0_uc50.mat';

% Let's load this data and see what it looks like:
Bz_data = load([path, Bz_file]);
QDM_figure(Bz_data.newLED,'led','true')
QDM_figure(Bz_data.Bz, 'colormap','turbo')

% And the upward continued data:
Bz_uc_data = load([path, Bz_uc_file]);
QDM_figure(Bz_uc_data.Bz, 'colormap','turbo')

%% MAGIC CONVERSION FOR ONE FILE
% The first thing we need to supply is a contribution ID (typically 
% acquired when the contribution is initiated in MagIC--this one is fake)
contrib_ID = 55555;

% Next, we set the name of the experiment to use and the name of the 
% specimen. In this case, each line describing steps of this experiment in 
% the MagIC contribution file will name the experiment name
% ‘QDM_5a_ARM-Demag-AF_##’ where ## is the step number within the
% experiment. 
experiment = 'QDM-ARM2G'; % name of the experiment
specimen = 'T1'; % specimen name
UC = 0; % the upward continuation amount in microns (the code will attempt to guess this if it isn't supplied) 

% We also make a structure that includes the details of the experiment
% we’re adding. This includes relevant method codes and lists describing
% the details of any treatments applied.
magic.method_codes = 'LP-QDM:LP-ARM'; % method codes for a QDM measurement and AF demag of ARM
magic.treat_ac_field = 3000/10000; % AF applied during step
magic.treat_dc_field = 2/10000; % bias field
magic.treat_dc_field_phi = 0; % bias field direction
magic.treat_dc_field_theta = 0; % bias field direction

% We also add the DOI of the person who performed the experiments and
% fitting for this data:
magic.analysts = '0000-0002-7149-2693'; % ORCID of the person who performed the analysis
magic.instrument = 'QDM-Antares'; % name of the QDM used to make the measurements

% It's easiest to convert this structure to a cell to supply keyword
% arguments :)
arg_array = namedargs2cell(magic);

% Now we can run Bfile2MagIC to get everything ready for a MagIC upload!
% Here we set the keyword arguments to do everything: write the fit data to
% the contribution file (writeToContrib), make the condensed MagIC files 
% for the Bz and B111 images (makeImgs), make and save the cropped Bz map
% (saveBzFitted), and make and save the context images (saveContext). 

Bfile2MagIC(path,[stem,'demo_MagIC' filesep],contrib_ID,...
    experiment,specimen,'UC',UC,'writeImg',true,'fileName',Bz_file,...
    'makeImgs',true,'saveContext',true,arg_array{:})

fclose('all');

%% FULL EXPERIMENT SETUP

stem = ['data_ims' filesep];
paths = {[stem,'AF0G' filesep],[stem, 'AF50G_rep1' filesep],[stem,'AF50G_rep2' filesep],...
    [stem,'AF75G' filesep],[stem,'AF100G' filesep]};
Bz_file = 'Bz_uc0_uc50.mat';

%% MAGIC CONVERSION

% Now, we want to put these images into MagIC. 

% The first thing we need to supply is a contribution ID (typically 
% acquired when the contribution is initiated in MagIC--this one is fake)
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
magic.analysts = '0000-0002-7149-2693'; % ORCID of the person who performed the analysis
magic.instrument = 'QDM-Antares'; % name of the QDM used to make the measurements

% It's easiest to convert this structure to a cell to supply keyword
% arguments :)
arg_array = namedargs2cell(magic);

% Now we can run results2MagIC to get everything ready for a MagIC upload!
% Here we set the keyword arguments to do everything: write the fit data to
% the contribution file (writeToContrib), make the condensed MagIC files 
% for the Bz and B111 images (makeImgs), make and save the cropped Bz map
% (saveBzFitted), and make and save the context images (saveContext). 

Bfile2MagIC(paths,[stem,'demo_MagIC' filesep],contrib_ID,...
    experiment,specimen,'UC',UC,'writeImg',true,'fileName',Bz_file,...
    'makeImgs',true,'saveContext',true,arg_array{:})

% A few notes about arguments here: in this demo, we are fitting data from
% a map that has already been upward-continued rather than cropping the
% source out and upward-continuing later. I therefore specified the UC
% value here (although the code would make a guess and check with me if I hadn't
% done so) and set 'doUC' to false (this prevents the code from duplicating 
% upward-continuation while plotting the BzFitted images).


