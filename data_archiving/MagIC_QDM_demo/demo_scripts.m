%% DEMO
% This demo will go through the whole process of fitting a source at
% each step of an experiment and then preparing that data to be added to
% MagIC. In this case, the data is a subset of an AF demagnetization
% of a 2G ARM applied to a slice of ALH 84001, and we fit the rightmost 
% source in the field of view.

%% SETUP
stem = 'data\';
paths = {[stem,'ARM2G\'],[stem, 'AF25G\'],[stem,'AF50G\'],[stem,'AF75G\'],[stem,'AF100G\']};
Bz_file = 'Bz_uc0.mat';

%% FITTING (OPTIONAL)
% Optionally, we can begin by producing the results file that we will
% eventually put into MagIC format in the next cell. 

% set transform reference file
tf_file = [stem,'demo_transform.mat'];

% set region of interest for fitting
nROI_1 = {[337, 163, 115, 105]};

% do fitting
results = fit_sources_series(paths, 'upCont',50,'transFormFile',tf_file,...
    'fileName',Bz_file,'nROI',nROI_1);

% save results
save([stem,'demo_results.mat'],'results')


%% MAGIC CONVERSION

% Now, we can take the fit data we just generated and get it ready to put
% into MagIC. 

% First, we make a list of the Bz files of interest (including the Bz file
% name). Maybe here we choose to only use some of the steps in the results
% file we made (in this case, numbering will completely ignore any skipped
% steps):
data_files = strcat(paths([1,2,4]),Bz_file);

% We also need to supply a contribution ID (typically acquired when the
% contribution is initiated in MagIC, but this one is fake)
contrib_ID = 55555;

% We can then load the results file (since typically we are not making and
% archiving the fit data simultaneously)
load([stem,'demo_results.mat'],'results')

% Next, we set the name of the experiment to use and the name of the 
% specimen. In this case, each line describing steps of this experiment in 
% the MagIC contribution file will name the experiment name
% ‘QDM_5a_ARM-Demag-AF_##’ where ## is the step number within the
% experiment. 
experiment = 'QDM-ARM-Demag-AF'; % name of the experiment
specimen = '5a'; % specimen name
magic.meas_orient_phi = 286.2; % orientation of transformation reference frame


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

results2MagIC(data_files,results,[stem,'demo_MagIC'],contrib_ID,...
    experiment,specimen,'useTform',true,'writeMeas',true,'writeImg',true,...
    'makeImgs',true,'saveBzFitted',true,'saveContext',true,arg_array{:})



