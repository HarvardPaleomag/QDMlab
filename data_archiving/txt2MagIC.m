function txt2MagIC(data_files,dest,contrib_id,experiment,specimen,txtkwargs,kwargs,magic)
% This function takes a results file (output of fit_source_series) and,
% optionally, a set of file paths and performs four (five) main tasks for 
% producing a MagIC dataset:
%   1) Writes fit information into a supplied MagIC measurements file
%   2) Converts Bz and B111 files into compact MagIC format and saves at a
%       specified location
%   3) Copies context (LED) image and saves at a specified location 
%   4) (optional) makes an image of the fitted(cropped) map and saves at 
%       a specified location 
%   5) 
%
% As a part of this, the script also makes headers for each B111/Bz image.
% The user should use the header keyword arguments to supply information
% for the file headers.
% 
% Parameters
% ------------
% data_files: cell array of strings
%   ordered list of file paths of folders containing experiment data
% results: struct
%   loaded results strucutre (output of dipole_fit_series) 
% dest: str
%   path to save MagIC files
% contrib_ID: integer
%   unique MagIC contribution identifier (obtained by initiating
%   contribution)
% experiment: str
%   experiment name
% specimen: str
%   unique identifier matching label in specimens table
% 
% Keyword arguments (kwargs)
% ------------
% makeImgs: bool [true]
%   whether to make new image files in MagIC format
% writeToMeas: bool [true]
%   whether to write measurement data into file
% writeToImages: bool [true]
%   whether to write image data into file
% saveContext: bool [true]
%   whether to make and save context (LED) images
% saveBzFitted: bool [false]
%   whether to make and save cropped Bz maps
% cax: ['max'] or float array
%   colorbar limits to use in plotting
% nPixels: int [4]
%   number of pixels used for binning when producing Bz map
% nameFunction: function [false]
%   function to produce measurement names from results info
% useTform: bool [false]
%   
% 
% MagIC file arguments (magic)
% ------------
%   NOTE: Unless otherwise defined, header arguments EXACTLY MATCH 
%   definitions in the MagIC data model (https://www2.earthref.org/MagIC/data-models/3.0).
%   Refer to the data model for definitions, units, and formatting requirements.
% ------------ TREATMENT PARAMETERS ------------
%   NOTE: for sequential experiments (e.g. AF demagnetization, anisotropy
%   measurements) the values of any properties relevant to the experiment 
%   should be supplied to the relevant keyword argument as a list of the 
%   same length as data_folders. 
% ------------
% treat_temp: list [false]
% treat_ac_field: list [false]
% treat_dc_field: list [false]
% treat_dc_field_phi: list [false]
% treat_dc_field_theta: list [false]
% ------------  GENERAL PARAMETERS  ------------
% GFFactor: float [0.25]
%   global fluorescence factor used to produce B111 file
% nameFunction: function
%   function to produce measurement names from results info
% standard: 's',['u']
% quality: ['g'],'b'
% method_codes: str ['LP-QDM']
%   method codes that describe measurement (https://www2.earthref.org/MagIC/method-codes)
% citations: str ['This study']
% citations_contrib: str ['This study']
%   use to include additional references for line in contribution file
% software_packages: str ['QDMLab']
% description: str ['bottom left is (0,0) and field of view is 226x141 microns']
% instrument_codes: str [false]
% meas_temp: float [291]
% analysts: str ['']
% dir_csd: float [0]
% magn_x_sigma: float [0]
% magn_y_sigma: float [0]
% magn_z_sigma: float [0]

arguments
    data_files;
    dest;
    contrib_id;
    experiment;
    specimen;
    
    % text kwargs
    txtkwargs.UC = 0;
    txtkwargs.fileName = 'Bz_uc0.mat';
    txtkwargs.tFormFile = false;
    
    % general kwargs
    kwargs.makeImgs = true;
    kwargs.writeMeas = true;
    kwargs.writeImg = true;
    kwargs.saveBzFitted = false;
    kwargs.saveContext = true;
    kwargs.cax = 'max';
    kwargs.nPixels = 4;
    kwargs.nameFunction = false;        % function to produce measurement names from results info
    kwargs.meas_n = false;
    kwargs.useTform = false;
    kwargs.doUC = true;
    kwargs.decAdjust = false;
    
    % MagIC kwargs
    magic.GFFactor = 0.25;          
    magic.standard = 'u';              
    magic.quality = 'g';
    magic.method_codes = 'LP-QDM';
    magic.citations = 'This study';
    magic.software_packages = 'QDMLab';
    magic.meas_orient_phi = false;
    magic.meas_duration = false;
    magic.description = 'Bottom left is (0,0) and field of view is 226x141 microns';
    magic.description_fit = 'Best-fit dipole';
    magic.citations_meas = 'This study';
    magic.instrument_codes = '';
    magic.meas_temp = 291;
    magic.analysts = '';
    
    % treatment variables
    magic.treat_temp = false;
    magic.treat_ac_field = false;
    magic.treat_dc_field = false;
    magic.treat_dc_field_phi = false;
    magic.treat_dc_field_theta = false;
    
    magic.dir_csd = 0;
    magic.magn_x_sigma = 0;
    magic.magn_y_sigma = 0;
    magic.magn_z_sigma = 0;
end

%% setup
% get pixel size
pixelsize = kwargs.nPixels/4 * 4.7*10^(-6);

% make file paths for measurements and images files
meas_fp = [dest,'/magic_measurements_',num2str(contrib_id),'.txt'];
img_fp = [dest,'/magic_images_',num2str(contrib_id),'.txt'];

dest = strrep(dest,'\','/');

%% build temporary results structure

% initialize arrays
preAllocatedArray = zeros([1,length(data_files)]);
moments = preAllocatedArray;
inclinations = preAllocatedArray;
declinations = preAllocatedArray;
heights = preAllocatedArray;
datas = cell(size(preAllocatedArray));
models = cell(size(preAllocatedArray));
residuals = preAllocatedArray;
dipolarity = preAllocatedArray;
xMin = preAllocatedArray;
xMax = preAllocatedArray;
yMin = preAllocatedArray;
yMax = preAllocatedArray;
xloc = preAllocatedArray;
yloc = preAllocatedArray;
fileResults = num2cell(preAllocatedArray);
% get results structure from DipoleInversions files
txtkw_cell = namedargs2cell(txtkwargs);
results = txt2results(data_files,specimen,txtkw_cell{:});

% call results2MagIC
kw_cell = namedargs2cell(kwargs);
mg_cell = namedargs2cell(magic);
results2MagIC(strcat(data_files,txtkwargs.fileName),results,dest,contrib_id,experiment,specimen,kw_cell{:},mg_cell{:})


%% 
% 
% % crawl through file paths and load dipole inversions data
% for i = 1:length(data_files)
%     % get file parts
%     filei = data_files(i);
%     [fpi,namei,exti] = fileparts(filei);
%     
%     % open DipoleInversions.txt
%     try
%         fID = fopen([fpi, '/DipoleInversions.txt'],'r');
%     catch
%         error('DipoleInversions:MissingFile',...
%                     'Error: %s not found.',[fpi, '/DipoleInversions.txt'])
%     end
%     
%     % reshape dipole inversions data
%     temp = textscan(fID,'%s','Delimiter','\t');
%     temp = temp{1};
%     temp = reshape(temp,8,[]); % FIX THE TOTAL # AFTER DEFINING FINAL TEXT OUTPUT FORMAT
%     
%     % search for desired specimen
%     try
%         spj = find(contains(temp(1,:),specimen));
%     catch
%         error('DipoleInversions:MissingSpecimen',...
%                     'Error: %s is missing from the supplied inversions file.',specimen)
%     end
%     
%     preAllocatedArray = zeros([1,length(data_files)]);
%     moments(i) = temp(2,spj);
%     inclinations(i) = temp(3,spj);
%     declinations(i) = temp(4,spj);
%     heights(i) = preAllocatedArray;
%     residuals(i) = preAllocatedArray;
%     dipolarity(i) = preAllocatedArray;
%     xMin(i) = temp(7,spj);
%     xMax(i) = temp(8,spj);
%     yMin(i) = temp(9,spj);
%     yMax(i) = temp(10,spj);
%     xloc(i) = temp(5,spj);
%     yloc(i) = temp(6,spj);
%     
%     
%     keyboard;
% end