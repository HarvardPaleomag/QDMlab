function results2MagIC(data_files,results,dest,contrib_id,experiment,specimen,kwargs,magic)
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
% 


arguments
    data_files;
    results;
    dest;
    contrib_id;
    experiment;
    specimen;
    
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

useTform = kwargs.useTform;

%% setup
% get pixel size
pixelsize = kwargs.nPixels/4 * 4.7*10^(-6);

% make file paths for measurements and images files
meas_fp = [dest,'/magic_measurements_',num2str(contrib_id),'.txt'];
img_fp = [dest,'/magic_images_',num2str(contrib_id),'.txt'];

dest = strrep(dest,'\',filesep);

%% prepare data

% fix files field name for legacy data
if find(contains(fields(results),'nFiles'))
    results.files = results.nFiles; 
end  

if isstruct(results.UC)
    results.UC = results.UC{1};
end


% check if supplied data_files are in results
for data_file = data_files
    if ~any(strcmp(results.files,data_file))
        error('DataInput:MissingFile',...
            'Error: %s is missing from the supplied results structure.',data_file)
    end
end

% check if orientation data is supplied
if isnumeric(magic.meas_orient_phi)
    % if results contains decAdjust, check that the two values are the same
    if any(contains(fields(results),'decAdjust'))
        if mod((360+results.decAdjust),360) ~= mod((360+magic.meas_orient_phi),360)
            error('DataInput:OrientationMismatch',...
            'Error: Conflicting orientation data supplied (results.decAdjust ~= meas_orient_phi).')
        end
    
    % if it doesn't, add meas_orient_phi temporarily to results structure
    else
        results.decAdjust = mod((360+magic.meas_orient_phi),360);
    end
% if absolutely no orientation data is supplied, set decAdjust to 0
elseif ~any(contains(fields(results),'decAdjust'))
    results.decAdjust = 0;
end

% check transform data
% see if tform data is saved in results structure
if useTform
    if ~any(contains(fields(results),'tforms'))
        % if it isn't, see if user supplied a transform file path
        if isbool(kwargs.tform_fp)
            % if not, ask user to supply transform file location
            [longfilename, pathname] = uigetfile('*.mat', 'Pick a transform file');
            kwargs.tform_fp=[pathname longfilename]; 
        end

        % try to load supplied transform file
        try
            tf = load(kwargs.tform_fp);
        catch
            error('Transform:TformFileNotFound','Error: Could not load transform file.')
        end

        % check if all data_files are in transform file
        for i = 1:length(data_files)
            data_file = data_files{i};
            if ~any(strcmp(keys(tf.nRefFrames),data_file))
                error('Transform:MissingFile',...
                    'Error: %s is missing from the supplied transform structure.',data_file)
            end
        end

    end
end


%% prepare output 
% fix dest if necessary
if ~strcmp(dest(end),'/') || ~strcmp(dest(end),'\')
    dest = [dest filesep];
end




meas_n = 0;
% make output file if it doesn't already exist
if ~isfolder(dest)
    mkdir(dest)
    mkdir([dest 'images' filesep])
elseif kwargs.meas_n > 0
    meas_n = kwargs.meas_n;
else
    % get next starting number for magic_measurements_*.## file naming
    dest_files = dir([dest  'magic_measurements_*']);
    meas_n = meas_n + length(dest_files);
end


% check whether we'll need to write column names for measurements table
writeColNames = 0;
if kwargs.writeMeas
    % open MagIC contribution file
    fileIDMeas = fopen(meas_fp,'a+');
    
    frewind(fileIDMeas)
    temp = textscan(fileIDMeas,'%s','Delimiter','\n');
    if isempty(temp{1})
        writeColNames = 1;
    elseif contains(temp{1}{end},'>>>>>>>>>>')
        writeColNames = 1;
    end
    clear temp;
end

% write column names if necessary
if writeColNames
    colNames = ['tab delimited\tmeasurements\n',...
        'measurement\texperiment\tspecimen\tfiles\tdir_csd\tdir_dec\t'...
        'dir_inc\ttreat_ac_field\ttreat_dc_field\ttreat_dc_field_phi\t',...
        'treat_dc_field_theta\ttreat_step_num\ttreat_temp\tcitations\t',...
        'instrument_codes\tmethod_codes\tquality\tstandard\t',...
        'meas_orient_phi\tmeas_temp\tmeas_pos_x\tmeas_pos_y\t',...
        'inversion_height\tinversion_residuals\tanalysts\tdescription\t',...
        'software_packages\tmagn_moment\tmagn_x_sigma\tmagn_y_sigma\tmagn_z_sigma\n'];
    
    fprintf(fileIDMeas,colNames);
end


% check whether we'll need to write column names for images table
writeColNames = 0;
if kwargs.writeImg
    % open MagIC contribution file
    fileIDImg = fopen(img_fp,'a+');
    
    frewind(fileIDMeas)
    temp = textscan(fileIDImg,'%s','Delimiter','\n');
    if isempty(temp{1})
        writeColNames = 1;
    elseif contains(temp{1}{end},'>>>>>>>>>>')
        writeColNames = 1;
    end
    clear temp;
end

% write column names if necessary
if writeColNames
    colNames = ['tab delimited\timages\n',...
        'specimen\tfile\ttype\ttitle\tkeywords\n'];
    
    fprintf(fileIDImg,colNames);
end


% get names for different context/crop files (Bz/B111 files follow standard
% naming format)
if ~kwargs.nameFunction
    names = cell(1,length(data_files));
    
    for s = 1:length(data_files)
        names(s) = {[specimen, '_',experiment, '_',num2str(s)]};
    end
else
    names = kwargs.nameFunction(results);
end


%% make images
for j = 1:length(data_files)
    file = data_files{j}; 
    
    % find position in array 
    i = find(strcmp(file,results.files));

    [stepPath,BzName,~] = fileparts(file);
    
    experiment_i = [experiment,'-',num2str(j)];

    % get image file names
    B111FName = [dest,'magic_measurements_',num2str(contrib_id),'.',num2str(meas_n+j),'.txt'];
    BzFName = [dest,'magic_measurements_',num2str(contrib_id),'.',num2str(meas_n+j),'.txt'];
    contextFName = [dest,'images' filesep 'QDM_context_image_',names{j},'.jpg'];
    fittedMapFName = [dest,'images' filesep 'QDM_fitted_map_', names{j},'.jpg'];
    
    
    % drop properties that are not used (and assemble header info while
    % we're at it)
    deetsStr = '';
    if ((length(magic.treat_dc_field)==1) && (magic.treat_dc_field == false)) || ((length(magic.treat_dc_field)>1) && (magic.treat_dc_field(i) == false))
        treat_dc_field_i = '';
        treat_dc_field_phi_i = '';
        treat_dc_field_theta_i = '';
    else
        treat_dc_field_i = num2str(magic.treat_dc_field(i),'%0.6e');
        treat_dc_field_phi_i = num2str(magic.treat_dc_field_phi(i));
        treat_dc_field_theta_i = num2str(magic.treat_dc_field_theta(i));
        
        deetsStr = [deetsStr,'\r\n* treat_dc_field\t',treat_dc_field_i,...
            '\r\n* treat_dc_field_phi\t',treat_dc_field_phi_i,...
            '\r\n* treat_dc_field_theta\t',treat_dc_field_theta_i];
    end
    
    if ((length(magic.treat_ac_field)==1) && (magic.treat_ac_field == false)) || ((length(magic.treat_ac_field)>1) && (magic.treat_ac_field(i) == false))
        treat_ac_field_i = '';
    else
        treat_ac_field_i = num2str(magic.treat_ac_field(i),'%0.6e');
        
        deetsStr = [deetsStr,'\r\n* treat_ac_field\t',treat_ac_field_i];
    end
    
    if ((length(magic.treat_temp)==1) && (magic.treat_temp == false)) || ((length(magic.treat_temp)>1) && (magic.treat_temp(i) == false))
        treat_temp_i = '';
    else
        treat_temp_i = num2str(magic.treat_temp(i));
        deetsStr = [deetsStr,'\r\n* treat_temp\t',treat_temp_i];
    end
    
    
    if kwargs.makeImgs
        
        
        
        %% b111
        % make derived_val string
        derivedVal = containers.Map({'QDM_GF_factor'},...  % keys
                            {[{num2str(magic.GFFactor,'%0.2f')},...     % values
                            {'10.1029/2020GC009147'}]});           % references

        derivedValStr = '';
        
        k=0;
        for keyi = keys(derivedVal)
            valRef = derivedVal(keyi{1});
            if k == 0
                derivedValStr = append(derivedValStr,[keyi{1},':',valRef{1},':',valRef{2}]);
            else 
                derivedValStr = append(derivedValStr,['; ',keyi{1},':',valRef{1},':',valRef{2}]);
            end
            k = k + 1;
        end
        
        % build header
        meas_header = ['tab\tmeasurements\r\n',...
            '* experiment\t',experiment_i,...
            '\r\n* specimen\t',specimen,...
            '\r\n* standard\t',magic.standard,...
            '\r\n* quality\t',magic.quality,...
            '\r\n* method_codes\t',	magic.method_codes,...
            '\r\n* citations\t',magic.citations,...
            '\r\n* description\t','B111 image of ' experiment_i,'; ',magic.description,... 
            '\r\n* derived_value\t',derivedValStr,...
            deetsStr,...
            '\r\nmeasurement\tmagn_z\tmeas_pos_x\tmeas_pos_y'
            ];

        % write header
        fileIDB111 = fopen(B111FName,'w');
        fprintf(fileIDB111, meas_header);
        
        % load B111 
        load([stepPath filesep 'B111dataToPlot.mat'],'B111ferro');
        
        % convert to A/m (B111 is in G initially)
        B111ferro = B111ferro * 79.5774715;
        
        % write lines to condensed MagIC file
        Bfile2MagIC(B111ferro,fileIDB111,num2str(meas_n+i),pixelsize)
        
        %% Bz
        derivedVal = containers.Map({'QDM_GF_factor'},...  % keys
                            {[{num2str(magic.GFFactor,'%0.2f')},...   % values
                            {'10.1029/2020GC009147'}]});      % references

        derivedValStr = '';
        k=0;
        for keyi = keys(derivedVal)
            valRef = derivedVal(keyi{1});
            if k == 0
                derivedValStr = append(derivedValStr,[keyi{1},',',valRef{1},',',valRef{2}]);
            else 
                derivedValStr = append(derivedValStr,['; ',keyi{1},',',valRef{1},',',valRef{2}]);
            end
            k = k + 1;
        end
        
        % build header
        meas_header = ['tab\tmeasurements\r\n',...
            '* experiment\t',experiment_i,...
            '\r\n* specimen\t',specimen,...
            '\r\n* standard\t',magic.standard,...
            '\r\n* quality\t',magic.quality,...
            '\r\n* method_codes\t',	magic.method_codes,...
            '\r\n* citations\t',magic.citations,...
            '\r\n* description\t','Bz image of ' experiment_i,'; ', magic.description,... 
            '\r\n* derived_value\t',derivedValStr,...
            deetsStr,...
            '\r\nmeasurement\tmagn_z\tmeas_pos_x\tmeas_pos_y'
            ];

        % write header
        fileIDBz = fopen(BzFName,'w');
        fprintf(fileIDBz, meas_header);
        
        % load Bz and write lines to file
        load(file,'Bz');
        
        % convert to A/m (Bz is in T initially)
        Bz = Bz * 795774.715;
        
        % write lines to condensed MagIC file
        Bfile2MagIC(Bz,fileIDBz,num2str(meas_n+i),pixelsize)
        
        %% LED and QDM_fitted_map
        % load Bz and write lines to file
        if ~exist([dest,'images'],'dir')
            mkdir([dest,'images'])
        end
        
        % load context (LED) image data
        try
            LED = load(file,'newLED').newLED;
        catch
            load(file,'LED').LED;
        end
        
        % make and save LED image
        if kwargs.saveContext
            figure();
            imshow(LED./255);
            axis xy, axis equal, axis tight, axis off
            saveas(gcf, contextFName)
            close all
        end

        % make and save crops
        if kwargs.saveBzFitted 
            % get transformed FOV
            if useTform
                try
                    BzTemp = tform_data(Bz, results.tforms{i}, results.refFrame);
                catch
                    BzTemp = tform_data(Bz,tf.nTransForms{file},tf.nRefFrames{file});
                end
            else
                BzTemp = Bz;
            end
            
            if kwargs.doUC
                BzTemp = UpCont(BzTemp, results.UC, 1/kwargs.nPixels);
            end
            
            BzTemp = BzTemp(results.yMin(i):results.yMax(i),...
                results.xMin(i):results.xMax(i));
            
            if strcmp(kwargs.cax,'max')
                kwargs.cax = max(abs(BzTemp),[],'all');
            end
            
            % make and save image of crop
            figure();
            imagesc(BzTemp);
            caxis([-1,1].*kwargs.cax);
            colormap(turbo)
            colorbar;
            axis xy, axis equal, axis tight, axis off
            
            saveas(gcf, fittedMapFName)
            close all
            
        end
    end
    
     % make files names and image strings for measurements/images tables
    if kwargs.writeMeas || kwargs.writeImg
        BzContextFName = ['QDM_context_image_',names{j},'.jpg'];
        fittedMapFName = ['QDM_fitted_map_',names{j},'.jpg'];

        if isfile([dest,'images/',fittedMapFName])
            imagesStr = ['qdm_context_image[',BzContextFName,']:',...
                    'qdm_b111_map[',B111FName(length(dest)+1:end),']:',...
                    'qdm_bz_map[',BzFName(length(dest)+1:end),']:',...
                    'qdm_fitted_map[',fittedMapFName,']'];
        else
            imagesStr = ['qdm_context_image[',BzContextFName,']:',...
                    'qdm_b111_map[',B111FName(length(dest)+1:end),']:',...
                    'qdm_bz_map[',BzFName(length(dest)+1:end),']'];
        end
        
        % replace any backslashes before writing to file
        imagesStr = strrep(imagesStr,'\','/');
    end

%% make lines in MagIC measurements table
    if kwargs.writeMeas
        % load Bz image
        load(file,'Bz');
        
        if any(strcmp(fields(results),'x'))
            x_tmp = results.x;
            y_tmp = results.y;
        else
            x_tmp = (results.xMax(i) + results.xMin(i))/2;
            y_tmp = (results.yMax(i) + results.yMin(i))/2;
        end
        
        % get transformed FOV
        if useTform
            try
                [xCorr,yCorr] = transformPointsInverse(results.tforms{i},x_tmp,y_tmp);
            catch
               [xCorr,yCorr] = transformPointsInverse(tf.nTransForms{file},x_tmp,y_tmp);
            end
        else
            xCorr = x_tmp;
            yCorr = y_tmp;
        end
        
        % write line in contrib file
        contribLine = [names{j},'\t',...
                        experiment_i,'\t',...
                        specimen,'\t',... 
                        imagesStr,'\t',...
                        num2str(magic.dir_csd),'\t',...
                        num2str(results.decs(i),'%0.1f'),'\t',...
                        num2str(results.incs(i),'%0.1f'),'\t',...
                        treat_ac_field_i,'\t',...            
                        treat_dc_field_i,'\t',...
                        treat_dc_field_phi_i,'\t',...
                        treat_dc_field_theta_i,'\t',...
                        num2str(j),'\t',...
                        treat_temp_i,'\t',...
                        magic.citations_meas,'\t',...
                        magic.instrument_codes,'\t',...
                        magic.method_codes,'\t',...
                        magic.quality,'\t',...
                        magic.standard,'\t',...
                        num2str(mod((360+results.decAdjust),360),'%0.1f'),'\t',...
                        num2str(magic.meas_temp),'\t',...
                        num2str(xCorr*pixelsize,'%0.2f'),'\t',...
                        num2str(yCorr*pixelsize,'%0.2f'),'\t',...
                        num2str(abs(results.heights(i)),'%0.14e'),'\t',...
                        num2str(results.res(i),'%0.2f'),'\t',...
                        magic.analysts,'\t',...
                        magic.description_fit,'\t',...
                        magic.software_packages,'\t',...
                        num2str(results.moments(i),'%0.14e'),'\t',...
                        num2str(magic.magn_x_sigma),'\t',...
                        num2str(magic.magn_y_sigma),'\t',...
                        num2str(magic.magn_z_sigma),'\r\n'];
        
        % print line to file
        fprintf(fileIDMeas, contribLine);
        
    end
    
    % write line to images table
    if kwargs.writeImg
        %specimen\tfile\ttype\ttitle\tkeywords\n fittedMapFName
        % write context image line to images table
        contribLine = [specimen,'\t',... 
                        'QDM_context_image_',names{j},'.jpg','\t',...
                        'Context image','\t',...
                        'QDM context image of ', names{j},'\t',...
                        'QDM:Visible light:Context','\n'];
        
        % print line to file
        fprintf(fileIDImg, contribLine);
        
        % write fitted map line to images table
        if isfile([dest 'images' filesep fittedMapFName])
            contribLine = [specimen,'\t',... 
                        'QDM_fitted_map_', names{j},'.jpg','\t',...
                        'Bz map','\t',...
                        'QDM Bz Map (fitted) of ', names{j},'\t',...
                        'QDM:Bz','\n'];
        
            % print line to file
            fprintf(fileIDImg, contribLine);
        end
    end
end

% cleanup
fclose(fileIDMeas);
close all

end

