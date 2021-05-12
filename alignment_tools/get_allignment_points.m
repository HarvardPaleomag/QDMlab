function allignment_points = get_allignment_points(nFolders, filename, savefile)
%{
These codes (1) register the maps and (2) fit the sources in each map%
parameters:
    folders:list
        list of absolute path to each data folder. First entry in the list
        is used as reference image.
    filename: str
        name of the .mat file to be used.
        default: 'Bz_uc0.mat'
    savefile: str
        absolute path of the file the shifts are in. If file does not exist
        the file will be created in the current directory,
        so that you dont have to do this every time.

%}
skip=0;

% check if filename was passed
if nargin <2
    filename ='Bz_uc0.mat';
end

% check if shift_file was passed
if nargin <3
    savefile = [pwd filesep 'shifts.mat'];
end

% loading of shift files in case you saved them
if isfile(savefile)
    disp('loading file')
    load(savefile, 'angles', 'dx', 'dy', 'all_alignment_points', '-mat');
else
    % check if the shift_file had not been created
    % i.e. allalignementpoints, anglesm dx, dy have to be 0
    all_aps = zeros(size(nFolders,1),4);
end

if all(all_aps==0) && all(angles==0) && all(dx==0) && all(dy==0)
    for i= 1 :size(nFolders,1)
        fname=[deblank(nFolders(i,:)) filesep filename];

        load(fname,'newLED')

        % use image 1 as reference
        if i == 1
            m1 = msgbox('select two points on the LED image as reference');
            image_aps = data_from_plot(newLED,2,true,false);
            set(gcf, 'Position', [1100,250,1000,1000])
            title('reference')
            hold off
            if exist('m1', 'var')
                delete(m1);
                clear('m1')
            end
        else
            % all other images
            pause(0.3)
            image_aps = data_from_plot(newLED,2,true,true);

        end

         % append to rest of the points
        all_aps(i,:) = image_aps;
    end
end
