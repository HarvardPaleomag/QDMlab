function [input] = automatic_input_ui__(input, kwargs)
% Select files or directories using a GUI. If the input argument is the
% string 'none', the function opens a GUI window to let the user select
% files or directories. Otherwise, the input argument is returned
% unchanged. 
%
% INPUTS:
% - input: A cell array containing the path(s) to the file(s) or
%   directory(ies) to be selected, or the string 'none' to trigger the GUI
%   window for selection.
% - kwargs: A struct containing optional keyword arguments for customizing
%   the GUI window. Fields include:
%     - 'type': A string specifying whether to select files or directories.
%       Default is 'dir'.
%     - 'title': A string specifying the title of the GUI window. Default
%       is 'select file/directory'.
%     - 'MultiSelect': A string specifying whether to allow selection of
%       multiple files or directories. Default is 'off'.
%     - 'single': A boolean specifying whether to return a single path as a
%       string, rather than a cell array containing a single string. Default
%       is false.
%     - 'filter': A string specifying the file extension filter to use for
%       selecting files. Default is '*.mat'.
%
% OUTPUTS:
% - input: A cell array containing the selected path(s) to the file(s) or
%   directory(ies), or a single string if the 'single' keyword argument is
%   true.

% If input is 'none', trigger GUI window for file or directory selection
if strcmp(input, 'none')
    input = {};
    logMsg('input', kwargs.title, 1, 0);

    if strcmp(kwargs.type, 'dir')
        % Select directory/directories
        if strcmp(kwargs.MultiSelect, 'on')
            % Allow multiple selections
            path = 1;
            while path ~= 0
                path = uigetdir('title', [kwargs.title ' cancel to stop selecting']);
                if path ~= 0
                    input{end+1} = path;
                end
            end
        else
            % Single selection
            input = {uigetdir('title', kwargs.title)};
        end
        
    end
    if strcmp(kwargs.type, 'file')
        % Select file/files
        [file, path] = uigetfile(kwargs.filter, kwargs.title, 'MultiSelect', kwargs.MultiSelect);
        if strcmp(kwargs.MultiSelect, 'on')
            % Allow multiple selections
            for f = file
                f = fullfile(path, f);
                input{end+1} = f{:};
            end
        else
            % Single selection
            input = {fullfile(path, file)};
        end
    end
end

% Ensure input is in cell array format
input = correct_cell_shape(input);

% Check if any files were selected
if size(input,1) == 1
    if strcmp(input{1}, '/')
        error('<>   ERROR: NO files/folders selected, please specify or pick files/folders.')
    end
end

% If 'single' is true, return a single string rather than a cell array with
% a single string
if kwargs.single
    input = input{:};
end

end
