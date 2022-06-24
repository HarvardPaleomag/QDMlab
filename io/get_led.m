function ledData = get_led(path, kwargs)
%[ledData] = get_led(path; 'data')
% loads/ returns the led (reflected light) image data from a folder or 
% data structure
%
% Parameters
% ----------
%   path: str
%       path to the folder that contains the laser data
%   data: struct (false)
%       allows to also pass a data structure, if that structure contains
%       the LED data it will be returned without checking for the LED
%       file
%
% Returns
% -------
%   laserData
%

arguments
    path
    kwargs.data = false
end

if ~isequal(kwargs.data, false)
    msg = sprintf('NOT IMPLEMENTED');
    logMsg('error',msg,1,0);
%     if isfield(kwargs.data, 'laser')
%         ledData = kwargs.data.led;
%         return
%     end
end


ledFiles = {};
n = 0;
path_ = path;

while numel(ledFiles) == 0
    if n == 2
        error('NO LED.* files found in << %s >> or the two directories above!', path)
    end
    
    ledFiles = dir(fullfile(path_,'LED.*'));

    if numel(ledFiles) == 0
        msg = sprintf('NO LED.* files found! Going up one directory', numel(ledFiles));
        logMsg('error',msg,1,0);
        splitPath = split(path_, filesep);
        path_ = fullfile(splitPath{1:end-1});
        n = n+1;
    end
end

msg = sprintf('found %i LED.* files', numel(ledFiles));
logMsg('debug',msg,1,0);

[~,~,ext] = fileparts(ledFiles(1).name);

msg = sprintf('loading: %s', fullfile(path, ledFiles(1).name));
logMsg('debug',msg,1,0);

if strcmp(ext, '.csv')
    ledData = load(fullfile(path, ledFiles(1).name));
end

if strcmp(ext, '.jpg')
    ledData = imread(fullfile(path, ledFiles(1).name));
    return
end
    

end