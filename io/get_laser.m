function laserData = get_laser(path, kwargs)
%[laserData] = get_laser(path; 'data')
% loads/ returns the laser image data from a folder or data structure
%
% Parameters
% ----------
%   path: str
%       path to the folder that contains the laser data
%   data: struct (false)
%       allows to also pass a data structure, if that structure contains
%       the laser data it will be returned without checking for the laser
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
    if isfield(kwargs.data, 'laser')
        laserData = kwargs.data.laser;
        return
    end
end

if isfile(path)
    msg = sprintf('File was passed instead of path.');
    logMsg('debug',msg,1,0);
    [path, ~, ~] = fileparts(path);
end

laserFiles = {};
n = 0;
path_ = path;

while numel(laserFiles) == 0
    if n == 2
        error('NO laser.* files found in << %s >> or the two directories above!', path)
    end
    
    laserFiles = dir(fullfile(path_,'laser.*'));

    if numel(laserFiles) == 0
        msg = sprintf('NO laser.* files found! Going up one directory');
        logMsg('error',msg,1,0);
        splitPath = split(path_, filesep);
        path_ = fullfile(splitPath{1:end-1});
        n = n+1;
    end
end

msg = sprintf('found %i laser.* files', numel(laserFiles));
logMsg('debug',msg,1,0);

[~,~,ext] = fileparts(laserFiles(1).name);

msg = sprintf('loading: %s', fullfile(path, laserFiles(1).name));
logMsg('debug',msg,1,0);

if strcmp(ext, '.csv')
    laserData = load(fullfile(path, laserFiles(1).name));
end

if strcmp(ext, '.jpg')
    laserData = imread(fullfile(path, laserFiles(1).name));
    return
end
    

end