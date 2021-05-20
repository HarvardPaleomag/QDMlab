function laserData = get_laser(path, kwargs)
%[laserData] = get_laser(path)

    arguments
        path
        kwargs.data = 'none'
    end
    
    if ~strcmp(kwargs.data, 'none')
        if exists_struct(kwargs.data, 'laser')
            laserData = kwargs.data.laser;
            return
        end
    end
    
    laserFiles = dir(fullfile(path,'laser.*'));
    [~,~,ext] = fileparts(laserFiles(1).name);
    
    if strcmp(ext, '.jpg')
        laserData = imread(fullfile(path, laserFiles(1).name));
        return
    end
    
    if strcmp(ext, '.csv')
        laserData = load(fullfile(path, laserFiles(1).name));
    end
end