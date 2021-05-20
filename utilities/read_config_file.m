function read_config_file()
    fname = '/Users/mike/Desktop/QDMlab_config.ini'; 
    fid = fopen(fname); 
    raw = fread(fid,inf); 
    str = char(raw'); 
    fclose(fid); 
    val = jsondecode(str);