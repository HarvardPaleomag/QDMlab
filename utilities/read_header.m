function header = read_header(file)

data = importdata(file);

for i=1:size(data.rowheaders)
    key = data.rowheaders{i};
    key = strrep(key, ' ', '_');
    key = strrep(key, '(', '');
    key = strrep(key, ')', '');
    key = strrep(key, '?', '');

    header.(key) = data.data(i);
end
end