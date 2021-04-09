function out = exists_struct(structIn, key)

if sum(strcmp(fieldnames(structIn), key)) == 1
    out = true;
else
    out = false;
end