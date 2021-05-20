function out = exists_struct(structIn, key)
%[out] = exists_struct(structIn, key)

if isfield(structIn, key)
    msg = sprintf('found key %s in the structure', key);
    logMsg('debug',msg,1,0);
    out = true;
else
    msg = sprintf('key %s does not exist in the structure', key);
    logMsg('debug',msg,1,0);
    out = false;
end