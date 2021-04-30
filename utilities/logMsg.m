function msg = logMsg(level,message,lineEnd,indent,kwargs)
% simple logging function
arguments
    level
    message
    lineEnd
    indent
    kwargs.returnOnly = false
end

levels = containers.Map({'DEBUG','INFO', 'WARN', 'ERROR', 'NONE'}, {0, 1 ,2 ,3, 4});

try 
    logLevel = evalin('base','logLevel');
catch
    logLevel = 'ERROR';
end

% if a non standard level is passed e.g. FIT
try
    levelIdx = levels(upper(level));
catch
    levelIdx = 4;
end

logLevelIdx = levels(upper(logLevel));

if logLevelIdx > levelIdx
    return
end

try
    caller = dbstack(1);
    callerName = caller.name; % get functionName
catch
    callerName = 'none';
end

indent = repmat('  ',1,indent);

msg = sprintf('<> %6s %s: %s:%s %s', datetime('now', 'Format', 'HH:mm:ss:SS'), ...
              pad(upper(level), 6,'left'),...
              ['QDMlab.' callerName], indent, message);
          
if ~kwargs.returnOnly
    if lineEnd
        fprintf('%s\n',msg);
    else
        fprintf('%s',msg);
    end
elseif kwargs.returnOnly && lineEnd
    msg = sprintf('%s\n',msg);
end

end