function msg = logMsg(level,message,lineEnd,indent,kwargs)
% simple logging function
arguments
    level
    message
    lineEnd
    indent
    kwargs.returnOnly = false
end

levels = containers.Map({'ALL', 'DEBUG','INFO', 'WARN', 'ERROR', 'NONE'}, {-1, 0, 1 ,2 ,3, 4});

try 
    logLevel = evalin('base','logLevel');
catch
    logLevel = 'INFO';
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
              pad(upper(level), 7,'left'),...
              ['QDMlab.' callerName], indent, message);

% while msg           
          
if ~kwargs.returnOnly
    if strcmpi(level, 'INPUT') | strcmpi(level, 'ERROR')
        if lineEnd
            fprintf(2, '%s\n',msg);
        else
            fprintf(2, '%s',msg);
        end
    else
        if lineEnd
            fprintf('%s\n',msg);
        else
            fprintf('%s',msg);
        end
    end
    
elseif kwargs.returnOnly && lineEnd
    msg = sprintf('%s\n',msg);
end

end


function newMsg = split_msg(textWidth, msg, message)
% function clips the text to a fixed length
% it does not work properly and can clip text if it becomes too large

    nChar = size(msg,2);
    nSplits = floor(nChar/textWidth);
    spaces = find(msg==' ');
    newMsg = '';
    
    if nSplits
        logChars = numel(strrep(msg,message, ''));
        indentation = repmat(' ',1,logChars);
        splits = [];
        for i = 1:nSplits
            n = i*textWidth;
            [~, spaceIdx] = min(abs(spaces - i*textWidth));
            [~, nextSpaceIdx] = min(abs(spaces - (i+1)*textWidth));

            lineString = msg(spaces(spaceIdx):spaces(nextSpaceIdx));
            if i == 1
                newMsg = [msg(1:spaces(spaceIdx)-1) '\n'  repmat(' ',1,logChars) lineString];
            else
                newMsg = [newmsg '\n' repmat(' ',1,logChars) lineString];
            end
        end
    end
end
    