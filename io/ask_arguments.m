function kwargs = ask_arguments(kwargs, defaults)
%[kwargs] = ask_arguments(defaults)
% ask for arguments
% get all default values in struct
%
% Parameters
% ----------
% 

fn = fieldnames(defaults);

for k=1:numel(fn)
    if strcmp(kwargs.(fn{k}), 'none')
        default = defaults.(fn{k});
        
        msg = sprintf('<< %s [%s] >>: ', fn{k}, string(default));
        prompt = logMsg('INPUT',msg,0,0,'returnOnly',true);
%         prompt = sprintf('<>   INPUT << %s [%s] >>: ', fn{k}, string(default));
        inP = input(prompt);
        
        if isempty(inP)
            inP = default;
        end
        
        kwargs.(fn{k}) = inP;
    end
end