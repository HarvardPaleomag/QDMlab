function kwargs = ask_arguments(kwargs, defaults)
% ask for arguments

% get all default values in struct
fn = fieldnames(defaults);

for k=1:numel(fn)
    if strcmp(kwargs.(fn{k}), 'none')
        default = defaults.(fn{k});
        
        prompt = sprintf('<>   INPUT << %s [%s] >>: ', fn{k}, string(default));
        inP = input(prompt);
        
        if ~ inP
            inP = default;
        end
        
        kwargs.(fn{k}) = inP;
    end
end