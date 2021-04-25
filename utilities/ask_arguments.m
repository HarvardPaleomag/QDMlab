function kwargs = ask_arguments(kwargs, defaults)

fn = fieldnames(kwargs);
for k=1:numel(fn)
    if strcmp(kwargs.(fn{k}), 'none')
        if exists_struct(defaults, fn{k})
            default = defaults.(fn{k});
        else
            default = 'not specified';
        end
        
        prompt = sprintf('<>   INPUT %s: ', fn{k});
        inP = input(prompt)
        
        if ~ inP
            inP = default
        end
        
        kwargs.(fn{k}) = inP;
    end
end