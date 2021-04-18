function out = is_B111(expData)
% out = is_B111(expData) checks the expData struct for the fieldnames
% 'B111ferro' or 'Bz'. Returns true if 'B111ferro' is found otherwise
% returns false.
arguments
    expData struct
end

if sum(strcmp(fieldnames(expData), 'B111ferro')) == 1
    out = true;
elseif sum(strcmp(fieldnames(expData), 'Bz')) == 1
    out = false;
end
