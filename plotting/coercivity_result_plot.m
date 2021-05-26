function coercivity_result_plot(results, kwargs)
%coercivity_result_plot(results; 'steps', 'stepUnit', 'led', 'mean')
% plots results from estimate_coercivity
% 
% Parameters
% ----------
%     positional
%     ==========
%         results: struct
%             results structure from 'estmate_coercivity'
%     keyword
%     =======
%         steps: bool (0)
%         led: bool (0)
%
arguments
    results struct
    kwargs.steps  (1,:) double = false
    kwargs.stepUnit  (1,:) = 'mT'
    kwargs.led  (1,1) {mustBeMember(kwargs.led, [1, 0])} = 0
    kwargs.mean (1,1) {mustBeBoolean(kwargs.mean)} = false
end
    

msg = sprintf('please use demag_behavior_plot. coercivity_result_plot will be removed in a later update');
logMsg('deprecated',msg,1,0);

kwargs = namedargs2cell(kwargs);
demag_behavior_plot(results, kwargs{:})