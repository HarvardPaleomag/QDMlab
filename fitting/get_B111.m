function [B111ferro, B111para] = get_B111(negRes, posRes)
% get B111 calcuates the B111 value from pos and negive resonances
%
% Parameters
% ----------
%   negRes: double
%        Values for the high and low MW resonances with neg. applied field
%   posRes: double
%        Values for the high and low MW resonances with pos. applied field

gamma = 0.0028;

negDiff = - real( (negRes.Resonance2-negRes.Resonance1)/2 / gamma );
posDiff =   real( (posRes.Resonance2-posRes.Resonance1)/2 / gamma );

B111ferro = (posDiff + negDiff)/2;
B111para = (posDiff - negDiff)/2;
end