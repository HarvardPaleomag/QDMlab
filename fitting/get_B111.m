function [B111ferro, B111para] = get_B111(negB111, posB111)
  gamma = 0.0028; 
  
  negDiff = - real( (negB111.Resonance2-negB111.Resonance1)/2 / gamma );
  posDiff =   real( (posB111.Resonance2-posB111.Resonance1)/2 / gamma );
  
  B111ferro = (posDiff + negDiff)/2;
  B111para = (posDiff - negDiff)/2;