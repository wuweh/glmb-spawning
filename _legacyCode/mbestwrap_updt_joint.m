function [assignments,costs]= murty_wrapper_update(P0,m)
    
  n1 = size(P0,1);
  n2 = size(P0,2);

  % Make costs non-negative (required by 'assignmentoptimal')
  x = min(min(P0));
  P0 = P0 - x;
  
  % Murty
  [assignments, costs] = murty_custom(P0,m);

  % Restore correct costs to assignments
  costs = costs + (x.*sum(assignments>0,2))';
      
    