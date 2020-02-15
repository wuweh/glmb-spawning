function Z= gen_observation(model,X)
% this observation generation function is for coordinate
% measurements

if isempty(X),
    Z= [];
else
    Z= model.C_posn*X+ model.D*randn(2,size(X,2)); %coordinate extraction
end;