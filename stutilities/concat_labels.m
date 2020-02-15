function labelSetOut  = concat_labels( labelSetIn,labelIn )

[setLen,setWid] = size(labelSetIn);
labLen = length(labelIn);
dLen = abs(setLen-labLen);

if setLen > labLen
    labelSetOut = [labelSetIn,[labelIn;zeros(dLen,1)]];
elseif setLen < labLen
    labelSetOut = [[labelSetIn;zeros(dLen,setWid)],labelIn];
else
    labelSetOut = [labelSetIn,labelIn];
end


end

