function zeroCrossExists = checkZeroCross(FTotal)

% print zero crossing
zeroCrossExists = false;
if sum(FTotal(2:end-1)>0) > 0 && sum(FTotal(2:end-1)<0) > 0
    zeroCrossExists = true;
end

