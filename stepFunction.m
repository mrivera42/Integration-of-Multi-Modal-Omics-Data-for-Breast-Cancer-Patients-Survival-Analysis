 function h = stepFunction(diff)
    if diff > 0
        h = 1;
    elseif diff == 0
        h = 0.5;
    else
        h = 0;
    end
 end