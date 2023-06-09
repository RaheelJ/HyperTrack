%Function to check if any parameter voilates the limits
function out = Check_Limit(in, upper_lim, lower_lim)
    if (in > upper_lim || in < lower_lim)
        out = 1;
    else
        out=0;
    end
end
