function [loc] = invCDF(w,seq,m)
% Input:
% w: the external weight to compare with
% seq: cumulative weights at each integer points
% m: maximum length
% uniform density between integer axis
% output: axis from 0 to m
% matlab accuracy is up to 16 digits
    if w < 0 || w >= 1+1e-6
        w
        error('w value out of bound');
    end
    if abs(w-1) < 1e-6
        loc = m;
        return;
    end
    if w <= seq(1) && seq(1) ~= 0
        loc = w / seq(1);
        return;
    end
    if w == seq(1) && seq(1) == 0
        loc = 0;
        return;
    end
    
    for id = 1:(m-1)
        if w == seq(id)
            loc = id;
            return;
        elseif w > seq(id) && w <= seq(id+1)
            loc = id + (w-seq(id)) / (seq(id+1)-seq(id));
            return;
        else
            continue;
        end
    end
    
    if id == m-1
       format long;
       disp(w);
       disp(seq);
       error('digit accuracy problem'); 
    end
    
end