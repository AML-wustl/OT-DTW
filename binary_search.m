function idx = binary_search(x, ordered_seq)
    len = length(ordered_seq);
    left_i = 1;
    right_i = len;
    if x <= ordered_seq(1)
        idx = 1;
        return;
    end
    while (right_i - left_i > 1)
        round_int = round((left_i+right_i)/2);
        if x > ordered_seq(round_int)
            left_i = round_int;
        elseif x < ordered_seq(round_int)
            right_i = round_int;
        else
            idx = round_int;
            return;
        end
    end
    idx = right_i;
end