%%% 
% This code is for peak-finding processing for Erik's data
% Change log: 
% 1. 03/10/17 remove fake peaks and double peaks problem by check loc
% diff length
%%%
function [pkphase] = getPhase_peakfinding_SCN(X)

    [T, N] = size(X);
    pkphase = zeros(T, N);
    peak_len = zeros(N, 1);
    time_seq = 1:T;

    for osc_id = 1 : N
        osc = X(time_seq, osc_id);
%         figure(1);
%         plot(osc); hold on;
        % filter to make the curve smooth
        [b, a] = butter(12, 0.2, 'low');
        osc = filtfilt(b, a, osc);
%         plot(osc); hold off;
        [pk, loc] = findpeaks(osc);
%         % remove fake peak by check the loc diff length
%         extra = find(diff(loc) < 18); % 18 for scn 2,3,4
%         if ~isempty(extra)
%             extra_id = extra(1) + 1;
%             loc(extra_id) = [];
%         end
        peak_len(osc_id) = length(loc);
%         scatter(loc, pk, 'x'); hold off;
%         legend('8');
        Y = (0 : (length(loc) - 1)) * 2 * pi;
        % linear interpolation with extrapolation
        p = interp1(loc, Y, time_seq, 'pchip', 'extrap');
        pkphase(:, osc_id) = p';
    end
    
    tabulate(peak_len);
    % choose peak_len number according to major post-TTX peaks from table
    % printed
%     chosen_id = find(peak_len == 9);
%     pkphase = pkphase(:, peak_len == 9); 
end