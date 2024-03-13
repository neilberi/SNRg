function [dphi] = phasediff(S2,S1)
    dphi = angle(S2) - angle(S1);
    for k = 1:numel(dphi)
        if(dphi(k) >= pi)
            dphi(k) = dphi(k) - 2*pi;
        elseif(dphi(k) <= -pi)
            dphi(k) = dphi(k) + 2*pi;
        end
    end
end