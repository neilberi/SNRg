function [phi2] = interp1Phase(x1, x3, phi1, phi3, x2)
    phi2 = interp1([x1, x3], [phi1, phi3], x2);
    dInds = phi3-phi1<-pi;
    uInds = phi3-phi1>pi;
    phi2 = phi2 + 2*pi*(x2-x1)/(x3-x1).*(dInds - uInds);
    phi2(phi2>pi) = phi2(phi2>pi) - 2*pi;
    phi2(phi2<-pi) = phi2(phi2<-pi) + 2*pi;
end