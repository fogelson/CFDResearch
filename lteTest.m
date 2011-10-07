for index = 1:maxIndex
    h = .25/(2*2^index); offsetFrac = .5; offset = h*offsetFrac;
    [x,y,v] = mxGetGridInfo(h,offset); covered = find(v == 0); v(covered) = 0; r = sqrt(x.^2 + y.^2); r(covered) = nan;
    
    g = cos(2*pi*r);
    Lg = mxApply(g,h,offset);
    f = -2*pi*sin(2*pi*r)./r - 4*(pi^2)*cos(2*pi*r);
    
    %g = (r.^2).*((1-r).^2);
    %Lg = mxApply(g,h,offset);
    %f = (2 - 6*r.^2 + 4*r.^3) + (2 - 12*r + 12*r.^2);
    
    
    tau = abs(f - Lg);
    [vI, vP, v3, v4, v5] = mxSplit(h,offset,v);
    [tauI, tauP, tau3, tau4, tau5] = mxSplit(h,offset,tau);
    lteI(index) = max(max(abs(tauI)));
    lteP(index) = max(max(abs(tauP)));
    lte3(index) = max(max(abs(tau3)));
    lte4(index) = max(max(abs(tau4)));
    lte5(index) = max(max(abs(tau5)));
    
    maxI(index) = max(max(vI));
    maxP(index) = max(max(vP));
    max3(index) = max(max(v3));
    max4(index) = max(max(v4));
    max5(index) = max(max(v5));
    
    minI(index) = min(min(vI));
    minP(index) = min(min(vP));
    min3(index) = min(min(v3));
    min4(index) = min(min(v4));
    min5(index) = min(min(v5));
end