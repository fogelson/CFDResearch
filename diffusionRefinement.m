function [norms, ratios] = diffusionRefinement(levels,hCoarsest)
    deltaT = .02;
    D = 1;
    H = 0;
    Q0 = 1;
    v1 = 2;
    v2 = 2;
    its = 3000;
    timesteps = 1;

    norms = zeros(levels,1);
    
    h = hCoarsest;
    
    for l = 1:levels
        h = hCoarsest/(2^(l-1));
        [fNumerical, fExact] = doSolve();
        diff = abs(fNumerical - fExact);
        norms(l) = max(max(diff));
    end
   
    ratios = norms(1:end-1)./norms(2:end);
%    ratios = norms(2:end)./norms(1:end-1);
    

    function [fNumerical, fExact] = doSolve()
        [x, y, ~] = mxGetGridInfo(h,h/2);
        f0 = 0*x;%cos(pi*x);%+cos(pi*y);
        rhs = (-pi^2)*(cos(pi*x));%+(-pi^2)*cos(pi*y);
        fNumerical = mxFENESolver(f0,rhs,h,deltaT,D,H,Q0,v1,v2,its,timesteps);
        t = deltaT*timesteps;
        fExact = cos(pi*x);%+cos(pi*y);
        %figure
        %subplot(1,2,1)
        %surf(x,y,fNumerical)
        %subplot(1,2,2)
        %surf(x,y,fExact)
        figure
        %graph = surf(x,y,abs(fExact - fNumerical));
        %graph = surf(x,y,rhs);
        graph = surf(x,y,abs(fNumerical-fExact));
        set(graph,'EdgeColor','None');
    end

end