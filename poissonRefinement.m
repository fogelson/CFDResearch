function [norms, ratios] = poissonRefinement(hCoarsest, levels)
    norms = zeros(1,levels);

    deltaT = 0.02;
    D = 1;
    H = 0;
    Q0 = 1;
    v1 = 3;
    v2 = 3;
    its = 30;
    timesteps = 1;
    offsetFrac = .3;
    
    for l=1:levels
        h = hCoarsest/(2^(l-1));
        [x,y,v] = mxGetGridInfo(h,h*offsetFrac);
        r = sqrt(x.^2 + y.^2);
        
        
        % Poisson problem settings
        f0 = 0*x;%sin(2*pi*x).*sin(2*pi*y);
        fE = cos(4*pi*x).*cos(12*pi*y);
        rhs = fE*(-16*pi^2 - 144*pi^2);
        
        % Radial poisson problem settings
        f0 = 0*x;
        fE = cos(2*pi*r);
        rhs = -(2*pi./r).*sin(2*pi*r) - 4*(pi^2)*cos(2*pi*r);
        
%         % Diffusion problem settings
%         t = deltaT*(timesteps);
%         f0 = cos(2*pi*x)+cos(2*pi*y);
%         fE = exp((-D*(2*pi)^2)*t)*f0;
%         rhs = 0*fE;
        
        fN = mxFENESolver(f0,rhs,h,deltaT,D,H,Q0,v1,v2,its,timesteps,offsetFrac);
        
         [n, ~] = size(fN);
         m = round(n/2);
         
         shift = fE(m,m) - fN(m,m);
         fE = fE - shift;
        
        
        %uncovered = find(~isnan(fN));
        
        %fNInt = sum(sum(fN(uncovered).*v(uncovered)*h^2));
        %fN = fN - fNInt;
        
        
        %[nX,~] = size(fN);
        %m = round(nX/2);
        
        err = abs(fN - fE);
        %err = err(2:end-1,2:end-1);
        
        norms(l) = max(max(err(2:end-1,2:end-1)));
        
        
        figure(l);
        subplot(2,3,1);
        mesh(x,y,fE)
        colorbar;
        subplot(2,3,2);
        mesh(x,y,fN)
        colorbar;
        subplot(2,3,3);
        mesh(x,y,err);
        colorbar;
        subplot(2,3,4)
        plot(y',err(1,:));
        subplot(2,3,5)
        plot(y',err(round(end/2),:))
        subplot(2,3,6)
        plot(y',err(end-1,:))
    end
    ratios = norms(1:end-1)./norms(2:end);
end