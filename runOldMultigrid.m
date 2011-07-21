function [x,y,u,its] = runOldMultigrid(h,v1,v2,tol,maxIts)
    cd ..
    cd Multigrid
    x = -1:h:1;
    y = -1:h:1;
    [x,y] = meshgrid(x,y);
    u0 = (x<0).*(y<0);% + (x>0).*(y>0);
    u0(2:end-1,2:end-1) = 1;
    f = 0*x;
    [u,its] = mxMultigrid(u0,f,v1,v2,tol,maxIts);
    cd ..
    cd CFDResearch
    u = u(2:end-1,2:end-1);
    x = x(2:end-1,2:end-1);
    y = y(2:end-1,2:end-1);
end