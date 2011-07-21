
u = u0;
for j = 1:1000
    u = mxSmootherTest(0,0,0,h,offset,u,f,20,0.0000005);
    graph = surf(x,y,u);
    set(graph, 'Edgecolor', 'None');
    drawnow;
end
