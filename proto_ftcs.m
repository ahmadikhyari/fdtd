M = 50;
dx = 0.05;  dt = 0.025;  vp = 1;
x = 0:M-1;
v(1:10) = x(1:10);
v(11:20) = 20-x(11:20);
v(21:M) = 0;
vs = v;
vftcs = v;
t = 0:20;
fh = figure(1);
for n=1:length(t)
    figure(1)
    plot(x,v,x,vftcs)
    getframe();
    v(2:M) = v(1:M-1);
    vftcs(2:M-1) = vs(2:M-1)-(vp*dt/2/dx).*(vs(3:M)-vs(1:M-2));
    vs = vftcs;
end