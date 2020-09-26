M = 50;
x = 0:M-1;
v(1:10) = x(1:10);
v(11:20) = 20-x(11:20);
v(21:M) = 0;
t = 0:20;
fh = figure(1);
for n=1:length(t)
    figure(1)
    plot(x,v)
    getframe();
    v(2:M) = v(1:M-1);
end