clear

mu0 = pi*4e-7; vp = 299792458; ep0 = 1/(vp^2*mu0);
lam = 1;
dx = lam / 20;
cfl = 1.00000;
dt = cfl * dx/abs(vp);

xx = 400;
width = (xx-1) * dx;
x = 0:dx:width;

Ezold = zeros(1,xx);
Hyold = zeros(1,xx);

tau = 10;
t0 = 3 * tau;

tstep = 500;

maxE = zeros(1,tstep);
maxexp = zeros(1,tstep);

myVideo = VideoWriter('fdtd_stab_1.00000');
open(myVideo);

for t = 1:tstep
    Ezold(1) = exp(-(t-t0).^2 / tau^2);
    
    Hynew(1:xx-1) = Hyold(1:xx-1) + dt/mu0/dx * (Ezold(2:xx) - Ezold(1:xx-1));
    Eznew(2:xx-1) = Ezold(2:xx-1) + dt/ep0/dx * (Hynew(2:xx-1) - Hynew(1:xx-2));
    % absorbing boundary condition (right end)
    Eznew(xx) = Ezold(xx-1) + (vp*dt-dx)/(vp*dt+dx) * (Eznew(xx-1) - Ezold(xx));
    
    figure(1);
    subplot(2,1,1);
    plot(x,Eznew);
    axis([0 width -0.1 1.1]);
    xlabel('x (meter)');
    ylabel('Ez');
    str = {['dt = ',num2str(dt)];['C = ',num2str(cfl,6)];['t = ',num2str(t*dt)]};
    text(15,0.8,str);
    
    maxE(t) = abs(max(Eznew));
    if t > 2*t0 && t < 2*t0+xx
        maxexp(t) = abs(max(Eznew(1:ceil((t-2*t0)*cfl))));
    else if t >= 2*t0+xx
        maxexp(t) = abs(max(Eznew));
        end
    end
    
    subplot(2,1,2);
    plot((1:tstep)*dt*1e9,maxE(1:tstep),(1:tstep)*dt*1e9,maxexp(1:tstep));
    axis tight;
    xlabel('time (ns)');
    ylabel('max(Ez)');
    
    frame = getframe(gcf);
    writeVideo(myVideo, frame);
    
    Ezold = Eznew;
    Hyold = Hynew;
end

close(myVideo);