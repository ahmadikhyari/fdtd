clear
clc

% inisiasi
xmin = 0;
xmax = 2.5;
t = 0;
tmax = 0.5;
dx = 0.05;
dt = 0.025;
vp = 1;
p = 6;      a = 0;       b = 1;

x = (xmin-dx):dx:(xmax+dx); %sumbu x
i = (xmax-xmin)/dx; % banyaknya kuantisasi i dari posisi x
% inisiasi kondisi awal dari sinyal v
v0 = (2*p/(b-a))*((x-a).*(x>=a & x<=(a+b)/2)+(b-x).*(x>(a+b)/2 & x<=b));
vs = v0;
vftcs = v0;
nstep = (tmax-0)/dt; % banyaknya kuantisasi n dari waktu t

% inisiasi variabel terkait video
myVideo = VideoWriter('ftcs'); %open video file
myVideo.FrameRate = 5;  % 5 - 10 works well
open(myVideo)

for n=0:nstep
    % variabel vs sebagai penampung nilai FTCS
    vs = vftcs;
    vtt(n+1,:) = vs(:);
    % rumus umum, sinyal segitiga merambat...
    vor = (2*p/(b-a))*((x-vp*t-a).*((x-vp*t)>=a & (x-vp*t)<=(a+b)/2)...
        +(b-(x-vp*t)).*((x-vp*t)>(a+b)/2 & (x-vp*t)<=b));
    % plot
    figure(1)
	plot(x,vor,x,vftcs,'.-')
    axis([0 xmax -2 8]) % batas sumbu-sumbu grafik
    str = {'t = ',t}; % tampilan waktu
    text(1.75,5,str) % tampilan waktu
    str1 = {'dt = ' dt, '','dx = ' dx, '', 'vp = ' vp}; % tampilan waktu
    text(0.25,6,str1) % tampilan waktu
	getframe();
    % increment of time
    t = t+dt;
    % perekaman gambar sebagai frame video
    frame = getframe(gcf); %get frame
    writeVideo(myVideo, frame);
    % rumus FTCS
    vftcs(2:i+2) = vs(2:i+2)-(vp*dt/2/dx)*(vs(3:i+3)-vs(1:i+1));
end

close(myVideo) % penulisan video berhenti.