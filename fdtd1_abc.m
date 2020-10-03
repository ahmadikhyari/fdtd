clear
clc

% inisiasi
% batas x di grafik
xmin = -10.8;  xmax = 10.8;
% batas t waktu pengamatan
t = 0;         tmax = 54e-9;
% variabel
dx = 0.054;  dt = 0.18e-9;
mu0 = pi*4e-7; vp = 299792458; ep0 = 1/(vp^2*mu0);
% variabel sinyal sinusoildal
f = 278e6;

xe = (xmin):dx:(xmax); % sumbu x untuk Ez
i = length(xe); % banyaknya kuantisasi i dari posisi xe
xsource = 0; % posisi sumber
xs = floor((xsource-xe(1)+dx)/dx);
xe = xe-xe(xs);
xh = xe(1:i-1)+dx/2; % sumbu x untuk Hy

nstep = ceil((tmax-t)/dt); % banyaknya kuantisasi n dari waktu t
% inisiasi kondisi awal dari sinyal v (sinyal sinusoidal)
Hy(1:2,1:i-1)=0.0; 
Ez(1:2,1:i)=0.0;
% *  #  *  #  *  # .... *  #  *
% E  H  E  H  E  H      E  H  E
% 1  1  2  2  3  3     i-1i-1 i

% inisiasi variabel terkait video
myVideo = VideoWriter('fdtd1'); %open video file
myVideo.FrameRate = 12;  % 5 - 10 works well
open(myVideo)

for n=0:nstep
    % nilai Ez di titik sumber (berubah terhadap n)
    Ez(1,xs) = sin(2*pi*f*n*dt); % sinusoidal
    Ez(2,xs) = exp(-(n-8)^2/(4^2)); % gaussian pulse
    %%% rumus leapfrog (TM)
    % *  #  *  #  *  # .... *  #  *
    % E  H  E  H  E  H      E  H  E
    % 1  1  2  2  3  3     i-1i-1 i
    Hy(:,1:i-1)=1.*Hy(:,1:i-1)+(dt/mu0/dx).*(Ez(:,2:i)-Ez(:,1:i-1));
    Ezbx(:,1) = Ez(:,1); Ezbx(:,2) = Ez(:,2); % variabel sementara
    Ezbx(:,3) = Ez(:,i-1); Ezbx(:,4) = Ez(:,i); % variabel sementara
    Ez(:,2:i-1)=1.*Ez(:,2:i-1)+(dt/ep0/dx).*(Hy(:,2:i-1)-Hy(:,1:i-2));
    % first-order mur boundary
    Ez(:,1) = Ezbx(:,2)+((vp*dt-dx)/(vp*dt+dx))*(Ez(:,2)-Ezbx(:,1));
    Ez(:,i) = Ezbx(:,3)+((vp*dt-dx)/(vp*dt+dx))*(Ez(:,i-1)-Ezbx(:,4));
    % plot
    figure(1)
    subplot(2,1,1);
    plot(xe,Ez(1,:),xh,Hy(1,:));
    axis([-12 12 -1.2 1.2]);
    ylabel('{\color{blue}E_z} and {\color{red}H_y}');
    title(['Sinusoidal Source at x = ',num2str(ceil(xe(xs))),...
        ', Snapshot at t = ',num2str(t)]);
    subplot(2,1,2);
    plot(xe,Ez(2,:),xh,Hy(2,:));
    axis([-12 12 -1.2 1.2]);
    ylabel('{\color{blue}E_z} and {\color{red}H_y}');
    title(['Gaussian Pulse Source at x = ',num2str(ceil(xe(xs))),...
        ', Snapshot at t = ',num2str(t)]);
    txt = xlabel('x (meter)');
    getframe();
    % increment of time
    t = t+dt;
    % perekaman gambar sebagai frame video
    frame = getframe(gcf); %get frame
    writeVideo(myVideo, frame);
end

close(myVideo)
