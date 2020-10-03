clear
clc

% inisiasi
% batas x di grafik
xmin = -10.8;  xmax = 10.8;
% batas t waktu pengamatan
t = 0;         tmax = 32.4e-9;
% variabel
dx = 0.054;  dt = 0.18e-9;
mu0 = pi*4e-7; vp = 299792458; ep0 = 1/(vp^2*mu0);
mur = 1;        epr = 1;
sigm = [0;0];       sige = [0;0.001];
% variabel sinyal sinusoildal
f = 278e6;

xe = (xmin):dx:(xmax); % sumbu x untuk Ez
i = length(xe); % banyaknya kuantisasi i dari posisi xe
xsource = 0; % posisi sumber
xs = floor((xsource-xe(1)+dx)/dx);
xe = xe-xe(xs);
xh = xe(1:i-1)+dx/2; % sumbu x untuk Hy

nstep = ceil((tmax-t)/dt); % banyaknya kuantisasi n dari waktu t

ca = (2*epr*ep0-sige*dt)./(2*epr*ep0+sige*dt);
cb = 2*dt./(2*epr*ep0+sige*dt);
cc = (2*mur*mu0-sigm*dt)./(2*mur*mu0+sigm*dt);
cd = 2*dt./(2*mur*mu0+sigm*dt);

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
    %Ez(1,xs) = sin(2*pi*f*n*dt); % sinusoidal
    Ez(1,xs) = exp(-(n-8)^2/(4^2)); % gaussian pulse
    Ez(2,xs) = exp(-(n-8)^2/(4^2)); % gaussian pulse
    %%% rumus leapfrog (TM)
    % *  #  *  #  *  # .... *  #  *
    % E  H  E  H  E  H      E  H  E
    % 1  1  2  2  3  3     i-1i-1 i
    Hy(:,1:i-1)=cc.*Hy(:,1:i-1)+cd/dx.*(Ez(:,2:i)-Ez(:,1:i-1));
    Ezbx(:,1) = Ez(:,1); Ezbx(:,2) = Ez(:,2); % variabel sementara
    Ezbx(:,3) = Ez(:,i-1); Ezbx(:,4) = Ez(:,i); % variabel sementara
    Ez(:,2:i-1)=ca.*Ez(:,2:i-1)+cb/dx.*(Hy(:,2:i-1)-Hy(:,1:i-2));
    % first-order mur boundary
    Ez(:,1) = Ezbx(:,2)+((vp*dt-dx)/(vp*dt+dx))*(Ez(:,2)-Ezbx(:,1));
    Ez(:,i) = Ezbx(:,3)+((vp*dt-dx)/(vp*dt+dx))*(Ez(:,i-1)-Ezbx(:,4));
    % plot
    figure(1)
    subplot(2,1,1);
    plot(xe,Ez(1,:),xh,Hy(1,:));
    axis([-12 12 -1.2 1.2]);
    ylabel('{\color{blue}E_z} and {\color{red}H_y}');
    title(['Snapshot at t = ',num2str(t),' \sigma = ',num2str(sige(1))]);
    subplot(2,1,2);
    plot(xe,Ez(2,:),xh,Hy(2,:));
    axis([-12 12 -1.2 1.2]);
    ylabel('{\color{blue}E_z} and {\color{red}H_y}');
    title(['Snapshot at t = ',num2str(t),' \sigma = ',num2str(sige(2))]);
    xlabel('x (meter)');
    getframe();
    % increment of time
    t = t+dt;
    % perekaman gambar sebagai frame video
    frame = getframe(gcf); %get frame
    writeVideo(myVideo, frame);
end

close(myVideo)
