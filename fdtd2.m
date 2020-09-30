clear
clc

% inisiasi
mu0 = pi*4e-7; vp = physconst('LightSpeed'); ep0 = 1/(vp^2*mu0);
mur = 1;        epr = 1;
sigm = 0;       sige = 0;
% batas x di grafik
dx = 0.054;   xmin = -5.4;  xmax = 5.4;
dy = 0.054;   ymin = -5.4;  ymax = 5.4;
% batas t waktu pengamatan
dt = dx/vp/sqrt(2);      t = 0;      tmax = 400*dt;
% variabel sinyal sinusoildal
f = 278e6;

xe = (xmin):dx:(xmax); % sumbu x untuk Ez
ye = (ymin):dy:(ymax); % sumbu y untuk Ez
ii = length(xe); % banyaknya kuantisasi i dari posisi xe
jj = length(ye); % banyaknya kuantisasi j dari posisi ye
xsource = 0; % posisi sumber
ysource = 0; % posisi sumber
xs = floor((xsource-xe(1)+dx)/dx);
xe = xe-xe(xs); % sumbu x untuk sumbu bilangan bulat
xh = xe(1:ii-1)+dx/2; % sumbu x untuk sumbu paruh
ys = floor((ysource-ye(1)+dy)/dy);
ye = ye-ye(ys); % sumbu y untuk sumbu bilangan bulat
yh = ye(1:jj-1)+dy/2; % sumbu y untuk sumbu paruh

cay = abs(xh.^(-0.5)./xh(xs)^(-0.5)); % plot 1/sqrt(2)

% koefisien terkait redaman
ca = (2*epr*ep0-sige*dt)/(2*epr*ep0+sige*dt);
cb = 2*dt/(2*epr*ep0+sige*dt);
cc = (2*mur*mu0-sigm*dt)/(2*mur*mu0+sigm*dt);
cd = 2*dt/(2*mur*mu0+sigm*dt);

% inisiasi kondisi awal dari sinyal v (sinyal sinusoidal)
Hxs(1:jj-1,1:ii)=0.0; 
Hys(1:jj,1:ii-1)=0.0; 
Ezs(1:jj,1:ii)=0.0;
Hxp(1:jj-1,1:ii)=0.0; 
Hyp(1:jj,1:ii-1)=0.0; 
Ezp(1:jj,1:ii)=0.0;
% *  y  *  y  *  y .... *  y  *             +---> x
% x     x     x ....    x     x             |
% *  y  *  y  *  y .... *  y  *             |
% x     x     x ....    x     x             
% *  y  *  y  *  y .... *  y  *             y

nstep = ceil((tmax-t)/dt); % banyaknya kuantisasi n dari waktu t

% inisiasi variabel terkait video
%myVideo = VideoWriter('fdtd2'); %open video file
%myVideo.FrameRate = 25;  % 5 - 10 works well
%open(myVideo)

for n=0:nstep
    % nilai Ez di titik sumber (berubah terhadap n)
    Ezs(ys,xs) = sin(2*pi*f*n*dt); % sinusoidal
    Ezp(ys,xs) = exp(-(n-8)^2/(4^2)); % gaussian pulse
    %%% rumus leapfrog (TM)
    % *  y  *  y  *  y .... *  y  *             +---> x
    % x     x     x ....    x     x             |
    % *  y  *  y  *  y .... *  y  *             |
    % x     x     x ....    x     x             
    % *  y  *  y  *  y .... *  y  *             y
    Hxs(1:jj-1,1:ii)=cc.*Hxs(1:jj-1,1:ii)-cd/dy...
        .*(Ezs(2:jj,1:ii)-Ezs(1:jj-1,1:ii));
    Hys(1:jj,1:ii-1)=cc.*Hys(1:jj,1:ii-1)+cd/dx...
        .*(Ezs(1:jj,2:ii)-Ezs(1:jj,1:ii-1));
    Ezsbx(:,1) = Ezs(:,2); Ezsbx(:,2) = Ezs(:,ii-1);
    Ezsby(1,:) = Ezs(2,:); Ezsby(2,:) = Ezs(jj-1,:);
    Ezs(2:jj-1,2:ii-1)=ca.*Ezs(2:jj-1,2:ii-1)+cb...
        .*(((Hys(2:jj-1,2:ii-1)-Hys(2:jj-1,1:ii-2)))/dx...
        -  ((Hxs(2:jj-1,2:ii-1)-Hxs(1:jj-2,2:ii-1))/dy));
    Ezs(:,1) = Ezsbx(:,1)+((vp*dt-dx)/(vp*dt+dx))*(Ezs(:,2)-Ezs(:,1));
    Ezs(:,ii) = Ezsbx(:,2)+((vp*dt-dx)/(vp*dt+dx))*(Ezs(:,ii-1)-Ezs(:,ii));
    Ezs(1,:) = Ezsby(1,:)+((vp*dt-dy)/(vp*dt+dy))*(Ezs(2,:)-Ezs(1,:));
    Ezs(jj,:) = Ezsby(2,:)+((vp*dt-dy)/(vp*dt+dy))*(Ezs(jj-1,:)-Ezs(jj,:));
    Hxp(1:jj-1,1:ii)=cc.*Hxp(1:jj-1,1:ii)-cd/dy...
        .*(Ezp(2:jj,1:ii)-Ezp(1:jj-1,1:ii));
    Hyp(1:jj,1:ii-1)=cc.*Hyp(1:jj,1:ii-1)+cd/dx...
        .*(Ezp(1:jj,2:ii)-Ezp(1:jj,1:ii-1));
    Ezpbx(:,1) = Ezp(:,2); Ezpbx(:,2) = Ezp(:,ii-1);
    Ezpby(1,:) = Ezp(2,:); Ezpby(2,:) = Ezp(jj-1,:);
    Ezp(2:jj-1,2:ii-1)=ca.*Ezp(2:jj-1,2:ii-1)+cb...
        .*(((Hyp(2:jj-1,2:ii-1)-Hyp(2:jj-1,1:ii-2)))/dx...
        -  ((Hxp(2:jj-1,2:ii-1)-Hxp(1:jj-2,2:ii-1))/dy));
    Ezp(:,1) = Ezpbx(:,1)+((vp*dt-dx)/(vp*dt+dx))*(Ezp(:,2)-Ezp(:,1));
    Ezp(:,ii) = Ezpbx(:,2)+((vp*dt-dx)/(vp*dt+dx))*(Ezp(:,ii-1)-Ezp(:,ii));
    Ezp(1,:) = Ezpby(1,:)+((vp*dt-dy)/(vp*dt+dy))*(Ezp(2,:)-Ezp(1,:));
    Ezp(jj,:) = Ezpby(2,:)+((vp*dt-dy)/(vp*dt+dy))*(Ezp(jj-1,:)-Ezp(jj,:));
    % plot
    figure(1)
    subplot(2,2,1);
    imagesc(xe,ye,abs(Ezs),[-1 1]);
    axis('square');
    subplot(2,2,2);
    imagesc(xe,ye,Ezp,[-1 1]);
    axis('square');
    colormap('gray');
    subplot(2,2,3);
    plot(xe,Ezs(xs,:),xh,cay,'-.');
    ylabel('Ez'); xlabel('x (meter)');
    axis([-5 5 -0.4 1]);
    title(['Snapshot at t = ',num2str(t/dt),' \Deltat']);
    subplot(2,2,4);
    plot(xe,Ezp(xs,:),xh,cay,'-.');
    ylabel('Ez'); xlabel('x (meter)');
    axis([-5 5 -0.2 0.25]);
    title(['Snapshot at t = ',num2str(t/dt),' \Deltat']);
    getframe();
    % increment of time
    t = t+dt;
    % perekaman gambar sebagai frame video
    %frame = getframe(gcf); %get frame
    %writeVideo(myVideo, frame);
end

%close(myVideo)
