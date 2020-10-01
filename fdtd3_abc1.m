clear
clc

mu0 = pi*4e-7; vp = 299792458; ep0 = 1/(vp^2*mu0);
mur = 1; epr = 1;
sigm = 0; sige = 0;
dx = 0.108;
dy = 0.108;
dz = 0.108;
dt = dx/vp/sqrt(3);

xmin = -5.4;  xmax = 5.4;
ymin = -5.4;  ymax = 5.4;
zmin = -5.4;  zmax = 5.4;
tmax = 150*dt;
t = -tmax:dt:tmax;

f = 278e6;

xs = 51;
ys = 51;
zs = 51;
x = (xmin:dx:xmax);
y = (ymin:dy:ymax);
z = (zmin:dz:zmax);
i = length(x);
j = length(y);
k = length(z);

ca = (2*epr*ep0-sige*dt)/(2*epr*ep0+sige*dt);
cb = 2*dt/(2*epr*ep0+sige*dt);
cc = (2*mur*mu0-sigm*dt)/(2*mur*mu0+sigm*dt);
cd = 2*dt/(2*mur*mu0+sigm*dt);

% pemanggilan (y,x,z) (baris,kolom,tinggi)
Hx(1:j-1,1:i,1:k-1)=0.0;
Hy(1:j,1:i-1,1:k-1)=0.0;
Hz(1:j-1,1:i-1,1:k)=0.0;
Ex(1:j,1:i-1,1:k)=0.0;
Ey(1:j-1,1:i,1:k)=0.0;
Ez(1:j,1:i,1:k-1)=0.0;

myVideo = VideoWriter('fdtd3');
open(myVideo);

for n=0:length(t)
    Ez(ys,xs,zs) = sin(2*pi*f*n*dt);
    Hx(1:j-1,1:i,1:k-1)=cc.*Hx(1:j-1,1:i,1:k-1)...
        +cd.*(((Ey(1:j-1,1:i,2:k)-Ey(1:j-1,1:i,1:k-1))/dz)...
        -((Ez(2:j,1:i,1:k-1)-Ez(1:j-1,1:i,1:k-1))/dy));
    Hy(1:j,1:i-1,1:k-1)=cc.*Hy(1:j,1:i-1,1:k-1)...
        +cd.*(((Ez(1:j,2:i,1:k-1)-Ez(1:j,1:i-1,1:k-1))/dx)...
        -((Ex(1:j,1:i-1,2:k)-Ex(1:j,1:i-1,1:k-1))/dz));
    Hz(1:j-1,1:i-1,1:k)=cc.*Hz(1:j-1,1:i-1,1:k)...
        +cd.*(((Ex(2:j,1:i-1,1:k)-Ex(1:j-1,1:i-1,1:k))/dy)...
        -((Ey(1:j-1,2:i,1:k)-Ey(1:j-1,1:i-1,1:k))/dx));
    % preparation for first-order Mur boundary calculation
    Ezbx(:,1,:) = Ez(:,2,:); Ezbx(:,2,:) = Ez(:,i-1,:);
    Ezby(1,:,:) = Ez(2,:,:); Ezby(2,:,:) = Ez(j-1,:,:);
    Ezbz(:,:,1) = Ez(:,:,2); Ezbz(:,:,2) = Ez(:,:,k-2);
    Ex(2:j-1,1:i-1,2:k-1)=ca.*Ex(2:j-1,1:i-1,2:k-1)...
        +cb.*(((Hz(2:j-1,1:i-1,2:k-1)-Hz(1:j-2,1:i-1,2:k-1))/dy)...
        -((Hy(2:j-1,1:i-1,2:k-1)-Hy(2:j-1,1:i-1,1:k-2))/dz));
    Ey(1:j-1,2:i-1,2:k-1)=ca.*Ey(1:j-1,2:i-1,2:k-1)...
        +cb.*(((Hx(1:j-1,2:i-1,2:k-1)-Hx(1:j-1,2:i-1,1:k-2))/dz)...
        -((Hz(1:j-1,2:i-1,2:k-1)-Hz(1:j-1,1:i-2,2:k-1))/dx));
    Ez(2:j-1,2:i-1,1:k-1)=ca.*Ez(2:j-1,2:i-1,1:k-1)...
        +cb.*(((Hy(2:j-1,2:i-1,1:k-1)-Hy(2:j-1,1:i-2,1:k-1))/dx)...
        -((Hx(2:j-1,2:i-1,1:k-1)-Hx(1:j-2,2:i-1,1:k-1))/dy));
    % First-order Mur boundary
    Ez(:,1,:) = Ezbx(:,1,:)+((vp*dt-dx)/(vp*dt+dx))*(Ez(:,2,:)-Ez(:,1,:));
    Ez(:,i,:) = Ezbx(:,2,:)+((vp*dt-dx)/(vp*dt+dx))*(Ez(:,i-1,:)-Ez(:,i,:));
    Ez(1,:,:) = Ezby(1,:,:)+((vp*dt-dy)/(vp*dt+dy))*(Ez(2,:,:)-Ez(1,:,:));
    Ez(j,:,:) = Ezby(2,:,:)+((vp*dt-dy)/(vp*dt+dy))*(Ez(j-1,:,:)-Ez(j,:,:));
    Ez(:,:,1) = Ezbz(:,:,1)+((vp*dt-dz)/(vp*dt+dz))*(Ez(:,:,2)-Ez(:,:,1));
    Ez(:,:,k-1) = Ezbz(:,:,2)+((vp*dt-dz)/(vp*dt+dz))*(Ez(:,:,k-2)-Ez(:,:,k-1));
    figure(1)
    subplot(2,2,1);
    slice(x,y,z(1:100),(abs(Ez)).^(0.5),x(xs),y(ys),z(zs));
    axis tight;
    xlabel('x (meter)');
    ylabel('y (meter)');
    zlabel('z (meter)');
    shading interp;
    colormap('hot');
    %colorbar;
    caxis([0 1]);
    subplot(2,2,2);
    imagesc(y,z(1:100),transpose(squeeze(abs(Ez(:,xs,:)).^(0.5))),[0 1]);
    axis('square');
    xlabel('y (meter)');
    ylabel('z (meter)');
    subplot(2,2,4);
    imagesc(x,y,abs((Ez(:,:,zs)).^(0.5)),[0 1]);
    axis('square');
    xlabel('x (meter)');
    ylabel('y (meter)');
    subplot(2,2,3);
    imagesc(x,z(1:100),transpose(squeeze(abs(Ez(ys,:,:)).^(0.5))),[0 1]);
    axis('square');
    xlabel('x (meter)');
    ylabel('z (meter)');
    frame = getframe(gcf);
    writeVideo(myVideo, frame);
end

close(myVideo);
