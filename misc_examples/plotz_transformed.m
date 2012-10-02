% Create a surface plot of a filter H(z)
clear all;

% Given coefficients of H(z)
Hnum = [ 1 0 ];
Hden = [ 1 -.5];

% x is real part of z, y is imaginary part of z
[x,y] = meshgrid(linspace(-1.25,1.25,1000),linspace(-1.25,1.25,1000)*1j);
z = x+y;
numOrder = length(Hnum) - 1;
denOrder = length(Hden) - 1;

wcNew = 0.1*pi;
wcOld = 0.9*pi;
wcNew2 = 0.8*pi;
ii = 0;
figure(1);
hold on;
% for fToTransformTo = linspace(wcOld, wcNew, 1)
% ii = ii + 1;
a = sin((wcOld - wcNew)/2)/sin((wcOld+wcNew)/2);

%a = cos( (wcNew2 + fToTransformTo)/2) / cos( (wcNew2 - fToTransformTo)/2 );
%B = cot( (wcNew2 - fToTransformTo)/2) * tan(wcOld/2);
%zTransformedNum = [ -1 2*a*B/(B+1) -(B-1)/(B+1) ];
%zTransformedDen = [(B-1)/(B+1) -2*a*B/(B+1) 1];
zTransformedNum = [ 1 -a];
zTransformedDen = [ -a 1];
zTransformedNumOrder = length(zTransformedNum) - 1;
zTransformedDenOrder = length(zTransformedDen) - 1;
zTransNumSum = 0;
zTransDenSum = 0;
wTransformNum = 0;
wTransformDen = 0;
w = linspace(0,2*pi,100);

dtftLinIndex = find((abs(z) > 0.999) & (abs(z) < 1.001));
dtftIndex = ind2sub(size(z),dtftLinIndex);

for nOrder = 0 : zTransformedNumOrder
   zTransNumSum = zTransNumSum + zTransformedNum(nOrder+1)*z.^(zTransformedNumOrder - nOrder);
   wTransformNum = wTransformNum + zTransformedNum(nOrder+1)*z(dtftIndex).^(zTransformedNumOrder - nOrder);
end
for dOrder = 0 : zTransformedDenOrder
   zTransDenSum = zTransDenSum + zTransformedDen(dOrder+1)*z.^(zTransformedDenOrder - dOrder);
   wTransformDen = wTransformDen + zTransformedDen(dOrder+1)*z(dtftIndex).^(zTransformedDenOrder - dOrder);
end
zTrans = zTransNumSum ./ zTransDenSum;
wTransform = wTransformNum./wTransformDen;
% plot(angle(z(dtftIndex)),angle(wTransform),'.');
plot(angle(z(dtftIndex)),angle(zTrans(dtftIndex)),'.');
return;
% Calculate H
numSum = 0;
denSum = 0;
for numPower = 0 : numOrder
    numSum = numSum + Hnum(numPower+1)*zTrans.^(numOrder - numPower); 
end

for denPower = 0 : denOrder
    denSum = denSum + Hden(denPower+1)*zTrans.^(denOrder - denPower);
end

H = numSum./denSum;

zeroIndeces = find(H == 0);

% Find coefficients of z that are at the unit circle
dtftLinIndex = find((abs(z) > 0.999) & (abs(z) < 1.001));
dtftIndex = ind2sub(size(z),dtftLinIndex);

dtft = H(dtftIndex);
w = angle(z(dtftIndex));
posFIndex = find(w>0);
[b,a] = invfreqz(dtft(posFIndex),w(posFIndex),numOrder+1,denOrder+1);
[zero, pole, gain] = nd2zpg(b, a);
%zplane(zero(ii),pole(ii));

plot(real(zero),imag(zero),'o','Color',[.1*ii 0 0]);
plot(real(pole),imag(zero),'x','Color',[.1*ii 0 0]);

% end
hold off;



%%
hDTFT = 20*log10(abs(H(dtftIndex))/max(abs(H(dtftIndex))));
cutOffIndex = find( min(abs(hDTFT(end/2:end) + 3) ) );

wcOld = angle(w(cutOffIndex));




% Plot the DTFT (2-dimensions) : |H(e^(jw))| vs w
% For some reason this skips pi * 0.1729 through pi * .8271
figure(2);
plot(angle(z(dtftIndex))/pi,20*log10(abs(H(dtftIndex))/max(abs(H(dtftIndex))) ),'.','MarkerSize',1)

% Trying to plot just the DTFT of the Z as a surface
% This does not work. How do you do it?
% figure(3);
% surface(x(dtftIndex),j*y(dtftIndex),abs(H(dtftIndex)));

% H = 0.5*(z + 1) ./ (z - 0.5);

figure(1);
surf(x,real(j*y),abs(H),'EdgeColor','interp','FaceLighting','none');
axis([-1.5 1.5 -1.5 1.5 -.5 5])
caxis([0 7])
camlight right


hold on;
% Plot the axes
[x,y] = meshgrid(zeros(1,2), linspace(-1.5,1.5,2) );
[y1,x1] = meshgrid(zeros(1,2), linspace(-1.5,1.5,2) );
z = meshgrid(linspace(0,2,2) );
alphaData = ones(size(z))*0;
surf(x,y,z,'EdgeColor','none','FaceAlpha',0.5,'FaceColor',[1,.8,.8]);
surf(x1,y1,z,'EdgeColor','none','FaceAlpha',0.5,'FaceColor',[1,.8,.8]);

[xCyl,yCyl,zCyl] = cylinder(1,50);
surf(xCyl,yCyl,zCyl*2,'EdgeColor','interp','EdgeColor','none','FaceAlpha',.5,'FaceColor',[1,.8,.8]);
hold off;
