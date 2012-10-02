% Create a surface plot of a filter H(z)
clear all;

% Given coefficients of H(z)
Hnum = [ 1 0 ];
Hden = [ 1 -0.5];

% x is real part of z, y is imaginary part of z
[x,y] = meshgrid(linspace(-1.25,1.25,1000),linspace(-1.25,1.25,1000)*1j);
z = x+y;
numOrder = length(Hnum) - 1;
denOrder = length(Hden) - 1;

% Calculate H
numSum = 0;
denSum = 0;
for numPower = 0 : numOrder
    numSum = numSum + Hnum(numPower+1)*z.^(numOrder - numPower); 
end

for denPower = 0 : denOrder
    denSum = denSum + Hden(denPower+1)*z.^(denOrder - denPower);
end

H = numSum./denSum;
 
% Find coefficients of z that are at the unit circle
dtftLinIndex = find( abs(abs(z) -1) < .001  );
dtftIndex = ind2sub(size(z),dtftLinIndex);
notDTFTIndex = find(abs(z) > 1.025 | abs(z) < .975);
% Plot the DTFT (2-dimensions) : |H(e^(jw))| vs w
% For some reason this skips pi * 0.1729 through pi * .8271
figure(2);
plot(angle(z(dtftIndex))/pi,20*log10(abs(H(dtftIndex))),'.','MarkerSize',1)


% Trying to plot just the DTFT of the Z as a surface
% This does not work. How do you do it?
R = zeros(size(z));
R(dtftIndex) = abs(H(dtftIndex));
R(notDTFTIndex) = NaN;
%R = NaN(size(z));

figure(3);
surf(x,real(j*y),R,'EdgeColor','interp');

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
