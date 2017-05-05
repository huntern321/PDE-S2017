%2D Heat Equation with Forced Boundry Conditions = 0
%Hunter Normandeau & Juan Carlos del Rio

res=100; % XY resolution
timeStep=100; %Time resolution
loops=20; %Number of fourier iterations
Xlength=2; % X spacial size
Ylength=2; % Y spacial size
endTime=3600; % Length of Simulation
runtime=5; %Animation runtime in seconds

%Calculate Constant
SpecificHeat=475; %J/kgK
ThermalConductivity=44.5; %W/(M*K)
Density=7.85; %g/cm^3
ThermalDiffusivity=ThermalConductivity/(SpecificHeat*Density);
c=ThermalDiffusivity/4; %constant

% Parameterizing dimensions of the frame 0-size
a=linspace(0,Xlength,res);
b=linspace(0,Ylength,res);
z=zeros(length(a),length(b));
lambda=zeros(loops,loops);
time = linspace(0,endTime,timeStep);
xbdy = linspace(0,Xlength,res);
ybdy = linspace(0,Ylength,res);
bnm=zeros(loops,loops);

%Initial Conditions
ic = @(x,y) 100;

for m=1:loops
    for n=1:loops

% Calculate Lambda
lambda(m,n)=(pi*c)*sqrt((m^2)/(a(end)^2)+(n^2)/b(end)^2);

% Calculate fourier coefficients
fun = @(x,y) (ic(x,y)).*sin((n*pi/b(end)).*y).*sin((m*pi/a(end)).*x);
q = integral2(fun,0,a(end),0,b(end));

bnm(n,m) = 4/(a(end)*b(end)).*q;
    end
end

Uxyt=zeros(length(a)+1,length(b)+1,length(time));

xind=0; %index for x loop
yind=0; %index for y loop
tind=0; %index for time loop

for t = 0:time(end)/timeStep:time(end)-(1/timeStep)
    tind=int16(t/(time(end)/timeStep)+1);
    
for x = 0:a(end)/res:a(end)
    xind=int16(x*res/a(end)+1);
    
for y = 0:b(end)/res:b(end)
    yind=int16(y*res/b(end)+1);
    
for m=1:loops
for n=1:loops
   
Uxyt(xind,yind,tind) = Uxyt(xind,yind,tind) + bnm(n,m).*sin((n*pi/b(end)).*y).*sin((m*pi/a(end)).*x).*exp(-lambda(n,m).^2.*t);

end
end

end
end

end
a(end+1) = a(end)+(x/res);
b(end+1) = a(end)+(y/res);

for t=1:timeStep
surf(a,b,Uxyt(:,:,t),'EdgeColor','none');
%urf(Uxyt(:,t,t),Uxyt(t,:,t),Uxyt(t,t,:),'EdgeColor','none');
colorbar
xlabel('X')
ylabel('Y')
zlabel('Z')
axis([0 x+(x/res) 0 y+(y/res) 0 110])
view(40,30)
pause(runtime/timeStep);

% for t=1:timeStep
% contour=contourf(a,b,Uxyt(:,:,t),'EdgeColor','none');
% colorbar
% xlabel('X')
% ylabel('Y')
% pause(runtime/timeStep);

end

