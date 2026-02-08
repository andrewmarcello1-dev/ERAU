close all; clear; clc

%% Load data for analysis
z = 1;
while z == 1
    z = 2;
    n = input('Enter 1 (flare) or 0 (arc): ');
    if n == 0
        data = readtable("arc.dat");
    elseif n == 1
        data = readtable("arc_flare_v2-1.dat");
    else
        fprintf('Enter 1 or 0 DUMBASS\n')
        z = 1;
    end
end
xu = table2array(data(:,1));
yu = table2array(data(:,2));

xl = xu(end:-1:1);
yl = -yu(end:-1:1);
Y = [yu;yl(2:end-1)];
X = [xu;xl(2:end-1)];

figure();
plot(X,Y)
xlabel('x (m)')
ylabel('y (m)')
xlim('auto')
ylim([-.5,.5])

Du = [xu,yu];
if n == 0
    t = [4.1,0];
    Du = [Du;t];
end
Du = Du(~isnan(Du(:,2)), :);
Du = unique(Du, 'rows');
if n == 1
    l = Du(57,:);
    m = Du(end,:);
    Du = Du(1:56,:);
    Du = [Du;m;l];
end

%% Analysis on upper airfoil

g = 1.4;
alpha = input('Angle of Attack (deg): ');

if n == 0
    chord = 4;
    cg = 2;
    R = 10;
    h = 0.202;
elseif n == 1
    chord = 2.4;
    cg = 1.7;
end

thetau = zeros(length(Du)-1,1);
thetal = zeros(length(Du)-1,1);
dely = Du(2,2)-Du(1,2);
delx = Du(2,1)-Du(1,1);
thetau(1) = atan(dely/delx)-deg2rad(alpha);
thetal(1) = atan(dely/delx)+deg2rad(alpha);
for r = 2:length(thetau)
    v1 = Du(r,:)-Du(r-1,:);
    v2 = Du(r+1,:)-Du(r,:);
    theta = acos(dot(v1,v2)/(norm(v1)*norm(v2)));
    cross_z = v1(1)*v2(2) - v1(2)*v2(1);
    if cross_z < 0
        theta = -theta;
    end
    thetau(r) = theta;
    thetal(r) = theta; % symmetry
end

M1 = 5;
M = zeros(length(Du(:,1))-1,1);
M(1) = M1;
Ku = zeros(length(thetau),1);
Ku(1) = M(1)*thetau(1);
P_ratio = zeros(length(Ku)-1,1);

for r = 2:length(Ku)+1
    if thetau(r-1) > 0
        B = sqrt(thetau(r-1)^2*((g+1)/4+sqrt(((g+1)/4)^2+1/(M(r-1)^2*thetau(r-1)^2)))^2);
        M1n = M(r-1)*sin(B);
        M2n = sqrt((M1n^2*(g-1)/2+1)/(g*M1n^2-(g-1)/2));
        M(r) = M2n/sin(B-thetau(r-1));
        if r < length(Ku)+1
            Ku(r) = M(r)*thetau(r);
        end
        P_ratio(r-1) = 1+g*(g+1)/4*Ku(r-1)^2+g*Ku(r-1)^2*sqrt(((g+1)/4)^2+1/Ku(r-1)^2);
    elseif thetau(r-1) < 0
        M(r) = M(r-1)/(1-(g-1)/2*M(r-1)*-thetau(r-1));
        if r < length(Ku)+1
            Ku(r) = M(r)*thetau(r);
        end
        P_ratio(r-1) = real((1-(g-1)/2*-Ku(r-1))^(2*g/(g-1)));
    end
end

Ml = zeros(length(Du(:,1))-1,1);
Ml(1) = M1;
Kl = zeros(length(thetal),1);
Kl(1) = Ml(1)*thetal(1);
P_ratiol = zeros(length(Kl)-1,1);
for r = 2:length(Kl)+1
    if thetal(r-1) > 0
        B = sqrt(thetal(r-1)^2*((g+1)/4+sqrt(((g+1)/4)^2+1/(Ml(r-1)^2*thetal(r-1)^2)))^2);
        M1n = Ml(r-1)*sin(B);
        M2n = sqrt((M1n^2*(g-1)/2+1)/(g*M1n^2-(g-1)/2));
        Ml(r) = M2n/sin(B-thetal(r-1));
        if r < length(Kl)+1
            Kl(r) = Ml(r)*thetal(r);
        end
        P_ratiol(r-1) = 1+g*(g+1)/4*Kl(r-1)^2+g*Kl(r-1)^2*sqrt(((g+1)/4)^2+1/Kl(r-1)^2);
    elseif thetal(r-1) < 0
        Ml(r) = Ml(r-1)/(1-(g-1)/2*Ml(r-1)*-thetal(r-1));
        if r < length(Kl)+1
            Kl(r) = Ml(r)*thetal(r);
        end
        P_ratiol(r-1) = real((1-(g-1)/2*-Kl(r-1))^(2*g/(g-1)));
    end
end

Pinf = 101325;
Pu = zeros(length(P_ratio),1);
Pu(1) = Pinf*P_ratio(1);
Pl = zeros(length(P_ratiol),1);
Pl(1) = Pinf*P_ratiol(1);

for r = 2:length(P_ratio)
    Pu(r) = P_ratio(r)*Pu(r-1);
    Pl(r) = P_ratiol(r)*Pl(r-1);
end

if real(Pu(end))<0 || real(Pl(end))<0
    Pu(end) = 1;
    fprintf('Negative Pressure, Real value enforced (1Pa)!\n')
end

Cpu = 2/(g*M1^2)*(Pu/Pinf-1);
Cpl = 2/(g*M1^2)*(Pl/Pinf-1); % symmetry

figure();
plot(Du(2:end-1,1),Cpu(1:end-1),'-k',Du(2:end-1,1),Cpl(1:end-1),'--',[0 4],[0 0],'--b')
ylabel('Coefficient of Pressure')
xlabel('Chord (m)')
title('Coefficient of Pressure Over Airfoil upper and Lower Surface')
grid on
legend('Upper','Lower','Zero')

theta_bar = zeros(length(thetau)-2,1);
cp_bar = zeros(length(Cpu)-2,1);
del = zeros(length(Du)-2,1);
theta_fs = zeros(length(xu),1);

theta_barl = zeros(length(thetal)-2,1);
cp_barl = zeros(length(Cpl)-2,1);
theta_fsl = zeros(length(xu),1);

Pressure_panel_u = zeros(length(Cpu)-2,1);
Pressure_panel_l = zeros(length(Cpu)-2,1);


for r = 1:length(Cpu)
    delx = Du(r+1,1)-Du(r,1);
    dely = Du(r+1,2)-Du(r,2);
    theta_fs(r) = atan(dely/delx)-deg2rad(alpha);
    theta_fsl(r) = atan(dely/delx)+deg2rad(alpha);
end

for r = 2:length(Cpu)
    cp_bar(r-1) = (Cpu(r-1)+Cpu(r))/2;
    theta_bar(r-1) = (theta_fs(r-1)+theta_fs(r))/2;

    cp_barl(r-1) = (Cpl(r-1)+Cpl(r))/2;
    theta_barl(r-1) = (theta_fsl(r-1)+theta_fsl(r))/2;

    delx = Du(r,1)-Du(r-1,1);
    dely = Du(r,2)-Du(r-1,2);
    del(r-1) = sqrt(delx^2+dely^2);

    Pressure_panel_u(r-1) = (Pu(r-1)+Pu(r))/2;
    Pressure_panel_l(r-1) = (Pl(r-1)+Pl(r))/2;
end


Cdu = 1/chord*sum(cp_bar.*sin(theta_bar).*del);
Cdl = 1/chord*sum(cp_barl.*sin(theta_barl).*del);
Cd = Cdu+Cdl;

Clu = 1/chord*sum(cp_bar.*cos(theta_bar).*del)
Cll = 1/chord*sum(cp_barl.*cos(theta_barl).*del);
Cl = Cll-Clu;

fprintf('Cd: %.4f\n',Cd)
fprintf('Cl: %.4f\n',Cl)

Du = Du(1:end-1,:);
panel_midpoints = (Du(1:end-1,1) + Du(2:end,1))/2;

Momsum = 1;
a = 0;
b = 4;
eps = 1e-8;
if alpha ~= 0
    vectors = Du(2:end,:)-Du(1:end-1,:);
    normals = [vectors(:,2) -vectors(:,1)];
    normals = normals./vecnorm(normals,2,2);
    midpoints = Du(1:end-1,:)+0.5*(Du(2:end,:)-Du(1:end-1,:));
    Fu = Pressure_panel_u.*normals;
    Fl = Pressure_panel_l.*normals;
    Fl(:,2) = -Fl(:,2);
    Fyu = Fu(:,2);
    Fyl = Fl(:,2);
    while abs(Momsum) > eps
        c = (a+b)/2;
        L = c - panel_midpoints;
        Momu = Fyu.*L;
        Moml = Fyl.*L;
        Momsum = sum(Momu+Moml);
        if Momsum > 0
            b = c;
        elseif Momsum < 0
            a = c;
        end
    end
    fprintf('Normalized Center of Pressure: %.3f\n',c)
    if c/chord > cg/chord
        fprintf('STABLE\n')
    elseif c/chord < cg/chord
        fprintf('UNSTABLE\n')
    end
end

