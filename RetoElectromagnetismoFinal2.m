%Programa de simulacion de frenado magnetico en caida libre creado por 
% Mauricio Perea González		        A01571406
% Luis Manuel González Martínez	A01722501
% Abdiel Fritsche Barajas		        A01234933

clc,clear, close all;
% Parámetros del sistema
mu0 = 4*pi*1e-7; % Permeabilidad magnética del vacío
sigma = 1e6; % Conductividad eléctrica del material conductor
R = 100; % Resistencia del circuito
m = 1000; % Masa del objeto
g = 9.81; % Aceleración gravitatoria
r = 2; % Radio de la espira
c = 0.9; %Constante de frenado
% Condiciones iniciales
z0 = 60; % Posición inicial de lanzamiento (en metros)
v0 = 0; % Velocidad inicial

teta = linspace(0,2*3.1416,360);
xcirculo = r*cos(teta);
ycirculo = r*sin(teta);
zcirculo = 0*teta;
ximan  = (r/1.25)*cos(teta);
yiman  = (r/1.25)*sin(teta);
xtorre  = (r/1.35)*cos(teta);
ytorre  = (r/1.35)*sin(teta);
ziman  = 60*(teta./teta);

xmax = 60; xmin = -1*xmax;
ymax = 60; ymin = -1*ymax;
zmax = 60; zmin = -1*zmax;

%Generacion de vectores de posicion
xmedida = linspace(xmax,xmin,60);
ymedida = linspace(ymax,ymin,60);
zmedida = linspace(zmax,zmin,60);
[CX,CY,CZ] = meshgrid(xmedida,ymedida,zmedida);


% Definición de las funciones
f_v = @(t, z, v) v;
f_a = @(t, z, v) (-m*g - ((9*sigma^2*mu0^2*r^4)/(4*R))*(z^2/(r^2+z^2))*v)/m;

% Parámetros de la simulación
dt = 0.01; % Paso de tiempo
tspan = 0:dt:10; % Intervalo de tiempo

% Inicialización de los vectores de resultados
z = zeros(size(tspan));
v = zeros(size(tspan));
a = zeros(size(tspan));

% Condiciones iniciales
z(1) = z0;
v(1) = v0;

% Implementación del método de Runge-Kutta de 4to orden
for i = 1:length(tspan)-1
    t = tspan(i);
    k1v = f_v(t, z(i), v(i));
    k1a = f_a(t, z(i), v(i));
    
    t = t + dt/2;
    k2v = f_v(t, z(i) + dt/2 * k1v, v(i) + dt/2 * k1a);
    k2a = f_a(t, z(i) + dt/2 * k1v, v(i) + dt/2 * k1a);
    
    t = t + dt/2;
    k3v = f_v(t, z(i) + dt/2 * k2v, v(i) + dt/2 * k2a);
    k3a = f_a(t, z(i) + dt/2 * k2v, v(i) + dt/2 * k2a);
    
    t = t + dt;
    k4v = f_v(t, z(i) + dt * k3v, v(i) + dt * k3a);
    k4a = f_a(t, z(i) + dt * k3v, v(i) + dt * k3a);
    
    z(i+1) = z(i) + dt/6 * (k1v + 2*k2v + 2*k3v + k4v);
    v(i+1) = v(i) + dt/6 * (k1a + 2*k2a + 2*k3a + k4a);
    
    % Frenado magnético cuando z se aproxima a cero
    if z(i+1) < 10
       v(i+1) = c * v(i+1); % Reducción del 10% en la velocidad
    end
end

a = gradient(v,tspan);
% Inicialización de gráficas
figure;
set(gcf,'Position', get(0,'Screensize'));
axis([0 11  -50 300 ])
hold on

pos_plot = animatedline('Color', 'b',LineWidth=5);
xlabel('Tiempo');
ylabel('Posición');
vel_plot = animatedline('Color', 'r',LineWidth=5);
xlabel('Tiempo');
ylabel('Velocidad');
acc_plot = animatedline('Color', 'g',LineWidth=5);
xlabel('Tiempo');
ylabel('Aceleración');
lgd = legend({'z(t)','v(t)','a(t)'},...
    'FontSize',30,'TextColor','blue');
Velocidad = text(0,100,['Velocidad  = ' num2str(v(1))]);
Posicion = text(0,120,['Posicion  = ' num2str(z(1))]);
Aceleracion = text(0,80,['Aceleracion  = ' num2str(a(1))]);
for i =1:length(tspan)
% Actualización de las gráficas animadas
    
    delete(Aceleracion); delete(Velocidad); delete(Posicion);

    addpoints(pos_plot, tspan(i), z(i));
    addpoints(vel_plot, tspan(i), v(i));
    addpoints(acc_plot, tspan(i), a(i));
    Velocidad = text(0,100,['Velocidad  = ' num2str(v(i))],"FontSize",20);
    Posicion = text(0,120,['Posicion  = ' num2str(z(i))],"FontSize",20);
    Aceleracion = text(0,80,['Aceleracion  = ' num2str(a(i))],"FontSize",20);   
    % Actualización y visualización de las gráficas
    drawnow;
    if i == length(z)
        delete(Aceleracion); delete(Velocidad); delete(Posicion);
        Velocidad = text(0,100,'Velocidad  =  0' ,"FontSize",20);
        Posicion = text(0,120,'Posicion  = 0',"FontSize",20);
        Aceleracion = text(0,80,'Aceleracion  = 0 ',"FontSize",20);
    end
    pause(0.025)
end
hold off

zvarilla = 0:1:75;
figure();
set(gcf,'Position', get(0,'Screensize'));
plot3(xcirculo,ycirculo,zcirculo,LineWidth=10,Color='black');
axis([-7 7  -7 7  0 100])
xlabel('x-axis');
ylabel('y-axis');
zlabel('z-axis');
hold on
plot3(xcirculo,ycirculo,zcirculo+75,LineWidth=10,Color='black')
p1 = plot3(xcirculo,ycirculo,zcirculo+62,LineWidth=10,Color='blue');
p2 = plot3(xcirculo,ycirculo,zcirculo+52,LineWidth=10,Color='red');
p3 = plot3(ximan,yiman,ziman,LineWidth=10,Color='red');
p4 = plot3(ximan,yiman,ziman-5,LineWidth=10,Color='blue');
for j = 0:pi/2:2*pi
    xvar = r*cos(j)*(zvarilla./zvarilla); yvar = r*sin(j)*(zvarilla./zvarilla);
    plot3(xvar,yvar,zvarilla,LineWidth=10,Color='black')  
end
for i = 1:length(z)
    delete(p1); delete(p2); delete(p3); delete(p4); delete(Aceleracion); delete(Velocidad); 
    delete(Posicion);
    p1 = plot3(xcirculo,ycirculo,zcirculo+z(i)+12.5,LineWidth=10,Color='blue');
    p2 = plot3(xcirculo,ycirculo,zcirculo+z(i)+2.5,LineWidth=10,Color='red');
    p3 = plot3(ximan,yiman,zcirculo+z(i)+7.5,LineWidth=10,Color='red');
    p4 = plot3(ximan,yiman,zcirculo+z(i)+3.5,LineWidth=10,Color='blue');
    Velocidad = text(6,0,0,['Velocidad  = ' num2str(v(i))],"FontSize",20);
    Posicion = text(6,0,10,['Posicion  = ' num2str(z(i))],"FontSize",20);
    Aceleracion = text(6,0,20,['Aceleracion  = ' num2str(a(i))],"FontSize",20);
    pause(0.005)
    if i == length(z)
        delete(Aceleracion); delete(Velocidad); delete(Posicion);
        Velocidad = text(6,0,0,'Velocidad  =  0' ,"FontSize",20);
        Posicion = text(6,0,10,'Posicion  = 0',"FontSize",20);
        Aceleracion = text(6,0,20,'Aceleracion  = 0 ',"FontSize",20);
    end
    drawnow
end
hold off
