%Programa de simulacion de campo magnetico en una espira creado por 
% Mauricio Perea González		        A01571406
% Luis Manuel González Martínez	A01722501
% Abdiel Fritsche Barajas		        A01234933

clc,clear, close all;
%Pedir Datos al usuario para la generacion de los vectores de posicion 
% asi como datos de corriente, radio del anillo etc.
fprintf('El rango de valores en los ejes X Y Z deben de ser iguales: \n');
xmax = input('Dame el valor de puntos en el Eje X: ');  xmin = xmax*-1;
ymax = input('Dame el valor de puntos en el Eje Y: ');  ymin = ymax *-1;
zmax = input('Dame el valor de puntos en el Eje Z: ');  zmin = zmax *-1;
n=input('Dame el numero de particiones del anillo que quieres: ');
dx = input('Dame el numero de particiones en los Ejes: ');
I = input('Dame el valor de la corriente que circula en el anillo: ');
R = input('Dame el valor del radio de la espira: ');

%Generacion de vectores de posicion
x = linspace(xmax,xmin,dx);
y = linspace(ymax,ymin,dx);
z = linspace(zmax,zmin,dx);
teta = linspace(0,2*3.1416,360);

%Generacion de coordenadas polares
xcirculo = R*cos(teta);
ycirculo = R*sin(teta);
zcirculo = 0*teta;
[CX,CY,CZ] = meshgrid(x,y,z);

%Creacion de matrices del campo magnetico del tamaño de los vectores de
%posicion
Bx = zeros(size(CX));
By = zeros(size(CY));
Bz = zeros(size(CZ));
Bxreal = zeros(size(CX));
Byreal = zeros(size(CY));
Bzreal = zeros(size(CZ));
%Ciclo que envia cada valor del Meshgrid para calcular su campo magnetico
%individualmente
for i = 1:1:numel(CX)
        campo_i = getCampoMagneticoEst("contrareloj",I,R,n,CX(i),CY(i),CZ(i));
        Bx(i) = campo_i(1);
        By(i) = campo_i(2);
        Bz(i) = campo_i(3);
end

for i = 1:1:numel(CX)
        campo_ireal = getCampoMagneticoReal("contrareloj",I,R,n,CY(i),CZ(i));
        Bxreal(i) = campo_ireal(1);
        Byreal(i) = campo_ireal(2);
        Bzreal(i) = campo_ireal(3);
end
%Normalizacion de las componentes
B = sqrt(Bx.^2+By.^2+Bz.^2);
Breal = sqrt(Bxreal.^2+Byreal.^2+Bzreal.^2);
a = 1;

fileID = fopen('DatosCampoMagnetico.txt','w');
fprintf(fileID,'Componente en X del Campo Magnetico %0.0f\n',a);
for i = 1:1:numel(CX)
    diferencia = abs(Bxreal(i)-Bx(i));
fprintf(fileID,['Componente X Real: %3.4eT  Componente X Estimada: %3.4eT ' ...
            ' En el punto  X: %3.1f Y: %3.1f Z: %3.1f  Con una diferencia de: %3.4eT\n'], ...
            Bxreal(i),Bx(i),CX(i),CY(i),CZ(i),diferencia);
end
fprintf(fileID,'\nComponente en Y del Campo Magnetico %0.0f\n\n',a);
for i = 1:1:numel(CX)
     diferencia = abs(Byreal(i)-By(i));
fprintf(fileID,['Componente Y Real: %3.4eT  Componente Y Estimada: %3.4eT ' ...
            ' En el punto  X: %3.1f Y: %3.1f Z: %3.1f Con una diferencia de: %3.4eT\n'], ...
            Byreal(i),By(i),CX(i),CY(i),CZ(i),diferencia);
end
fprintf(fileID,'\nComponente en Z del Campo Magnetico %0.0f\n\n',a);
for i = 1:1:numel(CX)
    diferencia = abs(Bzreal(i)-Bz(i));
fprintf(fileID,['Componente Z Real: %3.4eT  Componente Z Estimada: %3.4eT ' ...
            ' En el punto  X: %3.1f Y: %3.1f Z: %3.1f Con una diferencia de: %3.4eT\n'], ...
            Bzreal(i),Bz(i),CX(i),CY(i),CZ(i),diferencia);
end
fclose(fileID);
%Grafico del campo en 3D
figure()
quiver3(CX,CY,CZ,Bx./B,By./B,Bz./B);
hold on
plot3(xcirculo,ycirculo,zcirculo,LineWidth=30)
xlabel('x-axis');
ylabel('y-axis');
zlabel('z-axis');
grid on
title('Campo Magnetico Estimado 3D')
hold off

%Grafico del campo en 2D
figure()
quiver(CY,CZ,By./B,Bz./B)
hold on
plot([R,-1*R],[0,0],LineWidth=10)
hold off
title('Campo Magnetico Estimado 2D')

figure()
quiver3(CX,CY,CZ,Bxreal./Breal,Byreal./Breal,Bzreal./Breal);
hold on
plot3(xcirculo,ycirculo,zcirculo,LineWidth=30)
xlabel('x-axis');
ylabel('y-axis');
zlabel('z-axis');
grid on
title('Campo Magnetico Real 3D')
hold off

%Funcion que calcula el campo magnetico en cualquier punto para una espira
%que circula corriente
function campo = getCampoMagneticoEst(sentido,I, R,n,x, y, z)
 mu0 = 4*pi *10^(-7);
 puntoCruz= [0,0,0];

 %Ciclo iterativo que sirve para resolver de manera numerica la integral de
 %la ecuacion de Biot-Savart
 for i = 1:n 
         delta_teta = (2*pi)/n;
         tetan =(i-1)*delta_teta;
         r_rq = [x-R*cos(tetan),y-R*sin(tetan),z];
         r_rqMagnitud = norm(r_rq);
         r_rqMagnitudCubica = r_rqMagnitud^3;
    
         if sentido == "contrareloj"
                ds = [-sin(tetan), cos(tetan),0];
         elseif sentido == "reloj"
                ds = [sin(tetan),-cos(tetan),0];
         end

         puntoCruz = puntoCruz + cross(ds,(r_rq./(n*r_rqMagnitudCubica)).*delta_teta);
 end

 campo = ((mu0*R*I)/(4*pi))*puntoCruz;
end

function campo = getCampoMagneticoReal(sentido,I,R,n, y, z)
        mu0 = 4*pi *10^(-7);
        intBx=0;
        intBy=0;
        intBz=0;
        for i = 1:n 
             delta_teta = (2*pi)/n;
             tetan =(i-1)*delta_teta;
             D = ((R^2+y^2+z^2)-(2*R*y*sin(tetan)))^(3/2); 
             intBx = intBx +(cos(tetan)/D)*delta_teta;
             intBy = intBy +(sin(tetan)/D)*delta_teta;
             intBz = intBz +((R-(y*sin(tetan)))/D)*delta_teta;
        end
        if sentido =="contrareloj"
            B_x = ((mu0*R*I)/(4*pi))*z*intBx;
            B_y = ((mu0*R*I)/(4*pi))*z*intBy;
            B_z = ((mu0*R*I)/(4*pi))*intBz;
        elseif sentido =="reloj"
            B_x = -((mu0*R*I)/(4*pi))*z*intBx;
            B_y = -((mu0*R*I)/(4*pi))*z*intBy;
            B_z = -((mu0*R*I)/(4*pi))*intBz;
        end    
        campo = [B_x,B_y,B_z];
end
