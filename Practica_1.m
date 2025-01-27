function dXdt = Espacio_estados(t, z)
    %Parametros
    Xc = z(1);       
    dXc = z(2);      
    alpha = z(3);    
    dalpha = z(4);   
    Ip = 0.0079;     
    Mc = 0.7031;     
    lp = 0.3302;     
    Mp = 0.23;       
    Fc = 0;          
    Beq = 4.3;       
    g = 9.81;        
    Bp = 0.0024;     

    % Denominador comun
    Denominador_comun = (Mc + Mp) * Ip + Mc * Mp * lp^2 + Mp^2 * lp^2 * sin(alpha)^2;

    % Calculo de aceleraciones
    ddXc = ((Ip + Mp * lp^2) * Fc + Mp^2 * lp^2 * g * cos(alpha) * sin(alpha) - (Ip + Mp * lp^2) * Beq * dXc - (Ip * Mp * lp - Mp^2 * lp^3) * dalpha^2 * sin(alpha) - Mp * lp * dalpha * cos(alpha) * Bp) / Denominador_comun;

    ddalpha = ((Mc + Mp) * Mp * g * lp * sin(alpha) - (Mc + Mp) * Bp * dalpha + Fc * Mp * lp * cos(alpha) - Mp^2 * lp^2 * dalpha^2 * sin(alpha) * cos(alpha) - Beq * Mp * lp * dXc * cos(alpha)) / Denominador_comun;

    % Vector de salida
    dXdt = [dXc; ddXc; dalpha; ddalpha];
end

tspan = [0, 10];

% Condiciones iniciales
x0 = 0;             
dx0 = 0;            
alpha0 = deg2rad(1); 
dalpha0 = 0;        

X0 = [x0; dx0; alpha0; dalpha0]; 

[t, X] = ode45(@Espacio_estados, tspan, X0);

% Graficar resultados
figure;
subplot(2, 1, 1);
plot(t, X(:, 1), 'LineWidth', 1.5); 
xlabel('Tiempo [s]');
ylabel('Posición del carro [m]');
title('Evolución de la posición del carro');
grid on;

subplot(2, 1, 2);
plot(t, X(:, 2), 'LineWidth', 1.5); 
xlabel('Tiempo [s]');
ylabel('Ángulo del péndulo [rad]');
title('Evolución del ángulo del péndulo');
grid on;
