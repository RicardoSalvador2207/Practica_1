function Practica_1
    tspan = [0 10];
    x0 = [0, 0, deg2rad(1), 0];  

    % Resolver el sistema de ecuaciones diferenciales
    [t, x] = ode45(@pendulo, tspan, x0);

    % Graficar los resultados
    figure;
    plot(t, x);
    xlabel('Tiempo (s)');
    ylabel('Estados');
    legend('Xc', 'dXc', 'A', 'dA');
    title('Simulación del péndulo invertido');
    grid on;
end

function dx = pendulo(~, x)
    %Definición de parámetros 
    Ip = 0.0079;   % Momento de inercia del péndulo
    Mc = 0.7031;   % Masa del carrito
    lp = 0.3302;   % Longitud del péndulo
    Mp = 0.23;     % Masa del péndulo
    Fc = 0;        % Fuerza del motor
    Beq = 4.3;     % Coeficiente de amortiguamiento eq
    g = 9.81;      % Gravedad
    Bp = 0.0024;   % Coeficiente de amortiguamiento


    % Variables de estado
    x1 = x(1); % Xc (posición del carrito)
    x2 = x(2); % dXc (velocidad del carrito)
    x3 = x(3); % A (ángulo del péndulo)
    x4 = x(4); % dA (velocidad angular del péndulo)


    % Denominador común
    denominador_comun = (Mc + Mp) * Ip + Mc * Mp * lp^2 + Mp^2 * lp^2 * sin(x3)^2;

    % Sistema de ecuaciones diferenciales
    ddXc = ((Ip + Mp * lp^2) * Fc + Mp^2 * lp^2 * g * cos(x3) * sin(x3) - ...
            (Ip + Mp * lp^2) * Beq * x2 - (Ip * Mp * lp - Mp^2 * lp^3) * x4^2 * sin(x3) - ...
            Mp * lp * x4 * cos(x3) * Bp) / denominador_comun;

    dda =   ((Mc + Mp) * Mp * g * lp * sin(x3) - (Mc + Mp) * Bp * x4 + ...
            Fc * Mp * lp * cos(x3) - Mp^2 * lp^2 * x4^2 * sin(x3) * cos(x3) - ...
            Beq * Mp * lp * x2 * cos(x3)) / denominador_comun;

    dx = zeros(4,1);
    dx(1) = x2;
    dx(2) = ddXc;
    dx(3) = x4;
    dx(4) = dda;
end
