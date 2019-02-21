clear
close all
clc

global time enteredTime ball ball2 ball3 equation t x lineOfVelocity
global gravitationalAcceleration initialVelocity lineOfHeight
format long

x = 0:13;
y = [202.36 239.03 280.71 309.12 323.15 332.78 328.45 306.40 287.36 247.97 202.89 161.11 93.68 20.78];

n = length(x);

q = zeros(3);
q(3,3) = n;
p = zeros(1,3);

for i = 1:n
    q(1,1) = q(1,1) + x(i)^4;
    q(1,2) = q(1,2) + x(i)^3;
    q(1,3) = q(1,3) + x(i)^2;
    q(2,1) = q(2,1) + x(i)^3;
    q(2,2) = q(2,2) + x(i)^2;
    q(2,3) = q(2,3) + x(i)^1;
    q(3,1) = q(3,1) + x(i)^2;
    q(3,2) = q(3,2) + x(i)^1;
    
    p(1,1) = p(1,1) + x(i)^2 * y(i)^1;
    p(1,2) = p(1,2) + x(i)^1 * y(i)^1;
    p(1,3) = p(1,3) + y(i);
end

equation = p / q; % p*inv(q);
e = vpa(poly2sym(round(equation, 3)));
disp('Equation:'); 
pretty(e);


SSE = 0;
SST = 0;
for i = 1:n
    SSE = SSE + (y(i)-equation(1)*x(i)^2-equation(2)*x(i)-equation(3))^2;
    SST = SST + (y(i) - mean(y))^2;
end

R2 = 1 - (SSE / SST);
disp(['Relative Estimation Power: ', num2str(R2)]);

maximumHeightTime = -equation(2) / (2 * equation(1));
maximumHeight = polyval(equation, maximumHeightTime);
displacement = maximumHeight - polyval(equation, 0);
initialVelocity = displacement * 2 / maximumHeightTime;
disp(['Initial Velocity: ', num2str(initialVelocity), ' meter/seconds']);
gravitationalAcceleration = -initialVelocity / maximumHeightTime;
disp(['Gravitational Acceleration : ', num2str(gravitationalAcceleration), ' meter/seconds^2']);
fallTime = sqrt(2*maximumHeight/-gravitationalAcceleration) + maximumHeightTime;
finalVelocity = gravitationalAcceleration * fallTime + initialVelocity;
rootsOfEquation = roots(equation);

disp(['Ýnitial Height: ', num2str(polyval(equation, 0)), ' meter' ]);
disp(['Maximum Height: ', num2str(maximumHeight), ' meter']);
disp(['Final Velocity: ', num2str(finalVelocity), ' meter/seconds']);
disp(['Calculated Fall Time: ', num2str(fallTime), ' seconds']);

tic
subplot(2,2,[1,3]);
whitebg('k')
ball = animatedline('Color', 'w', 'Marker', 'o', 'MarkerSize', 10, 'MarkerFaceColor', 'w');
ball2 = animatedline('Color', 'k', 'Marker', 'o', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'LineStyle', 'none');
ball3 = animatedline('Color', 'k', 'Marker', 'o', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
axis([0, 10, 0, maximumHeight + 10]);
title('Ball');
ylabel('Height(m)');

subplot(2,2,2);
lineOfVelocity = animatedline('Color', 'w', 'LineWidth', 3);
axis([0 fallTime finalVelocity initialVelocity]);
ylabel('Height(m)');
xlabel('Time(s)');
title('Time-Velocity');

subplot(2,2,4);
lineOfHeight = animatedline('Color', 'w', 'LineWidth', 3);
axis([0 fallTime 0 maximumHeight]);
ylabel('Velocity(m/s)');
xlabel('Time(s)');
title('Time-Height');

time = 0;
enteredTime = 0;

t = timer;
t.TimerFcn = @timerCallback;
t.Period = 0.01;
t.TasksToExecute = round(1/t.Period * rootsOfEquation(1));
t.ExecutionMode = 'fixedRate';

addpoints(ball, 5, polyval(equation, time));
addpoints(lineOfVelocity, time, gravitationalAcceleration * time + initialVelocity);

start(t);
wait(t);
addpoints(ball, 5, 5);

toc

function timerCallback(~, ~)
    global time enteredTime ball ball2 ball3 equation t x lineOfVelocity
    global gravitationalAcceleration initialVelocity lineOfHeight
    
    time = time + t.Period;
    
    instantaneousHeight = polyval(equation, time);
    if instantaneousHeight < 5
        instantaneousHeight = 5;
    end
    
    for i = 1:14
        if time >= x(i) - 0.01 && time <= x(i) + 0.01
            height = polyval(equation, x(i));
            if time <= 5 
                addpoints(ball2, 7, height);
            else
                addpoints(ball3, 8, height);
            end
                enteredTime = time;
        end
    end
    
%     if enteredTime + 0.5 < time
%         clearpoints(ball2);
%         clearpoints(ball3);
%     end
    
    clearpoints(ball);
    addpoints(ball, 5, instantaneousHeight);
    
    addpoints(lineOfVelocity, time, gravitationalAcceleration * time + initialVelocity);
    
    addpoints(lineOfHeight, time, instantaneousHeight);
    
end % End of timerCallback




















