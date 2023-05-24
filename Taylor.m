clear;
%set up variables
t0 = 0; %initial time
tf = 1; %final time
h = 0.2; %time-step
Nit = (tf-t0)/h; %numerb of steps

%setting up variables
t = zeros(1,Nit+1); %time
x = zeros(1,Nit+1); %for TS(3)
y= zeros(1,Nit+1); %for TS(4)
  
%exact solution
exact = sin(1) + cos(1) + (1.^(3))/3 + 2;
%setting up initial conditions
t(1) = 0;
x(1) = 3;
y(1) = 3;
%Euler method
for i = 1 : Nit
   x(i+1) = x(i) + h*(cos(t(i)) - sin(t(i)) + t(i)^2) + ((h^2)/2)*(-sin(t(i)) - cos(t(i)) + 2*t(i)) + ((h^3)/6) * (-cos(t(i)) + sin(t(i)) + 2); %formula-for other methods you will change here
   y(i+1) = y(i) + h*(cos(t(i)) - sin(t(i)) + t(i)^2) + ((h^2)/2)*(-sin(t(i)) - cos(t(i)) + 2*t(i)) + ((h^3)/6) * (-cos(t(i)) + sin(t(i)) + 2) + ((h^4)/24) * (sin(t(i)) + cos(t(i))) ; %formula-for other methods you will change here
   t(i+1) = i*h; %time
end
%print of regular methods
fprintf('h\t|Answers:TS(3)\t TS(4) | GE:TS(3)\t    TS(4)\t    | GE for TS(3)/h^3 | GE for TS(4)/h^4\n');
fprintf('0.20|\t%d\t   %d  | %d\t%d| \t%d\t   | \t%d\n', x(end), y(end), x(end)-exact, y(end)-exact, (x(end)-exact)/(h^3), (y(end)-exact)/(h^4));
%graph
t1=[0 .2 .4 .6 .8 1];
u=@(t1)sin(t) + cos(t1) + ((t1).^(3))/3 + 2;

figure;
plot(t,x,'y-',"linewidth",8)
hold on
plot(t,y,'b--',"linewidth",6)
plot(t1,u(t1),"linewidth",4)
hold off


%set up variables

%now h=0.1!
h = 0.1; %time-step
Nit = (tf-t0)/h; %numerb of steps

%setting up variables
t = zeros(1,Nit+1); %time
x = zeros(1,Nit+1); %function x

%exact solution
exact = sin(1) + cos(1) + (1^(3))/3 + 2;

%setting up initial conditions
t(1) = 0;
x(1) = 3;
y(1) = 3;

%Euler method
for i = 1 : Nit
   x(i+1) = x(i) + h*(cos(t(i)) - sin(t(i)) + t(i)^2) + ((h^2)/2)*(-sin(t(i)) - cos(t(i)) + 2*t(i)) + ((h^3)/6) * (-cos(t(i)) + sin(t(i)) + 2); %formula-for other methods you will change here
   y(i+1) = y(i) + h*(cos(t(i)) - sin(t(i)) + t(i)^2) + ((h^2)/2)*(-sin(t(i)) - cos(t(i)) + 2*t(i)) + ((h^3)/6) * (-cos(t(i)) + sin(t(i)) + 2) + ((h^4)/24) * (sin(t(i)) + cos(t(i))); %formula-for other methods you will change here
   t(i+1) = i*h; %time
end
fprintf('0.10|\t%d\t   %d | %d    %d| \t%d\t   | \t%d\n', x(end), y(end), x(end)-exact, y(end)-exact, (x(end)-exact)/(h^3), (y(end)-exact)/(h^4));

t1=[0 .1 .2 .3 .4 .5 .6 .7 .8 .9 1];
u=@(t1)sin(t) + cos(t1) + ((t1).^(3))/3 + 2;

figure;
plot(t,x,'y-',"linewidth",8)
hold on
plot(t,y,'b--',"linewidth",6)
plot(t1,u(t1),"linewidth",4)
hold off


%now h=0.1!
h = 0.05; %time-step
Nit = (tf-t0)/h; %numerb of steps

%setting up variables
t = zeros(1,Nit+1); %time
x = zeros(1,Nit+1); %function x

%exact solution
exact = sin(1) + cos(1) + (1^(3))/3 + 2;

%setting up initial conditions
t(1) = 0;
x(1) = 3;
y(1) = 3;

%Euler method
for i = 1 : Nit
   x(i+1) = x(i) + h*(cos(t(i)) - sin(t(i)) + t(i)^2) + ((h^2)/2)*(-sin(t(i)) - cos(t(i)) + 2*t(i)) + ((h^3)/6) * (-cos(t(i)) + sin(t(i)) + 2); %formula-for other methods you will change here
   y(i+1) = y(i) + h*(cos(t(i)) - sin(t(i)) + t(i)^2) + ((h^2)/2)*(-sin(t(i)) - cos(t(i)) + 2*t(i)) + ((h^3)/6) * (-cos(t(i)) + sin(t(i)) + 2) + ((h^4)/24) * (sin(t(i)) + cos(t(i))); %formula-for other methods you will change here
   t(i+1) = i*h; %time
end
fprintf('0.05|\t%d\t   %d | %d   %d| \t%d\t   | \t%d', x(end), y(end), x(end)-exact, y(end)-exact, (x(end)-exact)/(h^3), (y(end)-exact)/(h^4));
t1=[0 .05 .1 .15 .2 .25 .3 .35 .4 .45 .5 .55 .6 .65 .7 .75 .8 .85 .9 .95 1];
u=@(t1)sin(t) + cos(t1) + ((t1).^(3))/3 + 2;

figure;
plot(t,x,'y-',"linewidth",8)
hold on
plot(t,y,'b--',"linewidth",6)
plot(t1,u(t1),"linewidth",4)
hold off
