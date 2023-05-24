clear
%set up variables
t0 = 0; %initial time
tf = 2; %final time
h = 0.25; %time-step
Nit = (tf-t0)/h; %numerb of steps

%setting up variables
t = zeros(1,Nit+1); %time
x = zeros(1,Nit+1); %for Method 1
y= zeros(1,Nit+1); %for Method 2

%exact solution
exact = exp(-40);

%setting up initial conditions
t(1) = 0;
x(1) = 1;
y(1) = 1;
for i = 1 : Nit

  %method 2
  k1 = -20 * x(i);
  k2 = -20 * (x(i) + 0.5*h*k1);
  k3 = -20 * (x(i) - (1*h*k1) + (2*h*k2));
  x(i+1) = x(i) + h*(k1*(1/6) + k2*(2/3) + k3*(1/6));
  %method 1
  k4 = -20 * y(i);
  k5 = -20 * (y(i) + 0.5*h*k4);
  k6 = -20 * (y(i) + 0.5*h*k5);  
  k7 = -20 * (y(i) + h*k6);
  y(i+1) = y(i) + h*( k4*(1/6) + k5*(1/3) + k6*(1/3) + k7*(1/6));
  t(i+1) = i*h; %time
end
fprintf('h\t  |Answers:RK(3)\t \RK(4)\t\t | GE:RK(3)\t \tRK(4)\n');
fprintf('0.25  |\t%d\t     %d | %d\t%d   \n', x(end), y(end), x(end)-exact, y(end)-exact);

t1=[0.00 0.25 0.50 0.75 1.00 1.25 1.50 1.75 2];
u=@(t1) 1*exp(-20*(t1));

figure;
plot(t,x,'y-',"linewidth",8)
hold on
plot(t1,u(t1),"linewidth",5)
plot(t,y,'b--',"linewidth",6)
hold off

%set up variables
t0 = 0; %initial time
tf = 2; %final time
h = 0.125; %time-step
Nit = (tf-t0)/h; %numerb of steps

%setting up variables
t = zeros(1,Nit+1); %time
x = zeros(1,Nit+1); %for Method 1
y= zeros(1,Nit+1); %for Method 2

%exact solution
exact = exp(-40);

%setting up initial conditions
t(1) = 0;
x(1) = 1;
y(1) = 1;
for i = 1 : Nit

  %method 2
  k1 = -20 * x(i);
  k2 = -20 * (x(i) + 0.5*h*k1);
  k3 = -20 * (x(i) - (1*h*k1) + (2*h*k2));
  x(i+1) = x(i) + h*(k1*(1/6) + k2*(2/3) + k3*(1/6));
  %method 1
  k4 = -20 * y(i);
  k5 = -20 * (y(i) + 0.5*h*k4);
  k6 = -20 * (y(i) + 0.5*h*k5);  
  k7 = -20 * (y(i) + h*k6);
  y(i+1) = y(i) + h*( k4*(1/6) + k5*(1/3) + k6*(1/3) + k7*(1/6));
  t(i+1) = i*h; %time
end

fprintf('0.125 |\t%d\t     %d  | %d\t    %d   \n', x(end), y(end), x(end)-exact, y(end)-exact);

t1=[0 0.125 0.25 0.375 0.5 0.625 0.75 0.875 1.0 1.125 1.25 1.375 1.5 1.625 1.75 1.875 2.0];
u=@(t1) 1*exp(-20*(t1));

figure;
plot(t,x,'y-',"linewidth",8)
hold on
plot(t1,u(t1),"linewidth",5)
plot(t,y,'b--',"linewidth",6)
hold off

%set up variables
t0 = 0; %initial time
tf = 2; %final time
h = 0.0625; %time-step
Nit = (tf-t0)/h; %numerb of steps

%setting up variables
t = zeros(1,Nit+1); %time
x = zeros(1,Nit+1); %for Method 1
y= zeros(1,Nit+1); %for Method 2

%exact solution
exact = exp(-40);

%setting up initial conditions
t(1) = 0;
x(1) = 1;
y(1) = 1;
for i = 1 : Nit

  %method 2
  k1 = -20 * x(i);
  k2 = -20 * (x(i) + 0.5*h*k1);
  k3 = -20 * (x(i) - (1*h*k1) + (2*h*k2));
  x(i+1) = x(i) + h*(k1*(1/6) + k2*(2/3) + k3*(1/6));
  %method 1
  k4 = -20 * y(i);
  k5 = -20 * (y(i) + 0.5*h*k4);
  k6 = -20 * (y(i) + 0.5*h*k5);  
  k7 = -20 * (y(i) + h*k6);
  y(i+1) = y(i) + h*( k4*(1/6) + k5*(1/3) + k6*(1/3) + k7*(1/6));
  t(i+1) = i*h; %time
end

fprintf('0.0625|\t%d\t     %d | %d\t%d   ', x(end), y(end), x(end)-exact, y(end)-exact);

t1=[0 0.0625 0.125 0.1875 0.25 0.3125 0.375 0.4375 0.5 0.5625 0.625 0.6875 0.75 0.8125 0.875 0.9375 1.0 1.0625 1.125 1.1875 1.25 1.3125 1.375 1.4375 1.5 1.5625 1.625 1.6875 1.75 1.8125 1.875 1.9375 2.0];
u=@(t1) 1*exp(-20*(t1));

figure;
plot(t,x,'y-',"linewidth",8)
hold on
plot(t1,u(t1),"linewidth",5)
plot(t,y,'b--',"linewidth",6)
hold off