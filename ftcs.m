clear
%ftcs
%bc and ic
ub = @(x) sin(pi*x) - sin(3*pi*x);
uleft  = @(t) 0;
uright = @(t) 0;

%exact sol
u = @(x,t) exp(-pi*pi*t).*sin(pi*x) - exp(-9*pi*pi*t)*sin(3*pi*x);

%parameters (r=k/h^2)
M = 10;
h = 1/M;
k = 1/600;
r = k/(h^2);
x = 0:h:1;
xi = x(2:end-1);
t = 0:k:1;

%setup matrices
D = (2*diag(ones(M-1,1)) - diag(ones(M-2,1),1) - diag(ones(M-2,1),-1));
I = eye(M-1);
A = I - r * D;
%U0
uin = ub(xi); 
u0 = uin';
usol = zeros(length(xi),length(t));
uex  = zeros(length(xi),length(t));
uex(:,1) = u0;
usol(:,1) = u0;
figure; 
hold on;
for ii = 2 : length(t)    
    u1 = A*u0;
    u0 = u1;
    usol(:,ii) = u1;    
    uex(:,ii) = u(xi,t(ii));
    if mod(ii,M) == 1 %&& t(ii)<0.05
        plot3(xi,t(ii)*ones(M-1,1),u1,'k');
    end
end
fprintf('Error=%d\n',max(abs(usol(:,end)-uex(:,end))))
hold off;

ub = @(x) sin(pi*x) - sin(3*pi*x);
uleft  = @(t) 0;
uright = @(t) 0;

%exact sol
u = @(x,t) exp(-pi*pi*t).*sin(pi*x) - exp(-9*pi*pi*t)*sin(3*pi*x);

%parameters (r=k/h^2)
M = 10;
h = 1/M;
k = 1/300;
r = k/(h^2);
x = 0:h:1;
xi = x(2:end-1);
t = 0:k:1;

%setup matrices
D = (2*diag(ones(M-1,1)) - diag(ones(M-2,1),1) - diag(ones(M-2,1),-1));
I = eye(M-1);
A = I - r * D;
%U0
uin = ub(xi); 
u0 = uin';
usol = zeros(length(xi),length(t));
uex  = zeros(length(xi),length(t));
uex(:,1) = u0;
usol(:,1) = u0;
figure; 
hold on;
for ii = 2 : length(t)    
    u1 = A*u0;
    u0 = u1;
    usol(:,ii) = u1;    
    uex(:,ii) = u(xi,t(ii));
    if mod(ii,M) == 1 %&& t(ii)<0.05
        plot3(xi,t(ii)*ones(M-1,1),u1,'k');
    end
end
fprintf('Error=%d\n',max(abs(usol(:,end)-uex(:,end))))
hold off;

ub = @(x) sin(pi*x) - sin(3*pi*x);
uleft  = @(t) 0;
uright = @(t) 0;

%exact sol
u = @(x,t) exp(-pi*pi*t).*sin(pi*x) - exp(-9*pi*pi*t)*sin(3*pi*x);

%parameters (r=k/h^2)
M = 10;
h = 1/M;
k = 1/100;
r = k/(h^2);
x = 0:h:1;
xi = x(2:end-1);
t = 0:k:1;

%setup matrices
D = (2*diag(ones(M-1,1)) - diag(ones(M-2,1),1) - diag(ones(M-2,1),-1));
I = eye(M-1);
A = I - r * D;
%U0
uin = ub(xi); 
u0 = uin';
usol = zeros(length(xi),length(t));
uex  = zeros(length(xi),length(t));
uex(:,1) = u0;
usol(:,1) = u0;
figure; 
hold on;
for ii = 2 : length(t)    
    u1 = A*u0;
    u0 = u1;
    usol(:,ii) = u1;    
    uex(:,ii) = u(xi,t(ii));
    if mod(ii,M) == 1 %&& t(ii)<0.05
        plot3(xi,t(ii)*ones(M-1,1),u1,'k');
    end
end
fprintf('Error=%d\n',max(abs(usol(:,end)-uex(:,end))))
hold off;

ub = @(x) sin(pi*x) - sin(3*pi*x);
uleft  = @(t) 0;
uright = @(t) 0;

%exact sol
u = @(x,t) exp(-pi*pi*t).*sin(pi*x) - exp(-9*pi*pi*t)*sin(3*pi*x);

%parameters (r=k/h^2)
M = 10;
h = 1/M;
k = 1/10;
r = k/(h^2);
x = 0:h:1;
xi = x(2:end-1);
t = 0:k:1;

%setup matrices
D = (2*diag(ones(M-1,1)) - diag(ones(M-2,1),1) - diag(ones(M-2,1),-1));
I = eye(M-1);
A = I - r * D;
%U0
uin = ub(xi); 
u0 = uin';
usol = zeros(length(xi),length(t));
uex  = zeros(length(xi),length(t));
uex(:,1) = u0;
usol(:,1) = u0;
figure; 
hold on;
for ii = 2 : length(t)    
    u1 = A*u0;
    u0 = u1;
    usol(:,ii) = u1;    
    uex(:,ii) = u(xi,t(ii));
    if mod(ii,M) == 1 %&& t(ii)<0.05
        plot3(xi,t(ii)*ones(M-1,1),u1,'k');
    end
end
fprintf('Error=%d\n',max(abs(usol(:,end)-uex(:,end))))
hold off;