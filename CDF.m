clear
%solving problem
%-0.05*u'' + r(x) * x' = 2
% u(0)=0, u(1)=1

%set values of the problem
alpha = 0;
beta  = 2;
%functions r and f
r = @(x) 2 + 0 * x;
f = @(x) 2 + 0 * x;
%exact solution
u = @(x) x - ( (exp(-40) - exp(-40*(1-x)) ) / (1 - exp(-40)) );
fprintf('M\t|Approximation|NGE\n');
Marray = [10 19 21 42];

imax = 4;
for ii =1: imax
%discretization parameters
  M = Marray(ii); %remember that the system will be (M-1)x(M-1)
  h = 1/M;
  x_m = 0:h:1;
  x_m = x_m(2:end-1);
  %setting up matrix for the system
  D = (2*.05 * h.^-2 *diag(ones(M-1,1)) );
  T1 = (diag(ones(M-2,1),-1))* (-1*(0.05 *h.^(-2) + h.^(-1)) );
  T2 = (diag(ones(M-2,1),1)) * (-1*(0.05 *h.^(-2) - h.^(-1)) );
  
  %change A to have right coefficents
  A = T2+ T1 + D;
  %setting up rhs
  b = f(x_m)';
  b(1)   = b(1)   + alpha/h/h;
  b(end) = b(end) + beta/h/h;
  %solving this system
  usol = A \ b;


  %calculating exact solution
  uex = u(x_m)';
  %calculating error
  GE = uex - usol;
  NGE(ii) = max(abs(GE));

  %printing table results
  fprintf('%d\t|%d\t  |%d  \n', M(1), usol(end),NGE(ii));
  
  %plotting ex and comp sols
  figure; 
  plot(x_m,uex,'-k');
  hold on;
  plot(x_m,usol,'-xr');
  hold off;
  
end
