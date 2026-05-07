% Newton-Iteration For Non-Linear, Three-Dimensional
% ODE-System Of Maxwell-Bloch Equations

% Step 1: Setup Of Initial Conditions

x_null = -0.5;
y_null = -0.5;
z_null = 0.5;

% Step 2: Initialization Of Newton Solution Procedure

h      = 0.1;
T      = 100;
t      = 0:h:T;

x      = zeros(length(t),1);
y      = zeros(length(t),1);
z      = zeros(length(t),1);
x(1)   = x_null;
y(1)   = y_null;
z(1)   = z_null;

F      = @(x_cur,y_cur,z_cur,x_prev,y_prev,z_prev,h_tilde) [(x_cur-x_prev)-(h_tilde/2)*(y_cur+y_prev); (y_cur-y_prev)-(h_tilde/4)*(x_prev+x_cur)*(z_prev+z_cur); (z_cur-z_prev)+(h_tilde/4)*(x_prev+x_cur)*(y_prev+y_cur)];
J      = @(x_cur,y_cur,z_cur,x_prev,y_prev,z_prev,h_tilde) [1, -(h_tilde/2), 0; -(h_tilde/4)*(z_prev+z_cur), 1, -(h_tilde/4)*(x_prev+x_cur); (h_tilde/4)*(y_prev+y_cur), (h_tilde/4)*(x_prev+x_cur), 1];

N_max  = 20; % Maximal Number of Newton-Iterations

for j = 1:1:(length(t)-1)
  x_help    = zeros((N_max+1),1);
  y_help    = zeros((N_max+1),1);
  z_help    = zeros((N_max+1),1);
  x_help(1) = x(j);
  y_help(1) = y(j);
  z_help(1) = z(j);
  Delta     = zeros(3,1);
  for k = 1:1:N_max
    Delta       = J(x_help(k),y_help(k),z_help(k),x(j),y(j),z(j),h)\(-F(x_help(k),y_help(k),z_help(k),x(j),y(j),z(j),h));
    x_help(k+1) = x_help(k) + Delta(1,1);
    y_help(k+1) = y_help(k) + Delta(2,1);
    z_help(k+1) = z_help(k) + Delta(3,1);
  endfor
  x(j+1) = x_help(end);
  y(j+1) = y_help(end);
  z(j+1) = z_help(end);
endfor

% Step 3: Plot Of Solutions

figure(4)
plot(t, x, 'color', 'red', 'linewidth', 0.75)
hold on
plot(t, y, 'color', 'blue', 'linewidth', 0.75)
hold on
plot(t, z, 'color', 'black', 'linewidth', 0.75)
hold off
xlabel('Time t', 'fontsize', 14)
ylabel('Solution Components', 'fontsize', 14)
title('Example 2: TDS of Second Iterative Solver', 'fontsize', 16)
legend({'x(t)', 'y(t)', 'z(t)'}, 'fontsize', 14, 'location', 'eastoutside')

figure(5)
plot(t, ((x.^2)/2+z), 'color', 'red', 'linewidth', 0.75)
hold on
plot(t, ((y.^2)/2)+((z.^2)/2), 'color', 'blue', 'linewidth', 0.75)
hold on
plot(t, ((x.^2)+(y.^2)+((z+1).^2)), 'color', 'black', 'linewidth', 0.75)
hold off
xlabel('Time t', 'fontsize', 14)
ylabel('Conservation Laws', 'fontsize', 14)
title('Example 2: TDCL of Second Iterative Solver', 'fontsize', 16)
legend({'Casimir', 'Hamiltonian', 'Sphere'}, 'fontsize', 14, 'location', 'eastoutside')

figure(6)
scatter3(x,y,z)

A_nit = ((x.^2)/2+z) - ((x_null^2)/2+z_null);
B_nit = (((y.^2)/2)+((z.^2)/2)) - ((y_null^2/2)+(z_null^2/2));
C_nit = ((x.^2)+(y.^2)+((z+1).^2)) - (x_null^2+y_null^2+(z_null+1)^2);
