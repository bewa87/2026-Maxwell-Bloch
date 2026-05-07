% Fixed-Point Iteration For Non-Linear, Three-Dimensional
% ODE-System Of Maxwell-Bloch Equations

% Step 1: Setup Of Initial Conditions

x_null = -0.5;
y_null = -0.5;
z_null = 0.5;
R      = sqrt(x_null^2+y_null^2+(z_null+1)^2);
t_max  = 1/(R+1);% (2/(sqrt(17)*(R+1))) (Frobenius norm)

% Step 2: Initialization Of Fixed-Point Solution Procedure

h      = 1;
T      = 100;
t      = 0:h:T;

x      = zeros(length(t),1);
y      = zeros(length(t),1);
z      = zeros(length(t),1);
x(1)   = x_null;
y(1)   = y_null;
z(1)   = z_null;

N_max  = 200; % Maximal Number of Iterations

for j = 1:1:(length(t)-1)
  x_help    = zeros((N_max+1),1);
  y_help    = zeros((N_max+1),1);
  z_help    = zeros((N_max+1),1);
  x_help(1) = x(j);
  y_help(1) = y(j);
  z_help(1) = z(j);
  for k = 1:1:N_max
    x_help(k+1) = x(j) + (h/2)*(y(j)+y_help(k));
    y_help(k+1) = y(j) + (h/4)*(x(j)+x_help(k))*(z(j)+z_help(k));
    z_help(k+1) = z(j) - (h/4)*(x(j)+x_help(k))*(y(j)+y_help(k));
  endfor
  x(j+1) = x_help(end);
  y(j+1) = y_help(end);
  z(j+1) = z_help(end);
endfor

% Step 3: Plotting Of Solutions

figure(1)
plot(t, x, 'color', 'red', 'linewidth', 0.75)
hold on
plot(t, y, 'color', 'blue', 'linewidth', 0.75)
hold on
plot(t, z, 'color', 'black', 'linewidth', 0.75)
hold off
xlabel('Time t', 'fontsize', 14)
ylabel('Solution Components', 'fontsize', 14)
title('Example 3: TDS of First Iterative Solver (h = 1)', 'fontsize', 16)
legend({'x(t)', 'y(t)', 'z(t)'}, 'fontsize', 14, 'location', 'eastoutside')

figure(2)
plot(t,((x.^2)/2+z), 'color', 'red', 'linewidth', 0.75)
hold on
plot(t,(((y.^2)/2)+((z.^2)/2)), 'color', 'blue', 'linewidth', 0.75)
hold on
plot(t,((x.^2)+(y.^2)+((z+1).^2)), 'color', 'black', 'linewidth', 0.75)
hold off
xlabel('Time t', 'fontsize', 14)
ylabel('Conservation Laws', 'fontsize', 14)
title('Example 3: TDCL of First Iterative Solver', 'fontsize', 16)
legend({'Casimir', 'Hamiltonian', 'Sphere'}, 'fontsize', 14, 'location', 'eastoutside')

A_fit = ((x.^2)/2+z) - ((x_null^2)/2+z_null);
B_fit = (((y.^2)/2)+((z.^2)/2)) - ((y_null^2/2)+(z_null^2/2));
C_fit = ((x.^2)+(y.^2)+((z+1).^2)) - (x_null^2+y_null^2+(z_null+1)^2);
