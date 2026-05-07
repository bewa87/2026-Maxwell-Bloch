% Fixed-Point Iteration For Non-Linear, Three-Dimensional
% ODE-System Of Maxwell-Bloch Equations

% Step 1: Setup Of Initial Conditions

x_null = -0.5;
y_null = -0.5;
z_null = 0.5;
R      = sqrt(x_null^2+y_null^2+(z_null+1)^2);
t_max  = 1/(R+1);% (2/(sqrt(17)*(R+1))) (Frobenius norm)

% Step 2: Initialization Of Fixed-Point Solution Procedure

h      = 0.0125;
T      = 1000;
t      = 0:h:T;

x      = zeros(length(t),1);
y      = zeros(length(t),1);
z      = zeros(length(t),1);
x(1)   = x_null;
y(1)   = y_null;
z(1)   = z_null;

N_max  = 20; % Maximal Number of Iterations

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

% Step 3: Comparison With Explicit Eulerian Method

x_ee      = zeros(length(t),1);
y_ee      = zeros(length(t),1);
z_ee      = zeros(length(t),1);
x_ee(1)   = x_null;
y_ee(1)   = y_null;
z_ee(1)   = z_null;

for j = 1:1:(length(t)-1)
  x_ee(j+1) = x_ee(j) + h*y_ee(j);
  y_ee(j+1) = y_ee(j) + h*x_ee(j)*z_ee(j);
  z_ee(j+1) = z_ee(j) - h*x_ee(j)*y_ee(j);
endfor

##x_ie      = zeros(length(t),1);
##y_ie      = zeros(length(t),1);
##z_ie      = zeros(length(t),1);
##x_ie(1)   = x_null;
##y_ie(1)   = y_null;
##z_ie(1)   = z_null;
##
##for j = 1:1:(length(t)-1)
##  x_ie_help    = zeros((N_max+1),1);
##  y_ie_help    = zeros((N_max+1),1);
##  z_ie_help    = zeros((N_max+1),1);
##  x_ie_help(1) = x(j);
##  y_ie_help(1) = y(j);
##  z_ie_help(1) = z(j);
##  for k = 1:1:N_max
##    Delta          = - inv([1 -h 0; -h*z_ie_help(k) 1 0; h*y_ie_help(k) h*x_ie_help(k) 1])*[x_ie_help(k) - x_ie(j) - h*y_ie_help(k); y_ie_help(k) - y_ie(j) - h*x_ie_help(k)*z_ie_help(k); z_ie_help(k) - z_ie(j) + h*x_ie_help(k)*y_ie_help(k)];
##    x_ie_help(k+1) = x_ie(j) + Delta(1);
##    y_ie_help(k+1) = y_ie(j) + Delta(2);
##    z_ie_help(k+1) = z_ie(j) + Delta(3);
##  endfor
##  x_ie(j+1) = x_ie_help(end);
##  y_ie(j+1) = y_ie_help(end);
##  z_ie(j+1) = z_ie_help(end);
##endfor

% Step 4: Comparison With Second-Order Runge-Kutta Method

x_rk2      = zeros(length(t),1);
y_rk2      = zeros(length(t),1);
z_rk2      = zeros(length(t),1);
x_rk2(1)   = x_null;
y_rk2(1)   = y_null;
z_rk2(1)   = z_null;

for j = 1:1:(length(t)-1)
  x_rk2_help = x_rk2(j) + h*y_rk2(j);
  y_rk2_help = y_rk2(j) + h*x_rk2(j)*z_rk2(j);
  z_rk2_help = z_rk2(j) - h*x_rk2(j)*y_rk2(j);
  x_rk2(j+1) = 0.5*(x_rk2(j)+x_rk2_help) + 0.5*h*y_rk2_help;
  y_rk2(j+1) = 0.5*(y_rk2(j)+y_rk2_help) + 0.5*h*x_rk2_help*z_rk2_help;
  z_rk2(j+1) = 0.5*(z_rk2(j)+z_rk2_help) - 0.5*h*x_rk2_help*y_rk2_help;
endfor

% Step 5: Plotting Of Solutions

figure(1)
plot(t, x, 'color', 'red', 'linewidth', 0.75)
hold on
plot(t, y, 'color', 'blue', 'linewidth', 0.75)
hold on
plot(t, z, 'color', 'black', 'linewidth', 0.75)
hold off
xlabel('Time t', 'fontsize', 14)
ylabel('Solution Components', 'fontsize', 14)
title('Example 1: TDS of First Iterative Solver', 'fontsize', 16)
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
title('Example 1: TDCL of First Iterative Solver', 'fontsize', 16)
legend({'Casimir', 'Hamiltonian', 'Sphere'}, 'fontsize', 14, 'location', 'eastoutside')

figure(3)
scatter3(x,y,z)
hold off
xlabel('x(t)', 'fontsize', 14)
ylabel('y(t)', 'fontsize', 14)
zlabel('z(t)', 'fontsize', 14)
title('Example 1: Spherical Surface for First Iterative Solver', 'fontsize', 16)

A_fit = ((x.^2)/2+z) - ((x_null^2)/2+z_null);
B_fit = (((y.^2)/2)+((z.^2)/2)) - ((y_null^2/2)+(z_null^2/2));
C_fit = ((x.^2)+(y.^2)+((z+1).^2)) - (x_null^2+y_null^2+(z_null+1)^2);

figure(4)
plot(t, x_ee, 'color', 'red', 'linewidth', 0.75)
hold on
plot(t, y_ee, 'color', 'blue', 'linewidth', 0.75)
hold on
plot(t, z_ee, 'color', 'black', 'linewidth', 0.75)
hold off
xlabel('Time t', 'fontsize', 14)
ylabel('Solution Components', 'fontsize', 14)
title('Example 1: TDS of Explicit Eulerian', 'fontsize', 16)
legend({'x(t)', 'y(t)', 'z(t)'}, 'fontsize', 14, 'location', 'eastoutside')

figure(5)
plot(t,((x_ee.^2)/2+z_ee), 'color', 'red', 'linewidth', 0.75)
hold on
plot(t,(((y_ee.^2)/2)+((z_ee.^2)/2)), 'color', 'blue', 'linewidth', 0.75)
hold on
plot(t,((x_ee.^2)+(y_ee.^2)+((z_ee+1).^2)), 'color', 'black', 'linewidth', 0.75)
hold off
xlabel('Time t', 'fontsize', 14)
ylabel('Conservation Laws', 'fontsize', 14)
title('Example 1: TDCL of Explicit Eulerian', 'fontsize', 16)
legend({'Casimir', 'Hamiltonian', 'Sphere'}, 'fontsize', 14, 'location', 'eastoutside')

figure(6)
plot(t, x_rk2, 'color', 'red', 'linewidth', 0.75)
hold on
plot(t, y_rk2, 'color', 'blue', 'linewidth', 0.75)
hold on
plot(t, z_rk2, 'color', 'black', 'linewidth', 0.75)
hold off
xlabel('Time t', 'fontsize', 14)
ylabel('Solution Components', 'fontsize', 14)
title('Example 1: TDS of Heun Method', 'fontsize', 16)
legend({'x(t)', 'y(t)', 'z(t)'}, 'fontsize', 14, 'location', 'eastoutside')

figure(7)
plot(t,((x_rk2.^2)/2+z_rk2), 'color', 'red', 'linewidth', 0.75)
hold on
plot(t,(((y_rk2.^2)/2)+((z_rk2.^2)/2)), 'color', 'blue', 'linewidth', 0.75)
hold on
plot(t,((x_rk2.^2)+(y_rk2.^2)+((z_rk2+1).^2)), 'color', 'black', 'linewidth', 0.75)
hold off
xlabel('Time t', 'fontsize', 14)
ylabel('Conservation Laws', 'fontsize', 14)
title('Example 1: TDCL of Heun Method', 'fontsize', 16)
legend({'Casimir', 'Hamiltonian', 'Sphere'}, 'fontsize', 14, 'location', 'eastoutside')

A_fit_rk2 = ((x_rk2.^2)/2+z_rk2) - ((x_null^2)/2+z_null);
B_fit_rk2 = (((y_rk2.^2)/2)+((z_rk2.^2)/2)) - ((y_null^2/2)+(z_null^2/2));
C_fit_rk2 = ((x_rk2.^2)+(y_rk2.^2)+((z_rk2+1).^2)) - (x_null^2+y_null^2+(z_null+1)^2);
