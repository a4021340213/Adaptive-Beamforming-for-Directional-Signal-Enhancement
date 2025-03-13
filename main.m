clc
clear all;
%% parameter set up
load('ASP_Final_Data.mat');
L = numel(theta_s_noisy);  % Time length
t = [1:L];                 % Time vector
N = numel(matX(:,1));      % Number of isotropic antennas
d = 0:(N-1);

%% Q2 denoising
theta_i_hat=AIC_SVD(theta_i_noisy,t);
theta_s_hat=AIC_SVD(theta_s_noisy,t);
figure;
plot(t,theta_s_noisy,t,theta_i_noisy);
legend('Signal $\tilde{\theta}_{s}(t)$','Interference $\tilde{\theta}_{i}(t)$','interpreter','latex','fontsize',14)
xlabel('Time')
ylabel('Degree')
ylim([-17,10]);
title('DOA Without Denoised')

figure;
plot(t,theta_s_hat,t,theta_i_hat);
legend('Signal $\hat{\theta}_{s}(t)$','Interference $\hat{\theta}_{i}(t)$','interpreter','latex','fontsize',14)
xlabel('Time')
ylabel('Degree')
ylim([-17,10]);
title('DOA With Denoised')

%% Q3 beamformer 

d=0:1:N-1;
theta=linspace(-pi/2,pi/2,2000);
theta_deg = rad2deg(theta);
theta_input = exp(j*pi*sin(theta).*d.');

%% uniform beamformer
[y_uniform w_uniform B_uniform] = uniformly_weighted(matX, theta_input, N, L);


figure;
image(t, theta_deg, abs(B_uniform)/max(max(abs(B_uniform)))*1000);
set(gca, "Ydir", "normal");
xlabel('Time')
ylabel('Degree')
title('Beampattern of uniformly weighted beamformer')
colorbar;

q_uniform=10.^(B_uniform/20);
q_uniform=q_uniform/max(max(q_uniform));
q_uniform=20*log(q_uniform);

figure;
surf(t, theta_deg, q_uniform);
shading interp;
xlabel('Time (snapshot index)');
ylabel('Angle (degrees)');
zlabel('Beam Pattern Magnitude(dB)');
zlim([-400 0]);
ylim([-90,90]);
title('3D Beampattern Over Time and Angle (uniform)');
colorbar;

figure;
subplot(2,1,1)
plot(t,real(y_uniform));
axis tight
xlabel('Time','fontsize',16)
ylabel('$Re\left\{\hat{s}(t)\right\}$','interpreter','latex','fontsize',16)
title('Real part of $\hat{s}(t)$ by uniformly weight beamformer','interpreter','latex','fontsize',18)
subplot(2,1,2)
plot(t,imag(y_uniform));
axis tight
xlabel('Time','fontsize',16)
ylabel('$Im\left\{\hat{s}(t)\right\}$','interpreter','latex','fontsize',16)
title('Imaginary part of $\hat{s}(t)$ by uniformly weight beamformer','interpreter','latex','fontsize',18)


%% Estimate source signal by array steering beamformer

[y_array_steering w_array_steering B_array_steering] = array_steering(matX, N, theta_s_hat, theta_input, d, L);


figure;
image(t, theta_deg, abs(B_array_steering)/max(max(abs(B_array_steering)))*1000);
set(gca, "Ydir", "normal");
xlabel('Time')
ylabel('Degree')
title('Beampattern of array steering beamformer')
colorbar;

q_steering=10.^(B_array_steering/20);
q_steering=q_steering/max(max(q_steering));
q_steering=20*log(q_steering);


figure;
surf(t, theta_deg, q_steering);
shading interp;
xlabel('Time (snapshot index)');
ylabel('Angle (degrees)');
zlabel('Beam Pattern Magnitude(dB)');
zlim([-400 0]);
ylim([-90,90]);
title('3D Beampattern Over Time and Angle(steering)');
colorbar;

figure;
subplot(2,1,1)
plot(t,real(y_array_steering));
axis tight
xlabel('Time','fontsize',16)
ylabel('$Re\left\{\hat{s}(t)\right\}$','interpreter','latex','fontsize',16)
title('Real part of $\hat{s}(t)$ by array steering beamformer','interpreter','latex','fontsize',18)
subplot(2,1,2)
plot(t,imag(y_array_steering));
axis tight
xlabel('Time','fontsize',16)
ylabel('$Im\left\{\hat{s}(t)\right\}$','interpreter','latex','fontsize',16)
title('Imaginary part of $\hat{s}(t)$ by array steering beamformer','interpreter','latex','fontsize',18)


%%  MVDR beamformer


delta = 0.01;
mu = 0.53;
k = zeros(N,1);
P = eye(N)/delta;

[y_MVDR w_MVDR B_MVDR ]= MVDR(matX, P, theta_s_hat, d, theta_input, 0.01, 0.999, N, L);


figure;
image(t, theta_deg, abs(B_MVDR)/max(max(abs(B_MVDR)))*1000);
set(gca, "Ydir", "normal");
xlabel('Time')
ylabel('magnitude')
title('Beampattern of MVDR beamformer')
colorbar;

q_MVDR=10.^(B_MVDR/20);
q_MVDR=q_MVDR/max(max(q_MVDR));
q_MVDR=20*log(q_MVDR);

figure;
surf(t, theta_deg, q_MVDR);
shading interp;
xlabel('Time (snapshot index)');
ylabel('Angle (degrees)');
zlabel('Beam Pattern Magnitude(dB)');
zlim([-400 0]);
ylim([-90,90]);
title('3D Beampattern Over Time and Angle(MVDR)');
colorbar;

figure;
plot(theta_deg,q_MVDR(:,1000));
axis tight
xlabel('Time','fontsize',16)
ylabel('$Re\left\{\hat{s}(t)\right\}$','interpreter','latex','fontsize',16)
title('Real part of $\hat{s}(t)$ by MVDR beamformer','interpreter','latex','fontsize',18)

figure;
subplot(2,1,1)
plot(t,real(y_MVDR));
axis tight
xlabel('Time','fontsize',16)
ylabel('$Re\left\{\hat{s}(t)\right\}$','interpreter','latex','fontsize',16)
title('Real part of $\hat{s}(t)$ by MVDR beamformer','interpreter','latex','fontsize',18)
subplot(2,1,2)
plot(t,imag(y_MVDR));
axis tight
xlabel('Time','fontsize',16)
ylabel('$Im\left\{\hat{s}(t)\right\}$','interpreter','latex','fontsize',16)
title('Imaginary part of $\hat{s}(t)$ by MVDR beamformer','interpreter','latex','fontsize',18)


%% LCMV beamformer

[y_LCMV w_LCMV B_LCMV ] = LCMV(matX, P, theta_s_hat, theta_i_hat, d,theta_input, 0.01, 0.999, N, L);


figure;
image(t, theta_deg, abs(B_LCMV)/max(max(abs(B_LCMV)))*1000);
set(gca, "Ydir", "normal");
xlabel('Time')
ylabel('Degree')
title('Beampattern of LCMV beamformer')
colorbar;


q_LCMV=10.^(B_LCMV/20);
q_LCMV=q_LCMV/max(max(q_LCMV));
q_LCMV=20*log(q_LCMV);

figure;
surf(t, theta_deg, q_LCMV);
shading interp;
xlabel('Time (snapshot index)');
ylabel('Angle (degrees)');
zlabel('Beam Pattern Magnitude(dB)');
zlim([-400 0]);
ylim([-90,90]);
title('3D Beampattern Over Time and Angle(LCMV)');
colorbar;


figure;
plot(theta_deg,q_LCMV(:,1000));
axis tight
xlabel('Time','fontsize',16)
ylabel('$Re\left\{\hat{s}(t)\right\}$','interpreter','latex','fontsize',16)
title('Real part of $\hat{s}(t)$ by LCMV beamformer','interpreter','latex','fontsize',18)

figure;
subplot(2,1,1)
plot(t,real(y_LCMV));
axis tight
xlabel('Time','fontsize',16)
ylabel('$Re\left\{\hat{s}(t)\right\}$','interpreter','latex','fontsize',16)
title('Real part of $\hat{s}(t)$ by LCMV beamformer','interpreter','latex','fontsize',18)
subplot(2,1,2)
plot(t,imag(y_LCMV));
axis tight
xlabel('Time','fontsize',16)
ylabel('$Im\left\{\hat{s}(t)\right\}$','interpreter','latex','fontsize',16)
title('Imaginary part of $\hat{s}(t)$ by LCMV beamformer','interpreter','latex','fontsize',18)

[s_t_hat w_proposed B_proposed] = RAIS_LCMV(matX,theta_i_hat,theta_s_hat,N,theta_input,L,d);

figure;
image(t, theta_deg, abs(B_proposed)/max(max(abs(B_proposed)))*1000);
set(gca, "Ydir", "normal");
xlabel('Time')
ylabel('Degree')
title('Beampattern of LCMV beamformer')
colorbar;


q_proposed=10.^(B_proposed/20);
q_proposed=q_proposed/max(max(q_proposed));
q_proposed=20*log(q_proposed);

figure;
surf(t, theta_deg, q_proposed);
shading interp;
xlabel('Time (snapshot index)');
ylabel('Angle (degrees)');
zlabel('Beam Pattern Magnitude(dB)');
zlim([-400 0]);
ylim([-90,90]);
title('3D Beampattern Over Time and Angle(proposed)');
colorbar;


figure;
plot(theta_deg,q_proposed(:,1000));
axis tight
xlabel('Time','fontsize',16)
ylabel('$Re\left\{\hat{s}(t)\right\}$','interpreter','latex','fontsize',16)
title('Real part of $\hat{s}(t)$ by proposed beamformer','interpreter','latex','fontsize',18)


figure;
subplot(2,1,1)
plot(t,real(s_t_hat));
axis tight
xlabel('Time','fontsize',16)
ylabel('$Re\left\{\hat{s}(t)\right\}$','interpreter','latex','fontsize',16)
title('Real part of $\hat{s}(t)$ by proposed beamformer','interpreter','latex','fontsize',18)
subplot(2,1,2)
plot(t,imag(s_t_hat));
axis tight
xlabel('Time','fontsize',16)
ylabel('$Im\left\{\hat{s}(t)\right\}$','interpreter','latex','fontsize',16)
title('Imaginary part of $\hat{s}(t)$ by proposed beamformer','interpreter','latex','fontsize',18)
save('ASP_Final_result.mat', 'theta_s_hat', 'theta_i_hat','s_t_hat');
save('weight.mat', 'w_uniform', 'w_array_steering','w_MVDR','w_LCMV','w_proposed');