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
