clc
clear all
close all


f_0 = 77e9;
c = 3e8;
lambda = c/f_0;
n_fft=59*20;

N = 59;
multipdx=1;
dx = multipdx*lambda/4;
f_sx = 1/dx;
L = N*dx;

tgt_range = 3;
tgt_angle = [-10];
antenna_locs = ((1-N)/2 : (N-1)/2)';
s_dem = zeros(N,length(tgt_angle));
s_dem_fft = zeros(n_fft,length(tgt_angle));
A_hat = exp(-1i*4*pi*tgt_range/lambda);                                  
f_x_axis = linspace(-f_sx/2,f_sx/2,n_fft);
theta_axis = asind((f_x_axis*lambda/2));


noise_inc = true;
P_N = 1;
 w = noise_inc*randn(N,2);
%w = noise_inc*wgn(N,2,P_N,'linear');
w(:,3) = w(:,1) + 1i*w(:,2);

for i = 1:length(tgt_angle)
    theta = tgt_angle(i);
    s_dem(:,i) = ((A_hat+w(:,3)).*exp(1i*2*pi*sind(theta)*2.*antenna_locs*dx/lambda));
    s_dem_fft(:,i) = fftshift(fft(s_dem(:,i),n_fft))*dx;
    plot(theta_axis,abs(s_dem_fft(:,i)),LineWidth=1)
    title("FFT of Received Signal")
    xlabel("Angle of the target[theta]")
    ylabel("Amplitude")
    hold on 

end

% Compare the estimated DoA with the true angular position of the target. Is it the same?
[a,b] = max(s_dem_fft);
estimated_theta = theta_axis(b);

if length(tgt_angle)~=1
    s_dem_multi_tgt_fft = fftshift(fft(sum(s_dem,2)+w(:,3),n_fft))*dx;
    figure
    plot(theta_axis,abs(s_dem_multi_tgt_fft))
    title("FFT of Received Signal When There Are More Than One Target")
    xlabel("Angle of the targets[theta]")
    ylabel("Amplitude")
end
