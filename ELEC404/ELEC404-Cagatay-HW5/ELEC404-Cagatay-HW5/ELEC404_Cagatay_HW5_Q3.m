clear
close all

[soundfile, Fs]=audioread('sa1.wav');

starting_sample=0.65*Fs; 
%starting_sample = 0.94*Fs;
L = 480; 
N = 512;
p = 12; %LPC order
hamming_window = hamming(L);

% windowed speech
window_frame = [(soundfile(starting_sample:starting_sample+L-1).*hamming(L))' zeros(1,N-L)];
windowed_speech = fft(window_frame, N);
freq = 0 : Fs/512 : Fs;
plot(freq(1:(N/2)), 20*log10(abs(windowed_speech(1:(N/2)))))
ylabel('Magnitude (dB)');
xlabel('Frequency (Hz)');
hold on

% autocorrelation 
autocorrelation_frame = soundfile(starting_sample : starting_sample+L-1);
autocorrelation_frame = autocorrelation_frame.*hamming_window;

R = zeros(p+1);
for i = 0:p
    R(i+1) = sum(autocorrelation_frame(1:L-i).*autocorrelation_frame(i+1:L));
end

eps = zeros(1,p);
eps(1)=R(1);

alpha = zeros(1,p);
alpha(1)=R(2)/eps(1);

autocor_matrice = zeros(p,p);
autocor_matrice(1,1) = alpha(1);
eps(2) = (1-alpha(1).^2)*eps(1);

for j = 2:p
    alpha(j) = (R(j+1)-sum(autocor_matrice(1:j-1,j-1)'.*R(j:-1:2)))/eps(j);
    autocor_matrice(j,j) = alpha(j);
    for j2 = 1:j-1
        autocor_matrice(j2,j) = autocor_matrice(j2,j-1)-alpha(j)*autocor_matrice(j-j2,j-1);
    end
    eps(j+1) = (1-alpha(j).^2)*eps(j);
end
autocor = autocor_matrice(1:p,p);
[autocor_sig, freq] = freqz(sqrt(eps(p+1)),[1 -autocor'], N, Fs);
plot(freq, 20*log10(abs(autocor_sig)), 'g');

% covariance
coveriance_frame = soundfile(starting_sample-p:starting_sample+L-1);

cov_matrice = zeros(p, p);
for i = 1:p
    for i2 = 1:p
        cov_matrice(i, i2) = sum(coveriance_frame(p+1-i:p+L-i).*coveriance_frame(p+1-i2:p+L-i2));
    end
end

v = zeros(1, p);
z = zeros(1, p);
for i=1:p
    v(i) = sum(coveriance_frame(p+1-i:p+L-i).*coveriance_frame(p+1:p+L));
    z(i) = sum(coveriance_frame(p+1:p+L).*coveriance_frame(p+1-i:p+L-i));
end

p0 = sum(coveriance_frame(p+1:p+L).^2);
coveriance = inv(cov_matrice)*v';
[coveriance_sig, freq] = freqz(sqrt(p0-sum(coveriance'.*z)),[1 -coveriance'], N, Fs);
plot(freq, 20*log10(abs(coveriance_sig)), 'r')
legend('windowed speech','autocorrelation','covariance');