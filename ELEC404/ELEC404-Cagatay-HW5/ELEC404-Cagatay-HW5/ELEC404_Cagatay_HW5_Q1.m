clear
close all

[soundfile, Fs] = audioread('sa1.wav');

soundfile_size = size(soundfile);
soundfile_size = soundfile_size(1);


L=480;
R=240;
hamming_window = hamming(L);


cepstrum = zeros(soundfile_size, L);
log_mag_spectrum = zeros(soundfile_size, L);

for i=1:R:soundfile_size    
    if (L+i-1>soundfile_size)
       break; 
    else
    log_mag_spectrum(i, 1:L) = log(abs(fft(hamming_window.*soundfile(i:L+i-1))));
    cepstrum(i, 1:L) = real(ifft(log(abs(fft(hamming_window.*soundfile(i:L+i-1))))));
    end
end

cepstrum(~any(cepstrum,2), :) = [];
log_mag_spectrum(~any(log_mag_spectrum,2), :) = []; 

cepstrum_size = size(cepstrum);
cepstrum_size1 = cepstrum_size(1);
cepstrum_size2 = cepstrum_size(2);

figure
subplot(4,1,1);
plot(soundfile(0.94*Fs : 0.97*Fs, 1))
xlabel('Time')
title('Sound Waveform of voiced signal')

subplot(4,1,2)
plot(log_mag_spectrum(49, :))
title('Log Magnitude Spectrum of voiced signal')
xlabel('Frequency')

subplot(4,1,3)
plot(cepstrum(49,:))
ylim([-0.75 0.75]);
title('Real Cepstrum of voiced signal')
xlabel('Time')

low_lifter=zeros(1,cepstrum_size2);
for i=1:cepstrum_size2
    if 4.5*i<cepstrum_size2
        low_lifter(1,i) = cepstrum(49,i);
    else
        low_lifter(1,i) = 0;
    end
end

low_lifter = abs(fft(low_lifter));
subplot(4,1,4)
plot(low_lifter(1,:))
title('Low Quefrency Liftered Log Magnitude Spectrum of voiced signal')
xlabel('Quefrency')


figure
subplot(4,1,1);
plot(soundfile(0.65*Fs:0.68*Fs,1))
xlabel('Time')
title('Sound Waveform of unvoiced signal')

subplot(4,1,2)
plot(log_mag_spectrum(34,:))
title('Log Magnitude Spectrum of unvoiced signal')
xlabel('frequency')

subplot(4,1,3)
plot(cepstrum(34,:))
ylim([-0.75 0.75]);
title('Real Cepstrum of unvoiced signal')
xlabel('Time')

for i=1:cepstrum_size2
    if 3.5*i<cepstrum_size2
        low_lifter(1,i)=cepstrum(34,i);
    else
        low_lifter(1,i)=0;
    end
end

low_lifter = abs(fft(low_lifter));
subplot(4,1,4)
plot(low_lifter(1,:))
title('Low Quefrency Liftered Log Magnitude Spectrum of voiced signal')
xlabel('Quefrency')