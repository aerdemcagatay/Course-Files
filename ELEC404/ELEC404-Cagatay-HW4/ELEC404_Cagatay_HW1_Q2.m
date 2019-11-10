clear
close all

[soundfile, Fs] = audioread('sx229.wav');

soundfile_size = size(soundfile);

t = [1:soundfile_size] / Fs;

figure

%Drawing Speech Waveform
subplot(4, 1, 1);
plot(t, soundfile);
title('Time Waveform');
xlabel('Time');
ylabel('Amplitude');


%Creating the Hamming Window
L = 401;
R = 20;
hamming_window = hamming(L);

%short-time energy
energy = zeros(soundfile_size(1),1);

for i=1:R:soundfile_size
    n_minus_m = L + i - 1 - soundfile_size;
    if (n_minus_m > 0)
        energy(i, 1)=sum((hamming_window(i:n_minus_m).*soundfile(i:n_minus_m)).^2);
    else
    energy(i, 1)=sum((hamming_window.*soundfile(i:(L+i-1))).^2);
    end
end

subplot(4,1,2)
plot(t, energy);
title('Short-Time Energy')
ylabel('Energy')
xlabel('Time')

%short-time magnitude
magnitude=zeros(soundfile_size(1),1);

for i=1:R:soundfile_size
    n_minus_m = L + i - 1 - soundfile_size;
    if (n_minus_m > 0)
        magnitude(i, 1)=sum(abs(hamming_window(i:n_minus_m).*soundfile(i:n_minus_m)));
    else
        magnitude(i, 1)=sum(abs(hamming_window.*soundfile(i:(L+i-1))));
    end
end

subplot(4,1,3)
plot(t, magnitude);
title('Short-Time Magnitude')
ylabel('Magnitude')
xlabel('Time')

%short-time zero crossing rate
zcrate=zeros(soundfile_size(1),1);

for i=1:R:soundfile_size
    n_minus_m = L + i - 1 - soundfile_size;
    if (n_minus_m > 0)
        for k = i:(n_minus_m - 1)
            zcrate(i, 1)=zcrate(i, 1) + abs(sign(soundfile(k)) - sign(soundfile(k+1)))/L;
        end
    else
        for j=i:(L-1+i)
            zcrate(i, 1) = zcrate(i, 1)+ abs(sign(soundfile(j)) - sign(soundfile(j+1)))/L;
        end
    end
end

subplot(4,1,4)
plot(t, zcrate);
title('Short-Time Zero Crossing Rate')
xlabel('Time')
ylabel('Zero Crossing Rate')