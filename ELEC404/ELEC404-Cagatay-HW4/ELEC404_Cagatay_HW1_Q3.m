clear
close all

[soundfile, Fs]=audioread('sx229.wav');
soundfile_size=size(soundfile);
soundfile_size=soundfile_size(1);

%(a)
figure 

L=1500;
K=1000;

autocorrelation = zeros(1,K);
mautocorrelation = zeros(1,K);

%Calculating autocorrelation and modified autocorrelation
r = 0.15;
for i=1:K
    autocorrelation(1,i)=sum(soundfile(r*Fs+i : r*Fs+L).*soundfile(r*Fs : r*Fs+L-i));
end
subplot(6,1,1)
plot(autocorrelation);
title("Short-Time Autocorrelation from " +r+ "th second");

for i=1:K
    mautocorrelation(1,i)=sum(soundfile(r*Fs+i : r*Fs+L+i-1).*soundfile(r*Fs : r*Fs+L-1));
end
subplot(6,1,2)
plot(mautocorrelation, 'r');
title("Modified Short-Time Autocorrelation from " +r+ "th second");

r = 0.50;
for i=1:K
    autocorrelation(1,i)=sum(soundfile(r*Fs+i : r*Fs+L).*soundfile(r*Fs : r*Fs+L-i));
end
subplot(6,1,3)
plot(autocorrelation);
title("Short-Time Autocorrelation from " +r+ "th second");

for i=1:K
    mautocorrelation(1,i)=sum(soundfile(r*Fs+i : r*Fs+L+i-1).*soundfile(r*Fs : r*Fs+L-1));
end
subplot(6,1,4)
plot(mautocorrelation, 'r');
title("Modified Short-Time Autocorrelation from " +r+ "th second");

r = 0.75;
for i=1:K
    autocorrelation(1,i)=sum(soundfile(r*Fs+i : r*Fs+L).*soundfile(r*Fs : r*Fs+L-i));
end
subplot(6,1,5)
plot(autocorrelation);
title("Short-Time Autocorrelation from " +r+ "th second");

for i=1:K
    mautocorrelation(1,i)=sum(soundfile(r*Fs+i : r*Fs+L+i-1).*soundfile(r*Fs : r*Fs+L-1));
end
subplot(6,1,6)
plot(mautocorrelation, 'r');
title("Modified Short-Time Autocorrelation from " +r+ "th second");

%(b)

rx=15;
window = ceil((soundfile_size-L-K)/rx); 
p = zeros(1, window);
th = 0.3;

for i = 0 : window-1
    
    for k = 1:K
        autocorrelation(1,k) = sum(soundfile(i*rx+1+k : i*rx+L+k).*soundfile(i*rx+1 : i*rx+L));
    end
    
    higher_than_th = autocorrelation > th;
    index = 2;
    
    if higher_than_th(1)~=0
        while index<length(higher_than_th)
            if higher_than_th(index) == 1 && higher_than_th(index-1) == 0
               p(1,i+1) = index;
               break
            end
            index=index+1;
        end

    else
        while index<length(higher_than_th)
            if higher_than_th(index) == 1 && higher_than_th(index-1) == 0
                p = index; 
                break
            end
            index = index+1;
        end
        while index<length(higher_than_th)
            if higher_than_th(index) == 1 && higher_than_th(index-1) == 0
                p(1,i+1) = index-p;
                break
            end
            index = index+1;
        end
    end
end

figure
plot((1:length(p))*rx/Fs, p*1000/Fs);
title("Pitch for threshold: " +th)
