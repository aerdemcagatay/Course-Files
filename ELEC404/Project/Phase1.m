clc
clear
close all

%%--------------------ENCODER PART---------------%%

prompt='Enter the audio file name(ex:s1010-010):  ';%Taking audio file name as an input
audioName=input(prompt,'s');
audioName=[audioName,'.wav'];%Add .wav extension to name
[y,fs]=audioread(audioName);%Read datas and sampling frequency from selected audio file
initiallength=length(y);
FrameSize=fs*0.02;%Converting frame size from seconds to samples
FrameShift=fs*0.01;%Convering frame shift from seconds to samples
NumberofFrames=ceil((length(y)-FrameSize)/FrameShift)+1;%Finding number of samples
paddedSize=(NumberofFrames-1)*FrameShift+FrameSize;%Adjusting padding size
y=[y;zeros(paddedSize-initiallength,1)];%Padding our signal with zeros so it size will be multiple of framesize
p=10;%LPC order

Parcor=zeros(NumberofFrames,p);%Initializing parcor coefficients array
Gains=zeros(NumberofFrames,1);%Initializing gain array
Excitations=zeros(NumberofFrames,FrameSize);%Initializing excitations array

win=hanning(FrameSize);%Creating hamming window

for n=1:(paddedSize-FrameSize)/FrameShift+1;%Used for incrementing frame number
    x=y((n-1)*FrameShift+1:(n-1)*FrameShift+FrameSize);%Selecting the frame
    x=x.*win;%Windowing the frame
    
    %Calculating autocorrelation vector
    R=zeros(1,p+1);
    for i=0:p
        R(1,i+1)=x(1:FrameSize-i)'*x(i+1:FrameSize);
    end
    
    %%Levinson Durbin Algorithm
    E=zeros(1,p+1);%%Short time averaged error
    k=zeros(1,p);%%Parcor Coefficients
    alpha=zeros(p,p);%%Initiliazing alphas array
    %In this part since in matlab there isn't index zero in matlab, i will
    %used index+1 values for every calculation.I implemented the levinson
    %durbin algorithm from lecture notes pseudocode.
    E(1)=R(1);%%E(0)=R(0);
    k(1)=R(2)/E(1);
    alpha(1,1)=k(1);
    E(2)=(1-k(1)^2)*E(1);
    
    for i=2:p
        term=0;
        for j=1:i-1
            term=term+alpha(j,i-1)*R(i-j+1);
        end
        k(i)=(R(i+1)-term)/E(i);
        alpha(i,i)=k(i);
        for j=1:i-1
            alpha(j,i)=alpha(j,i-1)-k(i)*alpha(i-j,i-1);
        end
        E(i+1)=(1-k(i)^2)*E(i);
    end
    GainforFrame=sqrt(E(p+1));%Calculating gain value
    Numerator=[1;-alpha(1:p,p)];%Creating transfer function denominator arrays
    excitationforFrame=filter(Numerator,1,x);%Creating excitaion signal
    
    %%OUTPUTS OF ENCODER:
    Parcor(n,1:p)=k(1:p);%adding parcor coefficients to main parcor matrix with frame number
    Gains(n)=GainforFrame;%Adding gain values to main  gain matrix with frame number
    Excitations(n,:)=excitationforFrame;%Adding excitiation signals to main excitation matrix with frame number
    
end

%%--------------DECODER PART------------%%
%%Since our decoder inputs are parcor coefficients,gain values and
%%excitation signals.Our decoder can only use these three datas.
OutputNumberofFrame=length(Gains);%Obtaining output total number of frames
OutputLength=(OutputNumberofFrame-1)*FrameShift+FrameSize;%Obtaining output file length
RecoveredSignal=zeros(OutputLength,1);%Initializing recovered signal
ExcitationSignal=zeros(OutputLength,1);%Initializing final combined excitation signal for plotting
for m=1:OutputNumberofFrame%For shifting frames
     parc=Parcor(m,:);%Taking parcor coefficients for present frame
    %Calculating alpha values then numerator of filter with Levinson Durbin
     AlphaSynthesis(1,1)=parc(1);
    for i=2:p
        AlphaSynthesis(i,i)=parc(i);
        for j=1:i-1
            AlphaSynthesis(j,i)=AlphaSynthesis(j,i-1)-parc(i)*AlphaSynthesis(i-j,i-1);
        end
    end
    Denominator=[1;-AlphaSynthesis(1:p,p)];%Calculating numerator matrix coeffients
    z=filter(1,Denominator,Excitations(m,:));%Filtering excitation signals
    %Creating recovering signal with overlap and add method
    RecoveredSignal((m-1)*FrameShift+1:(m-1)*FrameShift+FrameSize)=RecoveredSignal((m-1)*FrameShift+1:(m-1)*FrameShift+FrameSize)+z';
    %Creating combined excitation signal for plotting
    ExcitationSignal((m-1)*FrameShift+1:(m-1)*FrameShift+FrameSize)=Excitations(m,:);
    
end

%%PLOTTING
t=0:1/fs:(length(y)-1)/fs;

subplot(3,1,1);
plot(t,y);
title('Original Signal');
xlabel('Time(sec)');
ylabel('Amplitude');


subplot(3,1,2);
plot(t,ExcitationSignal);
title('Excitation Signal');
xlabel('Time(sec)');
ylabel('Amplitude');

subplot(3,1,3);
plot(t,RecoveredSignal);
title('Recovered Signal');
xlabel('Time(sec)');
ylabel('Amplitude');

%Calculating segSNR value with given formula
segSNR=0;
for n=1:NumberofFrames
    inputSigFrame=y((n-1)*FrameShift+1:(n-1)*FrameShift+FrameSize);
    outputSigFrame=RecoveredSignal((n-1)*FrameShift+1:(n-1)*FrameShift+FrameSize);
    a=sum(inputSigFrame.*inputSigFrame.*2);
    dif=inputSigFrame-outputSigFrame;
    b=sum(dif.*dif);
    segSNR=segSNR+10*log10(a/b);
end
segSNR=segSNR/NumberofFrames;%segSNR value
disp(['segSnr value for ', audioName,' is equals to ',num2str(segSNR)]);
sound(y,fs);
pause(5);
sound(RecoveredSignal,fs);