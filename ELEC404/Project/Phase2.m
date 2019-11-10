clc
clear
close all

%%--------------------ENCODER PART---------------%%

prompt='Enter the audio file name(ex:s1010-010):  ';%Taking audio file name as an input
audioName=input(prompt,'s');
audioName=[audioName,'.wav'];%Add .wav extension to name
[y,fs]=audioread(audioName);%Read datas and sampling frequency from selected audio file

prompt='Select codec with codec number(ex:1):   ';%%Selecting from 4 different codec
codecnumber=input(prompt);
codecs=[40 480
        40 240
        20 480
        20 240];
numbitforPARCOR=codecs(codecnumber,1);%Assigning number of bit allocated for parcor per frame
numbitforEXC=codecs(codecnumber,2);%Assigning number of bit allocated for excitation per frame

initiallength=length(y);
FrameSize=fs*0.02;%Converting frame size from seconds to samples
FrameShift=fs*0.01;%Convering frame shift from seconds to samples
NumberofFrames=ceil((length(y)-FrameSize)/FrameShift)+1;%Finding number of samples
paddedSize=(NumberofFrames-1)*FrameShift+FrameSize;%Adjusting padding size
y=[y;zeros(paddedSize-initiallength,1)];%Padding our signal with zeros so it size will be multiple of framesize
p=10;%LPC order


%%--------------------ENCODER PART---------------%%
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
    
    
    
    %%Outputs Of ENCODER(Not quantized yet):
    Parcor(n,1:p)=k(1:p);%adding parcor coefficients to main parcor matrix with frame numbers
    Gains(n)=GainforFrame;%Adding gain values to main gain matrix with frame numbers
    Excitations(n,:)=excitationforFrame;%Adding excitiation signals to main excitation matrix with frame numbers
   
end
%%--------------ENCODER QUANTÝZATÝON----------%%
%I used one kmaps for all values in excitation, therefore i have only one
%code book for all frames.
quantizedExcALL=[];%%Initializing excitation array for kmean
for i=1:NumberofFrames%%Adding all frames into one row for using kmap
    quantizedExcALL=[quantizedExcALL,Excitations(i,:)];
end
forcodebook=[];
%%Since we used vectoral quantization 2D, we split our excitation signal
%%into 2 dimension array with using X(2n) and X(2n+1 values)
for i=1:length(quantizedExcALL)/2
    forcodebook=[forcodebook;quantizedExcALL(2*i-1),quantizedExcALL(2*i)];
end

nbitEXC=2*numbitforEXC/FrameSize;%Bit per two sample calculated
[id_Exc,C_Exc]=kmeans(forcodebook,2^nbitEXC,'MaxIter',100000);%Creating clusters
quantizedExcitations=zeros(NumberofFrames,FrameShift);%Initializing array
id_Exc=id_Exc';
%Creating quantizied excitation arrays with frame numbers
for i=1:NumberofFrames
    quantizedExcitations(i,:)=id_Exc((i-1)*FrameShift+1:((i-1)*FrameShift+FrameShift));
end

forparcorcodebook=[];%%Initializing Parcor array for kmean
Parcor=asin(Parcor);
for i=1:NumberofFrames%%Adding all frames into one row for using kmap
    forparcorcodebook=[forparcorcodebook,Parcor(i,:)];
end
nbitPARCOR=numbitforPARCOR/p;%Bit per parcor coefficients calculated
[id_Par,C_Par]=kmeans(forparcorcodebook',2^nbitPARCOR,'MaxIter',1000);%Creating clusters
quantizedParcors=zeros(NumberofFrames,p);%Initializing array
id_Par=id_Par';
%Creating quantizied parcor arrays with frame numbers
for i=1:NumberofFrames
    quantizedParcors(i,:)=id_Par((i-1)*p+1:((i-1)*p+p));
end

%%--------------DECODER PART------------%%
%%Since our decoder inputs are quantizied parcor coefficients,gain values and
%%quantizied excitation signals.Our decoder can only use these three datas.
OutputNumberofFrame=length(Gains);%Obtaining output total number of frames
OutputLength=(OutputNumberofFrame-1)*FrameShift+FrameSize;%Obtaining output file length
RecoveredSignal=zeros(OutputLength,1);%Initializing recovered signal
ExcitationSignal=zeros(OutputLength,1);%Initializing final combined excitation signal for plotting


for m=1:OutputNumberofFrame

    SynthesQuantExc=quantizedExcitations(m,:);%Obtaining quantizied excitation signals for this frame number
    SynthesisExcitation=zeros(1,FrameSize);%Initializing excitation signal for this frame number
    %Obtaining real values of Excitation signals
    
    for k=1:FrameShift
        SynthesisExcitation((2*k-1):(2*k))=C_Exc(SynthesQuantExc(k));%%Obtaining quantized values from binary values to centroid values
    end
    
    SynthesQuantParcor=quantizedParcors(m,:);%Obtaining quantizied parcor coefficients for this frame number
    SynthesParcor=zeros(1,p);%Initializing parcor coefficient for this frame number
    %Obtaining quantized values of Parcor Coefficients
    for k=1:p
        SynthesParcor(k)=C_Par(SynthesQuantParcor(k));
    end
    parc=SynthesParcor;%Taking parcor coefficients for present frame
    parc=sin(parc);
    %Calculating alpha values then numerator of filter with Levinson Durbin
    AlphaSyn(1,1)=parc(1);
    for i=2:p
        AlphaSyn(i,i)=parc(i);
        for j=1:i-1
            AlphaSyn(j,i)=AlphaSyn(j,i-1)-parc(i)*AlphaSyn(i-j,i-1);
        end
    end
    Denominator=[1;-AlphaSyn(1:p,p)];%Calculating numerator matrix coeffients   
    z=filter(1,Denominator,SynthesisExcitation);%Filtering excitation signals
    
    %Creating recovering signal with overlap and add method
    RecoveredSignal((m-1)*FrameShift+1:(m-1)*FrameShift+FrameSize)=RecoveredSignal((m-1)*FrameShift+1:(m-1)*FrameShift+FrameSize)+z';
    ExcitationSignal((m-1)*FrameShift+1:(m-1)*FrameShift+FrameSize)=SynthesisExcitation;
end

%%Plotting results
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
disp(['segSnr value for ', audioName,' with codec ',num2str(codecnumber),' is equals to ',num2str(segSNR)]);
sound(y,fs);%Playing the original sound
pause(5);
sound(RecoveredSignal,fs);%Playing the recovered sound