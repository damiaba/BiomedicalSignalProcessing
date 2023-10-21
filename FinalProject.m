%Final Project
load('ecgproject.mat')
f = 200; %ECG Signal frequency
T = 1/f;  %Period of signal
x = ecg1*(1000/500); %ECG signal
s = length(x)*T-0.001; %signal end time
n = [0:T:s]; %signal time

figure
plot(n,x)
title('ECG Signal vs. Time')
xlabel('Time(s)')
ylabel('Voltage(mV)') %Plots ECG signal vs time

N = length(x); %Points for FFT
X = (fft(x,N));
X2 = abs(fftshift(X)); %Creates spectrum plot
F = f*[-N/2:N/2-1]/N; %frequency axis 
figure
plot(F,X2) %Plot for spectrum analysis
title('Spectrum of ECG signal')
xlabel('Frequency(Hz)')
ylabel('ECG signal phase (A.U)')
xlim([2 100])

%From inspection can see that the signal has a useful frequency content
%between 0.05-15 Hz
nf = f/2; %Nyquist freq
wo = 50/nf;%removes interference at 50 Hz
q = 35; %Q factor
bw = wo/q; %bandwith
[b,a] = iirnotch(wo,bw); %creates notch filter
s1 = filtfilt(b,a,x); %filters signal
Wp = [3.8/nf 15/nf]; %passband freq
Ws = [0.00001 0.999999]; %stopband freq
Rp = 0.5; %Ripple passband
Rs = 60; %Lowpass attenuation
[n1,Wc] = ellipord(Wp,Ws,Rp,Rs); %creates bandpass elliptic order filter
[b1, a1] = ellip(n1,Rp,Rs,Wc);
ecg2 = filtfilt(b1,a1,s1);  %filters notch


figure
plot(n,ecg2)
title('Filtered ECG Signal vs. Time')
xlabel('Time(s)')
ylabel('Voltage(mV)') %Plots filtered ECG signal vs time

%ECG Manipulation to find and segment data
absv = abs(ecg2); %rectification of signal
Wn = 5/nf;
[B,A] = butter(3, Wn); %3rd order Butterworth filter
esignal1 = filtfilt(B,A,absv);

l1 = 0.25;%time of noise before real signal activation
esignal = esignal1;
thresholdv = mean(esignal(n<l1)) + 3*std(esignal(n<l1));%threshold value 
actsignal = esignal > thresholdv;
actsignal(n < l1) = 0; %creates activation signal
me = esignal/max(esignal); %normalizes the enveloped signal 
acttime = n(find(diff(actsignal)== 1)); %points of activation
%%ti = 0.47;

%Plots enveloped signal
figure
plot(n,me)%normalized envelope signal
hold on
plot(n,actsignal)
xlabel('Time(s)')
ylabel('Voltage(mV)')
title('Envelope and Activation Signal vs Time')
legend('Envelope Signal','Activation Signal')

%plots ECG signal with segmented R waves
figure
plot(n,ecg2)
hold on
plot(n,actsignal.*max(ecg2))
xlabel('Time(s)')
ylabel('Voltage(mV)')
title('ECG and Activation Signal vs Time')
legend('ECG Signal','Activation Signal')

%2 features
[peaks,locs]=findpeaks(ecg2,f,'MinPeakHeight',2.325*thresholdv);
[peaks1,locs1]=findpeaks(ecg2,f,'MinPeakDistance',0.47);

%time variabilty feature
figure
plot(n,ecg2)
hold on
plot(n,actsignal.*max(ecg2))
hold on
plot(locs1,peaks1,'o')
xlabel('Time(s)')
ylabel('Voltage(mV)')
title('ECG and Activation ECG Signal vs Time')
legend('ECG Signal','Activation Signal','Beats')

%peak-height feature
figure
plot(n,ecg2)
hold on
plot(n,actsignal.*max(ecg2))
hold on
plot(locs,peaks,'o')
xlabel('Time(s)')
ylabel('Voltage(mV)')
title('Envelope and Activation ECG Signal vs Time')
legend('Enevelope Signal','Activation Signal','Beats')

%Classification Boundary Line
clb = ones(400)*3*thresholdv;


figure
plot(n,ecg2)
hold on
plot(n,actsignal.*max(ecg2))
hold on
plot(locs1,peaks1, 'o')
hold on
plot(locs1,clb,'--')
xlabel('Time(s)')
ylabel('Voltage(mV)')
title('ECG and Activation ECG Signal vs Time')
legend('Enevelope Signal','Activation Signal','Beats', 'Classification boundary line')


figure
plot(n,ecg2)
hold on
plot(n,actsignal.*max(ecg2))
hold on
plot(locs,peaks, 'o')
hold on
plot(locs,clb,'--')
xlabel('Time(s)')
ylabel('Voltage(mV)')
title('ECG and Activation ECG Signal vs Time')
legend('Enevelope Signal','Activation Signal','Beats', 'Classification boundary line')

%Classification Graph
figure
plot(peaks,locs,'o')
hold on
plot(clb,locs,'--')
xlabel('Location(A.U)')
ylabel('Beats(A.U)')
title('Classification Graph')
legend('Beats','Classification boundary line')

figure
plot(locs,peaks,'o')
hold on
plot(locs,clb,'--')
xlabel('Time(s)')
ylabel('Beat Amplitude(mV)')
title('Classification Graph')
legend('Beats','Classification boundary line')

%for feature 1
pvcp1 = peaks1 > 3*thresholdv; %Most PVC beats would be above this.
pvcp1 = pvcp1';

%for feature 2
pvcp = peaks > 3*thresholdv; %Most PVC beats would be above this.
pvcp = pvcp';

match = 0;
match1 = 0;

TP = 0;
FP = 0;
TN = 0;
FN = 0;

lengthy = [1:400];

for i = lengthy
    if pvcp(i) == isPVC(i)
        match = match +1;
    end  
end %Checks performance of test set

accuracy = (match/400)*100;
%printf('The classification model performance of feature 2 is %1.3f % accurate.',accuracy)
disp(accuracy)

nonPVC = peaks < 3*thresholdv;
nonPVC = nonPVC';

for i = lengthy
    if pvcp(i) == isPVC(i)
        match1 = match1 +1;
    end  
end %Checks performance of test sets

accuracy1 = (match1/400)*100;
disp(accuracy1)

%performances were the same so going with classification method using
%height
for i = lengthy
    if (pvcp(i) == 1 && isPVC(i) == 1)
        TP = TP+1;
    elseif (pvcp(i) == 1 && isPVC(i) == 0)
        FP = FP+1;
    elseif (pvcp(i) == 0 && isPVC(i) == 0)
        TN = TN+1;
    elseif (pvcp(i) == 0 && isPVC(i) == 1)
        FN = FN+1;
    end
end
TP
FP
TN
FN

sensitivity = TP/(FP+FN)%chance of identifying an actual PVC
specifity = TN/(FP+TN)%chance of identifying a non PVC
accuracy = (TP+TN)/(TN+FP+FN+FP)%The accuracy of the classification method