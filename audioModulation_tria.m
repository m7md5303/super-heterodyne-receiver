
%getting sampling rates & bandwidths:
name = ["QuranPalestine","BBCArabic2"];
audio_name_stereo = genvarname({'audio_stereo', 'audio_stereo'});
audio_name = genvarname({'audio', 'audio'});
for k=1:2
    filename = sprintf('Short_%s.wav',name(k));
    [audio_name_stereo{k},fs(k)] = audioread(filename);
    used = audio_name_stereo{k};
    audio_name{k} = used(:,1) + used(:,2);
    NN = -fs(k)/2:fs(k)/2-1;
    ploted_wave = fft(audio_name{k},fs(k));
    subplot(7,3,k)
    plot(NN,fftshift(abs(ploted_wave)));
    title(sprintf('Audio%d',k));
end
%sound( audio_name{2},fs(2));
%sound(audio_name{2} , fs(1));
%padding the short signal with zeros:
for i = 1:2
if(length(audio_name{i})<length(audio_name{2/i}))
    y = zeros(length(audio_name{2/i}),1);
    y(1:length(audio_name{i}),1) = audio_name{i};
    audio_name{i}=y;
end   
end

%increasing the audio sampling frequency:
for k = 1:2
    x = audio_name{k};
         z(:,1) = interp(x(:,1),10);
    audio_name{k} = z;
end

%getting the new bandwidths:
NN = -10*fs(k)/2:10*fs(k)/2-1;
for k = 1:2
    ploted_wave = fft(audio_name{k},10*fs(k));
    subplot(7,3,k+2)
    plot(NN,fftshift(abs(ploted_wave)));
    title(sprintf('Audio%d resampled',k));
end
l = linspace(0,1,10*fs(1));

for i=0:1
%getting the carrier signal:
    for n = 1:round((length(audio_name{1})/(10*fs(1)))+0.5)
        carrier((n-1)*10*fs(1)+1:n*10*fs(1),1) = cos(2*pi*n*(1e5+i*55e3)*l);
    end
%getting a part of the carrier signal of a length equal to the modulating
%signal length:
    final_carrier = zeros(length(audio_name{1}),1);
    final_carrier(:,1) = carrier(1:length(audio_name{1}),1);
    audio_name{i+1} = audio_name{i+1}.*final_carrier;
    ploted_wave = fft(audio_name{i+1},10*fs(1));
    subplot(7,3,5+i)
    plot(NN,fftshift(abs(ploted_wave)));
    title(sprintf('Audio%d after modulating',i+1));
end
modulated_signal =audio_name{1}+audio_name{2};
ploted_wave = fft(modulated_signal,10*fs(k));
subplot(7,3,7)
plot(NN,fftshift(abs(ploted_wave)));
title("The transmitted signal");
%The RF Stage
%defining the BPF parameters
bp1=[88600   111400];
bp2=[146500  163500];
%filtering the signals each at its receiving station and with the
%corresponding filter range
filtered_signal1=bandpass(modulated_signal,bp1,10*fs(1));
filtered_signal1_fft = fft(filtered_signal1,10*fs(1));
subplot(7,3,9)
plot(NN,fftshift(abs(filtered_signal1_fft)));
title("signal 1 at the receiver after the bandpass stage to remove any probable image interference")
filtered_signal2=bandpass(modulated_signal,bp2,10*fs(2));
filtered_signal2_fft = fft(filtered_signal2,10*fs(2));
subplot(7,3,10)
plot(NN,fftshift(abs(filtered_signal2_fft)));
title("signal 2 at the receiver after the bandpass stage to remove any probable image interference")
%The Oscillator Stage
n=linspace(0,length(filtered_signal1)-1,(length(filtered_signal1)));
%defining the carrier for each station
%getting the carrier signal:
    for n = 1:round((length(filtered_signal1)/(10*fs(1)))+0.5)
        carrier_osc1((n-1)*10*fs(1)+1:n*10*fs(1),1) = cos(2*pi*n*(127500)*l);
    end
%getting a part of the carrier signal of a length equal to the modulating
%signal length:
    final_carrier_osc1 = zeros(length(filtered_signal1),1);
    final_carrier_osc1(:,1) = carrier_osc1(1:length(filtered_signal1),1);
    filtered_signal1 = filtered_signal1.*final_carrier_osc1;
    
%//////////////////////////////////////////////////////////////////
 for n = 1:round((length(filtered_signal2)/(10*fs(2)))+0.5)
        carrier_osc2((n-1)*10*fs(1)+1:n*10*fs(2),1) = cos(2*pi*n*(182500)*l);
    end
%getting a part of the carrier signal of a length equal to the modulating
%signal length:
    final_carrier_osc2 = zeros(length(filtered_signal2),1);
    final_carrier_osc2(:,1) = carrier_osc2(1:length(filtered_signal2),1);
    filtered_signal2 = filtered_signal2.*final_carrier_osc2;
filtered_signal1_fftosc = fft(filtered_signal1,10*fs(1));
subplot(7,3,12)
plot(NN,fftshift(abs(filtered_signal1_fftosc)));
title("signal 1 at the receiver after the oscillator stage")
filtered_signal2_fftosc = fft(filtered_signal2,10*fs(2));
subplot(7,3,13)
plot(NN,fftshift(abs(filtered_signal2_fftosc)));
title("signal 2 at the receiver after the oscillator stage")
%The IF Stage
%Constructing the IF filters parameters each for each receiving station
bp1_if=[20190  34810];
bp2_if=[18210  36790];
%applying the filters on the signals
filtered_signal1=bandpass(filtered_signal1,bp1_if,10*fs(1));
filtered_signal1_fft = fft(filtered_signal1,10*fs(1));
subplot(7,3,14)
plot(NN,fftshift(abs(filtered_signal1_fft)));
title("signal 1 after being mixed in the IF stage")
filtered_signal2=bandpass(filtered_signal2,bp2_if,10*fs(2));
filtered_signal2_fft = fft(filtered_signal2,10*fs(2));
subplot(7,3,15)
plot(NN,fftshift(abs(filtered_signal2_fft)));
title("signal 2 after being mixed in the IF stage")
%The Base Band detection Stage
%Constructing the mixers signals properly for each station at frequency
%W_IF
%getting the carrier signal:
    for n = 1:round((length(filtered_signal1)/(10*fs(1)))+0.5)
        carrier_osc1((n-1)*10*fs(1)+1:n*10*fs(1),1) = cos(2*pi*n*(27500)*l);
    end
%getting a part of the carrier signal of a length equal to the modulating
%signal length:
    final_carrier_osc1 = zeros(length(filtered_signal1),1);
    final_carrier_osc1(:,1) = carrier_osc1(1:length(filtered_signal1),1);
    filtered_signal1 = filtered_signal1.*final_carrier_osc1;
    
%//////////////////////////////////////////////////////////////////
 for n = 1:round((length(filtered_signal2)/(10*fs(2)))+0.5)
        carrier_osc2((n-1)*10*fs(1)+1:n*10*fs(2),1) = cos(2*pi*n*(27500)*l);
    end
%getting a part of the carrier signal of a length equal to the modulating
%signal length:
    final_carrier_osc2 = zeros(length(filtered_signal2),1);
    final_carrier_osc2(:,1) = carrier_osc2(1:length(filtered_signal2),1);
    filtered_signal2 = filtered_signal2.*final_carrier_osc2;
%//////////////////////////////////////////////////////////////////
filtered_signal1_fft = fft(filtered_signal1,10*fs(1));
subplot(7,3,16)
plot(NN,fftshift(abs(filtered_signal1_fft)));
title("signal 1 after being mixed in the base band detection phase")
filtered_signal2_fft = fft(filtered_signal2,10*fs(2));
subplot(7,3,17)
plot(NN,fftshift(abs(filtered_signal2_fft)));
title("signal 2 after being mixed in the base band detection phase")
%Determining the LPF parameters fore each station
lp1=15000;
lp2=18130;
%applying LPF for each station to retrieve the received signal
filtered_signal1=lowpass(filtered_signal1,lp1,10*fs(1));
filtered_signal1_fft = fft(filtered_signal1,10*fs(1));
subplot(7,3,18)
plot(NN,fftshift(abs(filtered_signal1_fft)));
title("retrieved signal 1")
filtered_signal2=lowpass(filtered_signal2,lp2,10*fs(2));
filtered_signal2_fft = fft(filtered_signal2,10*fs(2));
subplot(7,3,19)
plot(NN,fftshift(abs(filtered_signal2_fft)));
title("retrieved signal 2");
[numer , denom] = rat (1/10);
down_sampled1 = resample(filtered_signal1,numer,denom);
down_sampled2 = downsample(filtered_signal2,10);
NN = -fs(1)/2:fs(1)/2-1;
downsampled_plotted_wave = fft(down_sampled1,fs(1));
subplot(7,3,20)
plot(NN,fftshift(abs(downsampled_plotted_wave)));
title("Audio1 after downsampling");
downsampled_plotted_wave = fft(down_sampled2,fs(1));
subplot(7,3,21)
plot(NN,fftshift(abs(downsampled_plotted_wave)));
title("Audio2 after downsampling");
sound(down_sampled1 , fs(1));