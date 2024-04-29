clc
close all
clear all

%getting sampling rates & bandwidths:
name = ["QuranPalestine","BBCArabic2"];
audio_name_stereo = genvarname({'audio_stereo', 'audio_stereo'});
audio_name = genvarname({'audio', 'audio'});
for k=1:2
    filename = sprintf('Short_%s.wav',name(k));
    [audio_name_stereo{k},fs(k)] = audioread(filename);
    used = audio_name_stereo{k};
    audio_name{k} = used(:,1) + used(:,2);
    plot_freq_range = -fs(k)/2:fs(k)/2-1;
    ploted_wave = fft(audio_name{k},fs(k));
    subplot(7,3,k)
    plot(plot_freq_range,fftshift(abs(ploted_wave)));
    title(sprintf('Audio%d',k));
end

%padding the short signal with zeros:
for i = 1:2
if(length(audio_name{i})<length(audio_name{2/i}))
    padding_array = zeros(length(audio_name{2/i}),1);
    padding_array(1:length(audio_name{i}),1) = audio_name{i};
    audio_name{i}=padding_array;
end   
end

%increasing the audio sampling frequency:
for k = 1:2
    audio_name{k}= interp(audio_name{k},10);
end

%getting the new bandwidths:
plot_freq_range = -10*fs(k)/2:10*fs(k)/2-1;
for k = 1:2
    ploted_wave = fft(audio_name{k},10*fs(k));
    subplot(7,3,k+2)
    plot(plot_freq_range,fftshift(abs(ploted_wave)));
    title(sprintf('Audio%d resampled',k));
end


for i=0:1
%getting the carrier signal:

time_range=linspace(0,length(audio_name{i+1}),(length(audio_name{i+1})));
final_carrier=cos(2*pi*(1e5+i*55e3)*time_range*(1/(10*fs(1))));
final_carrier=transpose(final_carrier);
    audio_name{i+1} = audio_name{i+1}.*final_carrier;
    ploted_wave = fft(audio_name{i+1},10*fs(1));
    subplot(7,3,5+i)
    plot(plot_freq_range,fftshift(abs(ploted_wave)));
    title(sprintf('Audio%d after modulating',i+1));
end
modulated_signal =audio_name{1}+audio_name{2};
ploted_wave = fft(modulated_signal,10*fs(k));
subplot(7,3,7)
plot(plot_freq_range,fftshift(abs(ploted_wave)));
title("The transmitted signal");

%sound( yy2,fs(1));
%The RF Stage
%defining the BPF parameters
bp1=[88600   111400];
bp2=[146500  163500];
%filtering the signals each at its receiving station and with the
%corresponding filter range

filtered_signal1=bandpass(modulated_signal,bp1,10*fs(2));
filtered_signal1_fft = fft(filtered_signal1,10*fs(1));
subplot(7,3,9)
plot(plot_freq_range,fftshift(abs(filtered_signal1_fft)));
title("signal 1 at the receiver after the bandpass stage to remove any probable image interference")
filtered_signal2=bandpass(modulated_signal,bp2,10*fs(2));
filtered_signal2_fft = fft(filtered_signal2,10*fs(2));
subplot(7,3,10)
plot(plot_freq_range,fftshift(abs(filtered_signal2_fft)));
title("signal 2 at the receiver after the bandpass stage to remove any probable image interference")
%The Oscillator Stage
time_range=linspace(0,length(filtered_signal1),(length(filtered_signal1)));
%defining the carrier for each station
osc_carrier1=cos(2*pi*127500*time_range*(1/(10*fs(1))));
osc_carrier1=transpose(osc_carrier1);
osc_carrier2=cos(2*pi*182500*time_range*(1/(10*fs(2))));
osc_carrier2=transpose(osc_carrier2);
%modulating the received signals
filtered_signal1=filtered_signal1.*osc_carrier1;
filtered_signal2=filtered_signal2.*osc_carrier2;
filtered_signal1_fft = fft(filtered_signal1,10*fs(1));
subplot(7,3,12)
plot(plot_freq_range,fftshift(abs(filtered_signal1_fft)));
title("signal 1 at the receiver after the oscillator stage")
filtered_signal2_fft = fft(filtered_signal2,10*fs(2));
subplot(7,3,13)
plot(plot_freq_range,fftshift(abs(filtered_signal2_fft)));
title("signal 2 at the receiver after the oscillator stage")
%The IF Stage
%Constructing the IF filters parameters each for each receiving station
bp1_if=[20190  34810];
bp2_if=[18210  36790];
%applying the filters on the signals
filtered_signal1=bandpass(filtered_signal1,bp1_if,10*fs(1));
filtered_signal1_fft = fft(filtered_signal1,10*fs(1));
subplot(7,3,14)
plot(plot_freq_range,fftshift(abs(filtered_signal1_fft)));
title("signal 1 after being mixed in the IF stage")
filtered_signal2=bandpass(filtered_signal2,bp2_if,10*fs(2));
filtered_signal2_fft = fft(filtered_signal2,10*fs(2));
subplot(7,3,15)
plot(plot_freq_range,fftshift(abs(filtered_signal2_fft)));
title("signal 2 after being mixed in the IF stage")
%The Base Band detection Stage
%Constructing the mixers signals properly for each station at frequency
%W_IF
osc_carrier1=cos(2*pi*27500*time_range*(1/(10*fs(1))));
osc_carrier1=transpose(osc_carrier1);
osc_carrier2=cos(2*pi*27500*time_range*(1/(10*fs(2))));
osc_carrier2=transpose(osc_carrier2);
%mixing the received signals
filtered_signal1=filtered_signal1.*osc_carrier1;
filtered_signal2=filtered_signal2.*osc_carrier2;
filtered_signal1_fft = fft(filtered_signal1,10*fs(1));
subplot(7,3,16)
plot(plot_freq_range,fftshift(abs(filtered_signal1_fft)));
title("signal 1 after being mixed in the base band detection phase")
filtered_signal2_fft = fft(filtered_signal2,10*fs(2));
subplot(7,3,17)
plot(plot_freq_range,fftshift(abs(filtered_signal2_fft)));
title("signal 2 after being mixed in the base band detection phase")
%Determining the LPF parameters fore each station
lp1=18260;
lp2=16080;
%applying LPF for each station to retrieve the received signal
filtered_signal1=lowpass(filtered_signal1,lp1,10*fs(1));
filtered_signal1_fft = fft(filtered_signal1,10*fs(1));
subplot(7,3,18)
plot(plot_freq_range,fftshift(abs(filtered_signal1_fft)));
title("retrieved signal 1")
filtered_signal2=lowpass(filtered_signal2,lp2,10*fs(2));
filtered_signal2_fft = fft(filtered_signal2,10*fs(2));
subplot(7,3,19)
plot(plot_freq_range,fftshift(abs(filtered_signal2_fft)));
title("retrieved signal 2");

%downsampling and amplifying the signals:
down_sampled1 = 4*downsample(filtered_signal1,10);
down_sampled2 = 4*downsample(filtered_signal2,10);
plot_freq_range = -fs(1)/2:fs(1)/2-1;
downsampled_plotted_wave = fft(down_sampled1,fs(1));
subplot(7,3,20)
plot(plot_freq_range,fftshift(abs(downsampled_plotted_wave)));
title("Audio1 after downsampling");
downsampled_plotted_wave = fft(down_sampled2,fs(1));
subplot(7,3,21)
plot(plot_freq_range,fftshift(abs(downsampled_plotted_wave)));
title("Audio2 after downsampling");
sound(down_sampled1 , fs(1));