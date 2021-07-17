N = 197;                 % order of the filter 
aph = (N-1) / 2;         % the shift factor 
Fs = 13000;              % sampling frequency
fc1 = 2000;              % low cutoff frequency 
fc2 = 6000;              % higher cutoff frequency 

n = [0:N-1];
wc1 = fc1*2*pi;
wc2 = fc2*2*pi;
w = 0:.01:pi;

% loading the audio signal
[x, Fso] = audioread("bird_chirp.wav");
xf = linspace(0, Fso, Fs);
mf = abs(fft(x, Fs));


figure(1)
plot(xf(1:Fs/2), mf(1:Fs/2));
title("plot of audio signal in frequency domain");
xlabel("frequency")
ylabel("magnitude of X(k)")


% filterDesigner  to design the filter with blackman window with 197 filter
% order

hd = (sin(wc2*(n-aph)-sin(wc1*(n-aph)))) ./ (pi*(n-aph));   %hd(n)
wn = blackman(N);  % w(n)
hn = hd .* wn';   % h(n)

figure(2);

subplot(2, 1, 1)
plot(n, hn)
title("plot of hn in time domain")
xlabel("n")
ylabel("h(n)")

subplot(2, 1, 2)
h = freqz(hn,1,w);
plot(w/pi, abs(h));
title("plot of hn in frequency domain")
xlabel("n")
ylabel("H")





% Create the window vector for the design algorithm.
win = blackman(N+1);
% Calculate the coefficients using the FIR1 function.
b  = fir1(N, [fc1 fc2]/(Fs/2), 'bandpass', win, 'scale');
y = filter(b, 1, x);

audiowrite('filtered.wav', y, Fso);


figure(3)
subplot(2, 1, 1)
plot(x)
title('time domain plot of the original signal')
subplot(2, 1, 2)
plot(y)

title('time domain plot of filtered signal')



figure(4)
subplot(2, 1, 1)
plot(w/pi, abs(freqz(x, 1, w)))
ylim([0 500])
title('frequency domain plot of original signal')
subplot(2, 1, 2)
plot(w/pi, abs(freqz(y, 1, w)))
ylim([0 500])
title('frequency domain plot of filtered signal')