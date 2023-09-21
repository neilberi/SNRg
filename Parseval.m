% Script to verify Parseval expectation for DFT coefficient magnitudes
% using the Matlab fft command conventions

% Expectation: for Gaussian white noise in the time domain with 
% a mean of zero and a standard deviation sigma, the real and imaginary
% components of the 2-sided DFT coefficients should be sigma *sqrt(Tcoh*fsamp/2)

% Show that this relation holds for various combinations of
% coherence time, sampling frequency and sigma

clear;
clc;
close all;

% Create sets of Tcoh, fsamp and sigma values
Tcohvals =  [100 1000 10000];
fsampvals = [512 1024 2048];
sigmavals = [0.5 2.0 8.0];
nvals = length(Tcohvals);

% Loop over sets
for i = 1:nvals

   % Get parameters and derived quantities for this set
   Tcoh = Tcohvals(i);
   fsamp = fsampvals(i);
   sigma = sigmavals(i);
   df = 1./Tcoh;
   deltat = 1./fsamp;
   t = [0:deltat:Tcoh];
   Nsample = length(t);
   fprintf('Set %d: (Tcoh = %f, fsamp = %f, sigma = %f)\n',i,Tcoh,fsamp,sigma);

   % Create time-domain data
   noise = sigma*random('norm',0.,1.,1,Nsample);

   % Compute DFT coefficients and compute powers
   rawfft = fft(noise,Nsample);
   power_2sided = abs(rawfft).^2;
   power_1sided = 2*power_2sided;

   % Plot powers vs frequency
   figure;
   freq = 0:df:fsamp/2.;
   plot(freq,power_1sided(1:length(freq)),'color','red');
   hold on

   % Get predicted and measured 2-sided real, imag parts and 1-sided powers (note: the 2-sided powers are half of 1-sided powers)
   pred_std = sqrt(0.5*Tcoh*fsamp)*sigma;
   real_std = std(real(rawfft(1:length(freq))));
   imag_std = std(imag(rawfft(1:length(freq))));
   power_mean_1sided = mean(power_1sided(1:length(freq)));
   fprintf('Predicted std dev of 2-sided real/imag coeff = %f, measured real coeff = %f, imag coeff = %f\n',pred_std,real_std,imag_std);
   power_predict_2sided = 2.0 * pred_std^2; 
   power_predict_1sided = 2.0 * power_predict_2sided;
   fprintf('Predicted mean 1-sided power = %f, Measured mean 1-sided power = %f\n',power_predict_1sided,power_mean_1sided);

   % Overlay predicted and measured powers
   
   plot([freq(1) freq(length(freq))],[power_predict_1sided power_predict_1sided],'color','green')
   plot([freq(1) freq(length(freq))],[power_mean_1sided power_mean_1sided],'color','black')
   legend('Measured 1-sided power','Parseval prediction','Mean measured');
   xlabel('Freq (Hz)');
   ylabel('DFT coefficients squared times two');
   titlestr = sprintf('Comparison of predicted to measured DFT power for Tcoh=%d s, fsamp=%d Hz, sigma = %f',Tcoh,fsamp,sigma);
   title(titlestr);

   % Plot histograms of real, imag parts and powers
   figure
   Nbin = 40;
   subplot(3,1,1);
   histogram(real(rawfft(1:length(freq))),Nbin);
   titlestr = sprintf('Real part of DFT, predict std=%f, meas=%f',pred_std,real_std);
   title(titlestr);
   subplot(3,1,2);
   histogram(imag(rawfft(1:length(freq))),Nbin);
   titlestr = sprintf('Imag part of DFT, predict std=%f, meas=%f',pred_std,imag_std);
   title(titlestr);
   subplot(3,1,3);
   histogram(power_1sided(1:length(freq)),Nbin);
   titlestr = sprintf('DFT power (1-sided), predict mean=%e, meas=%e',power_predict_1sided,power_mean_1sided);
   title(titlestr);
end



   
