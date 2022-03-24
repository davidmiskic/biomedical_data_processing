% accept the name of a record (.dat and .hea files of the database) 
% calculate a spectrogram (short-time Fourier transform) for each EHG signal of this record for chosen preprocessing filter
% calculate two estimators for each spectrogram: median frequency of the power spectrum, and peak frequency of the power spectrum
% plot the time course of each estimator along the spectrogram for all
% three EHG signals for chosen preprocessing filter - for each column

% sample entropy : measure of regularity - possibility that is not random
function Uterine(record, nF, preproccessingMode)
% => nF: high - better frequency resolution, low - better time resolution
% manjša variabilnost frekvence -> verjetno preterm
Fs = 20;
S = load(record);
[nsig, siglen] = size(S.val);
% 1, 5, 9 : unfiltered
% 2, 6, 10 : butter 0.08-4
% 3, 7, 11 : butter 0.3-3
% 4, 8, 12 : butter 0.3-4
sig0 = S.val(1 + preproccessingMode,:);
sig1 = S.val(5 + preproccessingMode,:);
sig2 = S.val(9 + preproccessingMode,:);

% power + spectrogram
% pspectrum(S.val(1,:), Fs, 'spectrogram')
% power spectra
% pspectrum(S.val(1,:), Fs, 'power')
% median of power spectrum
% medfreq(x, Fs)

% [pow, freq] = pspectrum(S.val(1,:), Fs)
% [maxval,index] = max(pow)
% maxfreqband = [freq(index), freq(index+1)]


Nr = siglen;	        % Length of signal in samples
f = ((0:Nr-1)/Nr)*Fs;	% Vector of frequencies
overlap = nF/2;
window = hamming(nF);   % could be window = ones(nF,1); or window = triang(nF); or window = hamming(nF); etc.
[B1,fr1,tm1] = spectrogram(sig0,window,overlap,nF,Fs);
[B2,fr2,tm2] = spectrogram(sig1,window,overlap,nF,Fs);
[B3,fr3,tm3] = spectrogram(sig2,window,overlap,nF,Fs);

% B is the fft, rows are frequencies starting from 0, columns are windows
B1pow = abs(B1).^2;
B2pow = abs(B2).^2;
B3pow = abs(B3).^2;
% for each row
[rows, cols] = size(B1pow);
medianfreq1 = [];
maxfreq1 = [];
medianfreq2 = [];
maxfreq2 = [];
medianfreq3 = [];
maxfreq3 = [];
for i=1:cols
    [y, idx] = min(abs(B1pow(:,i) - median(B1pow(:,i))));
    [y, idm] = max(B1pow(:,i));
    medianfreq1 = [medianfreq1 fr1(idx)];
    maxfreq1 = [maxfreq1 fr1(idm)];
    [y, idx] = min(abs(B2pow(:,i) - median(B2pow(:,i))));
    [y, idm] = max(B2pow(:,i));
    medianfreq2 = [medianfreq2 fr2(idx)];
    maxfreq2 = [maxfreq2 fr2(idm)];
    [y, idx] = min(abs(B3pow(:,i) - median(B3pow(:,i))));
    [y, idm] = max(B3pow(:,i));
    medianfreq3 = [medianfreq3 fr3(idx)];
    maxfreq3 = [maxfreq3 fr3(idm)];
end
% PLOT estimators
% figure(1);
% hold on;
% xlabel('Time');
% ylabel('Median frequency');
% scatter(tm1, medianfreq1, 5, 'MarkerFaceColor',[0 .7 .7], 'LineWidth',0.01);
% scatter(tm2, medianfreq2, 5, 'MarkerFaceColor',[.7 0 .7], 'LineWidth',0.01);
% scatter(tm3, medianfreq3, 5, 'MarkerFaceColor',[.7 .7 0], 'LineWidth',0.01);
% hold off;
% 
% figure(2);
% hold on;
% xlabel('Time');
% ylabel('Peak frequency');
% scatter(tm1, maxfreq1, 5, 'MarkerFaceColor',[0 .7 .7], 'LineWidth',0.01);
% scatter(tm2, maxfreq2, 5, 'MarkerFaceColor',[.7 0 .7], 'LineWidth',0.01);
% scatter(tm3, maxfreq3, 5, 'MarkerFaceColor',[.7 .7 0], 'LineWidth',0.01);
% hold off;

figure(3);
t = tiledlayout(2,2);
title(t, 'sig1');
p = nexttile;
scatter(tm1, maxfreq1, 5, 'MarkerFaceColor',[0 .7 .7], 'LineWidth',0.01); title(p, 'Peak');
p = nexttile;
scatter(tm1, medianfreq1, 5, 'MarkerFaceColor',[0 .7 .7], 'LineWidth',0.01); title(p, 'Median');
nexttile([1 2]);
imagesc(tm1,fr1,20*log10(abs(B1)/nF));
xlabel(t, 'Time'); ylabel(t, 'Frequency');
saveas(gcf, record + '-' +  nF + '_' + preproccessingMode + '-sig1' +  '.png')

figure(4);
t = tiledlayout(2,2);
title(t, 'sig2');
p = nexttile;
scatter(tm1, maxfreq2, 5, 'MarkerFaceColor',[0 .7 .7], 'LineWidth',0.01); title(p, 'Peak');
p = nexttile;
scatter(tm1, medianfreq2, 5, 'MarkerFaceColor',[0 .7 .7], 'LineWidth',0.01); title(p, 'Median');
nexttile([1 2]);
imagesc(tm2,fr2,20*log10(abs(B2)/nF));
xlabel(t, 'Time'); ylabel(t, 'Frequency');
saveas(gcf, record + '-' +  nF + '_' + preproccessingMode + '-sig2' +  '.png')

figure(5);
t = tiledlayout(2,2);
title(t, 'sig3');
p = nexttile;
scatter(tm1, maxfreq3, 5, 'MarkerFaceColor',[0 .7 .7], 'LineWidth',0.01); title(p, 'Peak');
p = nexttile;
scatter(tm1, medianfreq3, 5, 'MarkerFaceColor',[0 .7 .7], 'LineWidth',0.01); title(p, 'Median');
nexttile([1 2]);
imagesc(tm3,fr3,20*log10(abs(B3)/nF));
xlabel(t, 'Time'); ylabel(t, 'Frequency');
saveas(gcf, record + '-' +  nF + '_' + preproccessingMode + '-sig3' +  '.png')

end
