function rms = FindRMS (X)
% X is a passed in variable containg a emg channel
% Fs=1200; % make sure to check the sampling frequency
%- Compute MF for each column of X:
% plot(X)
%   NFFT=size(X,1);
%   % Compute power spectrum of fake signal:
%   [Pxx2,F2] = psd(X, NFFT, Fs);
%   Pxx=Pxx2(3:size(Pxx2,1),1);
%   F=F2(3:size(F2,1),1);
%   figure, plot(Pxx2)
%   % For Median
%   % Compute Normalized Cumulative Power Spectrum:
%   NCPxx = cumSum(Pxx) ./ max(cumSum(Pxx));
%   % Define indices of points just before and after NCPxx crosses 50% line:
%   I1 = max(find(NCPxx < 0.50));
%   I2 = I1 + 1;
%   %- Linearly interpolate to compute "Median Frequency":
%   %- linear interpolation between (x1,y1) and (x2,y2) to find x*, given y*
%   %-  x* = x1 + ((x2-x1)/(y2-y1)) * (y*-y1)
%   MedFr = F(I1,1) + ((F(I2,1)-F(I1,1)) ./ (NCPxx(I2,1)-NCPxx(I1,1)) .* (0.50-NCPxx(I1,1)));
%   
%   % For Mean
%   [Pxxa,iP]=sort(Pxx);
%   for k=1:size(iP,1)
%       Fa(k,1)=F(iP(k,1),1);
%   end
%   cumSumF(1,1)=(((Pxxa(2,1)-Pxxa(1,1))*(Fa(2,1)-Fa(1,1)))/2);
%   for j=2:size(Pxxa,1)-1
%       if (Fa(j+1,1)-Fa(j,1))<0
%           cumSumF(j,1)=cumSumF(j-1,1)+(Fa(j,1)*(Pxxa(j+1,1)-Pxxa(j,1)))-(((Pxxa(j+1,1)-Pxxa(j,1))*(Fa(j+1,1)-Fa(j,1)))/2);
%       else if (Fa(j+1,1)-Fa(j,1))>0
%               cumSumF(j,1)=cumSumF(j-1,1)+(Fa(j,1)*(Pxxa(j+1,1)-Pxxa(j,1)))+(((Pxxa(j+1,1)-Pxxa(j,1))*(Fa(j+1,1)-Fa(j,1)))/2);
%           else if (Fa(j+1,1)-Fa(j,1))==0
%                   cumSumF(j,1)=cumSumF(j-1,1)+(Fa(j,1)*(Pxxa(j+1,1)-Pxxa(j,1)));
%               end
%           end
%       end
%   end
%   cumSumP(1,1)=(Pxxa(size(Pxxa,1),1)-Pxxa(1,1));
%   
%   MeanFr=cumSumF(size(cumSumF,1),1)/cumSumP(1,1);

  % Demean
  X = detrend(X, 'constant');
  % 4-th order butterworth filter
  [b,a] = butter(4,[20/600 500/600],'bandpass');
  FilteredX = filtfilt(b,a,X);
  % For RMS
  rms=sqrt(mean(FilteredX.^2));
  
end
  
%   disp('IMP:First column in Results is Median Frequency, second is mean, thrid is RMS');
% Results=[MedFr MeanFr rms];
% dlmwrite('Results.CSV',Results);
