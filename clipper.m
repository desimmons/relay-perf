% This function takes a freq domain input signal, performs time domain
% clipping, and converts back to the freq domain. The second/third output
% gives the numerically calculated Bussgang scaling parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% D.E.Simmons

function [ T , zeta , epsd] = clipper( R , Pmax )
  R_hat = sqrt(length(R))*ifft(R); % Signal at relay (time dom)
  P_thresh = (abs(R_hat) > sqrt(Pmax)); % Saturation check
  T_hat_temp = exp(1i*angle(R_hat .* P_thresh)).* P_thresh; % find argument of saturated points
  T_hat = R_hat .* ~P_thresh + T_hat_temp * sqrt(Pmax);% replaced saturated points with arg*sqrt{Pmax)
  T = 1/sqrt(length(R))*fft(T_hat); % trans signal at relay (freq dom)
  zeta = mean(abs(T_hat.*conj(R_hat))) / mean(abs(R_hat).^2); 
  epsd = mean(abs(T_hat).^2) - mean(abs(R_hat).^2) *zeta^2;
end

