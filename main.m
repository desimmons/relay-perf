clear all
% This script calls the functions Po_analytical and Po_numerical
% that calculate the analytical and numerical outage probabilities
% of a two-hop two-way OFDM based amplify-and-forward relay network
% implementing either fixed-gain (FG) or variable-gain (VG) at the
% relay
%
% D.E.Simmons

s2A = 500; % input power at node A's amplifier
s2B= 300; % input power at node B's amplifier
s2R = 500; % input power at relay's amplifier
muA = 0.1;   muB = 1; % average subcarrier response
pmaxA = 1000; % node A maximum Tx power
pmaxB = 600; % node A maximum Tx power
pmaxR = 1000; % relay maximum Tx power
SNRth = 1; % outage SNR threshold
N0 = 1; % noise power sepctral density

N = 128; % number of subcarriers
trial = 1000; % number of monte carlo trials
taps = 32; % number of channel taps in each hop

s2range = logspace(3,6,15); % range of input powers at relay's amplifier
s2Arange = logspace(0,7,15); % range of Tx power at node A
count = 0;

for s2A = s2Arange
  count = count + 1;

  x = [s2A s2B s2R muA muB pmaxA pmaxB pmaxR SNRth N0];

  Po_FG(count) = Po_analytic(x, 'FG'); 
  Po_VG(count) = Po_analytic(x, 'VG'); 
  Po_num_FG(count) = Po_numerical([x N trial taps], 'FG');
  Po_num_VG(count) = Po_numerical([x N trial taps], 'VG');
end

figure(3)
semilogy(10*log10(s2Arange),Po_FG,'r','Linewidth',3)
hold on
semilogy(10*log10(s2Arange),Po_VG,'k','Linewidth',3)
semilogy(10*log10(s2Arange),Po_num_FG,'ro','Linewidth',3)
hold on
semilogy(10*log10(s2Arange),Po_num_VG,'ko','Linewidth',3)