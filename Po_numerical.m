% This function calculates the numerical outage probability for an FG or VG 
% relay network:
%
% (Node A) <----> (Node R) <----> (Node B)
%
% at Node B for a two-way system when the all nodes possess a soft envelope 
% limiter. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function takes vectorised first input x, where
%
% x(1) is Transmit power at node A
% x(2) is Transmit power at node B
% x(3) is Transmit power at node B
% x(4) and x(5)are verage channel gain for hops A and B
% x(6) is Relays maximum transmit power
% x(7) is Outage probability threshold
% x(8) is Noise power spectral density
% x(9) is Number of OFDM sybcarriers
% x(10) is number of monte carlo simulations
%
% and second string input 'type' = 'FG' or 'VG', specifying fixed-gain or
% variable-gain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% D.E.Simmons

function [ Po ] = Po_numerical( x , type)


    s2A=x(1); % Transmit power at node A
    s2B=x(2); % Transmit power at node B
    s2=x(3); % power into relay's clipper
    muA=x(4);muB=x(5); % Average channel gain for hops A and B
    PmaxA=x(6); % node A maximum transmit power
    PmaxB=x(7); % node B maximum transmit power
    Pmax=x(8); % Relay's maximum transmit power
    SNRth=x(9); % Outage probability threshold
    N0=x(10); % Noise power spectral density
    N = x(11); % Number of OFDM sybcarriers
    trial = x(12) ; % number of monte carlo simulations
    l = x(13); % number of channel taps

    Po = 0; % declare output
    for monte = 1:trial
    %%%%%%%%%%%%%%%%%%%%% node A and B Tx symbols %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      XA = exp(1i* 2 * pi * unidrnd(8,1,N)/8); % generate PSK symbols Tx at A
      XB = exp(1i* 2 * pi * unidrnd(8,1,N)/8); % """"""""""""""""""""""""""""
    %%%%%%%%%%%%%%%%%%%%% node A Clipping %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      TxA = sqrt(s2A).*XA;  % Signal at node A before amplification (freq dom)
      [TxA , zetaA , epsdA] = clipper(TxA , PmaxA); % Signal at node A after ...
      EA=s2A*(1-exp(-PmaxA/s2A)); % Tx power at node A         
    %%%%%%%%%%%%%%%%%%%%% node B Clipping %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
      if s2B > 0 % Stop division by zero if B is not transmitting
        TxB = sqrt(s2B).*XB; % Signal at node B before amplification (freq dom)
        [TxB , zetaB , epsdB] = clipper(TxB , PmaxB); % Signal at node B after ... 
        EB=s2B*(1-exp(-PmaxB/s2B)); % Tx power at node B    
      else
        zetaB=1;epsdB=0;
        EB = 0;TxB=0;
      end
    %%%%%%%%%%%%%% Calculate channel coefficients and noise %%%%%%%%%%%%%%%%%%%
      hATD = normrnd(0,1/sqrt(2),1,l) + 1i.*normrnd(0,1/sqrt(2),1,l); % hop A channel 
      hA = sqrt(muA/l)*fft(hATD,N); % freq domain channel coefficients

      hBTD = normrnd(0,1/sqrt(2),1,l) + 1i.*normrnd(0,1/sqrt(2),1,l); % hop B channel 
      hB = sqrt(muB/l)*fft(hBTD,N); % freq domain channel coefficients

      n1 = sqrt( N0/2)*(normrnd(0,1,1,N) + 1i.*normrnd(0,1,1,N));  %  "  1 noise 
    %%%%%%%%%%%%%% Calculate relay's gain %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if strcmp(type,'FG')
        G = sqrt(s2 ./ (EA * muA + EB*muB + N0)); % Fixed-gain (FG)
      elseif strcmp(type,'VG')
        G = sqrt(s2 ./ (EA*abs(hA).^2 + EB*abs(hB).^2 + N0)); % Variable-gain (VG)
      else
        error('Incorrect relaying type input. Please choose "FG" or "VG".')
      end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Relay Clipping %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      R = G .* (hA .*  TxA + hB .*  TxB + n1);  % Signal at relay (freq dom)
      [T , zeta , epsd] = clipper(R , Pmax);  
    %%%%%%%%%%%%%%%%%%%%% Calculate SNR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      Noise =  abs(hB).^2 .* (G.^2 .*zeta^2 .*(N0+abs(hA).^2 .*epsdA )+epsd)+...
               G.^2 .*abs(hB).^4 .*(s2B*zetaB^2*abs(zeta-1)^2+epsdB*zeta^2) ...
               +N0; % calc noise at B
      Sig = abs(G .* zeta .* hB .* hA ).^2 * s2A*zetaA^2; % calc signal at B
    %%%%%%%%%%%%%%%%%%%%% Calculate Outage Probability %%%%%%%%%%%%%%%%%%%%%%%%
      if SNRth > 0
        Po = sum( Sig./Noise < SNRth) + Po; % Count number of SNR<SNRth events
      else
        Po = N+Po;
      end

    end  
    Po = Po / (trial*N); % convert to probability     
end
