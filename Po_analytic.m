% This function calculates the analytic outage probability for an FG or VG 
% relay network
%
% (Node A) <----> (Node R) <----> (Node B)
%
% at Node B for a two-way system when all nodes possess a soft envelope 
% limiter. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function takes vectorised first input x, where
%
% x(1) is Transmit power at node A
% x(2) is Transmit power at node B
% x(3) is Transmit power at node B
% x(4) and x(5) are average channel gain for hops A and B
% x(6) is Relays maximum transmit power
% x(7) is Outage probability threshold
% x(8) is Noise power spectral density
%
% and second string input type = 'FG' or 'VG', specifying fixed-gain or
% variable-gain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% D.E.Simmons
function [ PoB ] = Po_analytic(x, type)

  s2A=x(1); % Transmit power at node A
  s2B=x(2); % Transmit power at node B
  s2R=x(3); % Transmit power at node B
  muA=x(4);muB=x(5); % Average channel gain for hops A and B
  PmaxA=x(6); % Relays maximum transmit power
  PmaxB=x(7); % Relays maximum transmit power
  Pmax=x(8); % Relays maximum transmit power
  SNRth=x(9); % Outage probability threshold
  N0=x(10); % Noise power spectral density

  EA=s2A*(1-exp(-PmaxA/s2A));
  if s2B > 0
    EB=s2B*(1-exp(-PmaxB/s2B));
    zetaB = 1 - exp(-PmaxB / s2B) + sqrt(pi*PmaxB)/(2*sqrt(s2B)) * erfc(sqrt(PmaxB/s2B));
    epsdB = EB - s2B*zetaB.^2;
  else
    zetaB=1;epsdB=0;
    EB = 0;
  end
  ER = s2R*(1-exp(-Pmax/s2R));
  zeta = 1 - exp(-Pmax / s2R) + sqrt(pi*Pmax)/(2*sqrt(s2R)) * erfc(sqrt(Pmax/s2R));
  epsd = ER - s2R*zeta.^2;
  if epsd> 0
  else
    epsd=0;
  end
  zetaA = 1 - exp(-PmaxA / s2A) + sqrt(pi*PmaxA)/(2*sqrt(s2A)) * erfc(sqrt(PmaxA/s2A));
  epsdA = EA - s2A*zetaA.^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  w = zetaB^2*s2B*(1-zeta)^2+epsdB*zeta^2;
  SNRth = EA*SNRth / (s2A*zetaA^2 - SNRth*epsdA);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Analytic calculation if relay applies a fixed-gain (FG)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  if strcmp(type,'FG')
    G2 = (s2R/(EA*muA + EB*muB +N0));
    SIG1B  = EA * muA/N0;
    SIG2B  = ER  * muB/N0;  

    u1 = 1 + epsd /((zeta^2 *G2)*N0);
    u2 = w/(zeta^2*ER);
    C =  ER / (N0*zeta^2*G2);

    if SNRth > 0     
      PoB = 1-2*sqrt((SNRth*C)/(SIG2B *(u2*SIG2B *SNRth+SIG1B )))...
      *exp(-u1 *SNRth/SIG1B )*besselk(1,2*sqrt((C)*SNRth*(SIG1B ...
      +u2*SIG2B *SNRth)/SIG2B )/SIG1B );
    else
      PoB = 1;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Analytic calculation if relay applies a variable-gain (VG)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  elseif strcmp(type,'VG')
    SIG_AR = EA*muA/N0;
    if s2B ~= 0 && PmaxB ~= 0
      SIG_BR = EB*muB/N0;

      a = w*s2R/(N0*EB) + epsd/N0; 
      b = (ER+EB)/N0; 
      c = EB/N0;
      q1 = zeta^2*s2R/N0 - SNRth*epsd/N0; 
      q2 = SNRth*a/(q1^2*SIG_AR)+1/(q1*SIG_BR); 
      q3 = 1/(SIG_AR)*(c*SNRth+b*c*SNRth^2/q1+a*c^2*SNRth^3/q1^2);
      q4 = c*SNRth/(q1*SIG_BR) + 1/SIG_AR*(b*SNRth/q1 + 2*a*c*SNRth^2/q1^2);

      if q1 > 0
        PoB = 1-2*q3/(q1*SIG_BR*sqrt(q3*q2))*exp(-q4)*besselk(1,2*sqrt(q3*q2));
      else
        PoB = 1;
      end
    else
      q1 = zeta^2*s2R/N0 - SNRth*epsd/N0; 
      q2q3 = SNRth/(SIG_AR*q1*muB)*( 1  + ER*SNRth / (q1*N0) );
      q4 = SNRth*(1/(q1*muB)  +  ER/(N0*SIG_AR*q1));

      if q1 > 0
        PoB = 1-2*sqrt(q2q3)*exp(-q4)*besselk(1,2*sqrt(q2q3));
      else
        PoB = 1;
      end

    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Catch error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  else
    error('Incorrect relaying type input. Please choose FG or VG.')
  end

end