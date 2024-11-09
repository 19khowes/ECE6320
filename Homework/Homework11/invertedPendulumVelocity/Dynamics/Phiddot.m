function phiddot = Phiddot(omega,phi,phidot,u1)
%PHIDDOT
%    PHIDDOT = PHIDDOT(OMEGA,PHI,PHIDOT,U1)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    04-Nov-2019 13:55:40

t2 = cos(phi);
t3 = sin(phi);
phiddot = ((t3.*-8.9237579508e10+u1.*2.110284e10+t2.*u1.*8.63e9+omega.^2.*t2.*t3.*2.71840685e8-phidot.^2.*t2.*t3.*6.38746861e8).*(5.0./7.3e1))./(t2.^2.*1.8619225e7-1.02941388e8);