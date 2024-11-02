function omegadot = Omegadot(omega,phi,phidot,u2)
%OMEGADOT
%    OMEGADOT = OMEGADOT(OMEGA,PHI,PHIDOT,U2)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    04-Nov-2019 13:55:40

t2 = sin(phi);
omegadot = ((u2.*2.0e6-omega.*phidot.*t2.*cos(phi).*6.2999e4).*2.46772582321671e14)./(t2.^2.*1.554642591368295e19+1.701040730077558e22);
