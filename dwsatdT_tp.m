function dwsdT=dwsatdT_tp(Tc,p)

dT = 0.001;

ws = wsat_tp([Tc-dT Tc+dT],p);
dwsdT = (ws(2)-ws(1))/(2*dT);