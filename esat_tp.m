function es=esat_tp(Tc,p,varargin)
%This function calculates saturated vapor mixing ratio from temperature (K)  
%and pressure (mb).  For saturated vapor mixing ratio over ice, pass an
%additional third argument (it can be literally anything).
%Created by Adele Igel - updated 6/2013


c0=0.6105851e3;
c1=0.4440316e2;
c2=0.1430341e1;
c3=0.2641412e-1;
c4=0.2995057e-3;
c5=0.2031998e-5;
c6=0.6936113e-8;
c7=0.2564861e-11;
c8=-0.3704404e-13;

x=Tc;
es=c0+x.*(c1+x.*(c2+x.*(c3+x.*(c4+x.*(c5+x.*(c6+x.*(c7+x.*c8)))))));
% x=T-273.15;
% es = 2.53e11*exp(-5.42e3./(x+273.15));
ws=0.622*es./(p*100-es);