VARIABLE LEGEND
01. s_t (s) 		elapsed Time in seconds since 0 UTC of flight (data file) start day
02. s_ap (m)		Pressure altitude (adjusted to radar altitude) 
03. s_hr (m)		Radar altitude
04. s_lat (degN)	LATitude from UCI's C-MIGITS III
05. s_lon (degE)	LONgitude from UCI's C-MIGITS III
06. s_hdg (deg)		true HeaDinG from UCI's C-MIGITS III range [0 360] deg 
07. s_wx (m/s)		Wind component in the east direction (X-axis)
08. s_wy (m/s)		Wind component in the north direction (Y-axis)
09. s_wz (m/s)		Wind component in the vertical direction (Z-axis)	[No data during MASE]
10. s_ah (g/m^3)	Absolute Humidity 		[POST TO01-TO04: Chilled mirror; TO05-TO17 LI-COR 7500]
11. s_ta (deg C)	static Ambient Temperature from UCI's Rosemount fast-response sensor
12. s_td (deg C)	ambient Dewpoint Temperature from CIRPAS's Edgtech Chilled mirror sensor
13. s_ts (deg C)	Sea surface Temperature from CIRPAS's downlooking Heiman KT 19.85 IR sensor 
14. s_ps (hPa)		Static atmospheric Pressure from fuselage flush ports and Setra 270 transducer 		 	 	  
15. s_tas (m/s)		True Air Speed (Dry Air)
16. s_rhoa (kg/m^3)	Moist Air density 	[No data during MASE]
17. s_mr (g/kg)		Mixing Ratio from UCI's LI-COR 7500 	[POST TO01-TO04: Chilled mirror; TO05-TO17 LI-COR 7500]
18. s_thet (K)		potential temperature (theta)
19. s_tvir (deg C)	VIRtual Temperature	[No data during MASE]
20. s_thete (K)		Equivalent potential temperature (thetae)
21. s_tirup (deg C)	Temperature from UCI's IR UPward-looking temperature sensor	[No data during MASE]
22. s_tdl (deg C)	Dewpoint Temperature from UCI's LI-COR 7500 	[No data from MASE; POST flights TO01-TO04: NaN]
23. s_lwc_xg (g/m^3)	Liquid water content from Gerber PVM-100A probe

24. s_conc_pdi (cm^-3)	PDI binned number concentration in units of dN/dlogDp
25. s_ntot_pdi (cm^-3)	PDI total drop number concentration
26. s_lwc_pdi (g/m^3)	PDI liquid water content
27. s_reff_pdi (um)	PDI effective radius
28. s_std_pdi (um)	Standard deviation of 1 Hz PDI DSDs
29. s_disp_pdi (-)	Relative dispersion of 1 Hz PDI DSDs
30. s_R_pdi (mm/d)	PDI rain rate/sedimentation flux

31. s_conc_cip (cm^-3)	CIP binned number concentration in units of dN/dlogDp
32. s_ntot_cip (cm^-3)	CIP total number concentration
33. s_lwc_cip (g/m^3)	CIP liquid water content
34. s_R_cip (mm/d)	CIP rain rate/sedimentation flux
35. s_R_merge (mm/d)	Combined PDI+CIP rain rate/sedimentation flux (crossover diameter 63 um)

36. s_ahkc (g/m3)	VOCALS ONLY: corrected absolute humidity from Krypton KH2O hygrometer


If you want to construct merged DSDs, use PDI bins 1-96 and CIP bins 3-62. Plot merged DSDs with the diameter grid given in dp_merge.mat.


