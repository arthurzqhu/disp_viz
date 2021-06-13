clear
opts = detectImportOptions('post_flight_details.csv');
tab = readtable('post_flight_details.csv', opts);

for iday = 1:size(tab,1)
    
    postpdi_flight_basics(iday).ti = strPDTto00Z(tab.StartTime(iday));
    postpdi_flight_basics(iday).tf = strPDTto00Z(tab.EndTime(iday));
    postpdi_flight_basics(iday).z_CB = tab.CldBase_m_(iday);
    postpdi_flight_basics(iday).z_CT = tab.CldTop_m_(iday);
    postpdi_flight_basics(iday).d_CT = tab.CTBumpiness(iday);
    
end

save('postpdi_flight_basics.mat','postpdi_flight_basics')

function time00Z = strPDTto00Z(strTime)

time00Z = ((datenum(strTime) - floor(datenum(strTime)))*24 + 8)*3600;

if time00Z >= 86400
    time00Z = time00Z - 86400;
end

end