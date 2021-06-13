clear

wind_files = dir('ORACLES_raw/WINDS_P3/*.ict');

%%
for iday = 1:length(wind_files)
    wind{iday,1}(:,:) = readmatrix([wind_files(iday).folder '/' wind_files(iday).name],...
        'FileType','text','NumHeaderLines',49);

    idx_1Hz = find(floor(wind{iday,1}(:,1))==wind{iday,1}(:,1));
    ORACLES_flight_basics(iday).date_str = wind_files(iday).name(12:17);
    ORACLES_flight_basics(iday).t_1Hz = wind{iday,1}(idx_1Hz, 1);
    ORACLES_flight_basics(iday).z_1Hz = wind{iday,1}(idx_1Hz, 4);
    ORACLES_flight_basics(iday).T_1Hz = wind{iday,1}(idx_1Hz, 7);
    
end

save('ORACLES_flight_basics.mat','ORACLES_flight_basics')