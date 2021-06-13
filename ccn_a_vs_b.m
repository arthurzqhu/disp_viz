c = 1;
campaigns = {'vocalspdi','postpdi','oraclespdi','gomaccspdi'};
camp=campaigns{c};

close all
for iday = 1:length(clouds.(camp))
    figure
    hold on
    plot(clouds.(camp)(iday).ccn_a)
    plot(clouds.(camp)(iday).ccn_b)
end