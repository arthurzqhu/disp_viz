clear

filepath = 'mase/050716_Ps.mat';
gen = load(filepath);

fields=fieldnames(gen);

for jfield=1:length(fields)
    val = gen.(fields{jfield});
    gen.(fields{jfield})(val==-9999) = nan;
end

save(filepath, '-struct','gen','s_*');

