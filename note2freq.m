function f = note2freq(scientific_name)

n_list = 1:88;
frequency_list_raw = 2^(1/12).^(n_list-49)*440;
frequency_list = [frequency_list_raw' frequency_list_raw'];

eq_temp_sharp = {'A','A#','B','C','C#','D','D#','E','F','F#','G','G#'};
eq_temp_flat = {'A','Bb','B','C','Db','D','Eb','E','F','Gb','G','Ab'};

for ikey = 1:88
    suff = num2str(floor((ikey+8)/12));
    
    inote_name = mod(ikey, 12);
    
    if inote_name==0
        inote_name=12;
    end
    
    pref_sharp = eq_temp_sharp{inote_name};
    pref_flat = eq_temp_flat{inote_name};
    
    key_list{ikey,1} = [pref_sharp suff];
    key_list{ikey,2} = [pref_flat suff];
end

if ischar(scientific_name) || isstring(scientific_name)
    scientific_name = {scientific_name};
elseif iscell(scientific_name)
    for inote = 1:length(scientific_name)
        ikey = find(ismember(key_list,scientific_name(inote)),1,'first');
        f(inote) = frequency_list(ikey);
    end
else
    error('please input its scientific name such as "A4" as in A4 = 440 Hz.')
end

end