function finishingTaskSound

smb3_extended_1up_notes = {'D6' 'E6' 'D6' 'Db6' 'C6' 'Bb5' 'C6' 'E6' 'G6' 'E7' 'C7' 'D7' 'G7'};
zelda_puzzle_solved_notes = {'G6' 'Gb6' 'E6' 'A5' 'Ab5' 'E6' 'G#6' 'C7'};

notes = smb3_extended_1up_notes;

dt = 0.01; % time res
freq = note2freq(notes); % converts notes to frequency
Sr = 44100; % sample rate
omega = findOmega(freq, Sr, dt); 
duration = repelem(0.13,length(notes)); % duration for each notes

wave_t = [];

for inote = 1:length(notes)
    wave_t = [wave_t 0.1*square(sin((1:dt:Sr*dt*duration)*omega(inote)))];
end

sound(wave_t,Sr)

end