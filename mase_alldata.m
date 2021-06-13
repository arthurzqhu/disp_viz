% VOCALS
%% Import data from text file.
filedir='MASE_raw/springston-summary/';
files=dir([filedir,'*.txt']);
ncol = 63;

for ff=1:length(files)
    %% Initialize variables.
    filename = [filedir,files(ff).name];
    delimiter = '\t';
    startRow = 40;
    
    %% Read columns of data as text:
    formatSpec = [repmat('%s',1,ncol),'%[^\n\r]'];
    
    %% Open the text file.
    fileID = fopen(filename,'r');
    
    %% Read columns of data according to the format.
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, ...
        'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    
    %% Close the text file.
    fclose(fileID);
    
    %% Convert the contents of columns containing numeric text to numbers.
    % Replace non-numeric text with NaN.
    raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
    for col=1:length(dataArray)-1
        raw(1:length(dataArray{col}),col) = dataArray{col};
    end
    numericData = NaN(size(dataArray{1},1),size(dataArray,2));
    
    for col=3:ncol
        % Converts text in the input cell array to numbers. Replaced non-numeric
        % text with NaN.
        rawData = dataArray{col};
        for row=1:size(rawData, 1)
            % Create a regular expression to detect and remove non-numeric prefixes and
            % suffixes.
            regexstr = ['(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}'...
                '\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}'...
                '\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)'];
            try
                result = regexp(rawData{row}, regexstr, 'names');
                numbers = result.numbers;
                
                % Detected commas in non-thousand locations.
                invalidThousandsSeparator = false;
                if any(numbers==',')
                    thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                    if isempty(regexp(numbers, thousandsRegExp, 'once'))
                        numbers = NaN;
                        invalidThousandsSeparator = true;
                    end
                end
                % Convert numeric text to numbers.
                if ~invalidThousandsSeparator
                    numbers = textscan(strrep(numbers, ',', ''), '%f');
                    numericData(row, col) = numbers{1};
                    raw{row, col} = numbers{1};
                end
            catch me
            end
        end
    end
    
    dateFormats = {'yyyy-MM-dd', 'HH:mm:ss'};
    dateFormatIndex = 1;
    blankDates = cell(1,size(raw,2));
    anyBlankDates = false(size(raw,1),1);
    invalidDates = cell(1,size(raw,2));
    anyInvalidDates = false(size(raw,1),1);
    for col=[1,2]% Convert the contents of columns with dates to MATLAB datetimes using the specified date format.
        try
            dates{col} = datetime(dataArray{col}, 'Format', ...
                dateFormats{col==[1,2]}, 'InputFormat', ...
                dateFormats{col==[1,2]}); %#ok<SAGROW>
        catch
            try
                % Handle dates surrounded by quotes
                dataArray{col} = cellfun(@(x) x(2:end-1), dataArray{col}, ...
                    'UniformOutput', false);
                dates{col} = datetime(dataArray{col}, 'Format', ...
                    dateFormats{col==[1,2]}, 'InputFormat', ...
                    dateFormats{col==[1,2]}); %%#ok<SAGROW>
            catch
                dates{col} = repmat(datetime([NaN NaN NaN]), size(dataArray{col})); %#ok<SAGROW>
            end
        end
        
        dateFormatIndex = dateFormatIndex + 1;
        blankDates{col} = cellfun(@isempty, dataArray{col});
        anyBlankDates = blankDates{col} | anyBlankDates;
        invalidDates{col} = isnan(dates{col}.Hour) - blankDates{col};
        anyInvalidDates = invalidDates{col} | anyInvalidDates;
    end
    dates = dates(:,[1,2]);
    blankDates = blankDates(:,[1,2]);
    invalidDates = invalidDates(:,[1,2]);
    
    %% Split data into numeric and cell columns.
    rawNumericColumns = raw(:,3:ncol);
    
    %% Replace non-numeric cells with NaN
    R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),rawNumericColumns); % Find non-numeric cells
    rawNumericColumns(R) = {NaN}; % Replace non-numeric cells
    
    %% Allocate imported array to column variable names
    %gen(ff).Date = dates{:, 1};
    %gen(ff).Time = dates{:, 2};
    gen(ff).Lat = cell2mat(rawNumericColumns(:, 1));
    gen(ff).Lon = cell2mat(rawNumericColumns(:, 2));
    gen(ff).Radar_Alt = cell2mat(rawNumericColumns(:, 3));
    gen(ff).Static_Press = cell2mat(rawNumericColumns(:, 4));
    %Paltitude1 = cell2mat(rawNumericColumns(:, 5));
    gen(ff).Ambient_Temp = cell2mat(rawNumericColumns(:, 6));
    %gen(ff).Theta1 = cell2mat(rawNumericColumns(:, 7));
    gen(ff).Dew_Point = cell2mat(rawNumericColumns(:, 8));
    %TIR
    %TRF
    gen(ff).RH_water = cell2mat(rawNumericColumns(:, 11));
    %H2OMR1 = cell2mat(rawNumericColumns(:, 12));
    gen(ff).Wind_speed = cell2mat(rawNumericColumns(:, 13));
    gen(ff).Wind_Dir = cell2mat(rawNumericColumns(:, 14));
    %UVsky1 = cell2mat(rawNumericColumns(:, 15));
    %UVgnd1 = cell2mat(rawNumericColumns(:, 16));
    gen(ff).CCNSSA = cell2mat(rawNumericColumns(:, 17));
    gen(ff).CCNNA = cell2mat(rawNumericColumns(:, 18));
    gen(ff).CCNSSC = cell2mat(rawNumericColumns(:, 19));
    gen(ff).CCNNC = cell2mat(rawNumericColumns(:, 20));
    gen(ff).CCNSSB = cell2mat(rawNumericColumns(:, 21));
    gen(ff).CCNNB = cell2mat(rawNumericColumns(:, 22));    
    %CPC3025 = cell2mat(rawNumericColumns(:, 23));
    %CPC3010 = cell2mat(rawNumericColumns(:, 24));
    %DMAtotN1 = cell2mat(rawNumericColumns(:, 25));
    %DMAtotA1 = cell2mat(rawNumericColumns(:, 26));
    %DMAtotV1 = cell2mat(rawNumericColumns(:, 27));
    %PCASPtotN1 = cell2mat(rawNumericColumns(:, 28));
    %PCASPtotA1 = cell2mat(rawNumericColumns(:, 29));
    %PCASPtotV1 = cell2mat(rawNumericColumns(:, 30));
    %gen(ff).CAStotN1 = cell2mat(rawNumericColumns(:, 31));
    %gen(ff).CAStotA1 = cell2mat(rawNumericColumns(:, 32));
    %gen(ff).CAStotV1 = cell2mat(rawNumericColumns(:, 33));
    
    gen(ff).CIP_RWC = cell2mat(rawNumericColumns(:, 34))*1e6;
    %CIPtotA1 = cell2mat(rawNumericColumns(:, 35));
    gen(ff).CIP_RWC = cell2mat(rawNumericColumns(:, 36))*1e6;
    LWCHotWire1 = cell2mat(rawNumericColumns(:, 37));
    LWCGerber1 = cell2mat(rawNumericColumns(:, 38)); 
    LWC=max([LWCHotWire1,LWCGerber1]')';
    gen(ff).Hgt_Abv_CB = gen(ff).Radar_Alt-max(0,min(gen(ff).Radar_Alt(LWC>=0.01)));
    %Bluebacks1 = cell2mat(rawNumericColumns(:, 39));
    %Greenbacks1 = cell2mat(rawNumericColumns(:, 40));
    %Redbacks1 = cell2mat(rawNumericColumns(:, 41));
    %Bluetots1 = cell2mat(rawNumericColumns(:, 42));
    %Greentots1 = cell2mat(rawNumericColumns(:, 43));
    %Redtots1 = cell2mat(rawNumericColumns(:, 44));
    %Blueabs1 = cell2mat(rawNumericColumns(:, 45));
    %Greenabs1 = cell2mat(rawNumericColumns(:, 46));
    %Redabs1 = cell2mat(rawNumericColumns(:, 47));   
    %O3 = cell2mat(rawNumericColumns(:, 48));
    %CO = cell2mat(rawNumericColumns(:, 49));
    %SO2 = cell2mat(rawNumericColumns(:, 50));
       
    % For code requiring serial dates (datenum) instead of datetime, uncomment
    % the following line(s) below to return the imported dates as datenum(s).
    
    % Date1=datenum(Date1);
    % Time1=datenum(Time1);
    
    
    %% Clear temporary variables
    clearvars filename delimiter startRow formatSpec fileID dataArray ans raw
    clearvars col numericData rawData row regexstr result numbers
    clearvars invalidThousandsSeparator thousandsRegExp me dateFormats
    clearvars dateFormatIndex dates blankDates anyBlankDates invalidDates
    clearvars anyInvalidDates rawNumericColumns R
end