function [data, timestamps, info] = load_open_ephys_multi_data(fileNames,block)

%
% [data, timestamps, info] = load_open_ephys_data(filename)
%
%   Loads continuous, event, or spike data files into Matlab.
%
%   Inputs:
%
%     filename: path to file
%
%
%   Outputs:
%
%     data: either an array continuous samples (in microvolts),
%           a matrix of spike waveforms (in microvolts),
%           or an array of event channels (integers)
%
%     timestamps: in seconds
%
%     info: structure with header and other information
%
%
%
%   DISCLAIMER:
%
%   Both the Open Ephys data format and this m-file are works in progress.
%   There's no guarantee that they will preserve the integrity of your
%   data. They will both be updated rather frequently, so try to use the
%   most recent version of this file, if possible.
%
%

%
%     ------------------------------------------------------------------
%
%     Copyright (C) 2014 Open Ephys
%
%     ------------------------------------------------------------------
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     <http://www.gnu.org/licenses/>.
%

if ~iscell(fileNames)
    fileNames={fileNames};
end

    % constants
    NUM_HEADER_BYTES = 1024;
    SAMPLES_PER_RECORD = 1024;
    RECORD_SIZE = 8 + 16 + SAMPLES_PER_RECORD*2 + 10; % size of each continuous record in bytes
    RECORD_MARKER = [0 1 2 3 4 5 6 7 8 255]';
    RECORD_MARKER_V0 = [0 0 0 0 0 0 0 0 0 255]';
    
    % constants for pre-allocating matrices:
    MAX_NUMBER_OF_SPIKES = 1e6;
    MAX_NUMBER_OF_RECORDS = 1e6;
    MAX_NUMBER_OF_CONTINUOUS_SAMPLES = 1e8;
    MAX_NUMBER_OF_EVENTS = 1e6;
    SPIKE_PREALLOC_INTERVAL = 1e6;
    
    for fileNum=1:size(fileNames,2)
        filetype = fileNames{fileNum}(max(strfind(fileNames{fileNum},'.'))+1:end); % parse filetype
        if ~strcmp(filetype, 'continuous')
            continue
        end
        
        fid(fileNum) = fopen(fileNames{fileNum});
        filesize(fileNum) = getfilesize(fid);
        
        switch nargin
            case 1
                block(fileNum)=filesize;
        end
    end
    
%% get continuous data 
   
    index = 0;
    
    hdr = fread(fid(1), NUM_HEADER_BYTES, 'char*1');
    eval(char(hdr'));
    info.header = header;
    
    if (isfield(info.header, 'version'))
        version = info.header.version;
    else
        version = 0.0;
    end
    
    % pre-allocate space for continuous data
    data{fileNum} = zeros(min([block MAX_NUMBER_OF_CONTINUOUS_SAMPLES]), 1,'int16');
%     info.ts = zeros(1, MAX_NUMBER_OF_RECORDS);
%     info.nsamples = zeros(1, MAX_NUMBER_OF_RECORDS);
%     
%     if version >= 0.2
%         info.recNum = zeros(1, MAX_NUMBER_OF_RECORDS);
%     end
    
    current_sample = 0;
    
    while ftell(fid(1)) + RECORD_SIZE < filesize % at least one record remains
        
        go_back_to_start_of_loop = 0;
        
        index = index + 1;
        
        if version >= 0.1
            timestamp = fread(fid, 1, 'int64', 0, 'l');
            nsamples = fread(fid, 1, 'uint16',0,'l');
            
            
            if version >= 0.2
                recNum = fread(fid, 1, 'uint16');
            end
            
        else
            timestamp = fread(fid, 1, 'uint64', 0, 'l');
            nsamples = fread(fid, 1, 'int16',0,'l');
        end
        
        
        if nsamples ~= SAMPLES_PER_RECORD && version >= 0.1
            
            disp(['  Found corrupted record...searching for record marker.']);
            
            % switch to searching for record markers
            
            last_ten_bytes = zeros(size(RECORD_MARKER));
            
            for bytenum = 1:RECORD_SIZE*5
                
                byte = fread(fid, 1, 'uint8');
                
                last_ten_bytes = circshift(last_ten_bytes,-1);
                
                last_ten_bytes(10) = double(byte);
                
                if last_ten_bytes(10) == RECORD_MARKER(end);
                    
                    sq_err = sum((last_ten_bytes - RECORD_MARKER).^2);
                    
                    if (sq_err == 0)
                        disp(['   Found a record marker after ' int2str(bytenum) ' bytes!']);
                        go_back_to_start_of_loop = 1;
                        break; % from 'for' loop
                    end
                end
            end
            
            % if we made it through the approximate length of 5 records without
            % finding a marker, abandon ship.
            if bytenum == RECORD_SIZE*5
                
                disp(['Loading failed at block number ' int2str(index) '. Found ' ...
                    int2str(nsamples) ' samples.'])
                
                break; % from 'while' loop
                
            end
            
            
        end
        
        if ~go_back_to_start_of_loop
            
            block = fread(fid, nsamples, 'int16', 0, 'b'); % read in data
            
            fread(fid, 10, 'char*1'); % read in record marker and discard
            
            data{fileNum}(current_sample+1:current_sample+nsamples) = block;
            
            current_sample = current_sample + nsamples;
            
            info.ts(index) = timestamp;
            info.nsamples(index) = nsamples;
            
            if version >= 0.2
                info.recNum(index) = recNum;
            end
            
        end
        
    end
    
    % crop data to the correct size
    data{fileNum} = data{fileNum}(1:current_sample);
    info.ts = info.ts(1:index);
    info.nsamples = info.nsamples(1:index);
    
    if version >= 0.2
        info.recNum = info.recNum(1:index);
    end
    
    % convert to microvolts
    data{fileNum} = data{fileNum}.*info.header.bitVolts;
    
    timestamps = nan(size(data{fileNum}));
    
    current_sample = 0;
    
    if version >= 0.1
        
        for record = 1:length(info.ts)
            
            ts_interp = info.ts(record):info.ts(record)+info.nsamples(record);
            
            timestamps(current_sample+1:current_sample+info.nsamples(record)) = ts_interp(1:end-1);
            
            current_sample = current_sample + info.nsamples(record);
        end
    else % v0.0; NOTE: the timestamps for the last record will not be interpolated
        
        for record = 1:length(info.ts)-1
            
            ts_interp = linspace(info.ts(record), info.ts(record+1), info.nsamples(record)+1);
            
            timestamps(current_sample+1:current_sample+info.nsamples(record)) = ts_interp(1:end-1);
            
            current_sample = current_sample + info.nsamples(record);
        end
        
    end
    
    fclose(fid); % close the file
    
%     if (isfield(info.header,'sampleRate'))
%         if ~ischar(info.header.sampleRate)
%             timestamps = timestamps./info.header.sampleRate; % convert to seconds
%         end
%     end

clearvars -except data timestamps info
data=vertcat(data{:});

end


function filesize = getfilesize(fid)

fseek(fid,0,'eof');
filesize = ftell(fid);
fseek(fid,0,'bof');

end
