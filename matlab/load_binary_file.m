function out = load_binary_file(filename,nchannels,samples,dataformat,memorymap,writable)
%LOAD_BINARY_FILE Load a binary file into memory or map it using memmapfile 
%
% out = LOAD_BINARY_FILE(filename,nchannels,samples,dataformat,memorymap,writable)
%
% This function is supposed to be used for dat files commonly used in
% Klusters suite but can be used for other things. The number of samples is
% computed from the filesize divided by the product of whatever values are
% in nchannels (typically the number of channels).
%
% Important note: When using 'memmap' always make sure to convert to double
% before making computations with the data.
% Why using memmap??
% When you use memmap the data is not loaded to memory so if you need only
% for instance one channel, it might be worthed to use it like this.
%
% Note: Use the memmap option for very large files)
% If in memmap mode, the data is a memmap file object; you can access the
% data in out.Data.data
%
% This function takes 2.4min to load 20 min recording in mcnanalysis01 when 
%used with memorymap = 0. 
% With memorymap = 1, it returns in 0.1s but it
%returns a memmap object.
%
% See also MEMMAPFILE FOPEN FREAD
%
% Joao Couto 1 April 2015

if ~exist('filename','var') || ~exist('nchannels','var')
    error('Required inputs: load_binary_file(filename,nchannels)')
end

if ~exist('dataformat','var')
    dataformat =[];
end
if isempty(dataformat)
    dataformat = 'int16';
end
if ~exist('memorymap','var')
    memorymap = false;
end

if ~exist(filename,'file')
    error('File does not exist [%s]',filename)
end
if ~exist('samples','var')
    samples = [];
end
if ~exist('writable','var')
    writable = false;
end

file = dir(filename);

tmp = eval([dataformat,'(1)']);
s=whos('tmp');
tmps = s.bytes;
nsamples = file.bytes/(tmps*prod(nchannels));


if memorymap
    if ~isempty(samples)
        warning('Dont specify samples when using memmap mode.')
    end
    dataformat =  {dataformat,floor([nchannels,nsamples]),'data'};

    out = memmapfile(filename,'Format',dataformat,'Writable',writable);
    %disp('Data memory mapped, to access index out.Data.data(channel,samples)')
else
    if isempty(samples)
        samples = [0,nsamples];
    end
    fid = fopen(filename,'r');
    fseek(fid,((samples(1)-1)*nchannels),'bof');
    out = fread(fid,(diff(samples)*nchannels),[dataformat,'=>double']);
    out = reshape(out,[nchannels,diff(samples)]);
    fclose(fid);

end
