function q = binread (fname,prec,varargin)
% q = binread (fname,varargin)
%
% read a binary file (with big-endian architecture and double precision
% data types) into a matlab array with specified dimensions
% 
% fname: filename or path (string)
% varargin: a comma-separated list of dimension sizes. 
%           note (dimsize1) x (dimsize 2) x ... (dimsize N) must be equal
%           to the size of the file (in bytes) divided by 8!

if (size(varargin,2) < 1)
    
    error ('enter array dimensions');
    
else
    
    arrsize = [];
    
    for ind=1:size(varargin,2);
      arg=varargin{ind};
      arrsize = [arrsize arg];
    end

    
end

fid = fopen (fname, 'r', 'b');
if (prec==8)
    q=fread(fid,inf,'real*8');
elseif (prec==4)
    q=fread(fid,inf,'real*4');
else
    error('give precision of data');
end

if (size(varargin,2)>1)
 q=reshape(q,arrsize);
end

fclose(fid);

return

    
