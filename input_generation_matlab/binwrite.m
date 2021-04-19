function q = binwrite (fname,V)
% q = binwrite (fname,V)
%
% write a matlab array of arbitrary dimension to binary file
% uses big-endian architecture and double precision (8-byte) size
%
% fname: filename or path (string)
% V: array of arbitrary dimension (storage is independent of dimension
% sizes)


fid = fopen (fname, 'w', 'b');
q=fwrite(fid,V,'real*8');
fclose(fid);

return