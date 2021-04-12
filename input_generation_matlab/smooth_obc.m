function val = smooth_obc(val_in,num);

gaussFilter = fct_GaussianFilter([num num], 1, 0);
[val,im_conv,count,NaNcount] = fct_convNaN(val_in, gaussFilter, 'same', .5);
val(isnan(val_in))=nan;
val = val_in;

return
