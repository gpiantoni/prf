
function [BINS_CENTER, rfsize_mean, rfsize_sem] = compute_bins(ecc, rfsize)
bin_start = 0.5;
bin_size = 0.5;
n_bins = 29;

BINS_EDGE = colon(bin_start, bin_size, n_bins * bin_size + bin_start);
binned = discretize(ecc, BINS_EDGE);

rfsize_mean = nan(1, length(BINS_EDGE));
rfsize_sem = nan(1, length(BINS_EDGE));
for i = 1:length(BINS_EDGE)
    val = rfsize(binned == i);
    rfsize_mean(i) = mean(val);
    rfsize_sem(i) = std(val) / sqrt(length(val));
end

BINS_CENTER = BINS_EDGE + bin_size / 2;
% 
% BINS_CENTER = BINS_CENTER(~isnan(rfsize_mean));
% rfsize_mean = rfsize_mean(~isnan(rfsize_mean));

end