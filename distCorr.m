%% Distance Correlation Function

function C = distCorr(x,y)

% Distance Matrix
a = pdist2(x,x);
b = pdist2(y,y);

% Double Centering
A = bsxfun(@minus, bsxfun(@minus, a, nanmean(a,1)), nanmean(a,2)) + nanmean(nanmean(a,1),2);
B = bsxfun(@minus, bsxfun(@minus, b, nanmean(b,1)), nanmean(b,2)) + nanmean(nanmean(b,1),2);

% Sample Covariance
C = sqrt(nanmean(nanmean(A.*B)) / sqrt(nanmean(nanmean(A.*A)) * nanmean(nanmean(B.*B))));

end

