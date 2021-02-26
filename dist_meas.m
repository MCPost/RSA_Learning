%% Distance Measure Function

function d = dist_meas(x, y, measure)

na_x = isnan(x);
na_y = isnan(y);

if(sum(na_x) == length(x) || sum(na_y) == length(y))
    d = NaN;
else
    
    x(na_x&na_y) = [];
    y(na_y&na_x) = [];
    
    if(strcmp(measure, 'pearson'))
        d = fast_corr(x, y);
    elseif(strcmp(measure, 'spearman'))
        d = corr(x, y, 'type', 'spearman');
    elseif(strcmp(measure, 'euclidian'))
        d = sqrt((x - y)'*(x - y));
    elseif(strcmp(measure, 'mahalanobis'))
        d = (x - nanmean(y))'/(cov(y))*((x - nanmean(y)));
    elseif(strcmp(measure, 'cosine'))
        d = dot(x, y)/(norm(x)*norm(y));
    end
    
end

end