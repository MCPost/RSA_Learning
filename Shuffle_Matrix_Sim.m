

Hyp_Mat_full = triu(ones(16)).*~eye(16);

Shuff_Mat = zeros(size(Hyp_Mat_full));
Ind_Mat = triu(reshape(1:size(Hyp_Mat_full,1)^2,size(Hyp_Mat_full,1),size(Hyp_Mat_full,2))).*~eye(16);
Rndp_Mat = Ind_Mat; Rndp_Mat(Rndp_Mat > 0) = shuffle(Rndp_Mat(Rndp_Mat > 0));
rand_idx = randperm(size(Hyp_Mat_full,1));
for row = 1:size(Hyp_Mat_full,1)-1
    for col = (row+1):size(Hyp_Mat_full,1)
        if(Ind_Mat(rand_idx(row),rand_idx(col)) ~= 0)
            Shuff_Mat(row,col) = Ind_Mat(rand_idx(row),rand_idx(col));
        else
            Shuff_Mat(row,col) = Ind_Mat(rand_idx(col),rand_idx(row));
        end
    end
end

figure
subplot(2,3,1)
imagesc(Ind_Mat); axis square
subplot(2,3,2)
imagesc(Shuff_Mat); axis square
subplot(2,3,3)
imagesc(Rndp_Mat); axis square


Hyp_Mat_mDiag = triu(ones(16)).*~eye(16).*~flip(eye(16));

Shuff_Mat = zeros(size(Hyp_Mat_mDiag));
Ind_Mat = triu(reshape(1:size(Hyp_Mat_mDiag,1)^2,size(Hyp_Mat_mDiag,1),size(Hyp_Mat_mDiag,2))).*~eye(16).*~flip(eye(16));
Rndp_Mat = Ind_Mat; Rndp_Mat(Rndp_Mat > 0) = shuffle(Rndp_Mat(Rndp_Mat > 0));
rand_idx = randperm(size(Hyp_Mat_mDiag,1));
for row = 1:size(Hyp_Mat_mDiag,1)-1
    for col = (row+1):size(Hyp_Mat_mDiag,1)
        if(Ind_Mat(rand_idx(row),rand_idx(col)) ~= 0)
            Shuff_Mat(row,col) = Ind_Mat(rand_idx(row),rand_idx(col));
        else
            Shuff_Mat(row,col) = Ind_Mat(rand_idx(col),rand_idx(row));
        end
    end
end

subplot(2,3,4)
imagesc(Ind_Mat); axis square
subplot(2,3,5)
imagesc(Shuff_Mat); axis square
subplot(2,3,6)
imagesc(Rndp_Mat); axis square



Dat1 = [ones(1,8) 2*ones(1,8)];
Dat2 = [ones(1,8) 2*ones(1,8)];

Dat_mat = zeros(16);
for i = 1:15
    for j = (i+1):16
        if((17-i) ~= j)
            if(Dat1(i) == 1 && Dat2(j) == 1)
                Dat_mat(i,j) = 1;
            elseif(Dat1(i) == 2 && Dat2(j) == 2)
                Dat_mat(i,j) = 2;
            else
                Dat_mat(i,j) = -1;
            end
        end
    end
end

figure
subplot(1,3,1)
imagesc(Dat_mat); axis square

rnd_idx = randperm(16);
Dat1_shuff = Dat1(rnd_idx);
Dat2_shuff = Dat2(rnd_idx);

Dat_mat_shuff = zeros(16);
for i = 1:15
    for j = (i+1):16
        if((17-i) ~= j)
            if(Dat1_shuff(i) == 1 && Dat2_shuff(j) == 1)
                Dat_mat_shuff(i,j) = 1;
            elseif(Dat1_shuff(i) == 2 && Dat2_shuff(j) == 2)
                Dat_mat_shuff(i,j) = 2;
            else
                Dat_mat_shuff(i,j) = -1;
            end
        end
    end
end

subplot(1,3,2)
imagesc(Dat_mat_shuff); axis square

Shuff_Mat = zeros(size(Dat_mat));
Ind_Mat = triu(reshape(1:size(Dat_mat,1)^2,size(Dat_mat,1),size(Dat_mat,2))).*~eye(16);
Rndp_Mat = Ind_Mat; Rndp_Mat(Rndp_Mat > 0) = shuffle(Rndp_Mat(Rndp_Mat > 0));
rand_idx = randperm(size(Dat_mat,1));
for row = 1:size(Dat_mat,1)-1
    for col = (row+1):size(Dat_mat,1)
        if(Ind_Mat(rand_idx(row),rand_idx(col)) ~= 0)
            Shuff_Mat(row,col) = Ind_Mat(rand_idx(row),rand_idx(col));
        else
            Shuff_Mat(row,col) = Ind_Mat(rand_idx(col),rand_idx(row));
        end
    end
end
Dat_mat_shuffmat = Dat_mat;
Dat_mat_shuffmat(Shuff_Mat > 0) = Dat_mat(Shuff_Mat(Shuff_Mat > 0));

subplot(1,3,3)
imagesc(Dat_mat_shuffmat); axis square


