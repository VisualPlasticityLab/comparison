%% Preparation
% Load cell data
[plane1, plane2, plane3] = db2()


for i=1:5
%% Concatenate plane data
 se = concatenate_plane(plane1{i},plane2{i},plane3{i})

%% Movie frame preparation
[movie_frames_pre, movie_frames_post, movie_sect_length] = set_movie_frame(se)

%% PCA
plot_pca(pre, post, movie_frames_pre, movie_frames_post, movie_sect_length, se)

%% Plot Skewness (pre vs post)
[pre_skew, post_skew] = plot_skewness(pre, post)

%% Plot Skewness Change vs Popn correlation
plot_skewness_popn_corr(pre, post, movie_frames_pre,movie_frames_post, npair)


%% Plot Skewness change vs. distance of COM

end