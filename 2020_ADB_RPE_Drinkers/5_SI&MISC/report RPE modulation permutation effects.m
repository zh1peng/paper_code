% report RPE modulation permutation effects
%% Model1
tmp=1:length(t1)
frn_points=tmp(t1>up_t1);
disp(['conscutive points number: ' num2str(length(frn_points(1:end)))])
disp(['time window left: ' num2str(frn_points(1)./0.512), ' ms in idx space: ', num2str(frn_points(1)+102)])
disp(['time window right: ' num2str(frn_points(end)./0.512),' ms in idx space: ', num2str(frn_points(end)+102)])

p3_points=tmp(t1<bot_t1);
disp(['conscutive points number: ' num2str(length(p3_points(10:end)))])
disp(['time window left: ' num2str(p3_points(10)./0.512), ' ms in idx space: ', num2str(p3_points(10)+102)])
disp(['time window right: ' num2str(p3_points(end)./0.512),' ms in idx space: ', num2str(p3_points(end)+102)])


%% Model2 |RPE|
tmp=1:length(t2)
frn_points=tmp(t2>up_t2);
disp(['conscutive points number: ' num2str(length(frn_points(1:end)))])
disp(['time window left: ' num2str(frn_points(1)./0.512), ' ms in idx space: ', num2str(frn_points(1)+102)])
disp(['time window right: ' num2str(frn_points(end)./0.512),' ms in idx space: ', num2str(frn_points(end)+102)])

p3_points=tmp(t2<bot_t2);
disp(['conscutive points number: ' num2str(length(p3_points(1:end)))])
disp(['time window left: ' num2str(p3_points(1)./0.512), ' ms in idx space: ', num2str(p3_points(1)+102)])
disp(['time window right: ' num2str(p3_points(end)./0.512),' ms in idx space: ', num2str(p3_points(end)+102)])

% sign of RPE
tmp=1:length(t2_2)
frn_points=tmp(t2_2>up_t2_2);
disp(['conscutive points number: ' num2str(length(frn_points(1:end)))])
disp(['time window left: ' num2str(frn_points(1)./0.512), ' ms in idx space: ', num2str(frn_points(1)+102)])
disp(['time window right: ' num2str(frn_points(end)./0.512),' ms in idx space: ', num2str(frn_points(end)+102)])

p3_points=tmp(t2_2<bot_t2_2);
disp(['conscutive points number: ' num2str(length(p3_points(1:end)))])
disp(['time window left: ' num2str(p3_points(1)./0.512), ' ms in idx space: ', num2str(p3_points(1)+102)])
disp(['time window right: ' num2str(p3_points(end)./0.512),' ms in idx space: ', num2str(p3_points(end)+102)])