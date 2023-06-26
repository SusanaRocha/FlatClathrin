%% set file path: manually select folder. Load data.
fp = uigetdir;
cd(fp)
list=dir('*.mat'); 
list={list.name};

%% threshold cells
for f=1:size(list, 2)
    fn=list{1,f}
    load([fn]) 
    cod2 = round(pts(:,4:5)); % take x and y coordinates of the localizations (column 4 = X position, column 5 = Y position in the pts function)
    cod2(cod2(:,1)<1, :)=[];
    cod2(cod2(:,2)<1, :)=[];
    hst2 = accumarray(cod2(:, [2, 1]), 1, [512, 512]);
    hst2=single(hst2);
    bw = imbinarize(hst2, 1); 
    bw2 = bwmorph(bw, 'dilate', 2); 
    bw3 = imfill(bw2, 'holes'); 
    bw4 = bwmorph( bw3, 'erode', 2); 
    figure, imagesc(bw4), axis xy

%% select 1 cell = the object with the largest surface area, and only take localizations inside this area into account
    cc = bwconncomp(bw4); 
    props = regionprops(cc, 'area'); 
    idx = find([props.Area] == max([props.Area])); 
    BW = ismember(labelmatrix(cc), idx); 
    B = cell2mat(bwboundaries(BW)); 
    in = inpolygon(pts(:,4),pts(:,5),B(:,2),B(:,1)); 
%% plot all localizations inside the selected area as green circles
    plot(B(:,2),B(:,1), '-r')
    hold on
    selpts = pts(in, [4,5]);
    figure, plot(selpts(:,1), selpts(:,2), 'og'), axis xy
%% knnsearch: find k-nearest neighbours
    [~, D] = knnsearch(selpts, selpts, 'K', 2);
%% plot diff limited image 2048 x 2048
    cm=1/25:1/25:1;cm=cm';
    pe=[0 0 0];
    pe=cat(1, pe, [cm cm/2 zeros(25, 1)]);
    pe=cat(1, pe, [ones(25, 1) 0.5+cm/2 zeros(25, 1)]);
    pe=cat(1, pe, [ones(25, 1) ones(25, 1) cm]);
    colormap(pe)
    mx = 512/0.25;
    cod = floor(pts(:,4:5)/0.25)+1;
    cod(cod(:,1)<1, :)=[];
    cod(cod(:,2)<1, :)=[];
    hst = accumarray(cod(:, [2, 1]), 1, [mx, mx]);
    hst=single(hst);
    hsts=imgaussfilt(hst,1);
    imagesc(hsts, [0 7]); colormap(pe);
    axis xy equal tight %ij xy
    line([1812.49,2000],[50,50],'Color','w','LineWidth',2) % put scalebar of 5µm

%% DBScan analysis with minimum 17 neighbors (MinPts) and search radius (ε) 0.55
    DB = dbscan(selpts, 0.55, 17);
    max(DB); % max number of clusters
    selpts = selpts.*0.107; % convert pixel to µm
%% calculate for all the clusters: the number of the cluster, number of localizations (NrPts), area, density, eccentricity and center of a plotted ellipse
    prop_cl = table('Size', [max(DB), 6], 'VariableTypes',{'single','single','single','single','single','single'}, 'VariableNames', {'ID', 'nrpts', 'area', 'density', 'ecc', 'center'});
    prop_cl.ID = [1:max(DB)]';
    flg_plot = 0; 
    h = waitbar(0, 'Calculating properties of the clusters...');
    for i=1:max(DB)
        cl_each = selpts(find(DB == i), 1:2);
        boundary(cl_each(:,1), cl_each(:,2)); %find edges of the cluster
        bx = boundary(cl_each(:,1), cl_each(:,2));
        prop_cl.area(i) = polyarea(cl_each(bx,1), cl_each(bx,2)); %calculate area of cluster. in µm²
        prop_cl.nrpts(i) = length(cl_each);  %calculate density of the cluster. nrpts per um²
        prop_cl.density(i) = length(cl_each)/prop_cl.area(i);
    try   %calculate eccentricity using the MinVolEllipse function (fit with an ellipse)
        [A, c] = MinVolEllipse (cl_each', .01); 
        prop_cl.center(i, 1:2) = c';
        [~, Q, ~] = svd(A); r1 = 1/sqrt(Q(1,1)); r2 = 1/sqrt(Q(2,2));
        o = max(r1, r2); p = min(r1, r2); prop_cl.ecc(i) = sqrt (1-p^2/o^2);
        if flg_plot
            plot(cl_each(:,1), cl_each(:,2), '.k');
            hold on, plot(cl_each(bx,1), cl_each(bx,2), 'or'), title(i);
            Ellipse_plot(A, c)
            hold off
         pause
        end
    catch 
       disp([fn ' not analyzed'])
    end
        waitbar(i/max(DB), h)
    end
    close(h)

%% calculate cell-specific parameters
    cell.area = max([props.Area])*0.0114; %cell area in um²
    cell.density = max(DB)/cell.area; %number of clusters per um2
    cell.ratio = sum(prop_cl.area)/cell.area; %ratio of area occupied by clusters


%% only retain clusters with ≥ 50 localizations and calculate again cell-specific parameters
    prop_50 = prop_cl;
    prop_50(prop_50.nrpts<50,:)=[];
    cell50.area = max([props.Area])*0.0114;
    cell50.density = size(prop_50, 1)/cell.area; %number of clusters per um2
    cell50.ratio = sum(prop_50.area)/cell.area; %ratio of area occupied by clusters

%% calculate mean and median cluster specific parameters per cell (both for all clusters and clusters ≥ localizations)
    MM = table('Size', [2, 6], 'VariableTypes',{'single','single','single','single','single','single'}, 'VariableNames', {'area_mean', 'area_median', 'density_mean', 'density_median', 'ecc_mean', 'ecc_median'}, 'RowNames', {'all loc', '>50'});
    MM (1,:) = {mean(prop_cl.area), median(prop_cl.area), mean(prop_cl.density), median(prop_cl.density), mean(prop_cl.ecc), median(prop_cl.ecc)};
    MM (2,:) = {mean(prop_50.area), median(prop_50.area), mean(prop_50.density), median(prop_50.density), mean(prop_50.ecc), median(prop_50.ecc)};

%% save 
    file = [fp fn];
    filename = file(27:strfind([fp fn], '.mat')-1); % adjust the number 27 depending on the name of the file in your pc
    outfile = [filename '_DBscan.mat'];
    save(outfile, 'DB', 'selpts', 'file', 'cell', 'prop_cl', 'cell50', 'prop_50','MM');
    filename_prop_50 = [filename '_prop_50.txt'];
    filename_prop_cl = [filename '_prop_cl.txt'];
    filename_MM = [filename '_means_medians.txt'];
    writetable(prop_50, filename_prop_50);
    writetable(prop_cl, filename_prop_cl);
    writetable(MM, filename_MM, 'WriteRowNames',true);
% Close figures and clean variables
  close all
  clearvars -except f fn fp list pts
end