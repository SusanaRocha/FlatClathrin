function  [prop_cl cellarea] = get_prop_clusters(file, flg, th)
%% load in parameters for graphical representation
cm=1/25:1/25:1;cm=cm';
pe=[0 0 0];
pe=cat(1, pe, [cm cm/2 zeros(25, 1)]);
pe=cat(1, pe, [ones(25, 1) 0.5+cm/2 zeros(25, 1)]);
pe=cat(1, pe, [ones(25, 1) ones(25, 1) cm]);
colormap(pe)

%% load data. Only analyze clusters with ≥ 50 localizations
 load(file, 'prop_50', 'DB', 'selpts')  
 ind_DB50 = ismember(DB, prop_50.ID);
 pts_DB = selpts(ind_DB50,:);
 %% calculate area of the cell
 j = boundary(selpts(:,1), selpts(:,2),0.8);
 cellarea = polyarea(selpts(j,1), selpts(j,2));
 
 %% create empty prop_cl table
 prop_cl = table('Size', [max(DB), 10], 'VariableTypes',{'single','single','single','single','single','single', 'single', 'single', 'single', 'string'}, 'VariableNames', {'ID', 'nrpts','center', 'area', 'perimeter', 'density', 'ecc', 'distance', 'classification_model', 'type'});
 
 %% prepare to plot data
 if flg 
     sx = round(512*50*0.107); 
     sy = sx;
     BW = zeros(sy, sx);
     BWcp=BW; BWfcl=BW;
 end
 %% For every cluster with ≥ 50 localizations, calculate cluster properties and classify as a pit or a lattice.
 for i = 1:size(prop_50,1)
     prop_cl.ID(i) = prop_50.ID(i);
     pts_cl = selpts(find(DB == prop_50.ID(i)), :); 
     bb = boundary(pts_cl(:,1), pts_cl(:,2), 0.6); % find edges cluster

     % calculate cluster properties
     prop_cl.area(i) = polyarea(pts_cl(bb,1), pts_cl(bb,2)); 
     prop_cl.nrpts(i) = length(pts_cl);
     prop_cl.density(i) = length(pts_cl)/prop_cl.area(i);  
     [A, c] = MinVolEllipse (pts_cl', .01); 
     prop_cl.center(i, 1:2) = c'; 
     [~, Q, ~] = svd(A); r1 = 1/sqrt(Q(1,1)); r2 = 1/sqrt(Q(2,2)); 
     prop_cl.ecc(i) = axes2ecc(max(r2,r1), min(r2,r1));
     prop_cl.perimeter(i) = perimeter(polyshape(pts_cl(bb,1), pts_cl(bb,2))); 
     ind = ismember(pts_DB, pts_cl, 'rows');
     DD = pdist2(pts_cl(bb,:),pts_DB(~ind,:));
     prop_cl.distance(i) = min(DD(:));

     % apply cluster classification model on every cluster, descide if cluster is a classical CCS (pit) or alternative CCS (lattice)
     prop_cl.classification_model(i)=((prop_50.nrpts(i)./1000).*prop_cl.area(i).*prop_cl.perimeter(i).*prop_cl.ecc(i))./(min(DD(:)));
     if prop_cl.classification_model(i) < th; prop_cl.type(i) = 'pit'; else; prop_cl.type(i) = 'lattice'; end

     % for plotting the figure
     if flg
        ROI = poly2mask(50*pts_cl(bb,1), 50*pts_cl(bb,2), sx, sy);
        BW = BW + ROI.*prop_cl.classification_model(i);
     end
 end

 prop_cl(prop_cl.ID(:,1)==0,:)=[];
 disp('done')

 disp('BW image calculated!')
 %% plot figure with 2 pannels. Left side = classified image (white = pits, pink = FCLs). Right side = super-resolution image for comparison.
 if flg
     % calculate classified image
     BWcp(BW(:)<th & BW(:)>0)=1; BWfcl(BW(:)>=th)=1;
     BW_RGD=zeros(size(BW,1), size(BW,2), 3);
     BW_RGD(:,:,1)=BWcp+BWfcl; BW_RGD(:,:,2)=BWcp+0.0736*BWfcl; BW_RGD(:,:,3)=BWcp+0.6471*BWfcl;
     h=figure; set(gcf, 'position', [300 300 900 450], 'color', 'w');
     subplot(1,2,1);image(BW_RGD); axis tight equal off; title('white = pits; pink = FLCs');
     
     % calculate super-resolution image 
     mx = round((512/0.05)*.107);
     cod = floor(selpts(:,1:2)/0.05)+1;
     cod(cod(:,1)<1, :)=[];
     cod(cod(:,2)<1, :)=[];
     hst = accumarray(cod(:, [2, 1]), 1, [mx, mx]);
     hst=single(hst);
     hsts=imgaussfilt(hst,1);
     subplot(1,2,2);imagesc(hsts, [0 10]); colormap(pe)
     axis equal tight 
     drawnow
     pause (0.01)
     outfile = [file(1:strfind(file, '.mat')-1) '_cluster_classification.png'];
     print(outfile, '-dpng')
     close(h)
 end