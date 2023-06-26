%% set file path: manually select folder
fp = uigetdir;
cd(fp)
%% load in parameters for graphical representation
cm=1/25:1/25:1;cm=cm';
pe=[0 0 0];
pe=cat(1, pe, [cm cm/2 zeros(25, 1)]);
pe=cat(1, pe, [ones(25, 1) 0.5+cm/2 zeros(25, 1)]);
pe=cat(1, pe, [ones(25, 1) ones(25, 1) cm]);
colormap(pe)
%% create list of.HIS files in the selected folder
list=dir('*.HIS'); 
list={list.name};
h = waitbar(0, 'Finding points...');
%% for all .HIS files in the list, localize single molecules using the Localizer function with PSF standard deviation factor 1.8 and intensity selection sigma factor 25
for fn=1:size(list, 2)
    clear 'pts' 'im'
    file=list{fn};
    if isfile([file(1:strfind(file, '.HIS')-1) '_LocRes.mat']) 
        rep = questdlg('File already analysed.', 'File exists', 'Repeat Analysis', 'Use existing file', 'Use existing file');
        if strcmp(rep,'Repeat Analysis')
            try
                im = Localizer('readccdimages', 0, -1, file); %from frame 0 to all frames (-1)
                pts = Localizer('localize', 1.8, 'glrt', 25, '2DGauss', im); %Localizer function fits 2D Gaussian with PSF standard deviation factor 1.8 and intensity selection sigma factor 25
                save([file(1:strfind(file, '.HIS')-1) '_LocRes.mat'], 'pts') %save information of localized points (pts). Filename = core name of the .HIS file followed by _LocRes.mat
            catch
                disp([file ' was not analysed!'])
            end
        else
            load([file(1:strfind(file, '.HIS')-1) '_LocRes.mat']);
        end
    else
        try
            im = Localizer('readccdimages', 0, -1, file);
            pts = Localizer('localize', 1.8, 'glrt', 25, '2DGauss', im);
            save([file(1:strfind(file, '.HIS')-1) '_LocRes.mat'], 'pts')
        catch
            disp([file ' was not analysed!'])
        end
    end
    % Plot reconstructed image in a 4x enlarged figure of 2048x2048 pixels (512 pixels/0.25). In the 512x512 image, 1 pixel = 106.66 nm. In the 4x enlarged image, 1 pixel = 26.665 nm)  
    if exist('pts', 'var')
        mx = 512/0.25;
        cod = floor(pts(:,4:5)/0.25)+1;
        cod(cod(:,1)<1, :)=[];
        cod(cod(:,2)<1, :)=[];
        hst = accumarray(cod(:, [2, 1]), 1, [mx, mx]);
        hst=single(hst);
        hsts=imgaussfilt(hst,1);
        imagesc(hsts, [0 7]); colormap(pe); %for adjusting LUT scale, for visual purposes only
        axis xy equal tight
        line([1812.49,2000],[50,50],'Color','w','LineWidth',2) % put scalebar of 5µm. 1 pixel = 26.665 nm, so 5 µm = 187.51 pixels.
       
        title(file(1:strfind(file, '.HIS')-1))
        print([file(1:strfind(file, '.HIS')-1) '.png'], '-dpng'); % export reconstructed image as png
    end
    waitbar(fn/size(list, 2), h);
end

%% combine 2 images with the same core name

list=dir('*LocRes.mat'); 
list={list.name}; 

while ~isempty(list)
    file = list{1};
    core = file(1: strfind(file, '_mov')-1); 
    index = find(contains(list, core));
    fr = 0; pts_total =[];
    for i=1:length(index)
        load(list{i});
        pts(:,1) = pts(:,1) + fr;
        pts_total = [pts_total; pts];
        fr = pts(end, 1) + 1;
    end
    pts = pts_total;
    save([core '_Res_all_movies.mat'], 'pts'); %  %save combined information of localized points (pts). Filename = core name of the .HIS file followed by _Res_all_movies.mat 
    mx = 512/0.25;
        figure
        cod = floor(pts(:,4:5)/0.25)+1;
        cod(cod(:,1)<1, :)=[];
        cod(cod(:,2)<1, :)=[];
        hst = accumarray(cod(:, [2, 1]), 1, [mx, mx]);
        hst=single(hst);
        hsts=imgaussfilt(hst,1);
        imagesc(hsts, [0 5]); colormap(pe);
        axis xy equal tight 
        line([1812.49,2000],[50,50],'Color','w','LineWidth',2)
        title(core)
        print([core 'all.png'], '-dpng');
    list(index)=[];
end
    