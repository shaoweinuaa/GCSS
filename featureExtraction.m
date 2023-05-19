clear all
clc
close all

Root='.\\Image\\';
maskRoot='.\\Mask\\';
saveMatFile='.\\SAVEMAT\\';
sampleDir=dir(Root);
superPixelNumber=30;
for i=3:length(sampleDir)
    t=0;
    Node=[];
%     sampleDir(i).name='ZXF_ROI005_ROI_005';
    sampleRoot=strcat(Root,sampleDir(i).name);
    sampleMaskRoot=strcat(maskRoot,sampleDir(i).name);
    saveFolder=strcat('.\\Superpixel\\',sampleDir(i).name);
    if ~exist(saveFolder,'dir')
        mkdir(saveFolder)
    end
    fileDir=dir(sampleRoot);
    
    for j=3:length(fileDir)
        fileName=fileDir(j).name;
        cellType=str2num(fileName(1:length(fileName)-4))
         
        filePath=strcat(sampleRoot,'\\',fileName);  
        fileMaskPath=strcat(sampleMaskRoot,'\\',fileName(1:length(fileName)-4),'.mat');
        img=imread(filePath);
        mask=load(fileMaskPath);
        mask=mask.mask;
        A=myRegionprops(mask)
        [L,NumLabels] = superpixels(img,superPixelNumber);
        BW = boundarymask(L);

        for k=1:NumLabels
            cellFeature=[];
            idxLable=find(L==k);
            [yPos,xPos]=find(L==k);
            pixelPos=[xPos,yPos];
            for ll=1:length(A)
                 cellPixel=A(ll).PixelIdxList;
                 inter_sect=intersect(cellPixel,idxLable);
                 if length(inter_sect)>0
                    temp=[cellType,A(ll).Centroid,A(ll).Area,A(ll).MajorAxisLength,A(ll).MinorAxisLength,A(ll).MajorAxisLength/A(ll).MinorAxisLength,A(ll).Circularity,A(ll).Eccentricity];
                    if size(cellFeature,1)==0 || norm(cellFeature-temp,2)~=0
                       cellFeature=[cellFeature;temp];
                    end
                
                 end
            end
            
            if size(cellFeature,1)==1
               
                t=t+1;
                cellFeature
                Node(t,:)=[cellType,cellFeature(:,2:9),zeros(1,6),zeros(1,3),zeros(1,3),1];
%                 saveFileName=strcat(saveFolder,'\\',sampleDir(i).name,'_SID_',num2str(t),'.mat');
% % %                 save(saveFileName,'t','cellType','pixelPos');

            elseif size(cellFeature,1)==2
               
               t=t+1;
               Centroid=cellFeature(:,2:3);
               dist=norm((Centroid(1,:)-Centroid(2,:)),2)
               Node(t,:)=[cellType,median(cellFeature(:,2:9)),std(cellFeature(:,4:9)),ones(1,3)*dist,zeros(1,3),2];    
%                saveFileName=strcat(saveFolder,'\\',sampleDir(i).name,'_SID_',num2str(t),'.mat');
%                save(saveFileName,'t','cellType','pixelPos');
        
            elseif size(cellFeature,1)>2   
               t=t+1;
               Centroid=cellFeature(:,2:3);   
               DT = delaunayTriangulation(double(Centroid));
               E = edges(DT);
               for k = 1:size(Centroid, 1)
                  edgesCell = E(sum(E==k, 2)~=0, :);
                  dist = zeros(size(edgesCell, 1), 1);
                  for m = 1:numel(dist)
                      p1 = Centroid(edgesCell(m, 1), :);
                      p2 = Centroid(edgesCell(m, 2), :);
                      dist(m) = norm(p1-p2);
                  end
                feas_dist(k, :) = [mean(dist), max(dist), min(dist)];
               end
               Node(t,:)=[cellType,median(cellFeature(:,2:9)),std(cellFeature(:,4:9)),median(feas_dist),std(feas_dist),size(cellFeature,1)];            
%                saveFileName=strcat(saveFolder,'\\',sampleDir(i).name,'_SID_',num2str(t),'.mat');
%                save(saveFileName,'t','cellType','pixelPos');
        
            end

                
                
            
            
            
            
            
            end
        
        
        end
        
        

 Node(isinf(Node))=0; 
% % %     Node(:,4:13)=1./(1+exp(-(zscore(Node(:,4:13)))));
 saveFileName=strcat(saveMatFile,'\\',sampleDir(i).name,'.mat');
 save(saveFileName,'Node')
end



