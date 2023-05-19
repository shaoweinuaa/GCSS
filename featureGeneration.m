clear all
clc
close all
clusterTypeNumber=23;
load punctureInfo.mat
load cellType.mat

for i=1:length(sampleName)
    if i==1
       currSampleName=sampleName{i};
       currROIName=roiName{i};
       nucleiPath=strcat('.\\Mat\\',currSampleName,'_',currROIName,'.mat');
       load(nucleiPath)
       clusterID(i)=find(strcmp(strtrim(clusterName{i}),strtrim(cellType))==1);
       pixelID=find(segmentaion==cellID(i));
       for o=1:clusterTypeNumber
           eval(['cluster',num2str(o),'=','zeros(size(segmentaion,1),size(segmentaion,2))',';'])
       end
       eval(['cluster',num2str(clusterID(i)),'(pixelID)','=255;'])
       

       
       
       
    else
       if strcmp(currSampleName,sampleName{i})==0 || strcmp(currROIName,roiName{i})==0
           saveFilePath=strcat('.\\Image\\',currSampleName,'_',currROIName);
           saveFilePath_mask=strcat('.\\Mask\\',currSampleName,'_',currROIName);
           if ~exist(saveFilePath)
               mkdir(saveFilePath);
           end
           if ~exist(saveFilePath_mask)
               mkdir(saveFilePath_mask);
           end
           for o=1:clusterTypeNumber
               mask=segmentaion;
               eval(['mask(find(cluster',num2str(o),'==0))=0;'])
%                segmentaion(find(cluster1==0))=0
               saveFile=strcat(saveFilePath,'\\',num2str(o),'.png');
               saveFile_mask=strcat(saveFilePath_mask,'\\',num2str(o),'.mat');
               eval(['imwrite(cluster',num2str(o),',saveFile);'])
               save(saveFile_mask,'mask');

           end
           
           
          currSampleName=sampleName{i};
          currROIName=roiName{i};
          nucleiPath=strcat('.\\MAT\\',currSampleName,'_',currROIName,'.mat');
          load(nucleiPath)
          clusterID(i)=find(strcmp(strtrim(clusterName{i}),strtrim(cellType))==1);
          pixelID=find(segmentaion==cellID(i));
          for o=1:clusterTypeNumber
           eval(['cluster',num2str(o),'=','zeros(size(segmentaion,1),size(segmentaion,2))',';'])
          end
          eval(['cluster',num2str(clusterID(i)),'(pixelID)','=255;'])
       

       else
           clusterID(i)=find(strcmp(strtrim(clusterName{i}),strtrim(cellType))==1);
           pixelID=find(segmentaion==cellID(i));
           eval(['cluster',num2str(clusterID(i)),'(pixelID)','=255;'])
 
       end
       end
 if i==length(sampleName)
    saveFilePath=strcat('.\\Image\\',currSampleName,'_',currROIName);
    saveFilePath_mask=strcat('.\\Mask\\',currSampleName,'_',currROIName);

    if ~exist(saveFilePath)
        mkdir(saveFilePath);
    end 
    if ~exist(saveFilePath_mask)
        mkdir(saveFilePath_mask);
    end
    
    for o=1:clusterTypeNumber
        mask=segmentaion;
        eval(['mask(find(cluster',num2str(o),'==0))=0;'])
%                segmentaion(find(cluster1==0))=0
        saveFile=strcat(saveFilePath,'\\',num2str(o),'.png');
        saveFile_mask=strcat(saveFilePath_mask,'\\',num2str(o),'.mat');
        eval(['imwrite(cluster',num2str(o),',saveFile);'])
        save(saveFile_mask,'mask');

    end  
 end     
 end



