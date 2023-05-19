clear all
clc
close all

root='.\\Image\\';

channels=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23];
labels=[];
dirRoot=dir(root);
for i=3:length(dirRoot)
  roiName=dirRoot(i).name;
  roiRoot=strcat(root,roiName,'\\');
  A=[];
  for j=1:length(channels) 
      img=imread(strcat(roiRoot,num2str(channels(j)),'.png'));
      if j==1
         A=zeros(size(img,1),size(img,2),length(channels));
         A(:,:,j)=img;
      else
         A(:,:,j)=img; 
      end
  end
  IMC_Similarity(i-2,:)=sim_split(A,4);
  sampleName{i-2}=roiName;

end
save mean_IMC_Sim_1.mat IMC_Similarity sampleName