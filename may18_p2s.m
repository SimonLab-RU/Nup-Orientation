function [ firstPOS] = may18_p2s( parentfilename )

% by Joan Pulupa October 2nd 2017
% Automated pvs analysis of the p over s imaging.

%outputs: firstPOS, each column represents different cell, 50 NPCs tested per cell. 
%zeros are NPCs that do not meet quality control.

%% Strategy for Automating the Selection and Intensity Analysis of NPCs

%Step 1: (Function, fileManager): Open current folder, identify 
%files to analyze. Create subfolder for saving the steps of data analysis 
%if subfolders do not exist. 

%Step 2: (Function: selectNucleus): Open image file and allow user to 
%select nucleus manually. If the user has already selected a "map" the 
%previously selected map will be used. This image is saved as a TIF in 

%Step 3: (Function: detectNPCs): Perform blob detection. The program is 
%external (python), but is run through the command line in MATLAB.

%Step 4: (Function:crossCorrelateNPCs): Cross correlate the NPCs generates.

%Step 5: Divide p over s and save matrices of intensities and 

%%

%Step 1: fileManager

[fileNPCImages,strFolderNames]=fileManager(parentfilename);


%%
%Step 2: selectNucleus and create NPC Image map

[filteredP, filteredS ] = selectNucleusanBKFilter(fileNPCImages,strFolderNames,parentfilename);

%%
% %Step 3:detectNPCs
% 
[xcoords,ycoords]=detectNPCs(filteredS,parentfilename,strFolderNames);
% 
% 
% %% 
% %Step 4: crossCorrelateNPCs
% 
% %Align and measure NPC intensities

[fileNameFirstPOS] = crossCorrelateNPCsForFirstPOS(xcoords,ycoords,filteredP,filteredS,strFolderNames,parentfilename,fileNPCImages);
% 
% 
% %%
% %Step 5: gaussFitNPCs
%Perform test of gaussianness to determine if image should be measured and
%extract intensities from each image. 
[firstPOS] = gaussFitForFirstPOS(strFolderNames,parentfilename,fileNameFirstPOS);
%[pIntensitiesOverTime,sIntensitiesOverTime] = gaussFitOverTime(strFolderNames,parentfilename,idImagesOverTime,fileNameOverTimeNPCs);



%plotImage1(filteredP,xcoords,ycoords,firstPOS);

%beep

end


%%
function [fileNPCImages,strFolderNames]=fileManager(parentfilename)

cd(parentfilename);

%create list of images within folder to process:
fileNPCImages = dir('*.tif');


%name output files for record keeping:
strFolderNames = strings(1,length(fileNPCImages));


%find names of images and create folders for saving intermediates if
%analysis folders do not exist

for index=1:size(fileNPCImages)
    lengthFileName=length(fileNPCImages(index).name);
    
    if (lengthFileName >15)
        strFolderNames(index) = fileNPCImages(index).name(1:15);
    else
        strFolderNames(index) = fileNPCImages(index).name;
    end
    
    cd(parentfilename);
    
    if (not(exist(char(strFolderNames(index)), 'File')))
        mkdir(char(strFolderNames(index)));
    end

end

end

%%
function [filteredP, filteredS] = selectNucleusanBKFilter(fileNPCImages,strFolderNames,parentfilename)

%get image width and height: 
FileTif=fileNPCImages(1).name;
InfoImage=imfinfo(FileTif);
mImage=InfoImage(1).Width;
nImage=InfoImage(1).Height;

%create temporary variable to display image for user to select nucleus if
%user has not already selected nucleus
imgToDisplay=zeros(mImage,nImage,length(fileNPCImages));

    for image = 1: length(fileNPCImages)
        %determine if there is a map of image to determine what area to image
        %process (only want to image process nucleus:
        foldName=strcat(parentfilename,'/',strFolderNames(image));
        cd(char(foldName))
        if exist('map.tif','file')
        else

            %open first ten frames of TIFF Stack:
            currentImage=zeros(mImage,nImage);
            cd(parentfilename);
            FileTif=fileNPCImages(image).name;
            TifLink = Tiff(FileTif, 'r');

            for i=1:10
               TifLink.setDirectory(i);
               newFrame=TifLink.read();
               newFrame=im2double(newFrame);
               currentImage=currentImage+newFrame;
            end

            TifLink.close();


            %display z-stack of first 10 frames of TIFF stack, filtered
            se = strel('disk',7);
            imgToDisplay(:,:,image)=imtophat(currentImage,se);
            imshow(imgToDisplay(:,:,image))
            caxis ([0 max(max(imgToDisplay(:,:,image)))])

            %User draws a polygon around the nucleus:
            bw=roipoly;

            %Saves map
            npcMapFilename=strcat(parentfilename,'/',strFolderNames(image),'/','map.tif');
            imwrite(bw,char(npcMapFilename));

        end

    end


%create blank images for background subtracted parental images
filteredP=zeros(mImage,nImage,1,length(fileNPCImages));
filteredS=zeros(mImage,nImage,1,length(fileNPCImages));

%run loop over images to be analysed. In loop, background subtraction will
%be performed and then the user will define region of interest or 
for image = 1: length(fileNPCImages)
    cd(parentfilename);
    FileTif=fileNPCImages(image).name;
    [ filteredP(:,:,:,image), filteredS(:,:,:,image)] = process_even_10_20( FileTif);

end

end

%%
function [xcoords,ycoords] = detectNPCs(filteredS,parentfilename,strFolderNames)
xcoords=zeros(50,length(strFolderNames));
ycoords=zeros(50,length(strFolderNames));
    for index = 1: length(strFolderNames)
        
        npcMapFilename=strcat(parentfilename,'/',strFolderNames(index),'/','map.tif');
               
        nucleusMap=imread(char(npcMapFilename));
        
        thisSImage=filteredS(:,:,:,index);
        thisSImage(nucleusMap==0)=0;
        
        

        %run detection module
        se = strel('disk',7);
        imgtoSearch=imtophat(thisSImage,se);
        [points,x] = kp_log(imgtoSearch,50);
        xcoords(:,index) = points(:,2);
        ycoords(:,index) = points(:,1);

        %Save a file with all selected NPCs marked.

        imshow(filteredS(:,:,:,index));
        caxis ([0 max(max(filteredS(:,:,:,index)))]);
        hold on;
        scatter(xcoords(:,index),ycoords(:,index),1,'r','.');
        cd(parentfilename);cd(char(strFolderNames(index)));
        numLabeledImage=length(dir('marked*.*'))+1;
        nameofMarkedImage=strcat('markedNPCs',int2str(numLabeledImage),'.eps');
        
        saveas(gcf,nameofMarkedImage);
        close all;
    end

end

%%
function [fileNameFirstPOS] = crossCorrelateNPCsForFirstPOS(xcoords,ycoords,filteredP,filteredS,strFolderNames,parentfilename,fileNPCImages)

% by Joan May 2015
% cross-Correlate the cut out NPCs 
   %add this intensity group to the calculations
   

numNPC=size(xcoords,1);


%Create Matrix for Uncorrelated NPC Images
UnCorrNPCImages=[];
NPCImages=[];

%pick size of cut out NPC

finalsize=11;
startsize=21;
halfSideSize=(startsize-1)/2;

bordersize=(startsize-finalsize)/2;

%set a empty matrix for holding image of nucleus
UnCorrNPCImages=zeros(22,22,1);
NPCImages=zeros(11, 11, numNPC,2);

for imageCount=1:length(strFolderNames)
    
    cd(parentfilename);
    FileTif=fileNPCImages(imageCount).name;
    InfoImage=imfinfo(FileTif);
    numImages=length(InfoImage);
    cd(char(strFolderNames(imageCount)));
    
    %Save each iteration of analysis in new folder
    attemptNum=length(dir('AlignedNPCsFirstPos*.*'))+1;
    fileNameFirstPOS=strcat('AlignedNPCsFirstPos',int2str(attemptNum));
    
    
    mkdir(fileNameFirstPOS)
    cd(fileNameFirstPOS);
    
    count=0;

    for npcCount = 1:numNPC

            x=xcoords(npcCount,imageCount);
            y=ycoords(npcCount,imageCount);
            
            
            %cut individual NPCs from image
            for imgFrame = 1:2

                    if imgFrame == 2 %numberiseven
                        thisImage = filteredS(:,:,1,imageCount);
                    else %numberisodd
                        thisImage = filteredP(:,:,1,imageCount);
                    end
                    UnCorrNPCImages(:,:,imgFrame)=imcrop(thisImage,[x-halfSideSize y-halfSideSize halfSideSize*2+1  halfSideSize*2+1]);
                    
                    if imgFrame==1
                        NPCImages(:,:,npcCount,1)=UnCorrNPCImages(bordersize+1:startsize-bordersize, bordersize+1:startsize-bordersize,1);
                    elseif imgFrame==2
                        NPCImages(:,:,npcCount,2)=UnCorrNPCImages(bordersize+1:startsize-bordersize, bordersize+1:startsize-bordersize,2); 
                    end

            end

            count=count+1;
            
            %save cut out NPCs as images
            for pol=1:2

                    thisNPCTiff=NPCImages(:,:,npcCount,pol);
                    thisNPCTiff=uint16(thisNPCTiff);

                    if pol==1
                        fileName=strcat('alignedNPC',int2str(count),'.tif');
                        imwrite(thisNPCTiff,fileName); 
                    else
                        imwrite(thisNPCTiff,fileName,'WriteMode','append'); 
                    end

            end

    end

    clear NPCImages  
    
end
    
    

end

%%
function  [firstPOS]=gaussFitForFirstPOS(strFolderNames,parentfilename,fileNameFirstPOS)

pIntensities = zeros(50,length(strFolderNames));
sIntensities = zeros(50,length(strFolderNames));
firstPOS = zeros(50,length(strFolderNames));


for numImages = 1 : length(strFolderNames)
    cd(parentfilename)
    strFolderNames(numImages)
    cd(char(strFolderNames(numImages)))

    circle=[0,0,1,0,0;0,1,1,1,0;1,1,1,1,1;0,1,1,1,0;0,0,1,0,0];


    %crop out nonTIFF files
    cd(fileNameFirstPOS)
    fileNames = dir('*.tif');
    names=cell(length(fileNames),1);

    fit = zeros(2,7,length(fileNames));
    
    

    fn = fullfile(strcat(char(strFolderNames(numImages)),'MarkedNPCs'));
    mkdir(fn);

    for i=1:length(fileNames)
        name=fileNames(i).name;
        names{i,1}=name;
        thisFit=FitTiffGauss20(fileNames(i).name);
 
        
        fit(:,:,i)=thisFit;
        clear thisFit
        

        lowerTestX=4<fit(:,5,i);
        upperTestX=fit(:,5,i)<8;
        lowerTestY=4<fit(:,6,i);
        upperTestY=fit(:,6,i)<8;

  
        
        test1= [lowerTestX(1),upperTestX(1),lowerTestY(1),upperTestY(1)];
        alltest1=all(test1);
        
        test2= [lowerTestX(2),upperTestX(2),lowerTestY(2),upperTestY(2)];
        alltest2=all(test2);

        if ( or(alltest1, alltest2) )

            for index = 1:2
                name=fileNames(i).name;
                thisNPC = imread(fileNames(i).name,index);
                thisNPC = double(thisNPC);


                if mod(index,2)==0 %numberiseven

                    %do not assign new x y

                    intensityNPC=thisNPC.*mask;
                    sIntensities(i,numImages)=max(max(intensityNPC));

                else %numberisodd
                    
                    if (~alltest1)
                        x_coord = round(fit(2,5,i));
                        y_coord = round(fit(2,6,i));
                    else
                        x_coord = round(fit(1,5,i));
                        y_coord = round(fit(1,6,i));
                    end
                    
                    mask = zeros(size(thisNPC));
                    mask(x_coord-2:x_coord+2,y_coord-2:y_coord+2)=circle;
                    intensityNPC=thisNPC.*mask;
                    pIntensities(i,numImages)=max(max(intensityNPC));

                    fileName=strcat(fn,'/',name,'_',int2str(index),'.tif');

                    imshow(thisNPC);
                    caxis ([0 1000])
                    hold on;
                    scatter(x_coord,y_coord,1,'r','.');


                    if index==1
                        saveas(gcf,fileName); 
                        close all;
                    else
                        saveas(gcf,fileName); 
                        close all;
                    end

                end



            end



        for indexfirstPoS = 1:2

            if indexfirstPoS == 1
                thisNPC = imread(fileNames(i).name,1);
                thisNPC = double(thisNPC);
                
                 if (~alltest1)
                        x_coord = round(fit(2,5,i));
                        y_coord = round(fit(2,6,i));
                 else
                        x_coord = round(fit(1,5,i));
                        y_coord = round(fit(1,6,i));
                 end
               
                mask = zeros(size(thisNPC));
                mask(x_coord-2:x_coord+2,y_coord-2:y_coord+2)=circle;
                intensityNPC=thisNPC.*mask;
                firstPIntensity=max(max(intensityNPC));

            elseif indexfirstPoS == 2

                thisNPC = imread(fileNames(i).name,2);
                thisNPC = double(thisNPC);
                intensityNPC=thisNPC.*mask;
                firstSIntensity=max(max(intensityNPC));
                firstPOS(i,numImages)=firstPIntensity/firstSIntensity;
            end
        end

        end

    end
end

end


 
%%
%%   
function [ filteredP, filteredS] = process_even_10_20( fileName )
%dataProcessing.m By Joan Pulupa, January 2017 
%   This function takes raw imaging data and preforms a background
%   subtraction as well as summing of 

FileTif=fileName;
InfoImage=imfinfo(FileTif);
mImage=InfoImage(1).Width;
nImage=InfoImage(1).Height;
NumberImages=length(InfoImage);
originalImage=zeros(nImage,mImage,NumberImages,'uint16');
 
TifLink = Tiff(FileTif, 'r');
for i=1:NumberImages
   TifLink.setDirectory(i);
   originalImage(:,:,i)=TifLink.read();
end
TifLink.close();

%subtract background
originalImage=originalImage-100;

npcImage=cat(3,originalImage(:,:,1:20));


filteredP=zeros(nImage,mImage,1,'uint16');
filteredS=zeros(nImage,mImage,1,'uint16');

filteredP(:,:,1) = 1 * npcImage(:,:,2) + 1 * npcImage(:,:,4) + 1 * npcImage(:,:,6) + 1 * npcImage(:,:,8) + 1 * npcImage(:,:,10) + 1 * npcImage(:,:,12) + 1 * npcImage(:,:,14) + 1 * npcImage(:,:,16) + 1 * npcImage(:,:,18) + 1 * npcImage(:,:,20);
filteredS(:,:,1) = 1 * npcImage(:,:,1) + 1 * npcImage(:,:,3) + 1 * npcImage(:,:,5) + 1 * npcImage(:,:,7) + 1 * npcImage(:,:,9) + 1 * npcImage(:,:,11) + 1 * npcImage(:,:,13) + 1 * npcImage(:,:,15) + 1 * npcImage(:,:,17) + 1 * npcImage(:,:,19);




end

 
 
%%
function fit = FitTiffGauss20(imageName)

options=optimset('Display','off');

fit = zeros(2,7);

exitflag = 0;


for m = 2:-1:1
    X = RearrangeGaussArray(double(mat2gray((imread(imageName,m)))));
    
    if exitflag == 1
        testFit = fit(m+1,:);
        
    else
        testFit = [0 max(X(3,:)) 3 3 6 6 0];
        %amp was 200
    end
    
    [Estimates, fval, exitflag]=fminsearch(@GaussErrorTheta,testFit,options,X);
    
    if exitflag == 1
        fit(m, :) = Estimates;
        
    end
end
end

function [points,blob_candidate_value] = kp_log(img,o_nb_blobs)
    % Extract keypoints using Laplacian of Gaussian (LoG) algorithm
    %
    % Author :: Vincent Garcia
    % Date   :: 05/12/2007
    %
    % INPUT
    % =====
    % img        : the graylevel image
    % o_nb_blobs : (optional) number of blobs detected
    %
    % OUTPUT
    % ======
    % points : the interest points extracted
    %
    % REFERENCES
    % ==========
    % Lindeberg, T. Feature Detection with Automatic Scale Selection
    % IEEE Transactions Pattern Analysis Machine Intelligence, 1998, 30, 77-116kp_harris(im)
    %
    % EXAMPLE
    % =======
    % points = kp_log(img)
   
    
    % input image
    img = double(img(:,:,1));
        
    iMax = max(max(img));
iMin = min(min(img));

img = (img-iMin) / (iMax-iMin);
img = imcomplement(img);
    % number of blobs detected
    if nargin==1
        nb_blobs = 120;
    else
        nb_blobs = o_nb_blobs;
    end
    
    % Laplacian of Gaussian parameters
    sigma_begin = 2;
    sigma_end   = 15;
    sigma_step  = 1;
    sigma_array = sigma_begin:sigma_step:sigma_end;
    sigma_nb    = numel(sigma_array);
        
    % variable
    img_height  = size(img,1);
    img_width   = size(img,2);
        
    % calcul scale-normalized laplacian operator
    snlo = zeros(img_height,img_width,sigma_nb);
    for i=1:sigma_nb
        sigma       = sigma_array(i);
        snlo(:,:,i) = sigma*sigma*imfilter(img,fspecial('log', floor(6*sigma+1), sigma),'replicate');
    end
        
    % search of local maxima
    snlo_dil             = imdilate(snlo,ones(3,3,3));
    blob_candidate_index = find(snlo==snlo_dil);
    blob_candidate_value = snlo(blob_candidate_index);
    min(blob_candidate_value);
    [tmp,index]          = sort(blob_candidate_value,'descend');
    blob_index           = blob_candidate_index( index(1:min(nb_blobs,numel(index))) );
    [lig,col,sca]        = ind2sub([img_height,img_width,sigma_nb],blob_index);
    points               = [lig,col,3*reshape(sigma_array(sca),[size(lig,1),1])];
    
end

function GaussData = RearrangeGaussArray(Data);

dim = size(Data);
GaussData = zeros(dim(1)*dim(2),3);  % Preallocate matrix
for m = 1:dim(1)
    for n = 1:dim(2)
        GaussData((m-1)*dim(2)+n,1) = m;
        GaussData((m-1)*dim(2)+n,2) = n;
        GaussData((m-1)*dim(2)+n,3) = Data(m,n);
    end
end
end

function E=GaussErrorTheta(params,GaussData)
A=params(1);     %offset
B=params(2);     %amplitude
C1=params(3);    %x width
C2=params(4);    %y width
D1=params(5);    %x offset
D2=params(6);    %y offset
theta=params(7); %rotation angle

X = GaussData(:,1);
Y = GaussData(:,2);
Z = GaussData(:,3);

Fitted_Curve=A+(B*exp(-(((X-D1)*cos(theta)+(Y-D2)*sin(theta)).^2)/(2*C1^2)-((-(X-D1)*sin(theta)+(Y-D2)*cos(theta)).^2)/(2*C2^2)));
Error_Vector=Fitted_Curve - Z;

E=sum(Error_Vector.^2);

%Utilize fminsearch to perform Simplex on the error function.  Below is an
%example.  

%Estimates=fminsearch(@GaussError,[1000 3000 2 2 7 7],options,Qadj)

%options=optimset('Display','iter');
end