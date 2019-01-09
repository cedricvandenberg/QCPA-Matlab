function QCPA_Version_1_32 %(Comment this for debugging, together with 'end' at the
%end of code, before functions)

%%%%%%%%%%%%%%%% QCPA version 1.32 (30.03.2018)%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Release notes:

% Spatial actuity is deactivated in this version and will be implemented in
% a later release. It will be implemented in ImageJ as per JT. A
% preliminary fix can be activated using the buttons in the code.

% The magnetic lasso, the function interp and clipper (copyright (c)2015-17, Prof. Erik A. Johnson <JohnsonE@usc.edu>, 01/28/17) are sourced from
% Mathworks forums and need to be credited.

% The paintbrush function has been adapted from Dr M.J.How
% The pattern analysis has been adapted from J.A. Endler's code
% Image manipulation is based on a modified compilation of segments
% initially designed by J.A. Endler.

% The concept and compilation of QCPA is my own intellectual property has
% been started in Feb 2015.

% This version produces fairly basic/standard parameters. Later releases
% will contain more sophisticated image statistics and parameters. (See
% QCPA in MICA)

%Rotation and resizing will be done in ImageJ prior to JND clustering in
%future releases. For now it is still done in Matlab.

%The LiveWire magentic lasso tool will be replaced by a more suitable tool
%in future releases.

% Animal and Background only analysis do not have all the features the
% animal vs. background analysis has.

% This code will not run smoothly on 4k screens. In case of issues see 4k
% version.

%%%%%%%%%%%%%%!!!!!Every image needs to contain a size standard!!!!%%%%%%%%%%%

%Files needed in same directory as QCPA code:

%livewire.m               %Livewire is the magnetic lasso
%fLiveWireCalcP.mexw64
%please_mex_me.m
%fLiveWireGetCostFcn.m
%fLiveWireCalcP.cpp
%fLiveWireCalcP.m
%fLiveWireGetPath.m
%clipper.m                %Clipper does the polygon security margin
%clipper.mexw64
%clipper.hpp

%%%%%%%%%%%%%%%%Set Parameters%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

target=5;       %target scale in pixels/mm (For correct choice see Endler 2012)
pixsp=1;        %number of pixels between transects (For correct choice see Endler 2012)
acd=1.75;       %Spatial acuity of viewer in cycles/degree (For triggerfish see Champ et al 2014)
dist=200;       %distance between viewer and target in mm
blr=0;          %main switch for visual acuity blurring 1=yes 0=no
pure=0;         %number of times a 1x1 pixel median filter is applied 
                %to reduce salt and pepper noise in zone map (reduces
                %impact of fringing and other relicts of the clustering),
                %switch to 0 to turn off
cutoff=5;      %size threshold of clusters that are removed and interpolated from neighbours

digits(4);

%%%%running of required modules%%%%%%%

%please_mex_me %C++ compiler required for live wire plugin, make sure it's in the work directory

%If not installed:
%You can install the freely available MinGW-w64 C/C++ compiler; see Install MinGW-w64
%Compiler. For more options, see  http://www.mathworks.com/support/compilers/R2015b/win64.html.

% While preferabe, the mex file doesn't need to work, the magnetic lasso
% works fine without it. Just a tad slower.

%%%%%%%calculation of kernel size for visual acuity blurring%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
adc=1/acd; %transforming cycles/degree to degrees/cycle
a=2*(tan(((adc/2*((pi/180))))*dist)); %a is kernel size in mm
a=round(a*target); % a is kernel size in pixels
if mod(a,2)==0 
    a=a+1;
    
end;

domain=fspecial('disk',((a/2)-1))>0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

scr=groot; scr=scr.ScreenSize; %reads properties of screen and sets dispay mode accordingly
scrsz = get(groot,'ScreenSize');

scr(3)=scr(3)-40;
scr(4)=scr(4)-80;

scrsz(3)=scrsz(3)-50;
scrsz(4)=scrsz(4)-50;


%%%%%%Select colour png%%%%%%%%%%%%

[fname,dr]=uigetfile('*.png','Select a CLUSTERED .png image corresponding to the zone map');%Opens ui command to select zone map
fn=fname;fp1=fullfile(dr,fn);

zones_colour=imread(fp1);

sz=size(zones_colour);
nc=sz(1);
nr=sz(2);

%%%%%select zone map%%%%%%%%%%%%%%

[fname,dr]=uigetfile('*.txt','Select a zone map');%Opens ui command to select zone map
fn=fname;fp=fullfile(dr,fn);zones=dlmread(fp);

zones_rgb=uint8(zones); %grey scale image of zone map

%%%%%select zone map ID file%%%%%%

[fname,dr]=uigetfile('*.xls','Select a zone ID file');%Opens ui command to select zone map
fn=fname;fp=fullfile(dr,fn);zones_ID=tdfread(fp);zones_ID=struct2table(zones_ID);

%%%%%rename zones in both files%%%%%%

pass=max(zones_ID.('Pass')); %maximum number of passes in zone ID file
zones_ID=zones_ID(zones_ID.Pass==pass,:); %discard unwanted pass data
rownr=size(zones_ID);
rownr=rownr(1,1);
zones_ID.IDnew=[1:rownr]'; %add new zones IDs in ID file

for i=1:rownr
zones(zones==zones_ID.ID(i))=zones_ID.IDnew(i);
end

%%%%%%%%%%size match check between input files%%%%%%%%%%%%%

sz=size(zones_colour);
    szz=size(zones);
    nc=sz(1);
    ncc=szz(1);
    nr=sz(2);
    nrr=szz(2);
    
    if ncc~=nc; close all; fprintf(1,'image size mismatch zone map vs. .png'); return;  end; 



%%%%%%%%%%%%%%%%%%%%%%removal of salt & pepper noise%%%%%%%%%%%%%%%
%%%%%%%%%% This removes zones with a size of 1 or 2 pixels%%%%%%%%%

if pure~=0 
for i=1:pure
        zones=medfilt2(zones,[3 3],'symmetric');
        round(zones);
        
        zones_colour(:,:,1)=medfilt2(zones_colour(:,:,1),[3 3],'symmetric');
        zones_colour(:,:,2)=medfilt2(zones_colour(:,:,2),[3 3],'symmetric');
        zones_colour(:,:,3)=medfilt2(zones_colour(:,:,3),[3 3],'symmetric');
        round(zones_colour(:,:,1));
        round(zones_colour(:,:,2));
        round(zones_colour(:,:,3));       
end;
end;

imshow(zones_colour);%Opens the image and displays it as figure 1
set(gcf,'Position',scr);

%%%%%%%%%%Click ruler and scale image%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    title('Click on two points exactly 1cm apart along size standard','FontSize',20);%adds title to figure 1
    [x,y] = ginput(2); hold on; plot(x,y,'b','LineWidth',4); title(fn); drawnow;%draws the line
    pp3=sqrt((x(1)-x(2))^2 + (y(1)-y(2))^2); %calculates the distance between the two points
    pp3=pp3/10; %pixels per mm
    
    if target>=pp3; close all; fprintf(1,'********Original resolution too low for target scale***********'); return;  end;
    
    scale=target/pp3; fprintf(1,'Scaling factor %5.3f\n',scale); clf;%calculates rescaling factor
    if scale~=1 zones_colour=imresize(zones_colour,scale,'nearest'); end;  %rescale image to common pixel density until ~1
    if scale~=1 zones=imresize(zones, scale,'nearest');end; %rescale zone map to common pixel density until ~1
    
    sz=size(zones_colour);
    szz=size(zones);
    nc=sz(1);
    ncc=szz(1);
    nr=sz(2);
    nrr=szz(2);
    
    if ncc~=nc; close all; fprintf(1,'************size mismatch zone map vs. .png***************'); return;  end; 
    
    imshow(zones_colour); set(gcf,'Position',scr); hold on; drawnow;

%%%%Re-orientate the nudibranch/Background%%%%%%%%%%%%%%%%%%%

    title('Click in front of head then click-drag for object''s long axis','FontSize',20); 
    [x1,y1]=ginput(1); [x2,y2,~]=DrawFromXY(1,x1,y1); plot([x1 x2],[y1 y2],'w--'); 
    s1=x2-x1;
    h=y1-y2 ; 
    ang=atand(h/s1);
    
    if x1<x2 &&  y1<y2 
         
         ang=-(90-abs(ang));       
    end
    
   if x1>x2 && y1<y2
          ang=90-ang;
   end
    
   if x1<x2 && y1>y2
         ang=-(180-(90-ang));
   end  
    
   if x1>x2 &&  y1>y2 && h<s1
         ang=-(180-ang);
   end   
    
   if x1>x2 &&  y1>y2 && h>s1
         ang=-(90-ang);
    end  
   
    %makes sure the animal is heads up
    
    zones_colour1=zones_colour;
    zones_colour=rotateAround(zones_colour,nc/2,nr/2,ang,'nearest'); %rotates the image according to "ang"
    zones_colour_backup=rotateAround(zones_colour1,nc/2,nr/2,ang,'nearest'); %rotates the image according to "ang"
    
    zones1=zones;
    
    zones=rotateAround(zones,nc/2,nr/2,ang,'nearest'); %rotates the zone map according to "ang"
    
    zones_backup=rotateAround(zones1,nc/2,nr/2,ang,'nearest'); %rotates the zone map according to "ang"
    
    sz=size(zones_colour);
    nc=sz(1);
    nr=sz(2);
    
    
    close all;
    
%%%%%%%%%%%remove small zones%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
if cutoff~=0
    
too_small=find(any(zones_ID{:,12}<cutoff,2));
s=size(too_small);

for q=1:s
    culprit=too_small(q);
for i=1:nc
    for j=1:nr
        if zones(i,j)==culprit
            zones(i,j)=NaN;
        end
    end
end
end


end;

zones=n_interp(zones); %replaces all NaN values (if any were left after salt&pepper removal) with closest non-0 value
    
%%%%%%%%%%%%%%%Blur image and zone map%%%%%%%%%%%%%%%%%%%%%%

if blr==1
    
fprintf(1,'Applying spatial acuity to image and zone map, this can take a while\n')   

zones=ordfilt2(zones,round(0.5*numel(find(domain))),domain,'symmetric');

zones_colour(:,:,1)=ordfilt2(zones_colour(:,:,1),round(0.5*numel(find(domain))),domain,'symmetric');
zones_colour(:,:,2)=ordfilt2(zones_colour(:,:,2),round(0.5*numel(find(domain))),domain,'symmetric');
zones_colour(:,:,3)=ordfilt2(zones_colour(:,:,3),round(0.5*numel(find(domain))),domain,'symmetric');

zones_backup=round(zones);
zones=round(zones);

zones_colour_backup(:,:,1)=round(zones_colour(:,:,1));
zones_colour_backup(:,:,2)=round(zones_colour(:,:,2));
zones_colour_backup(:,:,3)=round(zones_colour(:,:,3));
zones_colour(:,:,1)=round(zones_colour(:,:,1));
zones_colour(:,:,2)=round(zones_colour(:,:,2));
zones_colour(:,:,3)=round(zones_colour(:,:,3));

end;

%%%%%%%%%%%%%%%%%%%remove salt&pepper noise after blurring%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if blr==1
    
for i=1:pure
        zones=medfilt2(zones,[3 3],'symmetric');
        round(zones);
        
        zones_colour(:,:,1)=medfilt2(zones_colour(:,:,1),[3 3],'symmetric');
        zones_colour(:,:,2)=medfilt2(zones_colour(:,:,2),[3 3],'symmetric');
        zones_colour(:,:,3)=medfilt2(zones_colour(:,:,3),[3 3],'symmetric');
        round(zones_colour(:,:,1));
        round(zones_colour(:,:,2));
        round(zones_colour(:,:,3));
        
end
end

        
%%%%%%%%%%%%%%%%%%%%%%%IMAGE SEGMENTATION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%identify object zones in zone map%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%% Removal of colour standard%%%%%%%%%%%%%%%%%%%%%%%%%%%%

agn=1; %Main switch for image zone exclusion

while agn==1
 
imshow(zones_colour);   %Opens the clustered RGB image and displays it as figure 1
set(gcf,'Position',scr);

ch=questdlg('Do you want to exclude a colour standard?','--','Yes','No','Yes');  %last one is default
switch ch
    case 'Yes'
        agn=1;
        alt=0;
        close all;
    case 'No'
        agn=0;
        alt=1;
        
        close all;
end;

bt=0;
if alt==1
    bt=1;
end;
while bt==0

figure(1); imshow(zones_colour);set(figure(1),'Position',[1 1 scrsz(3)/2 scrsz(4)]);title({'PNG FILE';'Use for comparison, draw on right image -->'},'Color','k');
figure(2);imshow(zones,[]);set(figure(2),'Position',[scrsz(3)/2 1 scrsz(3)/2 scrsz(4)]);title({'ZONE MAP (Actual Data)';'draw polygon, dbl click to finish'});

hold on;

crop1=roipoly(mat2gray(zones));
zones(crop1)=0;
crop1=repmat(crop1,[1,1,3]); %multiply mask to all RGB channels
zones_colour(crop1)=0;


close all;

figure(1);imshow(zones_colour);%Display cropped image
set(figure(1),'Position',[1 1 scrsz(3)/2 scrsz(4)]);title('Zones Colour','Color','k')
figure(2);imshow(zones,[]);set(figure(2),'Position',[scrsz(3)/2 1 scrsz(3)/2 scrsz(4)]);title('Zones','Color','r');

ch=questdlg('Is this correct?','--',...
			'Yes','Do again','Modify','Yes');  %last one is default
		switch ch
			case 'Yes'
				bt=1; %switches bt to 1
                agn=0;
				
                close all;
			case 'Do again'
				bt=0; %switches bt to 0
                
                zones=zones_backup; %restores original data
                zones_colour=zones_colour_backup; %restores original data
                
			case 'Modify'
                close all;
                
				zones1 = modifyzonemap2(zones); %starts modify function
                zones1=zones1(1:nc,1:nr);
                mask=zones1==0;
                mask=repmat(mask,[1,1,3]);
                zones_colour1=zones_colour;
                zones_colour1(mask)=0;          %create zones_colour1 using
                                                %mask from modify
                
                figure(1);imshow(zones1(1:nc,1:nr),[]);
                set(figure(1),'Position',[1 1 scrsz(3)/2 scrsz(4)]);title('Zones','Color','k')
                
                figure(2),imshow(zones_colour1);
                set(figure(2),'Position',[scrsz(3)/2 1 scrsz(3)/2 scrsz(4)]);title('Colour PNG','Color','r');
                
                
                nb=1;
                
                while nb==1
                    
				ch3=questdlg('Are you happy with the exclusion of the colour standard?','--',...
					'Yes','No','Yes');  %last one is default
				switch ch3
					case 'Yes'
						nb=0; 
                        agn=0;
                        bt=1;
                        
                        zones=zones1(1:nc,1:nr,:);
                        zones_colour=zones_colour1;
                        close all;

                    case 'No'
						nb=0;
                        bt=0;
                        
                    zones=zones_backup;
                    zones_colour=zones_colour_backup; 
                    
                end;
                end;
        end;
end;
end;
zones2=zones; 
close all;
cln=1;

%%%%%%%%%%%%%%%% Removal of other parts of the image %%%%%%%%%%%%%%%%%%%%%%%%%

while cln==1
imshow(zones2,[]);set(gcf,'Position',scr);

ch=questdlg('Do you Want to exclude other parts of the image?','--','Yes','No','Yes');  %last one is default

switch ch
    case 'Yes'
        cln=1;
        close all;
        imshow(zones2,[]);set(gcf,'Position',scr);
        title('Outline area to exlude. Press Del to reset anchor. Enter when done. Program will stop without input');
        hold on;
        zones_colour2=zones_colour;
        crop2=livewire(zones_colour2);
        crop2=repmat(crop2,[1,1,3]); %multiply mask to all RGB channels
        zones_colour2(crop2)=0;
        imshow(zones_colour2);set(gcf,'Position',scr); %Display cropped image

        ch=questdlg('Is this correct?','--',...
			'Yes','Do again','Modify','Yes');  %last one is default
       nb=0;
		switch ch
			case 'Yes'
				bt=1; %switches bt to 1
                cln=0;
				zones_colour=zones_colour2;
                close all;
			case 'Do again'
				bt=0; %switches bt to 0
                cln=1;
			case 'Modify'
                close all;
				zones_colour2 = modifyzonemap(zones_colour2); %starts modify function
               
                imshow(zones_colour2(1:nc,1:nr,:));set(gcf,'Position',scr);
                
                nb=1;
                
                while nb==1
                    
				ch3=questdlg('Are you happy with the exclusion of the selected area?','--',...
					'Yes','No','Yes');  %last one is default
				switch ch3
					case 'Yes'
						nb=0; %switches bt to 1
                        agn=0;
                        bt=1;
                        
                        zones_colour=zones_colour2(1:nc,1:nr,:);
                        close all;

                    case 'No'
						nb=0; %switches nb to 0
                end;
                end;
        end;
    case 'No'
        
        cln=0;
end;
end;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%select animal outline (If present)%%%%%%%%%%%%%%%%%%

imshow(zones_colour);set(gcf,'Position',scr);

abt1=0; %main switch

while abt1==0
    ch=questdlg('Do you want to outline an animal in the image?','--',...
					'Yes','No','Yes');  %last one is default
    switch ch
        case 'Yes'
         anp=1;
         title('SELECT RECTANGLE AROUND ANIMAL','Color','r');
         hold on;
         rect = getrect(1);rect=round(rect(1:4));
         ratio1=nr/rect(3);
         ratio2=nc/rect(4);
         
         if ratio2<ratio1
             zoom_factor=ratio2;
         else
             zoom_factor=ratio1;
         end;
         
         center_x=rect(1)+(rect(3)/2);
         center_y=rect(2)+(rect(4)/2);
         zoomcenter(center_x,center_y,zoom_factor);
         
         clr=figure(1); set(clr,'Position',[1 1 scrsz(3)/2 scrsz(4)]);title('Use for comparison, draw on right image -->','Color','k');
         crp=figure(2);imshow(zones);zoomcenter(center_x,center_y,zoom_factor);
         title('Outline area to exclude, be conservative! Press CTRL to turn off auto fit / Press Enter to finish', 'Color','r');
         set(crp,'Position',[scrsz(3)/2 1 scrsz(3)/2 scrsz(4)]);hold on;
         [crop, xi, yi]=livewire(zones);
         coords=horzcat(xi(1:8:end)',yi(1:8:end)');
         
        
         close all;
         imshow(zones,'InitialMagnification','fit','DisplayRange',[]);set(gcf,'Position',scr);
         title({'drag & drop until border is defined. A for new vertex';'Place border INSIDE animal. DBL Click on polygon line to finish'});
         hold on;
         fprintf(1,'Translating object outline into polygon, this can take a while\n')
         h=impoly(gca,coords);zoomcenter(center_x,center_y,zoom_factor);
         api = iptgetapi(h);
         api.setColor('red');
         accepted_pos = wait(h);
         nudi_area_pixel=polyarea(accepted_pos(:,1),accepted_pos(:,2));
         abt1=1; %set switch
         close all;
         
         
        case 'No'
            close all;
            abt1=1;
            anp=0;
     end


%%%%%%%%%%%%%%%%%%%%%%%%%Pattern Analysis%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 ch=questdlg('What kind of pattern analysis do you want to make?','--',...
					'Animal only','Background only','Animal-Background','Animal only');  %last one is default
    switch ch
        case 'Animal only'
%%%%%Animal only%%%%%%%%%%%%

animal_mask=poly2mask(accepted_pos(:,1),accepted_pos(:,2),nc,nr);%turn animal polygon into mask

%apply mask to zone map and image

zones_animal=zones;
zones_animal(~animal_mask)=0;
animal_mask=repmat(animal_mask,[1,1,3]); %multiply mask to all RGB channels
zones_animal_colour=zones_colour;
zones_animal_colour(~animal_mask)=0;

%%%display image to control and modify if necessary%%%%%%%%%%%%%%



imshow(zones_animal_colour);zoomcenter(center_x,center_y,zoom_factor);
clr=figure(1); set(clr,'Position',[1 1 scrsz(3)/2 scrsz(4)]);title('zones colour animal','Color','k');

crp=figure(2);imshow(zones_animal,[]);zoomcenter(center_x,center_y,zoom_factor);
title('zones animal', 'Color','r');
set(crp,'Position',[scrsz(3)/2 1 scrsz(3)/2 scrsz(4)]);

ch=questdlg('Need to modify image?','--',...
				'Yes','No','Yes');  %last one is default
 ao1=0; %switch
 
 while ao1==0
 switch ch
        case 'Yes'
        zones_animal1=modifyzonemap2(zones_animal);
        zones_animal1=zones_animal1(1:nc,1:nr);
                mask=zones_animal1==0;
                mask=repmat(mask,[1,1,3]);
                zones_animal_colour1=zones_animal_colour;
                zones_animal_colour1(mask)=0; 
                

clr=figure(1); imshow(zones_animal_colour1);zoomcenter(center_x,center_y,zoom_factor);
title('zones animal colour', 'Color','k');
set(clr,'Position',[1 1 scrsz(3)/2 scrsz(4)]);

crp=figure(2);imshow(zones_animal1,[]);zoomcenter(center_x,center_y,zoom_factor);
title('zones animal', 'Color','r');
set(crp,'Position',[scrsz(3)/2 1 scrsz(3)/2 scrsz(4)]);
        
        ch=questdlg('need to redo modifying?','--',...
					'Yes','No','Yes');  %last one is default
        switch ch
                case 'Yes'
                    
                case 'No'
                zones_animal=zones_animal1;
                zones_animal_colour=zones_animal_colour1;
                ao1=1;
        
            
        end
     case 'No'
     ao1=1;
     close all;
 end;
 end %Animal only zone map is accepted
 
 close all;
 
%%%%%%%%%%%%re-arrange zone IDs%%%%%%%%%%%

%This is necessary as the numbering of zones must be continuous without gaps for the Adjacency code to work

old_zone_ID=unique(zones_animal);
n_zones_animal=size(old_zone_ID);

for q=2:n_zones_animal(1)
    for i=1:nc
        for j=1:nr
            if zones_animal(i,j)==old_zone_ID(q)
                zones_animal(i,j)=q-1; %replaces zone value
            end
        end
    end
end

for i=1:nc
    for j=1:nr
if zones_animal(i,j)==0
    zones_animal(i,j)=n_zones_animal(1); %Makes the zeros the highest zone number (last one), requirement of doAdjacency
end
    end
end
 
%%%%%%%%%%%%Re-arrange zone names in zone_ID file (for pattern stats)%%%%%%

%New list of zones in zone map%

new_zone_ID=unique(zones_animal);

%Identify rows with the zones that are used

zones_ID=zones_ID(ismember(zones_ID{:,14},old_zone_ID),:);%creates subset and overwrites old (gets rid of unused ones)

for i=1:(n_zones_animal(1)-1)
    zones_ID{i,14}=i;
end


zlist=unique(zones_animal); nz=length(zlist);

nz=nz-1;

figure(1);imshow(zones_animal,[]);zoomcenter(center_x,center_y,zoom_factor); %Image which is displayed at the end
if blr==1; title(['blurred zone map, k=' num2str(nz) ,' Cycles/degree =' num2str(acd) ,' Distance(mm)=' num2str(dist);]);end;
if blr==0; title(['unblurred zone map, k=' num2str(nz)]);end;

figure(2);imshow(zones_animal_colour);zoomcenter(center_x,center_y,zoom_factor); %Image which is displayed at the end
if blr==1; title(['blurred colour image, k=' num2str(nz) ,' Cycles/degree =' num2str(acd) ,' Distance(mm)=' num2str(dist);]);end;
if blr==0; title(['unblurred colour image, k=' num2str(nz)]);end;

%%%%%%%%get pattern geometry stats%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%Do Adjaceny for Animal zone map%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~,cdiv,~,tdiv,cmplx,pxud,pxlr,asptr,inout,seq]=DoAdjacency(zones_animal,pixsp); 
adjacency=[cdiv tdiv cmplx pxud pxlr asptr inout nz];


%%%%%% get visual contrast stats%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



zf=zeros(nz,1);

for i=1:nz
zf(i)=sum(sum(zones_animal==zlist(i)));
end;

zf=zf/sum(zf);  %these are the fractions (relative areas) of each zone in the animal

%%%get hue, chroma and luminance of zones%%%%

figure;
hold on;
title('Maxwell Triangle');
hue=zeros(nz,1); chr=hue; lum=hue;

for cl=1:nz
cc=[zones_ID.swMean(zones_ID.IDnew==zlist(cl)) zones_ID.mwMean(zones_ID.IDnew==zlist(cl)) zones_ID.lwMean(zones_ID.IDnew==zlist(cl))];
lum(cl)=zones_ID.dblMean(zones_ID.IDnew==zlist(cl));

gy=[1/3 1/3 1/3]; [ox,oy]=ToTriangle(gy);
[tx,ty]=ToTriangle(cc); plot(tx,ty,'x','LineWidth',2);
legend('labels',{num2str(zlist)},'Interpreter','none','FontSize',5,'FontWeight','bold');
legend('boxoff');


tx=tx-ox; ty=ty-oy;
[hua,chrv]=cart2pol(tx,ty); 
if chrv<0.05 hua=NaN; end; %hue angle for chr>0.05
if hua<0 hua=hua+pi; end; %convert -180 to 180 to 0 to 360 in radians
hua=hua/(2*pi); % convert 0 to 360 into 0 to 1
     hue(cl)=hua; chr(cl)=chrv;
     clear rfl cc ccd gy ox oy tx ty hua crv;   
end;

DrawTri;
hold off;
%take weighted means and sd of hue chr lum, weighted by zone fractions

[mnhu,sd]=WeightedMnSD(hue,zf); cvhu=sd/mnhu;
[mnch,sd]=WeightedMnSD(chr,zf); cvch=sd/mnch;
[mnlm,sd]=WeightedMnSD(lum,zf); cvlm=sd/mnlm;

contrasts=[mnhu mnch mnlm cvhu cvch cvlm];

%%%%%%Safe results%%%%%%%%%

ohd='MnHue,MnChr,MnLum,CVhue,CVChr,CVLum,CDiv,Tdiv,Cmplx,PxUD,PxLR,asptr,k';


if blr==1
    
    oname1=[num2str(dist) '_mm' '_PatternStatistics_animal_only_blurred.csv']; fid1=fopen(oname1,'wt'); %Gives name to files
    oname2=[num2str(dist) '_mm' '_transitionmatrix_animal_only_blurred.csv']; fid2=fopen(oname2,'wt');
    oname3=[num2str(dist) '_mm' '_colour_pattern_details_animal_only_blurred.csv']; fid3=fopen(oname3,'wt');
    
    saveas(figure(1),[num2str(dist) '_mm' '_animal_blurred_zones.png']);  %saves zone map as .png file
    [num2str(dist) '_mm' 'animal_blurred__zones.png'];
    
    saveas(figure(2),[num2str(dist) 'mm' '_animal_blurred_zones_maxwell.png']);  %saves png as .png file
    [num2str(dist) '_mm' 'animal_blurred__colour.png'];
    
    saveas(figure(3),[num2str(dist) 'mm' '_animal_blurred_zones_maxwell.png']);  %saves zone map as .png file
    [num2str(dist) '_mm' 'animal_blurred__zones_maxwell.png'];
    
    save([num2str(dist) 'mm' '_animal_blurred_zones.mat'],'zones_animal'); %saves zone map as matrix
    [num2str(dist) '_mm''animal_blurred_zones.mat'];
    
end;
   
if blr==0
    
     oname1='PatternStatistics_animal_only_unblurred.csv'; fid1=fopen(oname1,'wt'); %gives name to files
     oname2='Transitionmatrix_animal_only_unblurred.csv'; fid2=fopen(oname2,'wt');
     oname3='Colour_pattern_element_details_animal_only_unblurred.csv'; fid3=fopen(oname3,'wt');
     
    saveas(figure(1),['animal_unblurred_zones.png']);  %saves zone map as .png file
    ['animal_unblurred__zones.png'];
    
    saveas(figure(2),['animal_unblurred_zones_maxwell.png']);  %saves png as .png file
    ['animal_unblurred__colour.png'];
    
    saveas(figure(3),['animal_unblurred_zones_maxwell.png']);  %saves zone map as .png file
    ['animal_unblurred__zones_maxwell.png'];
    
    save(['animal_unblurred_zones.mat'],'zones_animal'); %saves zone map as matrix
    ['animal_unblurred_zones.mat'];
    
end;

fprintf(fid1,'%s\n',ohd);


for j=1 fprintf(fid1,'%6.4f',contrasts(j)); end; %selects the parameters from the adjacency and contrasts output 
for j=2:6 fprintf(fid1,',%6.4f',contrasts(j)); end;
for j=1:3 fprintf(fid1,',%6.4f',adjacency(j)); end;
for j=4:5 fprintf(fid1,',%6.2f',adjacency(j)); end;
for j=6 fprintf(fid1,',%6.4f',adjacency(j)); end;
for j=8 fprintf(fid1,',%6f',adjacency(j));end;
     
     
fprintf(fid1,'\n\n\n\n\n'); 
fprintf(fid1,'Meaning of columns\n');
fprintf(fid1,'MnHue,mean Hue angle (0 to 1=360degrees) weighted by patch relative area ignoring low chroma\n');
fprintf(fid1,'MnChr,mean Chroma, weighted by patch relative areas\n');
fprintf(fid1,'MnLum,mean Luminance, weighted by patch relative areas\n');
fprintf(fid1,'CVhue,CV (coefficient of variation) of Hue angle, weighted as above\n');
fprintf(fid1,'CVChr,CV of Chroma\n');
fprintf(fid1,'CVLum,CV of luminance\n');
fprintf(fid1,'CDiv,colour diversity\n');
fprintf(fid1,'Tdiv,transition diversity (diversity of off diagonals)\n');
fprintf(fid1,'Cmplx,transitions/sample (animal) a measure of pattern complexity\n');
fprintf(fid1,'PxUD,distance between transitions (pixels) up and down (vertical axis)\n');
fprintf(fid1,'PxLR,distance between transitions left and right (horizontal axis)\n');
fprintf(fid1,'asptr,transition aspect ratio (larger if longer in vertical axis)\n');

fclose(fid1);

fprintf(1,'\n');
fprintf(1,'Pattern stats in %s\n',oname1);
fprintf(1,'\n');
fprintf(1,'Transition Matrix in %s\n',oname2);
fprintf(1,'\n');
fprintf(1,'Pattern element details in %s\n',oname3);

%%%%print transition matrix of animal only%%%%%%%%%%%%%%%%%%

fprintf(fid2, 'Adjacency_Transition_Matrix_Animal_only\n\n');
fprintf(fid2,' ');
for j=1:nz fprintf(fid2,'\t%15.0f',zlist(j)); end; fprintf(fid2,'\n');
     for k=1:nz
       fprintf(fid2,'%15.0f',zlist(k));
       for j=1:nz fprintf(fid2,'\t%15.0f',seq(k,j)); end; fprintf(fid2,'\n');
     end;
     fprintf(fid2,'\n');
     
%%%print pattern element details of animal only%%%%

fprintf(fid3,'ZoneID');
fprintf(fid3,',HueAngle');
fprintf(fid3,',Chroma');
fprintf(fid3,',Luminance');
fprintf(fid3,',SWccq');
fprintf(fid3,',MWccq');
fprintf(fid3,',LWccq');
fprintf(fid3,',DBLccq');
fprintf(fid3,',RelativeSize');
fprintf(fid3,',TotalSize\n');

for i=1:nz
    fprintf(fid3,'%6f',i);
    fprintf(fid3,',%6.4f',hue(i));
    fprintf(fid3,',%6.4f',chr(i));
    fprintf(fid3,',%6.4f',lum(i));
    fprintf(fid3,',%6.4f',zones_ID{i,4});
    fprintf(fid3,',%6.4f',zones_ID{i,6});
    fprintf(fid3,',%6.4f',zones_ID{i,8});
    fprintf(fid3,',%6.4f',zones_ID{i,10});
    fprintf(fid3,',%6.4f',zf(i));
    fprintf(fid3,',%6.4f',zones_ID{i,12});
    fprintf(fid3,'\n');
end
   
%add other parameters here such as individual colour pattern complexity
%(number of transitions to all other patches/area)


%writetable(zones_ID,oname3,'WriteVariableNames',true,'Delimiter',',');
%%for MICA output per zones

%%%%%%%%%%%%%%%%%% finished Animal only %%%%%%%%%%%%%%

case 'Background only'
    
%%%%%Background only%%%%%%%%%%%%

if anp==1
background_mask=poly2mask(accepted_pos(:,1),accepted_pos(:,2),nc,nr);%turn animal polygon into mask

%apply mask to zone map and image

zones_background=zones;
zones_background(background_mask)=0;

background_mask=repmat(background_mask,[1,1,3]); %multiply mask to all RGB channels

zones_background_colour=zones_colour;
zones_background_colour(background_mask)=0;
end;

if anp==0
    zones_background_colour=zones_colour;
    zones_background=zones;
end;
%%%display image to control and modify if necessary%%%%%%%%%%%%%%

imshow(zones_background_colour);
clr=figure(1); set(clr,'Position',[1 1 scrsz(3)/2 scrsz(4)]);title('zones colour animal','Color','k');

crp=figure(2);imshow(zones_background,[]);
title('zones animal', 'Color','r');
set(crp,'Position',[scrsz(3)/2 1 scrsz(3)/2 scrsz(4)]);


ch=questdlg('Need to modify image?','--',...
				'Yes','No','Yes');  %last one is default
 ao1=0; %switch
 
 while ao1==0
 switch ch
        case 'Yes'
        zones_background1=modifyzonemap2(zones_background);
        zones_background1=zones_background1(1:nc,1:nr);
                mask=zones_background1==0;
                mask=repmat(mask,[1,1,3]);
                zones_background_colour1=zones_background_colour;
                zones_background_colour1(mask)=0; 
                

clr=figure(1); imshow(zones_background_colour1);
title('zones animal colour', 'Color','k');
set(clr,'Position',[1 1 scrsz(3)/2 scrsz(4)]);

crp=figure(2);imshow(zones_background1,[]);
title('zones animal', 'Color','r');
set(crp,'Position',[scrsz(3)/2 1 scrsz(3)/2 scrsz(4)]);
        
        
        ch=questdlg('need to redo modifying?','--',...
					'Yes','No','Yes');  %last one is default
        switch ch
                case 'Yes'
                    
                case 'No'
                    
                zones_background_colour=zones_background_colour1;
                zones_background=zones_background1;
                
                ao1=1;
        
            
        end
     case 'No'
     ao1=1;
     close all;
 end;
 end %Background only zone map is accepted
close all;

%%%%%%%%%%%%re-arrange zone IDs%%%%%%%%%%%

%This is necessary as the numbering of zones must be continuous without gaps for the Adjacency code to work

old_zone_ID=unique(zones_background);
n_zones_background=size(old_zone_ID);

for q=2:n_zones_background(1)
    for i=1:nc
        for j=1:nr
            if zones_background(i,j)==old_zone_ID(q)
                zones_background(i,j)=q-1; %replaces zone value
            end
        end
    end
end

for i=1:nc
    for j=1:nr
if zones_background(i,j)==0
    zones_background(i,j)=n_zones_background(1); %Makes the zeros the highest zone number (last one), requirement of doAdjacency
end
    end
end
 
%%%%%%%%%%%%Re-arrange zone names in zone_ID file (for pattern stats)%%%%%%

%New list of zones in zone map%

new_zone_ID=unique(zones_background);

%Identify rows with the zones that are used

zones_ID=zones_ID(ismember(zones_ID{:,14},old_zone_ID),:);%creates subset and overwrites old (gets rid of unused ones)

for i=1:(n_zones_background(1)-1)
    zones_ID{i,14}=i;
end


zlist=unique(zones_background); nz=length(zlist);

nz=nz-1;


figure(1);imshow(zones_background,[]); %Image which is displayed at the end
if blr==1; title(['blurred zone map, k=' num2str(nz) ,' Cycles/degree =' num2str(acd) ,' Distance(mm)=' num2str(dist);]);end;
if blr==0; title(['unblurred zone map, k=' num2str(nz)]);end;

figure(2);imshow(zones_background_colour); %Image which is displayed at the end
if blr==1; title(['blurred colour image, k=' num2str(nz) ,' Cycles/degree =' num2str(acd) ,' Distance(mm)=' num2str(dist);]);end;
if blr==0; title(['unblurred colour image, k=' num2str(nz)]);end;

%%%%%%%%get pattern geometry stats%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~,cdiv,~,tdiv,cmplx,pxud,pxlr,asptr,inout,seq]=DoAdjacency(zones_background,pixsp); 
adjacency=[cdiv tdiv cmplx pxud pxlr asptr inout nz];

%%%%%% get visual contrast stats%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

zf=zeros(nz,1);

for i=1:nz
zf(i)=sum(sum(zones_background==zlist(i)));
end;

zf=zf/sum(zf);  %these are the fractions (relative areas) of each zone in the animal

%%%get hue, chroma and luminance of zones%%%%

figure;hold on;
title('Maxwell Triangle');

hue=zeros(nz,1); chr=hue; lum=hue;

for cl=1:nz
    
cc=[zones_ID.swMean(zones_ID.IDnew==zlist(cl)) zones_ID.mwMean(zones_ID.IDnew==zlist(cl)) zones_ID.lwMean(zones_ID.IDnew==zlist(cl))];
lum(cl)=zones_ID.dblMean(zones_ID.IDnew==zlist(cl));

gy=[1/3 1/3 1/3]; [ox,oy]=ToTriangle(gy);
[tx,ty]=ToTriangle(cc); plot(tx,ty,'x','LineWidth',2)

legend('labels',{num2str(zlist)},'Interpreter','none','FontSize',5,'FontWeight','bold');
legend('boxoff');

tx=tx-ox; ty=ty-oy;
[hua,chrv]=cart2pol(tx,ty); 
if chrv<0.05 hua=NaN; end; %hue angle for chr>0.05
if hua<0 hua=hua+pi; end; %convert -180 to 180 to 0 to 360 in radians
hua=hua/(2*pi); % convert 0 to 360 into 0 to 1
     hue(cl)=hua; chr(cl)=chrv;
     clear rfl cc ccd gy ox oy tx ty hua crv;   
end;


DrawTri;hold off;

%take weighted means and sd of hue chr lum, weighted by zone fractions

[mnhu,sd]=WeightedMnSD(hue,zf); cvhu=sd/mnhu;
[mnch,sd]=WeightedMnSD(chr,zf); cvch=sd/mnch;
[mnlm,sd]=WeightedMnSD(lum,zf); cvlm=sd/mnlm;

contrasts=[mnhu mnch mnlm cvhu cvch cvlm];

%%%%%%Safe results%%%%%%%%%

ohd='MnHue,MnChr,MnLum,CVhue,CVChr,CVLum,CDiv,Tdiv,Cmplx,PxUD,PxLR,asptr,k';


if blr==1
    oname1=[num2str(dist) '_mm_PatternStatistics_Background_only_blurred.csv']; fid1=fopen(oname1,'wt');
    oname2=[num2str(dist) '_mm' '_transitionmatrix_background_only_blurred.csv']; fid2=fopen(oname2,'wt');
    oname3=[num2str(dist) '_mm' '_colour_pattern_details_background_only_blurred.csv']; fid3=fopen(oname3,'wt');
    
    saveas(figure(1),[num2str(dist) '_mm' '_background_blurred_zones.png']);  %saves zone map as .png file
    [num2str(dist) '_mm' 'background_blurred__zones.png'];
    
    saveas(figure(2),[num2str(dist) 'mm' '_background_blurred_zones_maxwell.png']);  %saves zone map as .png file
    [num2str(dist) '_mm' 'background_blurred__zones_maxwell.png'];
    
    saveas(figure(3),[num2str(dist) 'mm' 'background_blurred_zones_maxwell.png']);  %saves zone map as .png file
    [num2str(dist) '_mm' 'background_unblurred_zones_maxwell.png'];
    
    save([num2str(dist) 'mm' '_background_blurred_zones.mat'],'zones_background'); %saves zone map as matrix
    [num2str(dist) '_mm''background_blurred_zones.mat'];
    
end;
   
if blr==0
     oname1='PatternStatistics_Background_only_unblurred.csv'; fid1=fopen(oname1,'wt');
     oname2='Transitionmatrix_animal_only_unblurred.csv'; fid2=fopen(oname2,'wt');
     oname3='Colour_pattern_element_details_animal_only_unblurred.csv'; fid3=fopen(oname3,'wt');
    
     
    saveas(figure(1),['background_unblurred_zones.png']);  %saves zone map as .png file
    ['background_unblurred__zones.png'];
    
    saveas(figure(2),['background_unblurred_png.png']);  %saves zone map as .png file
    ['background_unblurred__zones_maxwell.png'];
    
    saveas(figure(3),['background_unblurred_zones_maxwell.png']);  %saves zone map as .png file
    ['background_unblurred_zones_maxwell.png'];
    
    save(['background_unblurred_zones.mat'],'zones_background'); %saves zone map as matrix
    ['background_unblurred_zones.mat'];
    
end;

fprintf(fid1,'%s\n',ohd);


for j=1 fprintf(fid1,'%6.4f',contrasts(j)); end;
for j=2:6 fprintf(fid1,',%6.4f',contrasts(j)); end;
for j=1:3 fprintf(fid1,',%6.4f',adjacency(j)); end;
for j=4:5 fprintf(fid1,',%6.2f',adjacency(j)); end;
for j=6 fprintf(fid1,',%6.4f',adjacency(j)); end;
for j=8 fprintf(fid1,',%6f',adjacency(j));end;     
     
fprintf(fid1,'\n\n\n\n\n'); 
fprintf(fid1,'Meaning of columns\n');
fprintf(fid1,'MnHue,mean Hue angle (0 to 1=360degrees) weighted by patch relative area ignoring low chroma\n');
fprintf(fid1,'MnChr,mean Chroma, weighted by patch relative areas\n');
fprintf(fid1,'MnLum,mean Luminance, weighted by patch relative areas\n');
fprintf(fid1,'CVhue,CV (coefficient of variation) of Hue angle, weighted as above\n');
fprintf(fid1,'CVChr,CV of Chroma\n');
fprintf(fid1,'CVLum,CV of luminance\n');
fprintf(fid1,'CDiv,colour diversity\n');
fprintf(fid1,'Tdiv,transition diversity (diversity of off diagonals)\n');
fprintf(fid1,'Cmplx,transitions/sample (animal) a measure of pattern complexity\n');
fprintf(fid1,'PxUD,distance between transitions (pixels) up and down (vertical axis)\n');
fprintf(fid1,'PxLR,distance between transitions left and right (horizontal axis)\n');
fprintf(fid1,'asptr,transition aspect ratio (larger if longer in vertical axis)\n');

fclose(fid1);

fprintf(1,'\n');
fprintf(1,'Pattern stats in %s\n',oname1);
fprintf(1,'\n');
fprintf(1,'Transition Matrix in %s\n',oname2);
fprintf(1,'\n');
fprintf(1,'Pattern element details in %s\n',oname3);

%%%%print transition matrix of animal only%%%%%%%%%%%%%%%%%%

fprintf(fid2, 'Adjacency_Transition_Matrix_background_only\n\n');
fprintf(fid2,' ');
for j=1:nz fprintf(fid2,'\t%15.0f',zlist(j)); end; fprintf(fid2,'\n');
     for k=1:nz
       fprintf(fid2,'%15.0f',zlist(k));
       for j=1:nz fprintf(fid2,'\t%15.0f',seq(k,j)); end; fprintf(fid2,'\n');
     end;
     fprintf(fid2,'\n');
     
%%%print pattern element details of animal only%%%%

fprintf(fid3,'ZoneID');
fprintf(fid3,',HueAngle');
fprintf(fid3,',Chroma');
fprintf(fid3,',Luminance');
fprintf(fid3,',SWccq');
fprintf(fid3,',MWccq');
fprintf(fid3,',LWccq');
fprintf(fid3,',DBLccq');
fprintf(fid3,',RelativeSize');
fprintf(fid3,',TotalSize\n');

for i=1:nz
    fprintf(fid3,'%6f',i);
    fprintf(fid3,',%6.4f',hue(i));
    fprintf(fid3,',%6.4f',chr(i));
    fprintf(fid3,',%6.4f',lum(i));
    fprintf(fid3,',%6.4f',zones_ID{i,4});
    fprintf(fid3,',%6.4f',zones_ID{i,6});
    fprintf(fid3,',%6.4f',zones_ID{i,8});
    fprintf(fid3,',%6.4f',zones_ID{i,10});
    fprintf(fid3,',%6.4f',zf(i));
    fprintf(fid3,',%6.4f',zones_ID{i,12});
    fprintf(fid3,'\n');
end

%%%%%%%%%%%%%%%%%%%Animal background comparison%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%Global comparisons%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%These parameters result from a simple substraction of the animal
%parameters from the background parameters and include:

% AspDiff : Background minus animal aspect ratio
% PxUDDiff : Background minus animal vertical average patch size
% PxLRDiff : Background minus animal horizontal average patch size
% TdivDiff : Background minus animal transition diversity
% CDivDiff : Background minus animal colour diversity
% CVLumDiff : Background minus animal coefficient of variance for luminance
% CVChrDiff : Background minus animal coefficient of variance for chroma
% CVHueDiff : Background minus animal coefficient of variance for hue angle
% MnLumDiff : Background minus animal mean luminance
% MnChrDiff : Background minus animal mean chroma
% MnHueDiff : Background minus animal mean hue angle
%animal rim
%Background rim
%total rim

case 'Animal-Background'

%%%%%%%%%%%Define 6 pixel security layer around animal outline%%%%%%%%%%%%%%%
    
%use please_mex_me2.m (to be written) in order to compile clipper
    
new_nudi1=polyout(accepted_pos(:,1),accepted_pos(:,2),(a/2),'r'); %creates a 3 pixel offset polygon to the outside (enlarged)
new_nudi2=polyout(accepted_pos(:,1),accepted_pos(:,2),-(a/2),'r'); %creates a 3 pixel offset polygon to the inside (shrunk)

x_inside=cell2mat(new_nudi2{1}); %transforms output of polyout
x_outside=cell2mat(new_nudi1{1});
y_inside=cell2mat(new_nudi2{2});
y_outside=cell2mat(new_nudi1{2});

coords_inside=horzcat(x_inside,y_inside); %fuses x&y output
coords_outside=horzcat(x_outside,y_outside);

close all;

         imshow(zones,'InitialMagnification','fit','DisplayRange',[]);set(gcf,'Position',scr);
         title({'drag & drop until internal border is defined. A for new vertex';'Place border INSIDE animal. DBL Click on polygon line to finish'});
         hold on;
         fprintf(1,'Translating object outline into polygon, this can take a while\n')
         h_inside=impoly(gca,coords_inside);zoomcenter(center_x,center_y,zoom_factor);
         api = iptgetapi(h_inside);
         api.setColor('red');
         inside_pos = wait(h_inside);
         
         hold on;
         
         fprintf(1,'Translating object outline into polygon, this can take a while\n')
         title({'drag & drop until external border is defined. A for new vertex';'Place border OUTSIDE animal. DBL Click on polygon line to finish'});
         h_outside=impoly(gca,coords_outside);
         api = iptgetapi(h_outside);
         api.setColor('blue');
         outside_pos = wait(h_outside);
         
         close all;
         
         zones_interp=zones; %create copy of zones
         
%create masks for values to be interpolated
                middle_mask=poly2mask(xi,yi,nc,nr);
                inside_mask=poly2mask(x_inside,y_inside,nc,nr);
                outside_mask=poly2mask(x_outside,y_outside,nc,nr);

                donut_mask=logical(outside_mask-inside_mask);
                inner_donut_mask=logical(middle_mask-inside_mask);
                outer_donut_mask=logical(outside_mask-middle_mask);
                
ch=questdlg('Are there significant shadows at the animal/background border?','--',...
			'Yes','No','Yes');  %last one is default
		switch ch
			case 'Yes'
                

                zones_interp(outer_donut_mask)=NaN;
                zones_interp(inner_donut_mask)=NaN;
                %punch animal hole
                % zones_interp(inside_mask)=0;

                %%%%%interpolate from background edge to middle edge%%%%%%%%

                fprintf(1,'Interpolating Buffer Area, this can take a while\n')
                zones_interp=n_interp(zones_interp);

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

                % zones_interp=zones; %create copy of zones
 
                %  zones_interp(~inside_mask)=NaN; %remove anything but animal
                %  zones_interp(~outside_mask)=0;

                % zones_interp2=n_interp(zones_interp);

                % zones(donut_mask)=0; %punch donut into original zones

                zones(inner_donut_mask)=zones_interp(inner_donut_mask);%replace interpolated values
                zones(outer_donut_mask)=zones_interp(outer_donut_mask);

                % for i=rect(2):(rect(2)+rect(4))
                %     for j=rect(1):(rect(1)+rect(3))
                %         if zones(i,j)==0
                %             zones(i,j)=NaN;
                %         end
                %     end
                % end

                %zones=n_interp(zones);				
				
                figure(1);
                imshow(zones,'InitialMagnification','fit','DisplayRange',[]);set(gcf,'Position',scr);
                title({'drag & drop until border is defined. A for new vertex';'Place border INSIDE animal. DBL Click on polygon line to finish'});
                hold on;
                fprintf(1,'Translating object outline into polygon, this can take a while\n')
                h=impoly(gca,coords);zoomcenter(center_x,center_y,zoom_factor);
                api = iptgetapi(h);
                api.setColor('red');
                accepted_pos = wait(h);
               
                abt1=1; %set switch
                close all;
         
case 'No'
				         
        end;

 nudi_area_pixel=polyarea(accepted_pos(:,1),accepted_pos(:,2));            


 
%%%%%%%%%%%%%%%%%%%%%End of image segmentation for Animal vs background%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%use coords to make masks

background_mask=poly2mask(accepted_pos(:,1),accepted_pos(:,2),nc,nr);%turn polygon into mask
animal_mask=poly2mask(accepted_pos(:,1),accepted_pos(:,2),nc,nr);

%apply mask to zone map and image

zones_background=zones; %Create background zones
zones_background(background_mask)=0;

background_mask_colour=repmat(background_mask,[1,1,3]); %multiply mask to all RGB channels

zones_background_colour=zones_colour;
zones_background_colour(background_mask_colour)=0;



zones_animal=zones;
zones_animal(~animal_mask)=0;

animal_mask_colour=repmat(animal_mask,[1,1,3]); %multiply mask to all RGB channels

zones_animal_colour=zones_colour;
zones_animal_colour(~animal_mask_colour)=0;

%create margin zone map

%Entire Donut

zones_donut=zones;

zones_donut(~donut_mask)=0;

donut_mask_colour=repmat(donut_mask,[1,1,3]); %multiply mask to all RGB channels

zones_donut_colour=zones_colour;
zones_donut_colour(~donut_mask_colour)=0;

%this can be further broken down into inner donut and outer donut

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%Run pattern analaysis for animal%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%re-arrange zone IDs%%%%%%%%%%%

%This is necessary as the numbering of zones must be continuous without gaps for the Adjacency code to work

old_zone_ID=unique(zones_animal);
n_zones_animal=size(old_zone_ID);

for q=2:n_zones_animal(1)
    for i=1:nc
        for j=1:nr
            if zones_animal(i,j)==old_zone_ID(q)
                zones_animal(i,j)=q-1; %replaces zone value
            end
        end
    end
end

for i=1:nc
    for j=1:nr
if zones_animal(i,j)==0
    zones_animal(i,j)=n_zones_animal(1); %Makes the zeros the highest zone number (last one), requirement of doAdjacency
end
    end
end
 
%%%%%%%%%%%%Re-arrange zone names in zone_ID file (for pattern stats)%%%%%%

%New list of zones in zone map%

new_zone_ID=unique(zones_animal);

%Identify rows with the zones that are used

zones_ID_animal=zones_ID(ismember(zones_ID{:,14},old_zone_ID),:);%creates subset and overwrites old (gets rid of unused ones)

for i=1:(n_zones_animal(1)-1)
    zones_ID_animal{i,14}=i;
end


zlist_animal=unique(zones_animal); nz_animal=length(zlist_animal);


nz_animal=nz_animal-1;


%%%%%%%%get pattern geometry stats%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%Do Adjaceny for Animal zone map%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~,cdiv_animal,~,tdiv_animal,cmplx_animal,pxud_animal,pxlr_animal,asptr_animal,inout_animal,seq_animal]=DoAdjacency(zones_animal,pixsp); 
adjacency_animal=[cdiv_animal tdiv_animal cmplx_animal pxud_animal pxlr_animal asptr_animal inout_animal nz_animal];
mat_animal=seq_animal;


%%%%%% get visual contrast stats%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%for animal

zf_animal=zeros(nz_animal,1);

for i=1:nz_animal
zf_animal(i)=sum(sum(zones_animal==zlist_animal(i)));
end;

zf_animal=zf_animal/sum(zf_animal);  %these are the fractions (relative areas) of each zone in the animal

%%%get hue, chroma and luminance of zones for animal%%%%

figure(1);hold on;
title('Maxwell Triangle Animal');

hue_animal=zeros(nz_animal,1); chr_animal=hue_animal; lum_animal=hue_animal;

tx_list=zeros(nz_animal,1);
ty_list=zeros(nz_animal,1);

for cl=1:nz_animal
    
cc=[zones_ID_animal.swMean(zones_ID_animal.IDnew==zlist_animal(cl)) zones_ID_animal.mwMean(zones_ID_animal.IDnew==zlist_animal(cl)) zones_ID_animal.lwMean(zones_ID_animal.IDnew==zlist_animal(cl))];
lum_animal(cl)=zones_ID.dblMean(zones_ID.IDnew==zlist_animal(cl));

gy=[1/3 1/3 1/3]; [ox,oy]=ToTriangle(gy);
[tx,ty]=ToTriangle(cc);



ty_list(cl)=ty;
tx_list(cl)=tx;

tx=tx-ox; ty=ty-oy;

[hua,chrv]=cart2pol(tx,ty); 
if chrv<0.05 hua=NaN; end; %hue angle for chr>0.05
if hua<0 hua=hua+pi; end; %convert -180 to 180 to 0 to 360 in radians
hua=hua/(2*pi); % convert 0 to 360 into 0 to 1
     hue_animal(cl)=hua; chr_animal(cl)=chrv;
     clear rfl cc ccd gy ox oy tx ty hua crv;   
end;

for i=1:nz_animal
    if zf_animal(i)<0.01 
        
           plot(tx_list(i),ty_list(i),'x','LineWidth',2,'MarkerSize',3)
    end
    
    if zf_animal(i)<0.01 && zf_animal(i)>0.05
        
            plot(tx_list(i),ty_list(i),'x','LineWidth',2 ,'MarkerSize', 4)
    end   

    if zf_animal(i)<0.10 && zf_animal(i)>0.05
        
            plot(tx_list(i),ty_list(i),'x','LineWidth',2 ,'MarkerSize', 6)
    end   
    if zf_animal(i)<0.30 && zf_animal(i)>0.10
        
            plot(tx_list(i),ty_list(i),'x','LineWidth', 2,'MarkerSize',8)
    end      
    if zf_animal(i)<0.50 && zf_animal(i)>0.30
        
            plot(tx_list(i),ty_list(i),'x','LineWidth',2,'MarkerSize',12)
    end        
    if zf_animal(i)<1 && zf_animal(i)>0.5
        
            plot(tx_list(i),ty_list(i),'x','LineWidth',2,'MarkerSize',14)             
    end
 
end


legend('labels',{num2str(zlist_animal(1:nz_animal))},'Interpreter','none','FontSize',7,'FontWeight','bold');
legend('boxoff');

DrawTri;hold off;

%take weighted means and sd of hue chr lum, weighted by zone fractions

[mnhu,sd]=WeightedMnSD(hue_animal,zf_animal); cvhu=sd/mnhu;
[mnch,sd]=WeightedMnSD(chr_animal,zf_animal); cvch=sd/mnch;
[mnlm,sd]=WeightedMnSD(lum_animal,zf_animal); cvlm=sd/mnlm;

contrasts_animal=[mnhu mnch mnlm cvhu cvch cvlm];

%%%%%%%%%%%%%%Run pattern analysis for donut%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%re-arrange zone IDs%%%%%%%%%%%

%This is necessary as the numbering of zones must be continuous without gaps for the Adjacency code to work

old_zone_ID=unique(zones_donut);
n_zones_donut=size(old_zone_ID);

for q=2:n_zones_donut(1)
    for i=1:nc
        for j=1:nr
            if zones_donut(i,j)==old_zone_ID(q)
                zones_donut(i,j)=q-1; %replaces zone value
            end
        end
    end
end

for i=1:nc
    for j=1:nr
if zones_donut(i,j)==0
    zones_donut(i,j)=n_zones_donut(1); %Makes the zeros the highest zone number (last one), requirement of doAdjacency
end
    end
end
 
%%%%%%%%%%%%Re-arrange zone names in zone_ID file (for pattern stats)%%%%%%

%New list of zones in zone map%

new_zone_ID=unique(zones_donut);

%Identify rows with the zones that are used

zones_ID_donut=zones_ID(ismember(zones_ID{:,14},old_zone_ID),:);%creates subset and overwrites old (gets rid of unused ones)

for i=1:(n_zones_donut(1)-1)
    zones_ID_donut{i,14}=i;
end


zlist_donut=unique(zones_donut); nz_donut=length(zlist_donut);

nz_donut=nz_donut-1;


%%%%%%%%get pattern geometry stats%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~,cdiv_donut,~,tdiv_donut,cmplx_donut,pxud_donut,pxlr_donut,asptr_donut,inout_donut,seq_donut]=DoAdjacency(zones_donut,pixsp); 
adjacency_donut=[cdiv_donut tdiv_donut cmplx_donut pxud_donut pxlr_donut asptr_donut inout_donut nz_donut];

normTDiv_donut=(adjacency_donut(2)/((adjacency_donut(8))*(adjacency_donut(8)-1)/2));
normCDiv_donut=(adjacency_donut(1)/(adjacency_donut(8)));
ClrCmplx_donut=(1-normTDiv_donut)+(1-normCDiv_donut);

%%%%%% get visual contrast stats%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

zf_donut=zeros(nz_donut,1);

for i=1:nz_donut
zf_donut(i)=sum(sum(zones_donut==zlist_donut(i)));
end;

zf_donut=zf_donut/sum(zf_donut);  %these are the fractions (relative areas) of each zone in the animal

%%%get hue, chroma and luminance of zones%%%%

figure(9);hold on;
title('Maxwell Triangle donut');

hue_donut=zeros(nz_donut,1); chr_donut=hue_donut; lum_donut=hue_donut;

tx_list=zeros(nz_donut,1);
ty_list=zeros(nz_donut,1);

for cl=1:nz_donut
    
cc=[zones_ID_donut.swMean(zones_ID_donut.IDnew==zlist_donut(cl)) zones_ID_donut.mwMean(zones_ID_donut.IDnew==zlist_donut(cl)) zones_ID_donut.lwMean(zones_ID_donut.IDnew==zlist_donut(cl))];
lum_donut(cl)=zones_ID.dblMean(zones_ID.IDnew==zlist_donut(cl));

gy=[1/3 1/3 1/3]; [ox,oy]=ToTriangle(gy);
[tx,ty]=ToTriangle(cc);



ty_list(cl)=ty;
tx_list(cl)=tx;

tx=tx-ox; ty=ty-oy;

[hua,chrv]=cart2pol(tx,ty); 
if chrv<0.05 hua=NaN; end; %hue angle for chr>0.05
if hua<0 hua=hua+pi; end; %convert -180 to 180 to 0 to 360 in radians
hua=hua/(2*pi); % convert 0 to 360 into 0 to 1
     hue_donut(cl)=hua; chr_donut(cl)=chrv;
     clear rfl cc ccd gy ox oy tx ty hua crv;   
end;

for i=1:nz_donut
    if zf_donut(i)<0.01 
        
           plot(tx_list(i),ty_list(i),'x','LineWidth',2,'MarkerSize',3)
    end
    
    if zf_donut(i)<0.01 && zf_donut(i)>0.05
        
            plot(tx_list(i),ty_list(i),'x','LineWidth',2 ,'MarkerSize', 4)
    end   

    if zf_donut(i)<0.10 && zf_donut(i)>0.05
        
            plot(tx_list(i),ty_list(i),'x','LineWidth',2 ,'MarkerSize', 6)
    end   
    if zf_donut(i)<0.30 && zf_donut(i)>0.10
        
            plot(tx_list(i),ty_list(i),'x','LineWidth', 2,'MarkerSize',8)
    end      
    if zf_donut(i)<0.50 && zf_donut(i)>0.30
        
            plot(tx_list(i),ty_list(i),'x','LineWidth',2,'MarkerSize',12)
    end        
    if zf_donut(i)<1 && zf_donut(i)>0.5
        
            plot(tx_list(i),ty_list(i),'x','LineWidth',2,'MarkerSize',14)             
    end
 
end


legend('labels',{num2str(zlist_donut(1:nz_donut))},'Interpreter','none','FontSize',7,'FontWeight','bold');
legend('boxoff');

DrawTri;hold off;

%take weighted means and sd of hue chr lum, weighted by zone fractions

[mnhu,sd]=WeightedMnSD(hue_donut,zf_donut); cvhu=sd/mnhu;
[mnch,sd]=WeightedMnSD(chr_donut,zf_donut); cvch=sd/mnch;
[mnlm,sd]=WeightedMnSD(lum_donut,zf_donut); cvlm=sd/mnlm;

contrasts_donut=[mnhu mnch mnlm cvhu cvch cvlm];

%%%%%%%%%%%%%%Run pattern analysis for background%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%re-arrange zone IDs%%%%%%%%%%%

%This is necessary as the numbering of zones must be continuous without gaps for the Adjacency code to work

old_zone_ID=unique(zones_background);
n_zones_background=size(old_zone_ID);

for q=2:n_zones_background(1)
    for i=1:nc
        for j=1:nr
            if zones_background(i,j)==old_zone_ID(q)
                zones_background(i,j)=q-1; %replaces zone value
            end
        end
    end
end

for i=1:nc
    for j=1:nr
if zones_background(i,j)==0
    zones_background(i,j)=n_zones_background(1); %Makes the zeros the highest zone number (last one), requirement of doAdjacency
end
    end
end
 
%%%%%%%%%%%%Re-arrange zone names in zone_ID file (for pattern stats)%%%%%%

%New list of zones in zone map%

new_zone_ID=unique(zones_background);

%Identify rows with the zones that are used

zones_ID_background=zones_ID(ismember(zones_ID{:,14},old_zone_ID),:);%creates subset and overwrites old (gets rid of unused ones)

for i=1:(n_zones_background(1)-1)
    zones_ID_background{i,14}=i;
end


zlist_background=unique(zones_background); nz_background=length(zlist_background);

nz_background=nz_background-1;


%%%%%%%%get pattern geometry stats%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~,cdiv_background,~,tdiv_background,cmplx_background,pxud_background,pxlr_background,asptr_background,inout_background,seq_background]=DoAdjacency(zones_background,pixsp); 
adjacency_background=[cdiv_background tdiv_background cmplx_background pxud_background pxlr_background asptr_background inout_background nz_background];



%%%%%% get visual contrast stats%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

zf_background=zeros(nz_background,1);

for i=1:nz_background
zf_background(i)=sum(sum(zones_background==zlist_background(i)));
end;

zf_background=zf_background/sum(zf_background);  %these are the fractions (relative areas) of each zone in the animal

%%%get hue, chroma and luminance of zones%%%%

figure(2);hold on;
title('Maxwell Triangle Background');

hue_background=zeros(nz_background,1); chr_background=hue_background; lum_background=hue_background;

tx_list=zeros(nz_background,1);
ty_list=zeros(nz_background,1);

for cl=1:nz_background
    
cc=[zones_ID_background.swMean(zones_ID_background.IDnew==zlist_background(cl)) zones_ID_background.mwMean(zones_ID_background.IDnew==zlist_background(cl)) zones_ID_background.lwMean(zones_ID_background.IDnew==zlist_background(cl))];
lum_background(cl)=zones_ID.dblMean(zones_ID.IDnew==zlist_background(cl));

gy=[1/3 1/3 1/3]; [ox,oy]=ToTriangle(gy);
[tx,ty]=ToTriangle(cc);



ty_list(cl)=ty;
tx_list(cl)=tx;

tx=tx-ox; ty=ty-oy;

[hua,chrv]=cart2pol(tx,ty); 
if chrv<0.05 hua=NaN; end; %hue angle for chr>0.05
if hua<0 hua=hua+pi; end; %convert -180 to 180 to 0 to 360 in radians
hua=hua/(2*pi); % convert 0 to 360 into 0 to 1
     hue_background(cl)=hua; chr_background(cl)=chrv;
     clear rfl cc ccd gy ox oy tx ty hua crv;   
end;

for i=1:nz_background
    if zf_background(i)<0.01 
        
           plot(tx_list(i),ty_list(i),'x','LineWidth',2,'MarkerSize',3)
    end
    
    if zf_background(i)<0.01 && zf_background(i)>0.05
        
            plot(tx_list(i),ty_list(i),'x','LineWidth',2 ,'MarkerSize', 4)
    end   

    if zf_background(i)<0.10 && zf_background(i)>0.05
        
            plot(tx_list(i),ty_list(i),'x','LineWidth',2 ,'MarkerSize', 6)
    end   
    if zf_background(i)<0.30 && zf_background(i)>0.10
        
            plot(tx_list(i),ty_list(i),'x','LineWidth', 2,'MarkerSize',8)
    end      
    if zf_background(i)<0.50 && zf_background(i)>0.30
        
            plot(tx_list(i),ty_list(i),'x','LineWidth',2,'MarkerSize',12)
    end        
    if zf_background(i)<1 && zf_background(i)>0.5
        
            plot(tx_list(i),ty_list(i),'x','LineWidth',2,'MarkerSize',14)             
    end
 
end


legend('labels',{num2str(zlist_background(1:nz_background))},'Interpreter','none','FontSize',7,'FontWeight','bold');
legend('boxoff');

DrawTri;hold off;

%take weighted means and sd of hue chr lum, weighted by zone fractions

[mnhu,sd]=WeightedMnSD(hue_background,zf_background); cvhu=sd/mnhu;
[mnch,sd]=WeightedMnSD(chr_background,zf_background); cvch=sd/mnch;
[mnlm,sd]=WeightedMnSD(lum_background,zf_background); cvlm=sd/mnlm;

contrasts_background=[mnhu mnch mnlm cvhu cvch cvlm];

%%%%%%%%%%%%%%Substract animal from background%%%%%%%%%%%%%%%%%%%%%%%

% AspDiff : Background minus animal aspect ratio

AspDiff=adjacency_background(6)-adjacency_animal(6);
rel_AspDiff=((adjacency_background(6)/adjacency_animal(6))*100)-100;

% PxUDDiff : Background minus animal vertical average patch size

PxUDDiff=adjacency_background(4)-adjacency_animal(4);
rel_PxUDDiff=((adjacency_background(4)/adjacency_animal(4))*100)-100;

% PxLRDiff : Background minus animal horizontal average patch size

PxLRDiff=adjacency_background(5)-adjacency_animal(5);
rel_PxLRDiff=((adjacency_background(5)/adjacency_animal(5))*100)-100;

% TdivDiff : Background minus animal transition diversity

TdivDiff=adjacency_background(2)-adjacency_animal(2);
normTDiv_background=(adjacency_background(2)/((adjacency_background(8))*(adjacency_background(8)-1)/2));
normTDiv_animal=(adjacency_animal(2)/((adjacency_animal(8))*(adjacency_animal(8)-1)/2));
normTdivDiff=normTDiv_background-normTDiv_animal;
rel_normTdivDiff=(normTDiv_background/normTDiv_animal*100)-100;

% CDivDiff : Background minus animal colour diversity

CDivDiff=adjacency_background(1)-adjacency_animal(1);
normCDiv_background=(adjacency_background(1)/(adjacency_background(8)));
normCDiv_animal=(adjacency_animal(1)/(adjacency_animal(8)));
normCDivDiff=normCDiv_background-normCDiv_animal;
rel_normCDivDiff=(normCDiv_background/normCDiv_animal*100)-100;

%CmplxDiff : Background minus animal complexity

CmplxDiff=adjacency_background(3)-adjacency_animal(3);
rel_CmplxDiff=((adjacency_background(3)/adjacency_animal(3))*100)-100;

%ClrCmplx : combined normalised transition diversity and colour diversity
ClrCmplx_animal=(1-normCDiv_background)+(1-normTDiv_background);
ClrCmplx_background=(1-normCDiv_animal)+(1-normTDiv_animal);
ClrCmplxDiff=ClrCmplx_background-ClrCmplx_animal;
rel_ClrCmplxDiff=(ClrCmplx_background/ClrCmplx_animal*100)-100;

% CVLumDiff : Background minus animal coefficient of variance for luminance

CVLumDiff=contrasts_background(6)-contrasts_animal(6);
rel_CVLumDiff=((contrasts_background(6)/contrasts_animal(6))*100)-100;

% CVChrDiff : Background minus animal coefficient of variance for chroma

CVChrDiff=contrasts_background(5)-contrasts_animal(5);
rel_CVChrDiff=((contrasts_background(5)/contrasts_animal(5))*100)-100;

% CVHueDiff : Background minus animal coefficient of variance for hue angle

CVHueDiff=contrasts_background(4)-contrasts_animal(4);
rel_CVHueDiff=((contrasts_background(4)/contrasts_animal(4))*100)-100;

% MnLumDiff : Background minus animal mean luminance

MnLumDiff=contrasts_background(3)-contrasts_animal(3);
rel_MnLumDiff=((contrasts_background(3)/contrasts_animal(3))*100)-100;

% MnChrDiff : Background minus animal mean chroma

MnChrDiff=contrasts_background(2)-contrasts_animal(2);
rel_MnChrDiff=((contrasts_background(2)/contrasts_animal(2))*100)-100;

% MnHueDiff : Background minus animal mean hue angle

MnHueDiff=contrasts_background(1)-contrasts_animal(1);
rel_MnHueDiff=((contrasts_background(1)/contrasts_animal(1))*100)-100;

% AvgPatch : Background minus animal average patch area

AvgPatchDiff=((adjacency_background(4)/target)*(adjacency_background(5)/target))-((adjacency_animal(4)/target)*(adjacency_animal(5)/target));
rel_AvgPatchDiff=(((adjacency_background(4)/target)*(adjacency_background(5)/target))/((adjacency_animal(4)/target)*(adjacency_animal(5)/target))*100)-100;


%%%%%%%%%%%%%%Save output%%%%%%%%%%%%%%%%%%%%%%%%

parms_bgminusanimal=array2table([rel_ClrCmplxDiff ClrCmplxDiff CmplxDiff rel_CmplxDiff AspDiff rel_AspDiff PxUDDiff rel_PxUDDiff PxLRDiff rel_PxLRDiff AvgPatchDiff rel_AvgPatchDiff TdivDiff normTdivDiff rel_normTdivDiff CDivDiff normCDivDiff rel_normCDivDiff CVLumDiff rel_CVLumDiff CVChrDiff rel_CVChrDiff CVHueDiff rel_CVHueDiff MnLumDiff rel_MnLumDiff MnChrDiff rel_MnChrDiff MnHueDiff rel_MnHueDiff]);
parms_bgminusanimal.Properties.VariableNames = {'rel_ClrCmplxDiff' 'ClrCmplxDiff' 'CmplxDiff' 'rel_CmplxDiff' 'AspDiff' 'rel_AspDiff' 'PxUDDiff' 'rel_PxUDDiff' 'PxLRDiff' 'rel_PxLRDiff' 'AvgPatchDiff_mm2' 'rel_AvgPatchDiff_mm' 'TdivDiff' 'normTdivDiff' 'rel_normTdivDiff' 'CDivDiff' 'normCDivDiff' 'rel_normCDivDiff' 'CVLumDiff' 'rel_CVLumDiff' 'CVChrDiff' 'rel_CVChrDiff' 'CVHueDiff' 'rel_CVHueDiff' 'MnLumDiff' 'rel_MnLumDiff' 'MnChrDiff' 'rel_MnChrDiff' 'MnHueDiff' 'rel_MnHueDiff'};

if blr==1
     writetable(parms_bgminusanimal,[num2str(dist) 'mm' '_blurred_animal_vs_background.csv'],'Delimiter',';'); %saves pattern parms output
    
end

if blr==0
    writetable(parms_bgminusanimal,'Pattern_stats_unblurred_animal_vs_background.csv','Delimiter',';'); %saves pattern parms output
end


%Display all the graphs (animal, background and donut)

figure(3);
imshow(zones_background,[]);
figure(4);
imshow(zones_animal,[]);
figure(5);
imshow(zones_animal_colour);
figure(6);
imshow(zones_background_colour);
figure(7);
imshow(zones_donut,[])
figure(8);
imshow(zones_donut_colour);

% safe background only outputs in animal vs. background analysis

ohd='MnHue;MnChr;MnLum;CVhue;CVChr;CVLum;CDiv;Tdiv;Cmplx;PxUD;PxLR;asptr;k;norm_TDiv;norm_CDiv;ClrCmplx';

if blr==1
    oname1=[num2str(dist) '_mm_PatternStatistics_Background_only_blurred.csv']; fid1=fopen(oname1,'wt');
    oname2=[num2str(dist) '_mm' '_transitionmatrix_background_only_blurred.csv']; fid2=fopen(oname2,'wt');
    oname3=[num2str(dist) '_mm' '_colour_pattern_details_background_only_blurred.csv']; fid3=fopen(oname3,'wt');
    
    saveas(figure(3),[num2str(dist) '_mm' '_background_blurred_zones.png']);  %saves zone map as .png file
    [num2str(dist) '_mm' 'background_blurred__zones.png'];
    
    saveas(figure(2),[num2str(dist) 'mm' '_background_blurred_zones_maxwell.png']);  %saves zone map as .png file
    [num2str(dist) '_mm' 'background_blurred__zones_maxwell.png'];
    
    saveas(figure(6),[num2str(dist) 'mm' 'background_blurred_zones_png.png']);  %saves zone map as .png file
    [num2str(dist) '_mm' 'background_unblurred_zones_maxwell.png'];
    
    save([num2str(dist) 'mm' '_background_blurred_zones.mat'],'zones_background'); %saves zone map as matrix
    [num2str(dist) '_mm''background_blurred_zones.mat'];
    
end;
   
if blr==0
     oname1='PatternStatistics_Background_only_unblurred.csv'; fid1=fopen(oname1,'wt');
     oname2='Transitionmatrix_Background_only_unblurred.csv'; fid2=fopen(oname2,'wt');
     oname3='Colour_pattern_element_details_Background_only_unblurred.csv'; fid3=fopen(oname3,'wt');
    
     
    saveas(figure(3),['background_unblurred_zones.png']);  %saves zone map as .png file
    ['background_unblurred_zones.png'];
    
    saveas(figure(6),['background_unblurred_png.png']);  %saves zone map as .png file
    ['background_unblurred_zones_maxwell.png'];
    
    saveas(figure(2),['background_unblurred_zones_maxwell.png']);  %saves zone map as .png file
    ['background_unblurred_zones_maxwell.png'];
    
    save(['background_unblurred_zones.mat'],'zones_background'); %saves zone map as matrix
    ['background_unblurred_zones.mat'];
    
end;

fprintf(fid1,'%s\n',ohd);


for j=1 fprintf(fid1,'%6.4f',contrasts_background(j)); end;
for j=2:6 fprintf(fid1,';%6.4f',contrasts_background(j)); end;
for j=1:3 fprintf(fid1,';%6.4f',adjacency_background(j)); end;
for j=4:5 fprintf(fid1,';%6.2f',adjacency_background(j)); end;
for j=6 fprintf(fid1,';%6.4f',adjacency_background(j)); end;
for j=8 fprintf(fid1,';%6f',adjacency_background(j));end;
for j=9 fprintf(fid1,';%6.4f',normTDiv_background); end;
for j=10 fprintf(fid1,';%6.4f',normCDiv_background); end;
for j=11 fprintf(fid1,';%6.4f',ClrCmplx_background); end;

fprintf(fid1,'\n\n\n\n\n'); 
fprintf(fid1,'Meaning of columns\n');
fprintf(fid1,'MnHue,mean Hue angle (0 to 1=360degrees) weighted by patch relative area ignoring low chroma\n');
fprintf(fid1,'MnChr,mean Chroma, weighted by patch relative areas\n');
fprintf(fid1,'MnLum,mean Luminance, weighted by patch relative areas\n');
fprintf(fid1,'CVhue,CV (coefficient of variation) of Hue angle, weighted as above\n');
fprintf(fid1,'CVChr,CV of Chroma\n');
fprintf(fid1,'CVLum,CV of luminance\n');
fprintf(fid1,'CDiv,colour diversity\n');
fprintf(fid1,'Tdiv,transition diversity (diversity of off diagonals)\n');
fprintf(fid1,'Cmplx,transitions/sample (animal) a measure of pattern complexity\n');
fprintf(fid1,'PxUD,distance between transitions (pixels) up and down (vertical axis)\n');
fprintf(fid1,'PxLR,distance between transitions left and right (horizontal axis)\n');
fprintf(fid1,'asptr,transition aspect ratio (larger if longer in vertical axis)\n');

fclose(fid1);

fprintf(1,'\n');
fprintf(1,'Pattern stats in %s\n',oname1);
fprintf(1,'\n');
fprintf(1,'Transition Matrix in %s\n',oname2);
fprintf(1,'\n');
fprintf(1,'Pattern element details in %s\n',oname3);

%%%%print transition matrix of background only%%%%%%%%%%%%%%%%%%

fprintf(fid2, 'Adjacency_Transition_Matrix_background_only\n\n');
fprintf(fid2,' ');
for j=1:nz_background fprintf(fid2,'\t%15.0f',zlist_background(j)); end; fprintf(fid2,'\n');
     for k=1:nz_background
       fprintf(fid2,'%15.0f',zlist_background(k));
       for j=1:nz_background fprintf(fid2,'\t%15.0f',seq_background(k,j)); end; fprintf(fid2,'\n');
     end;
     fprintf(fid2,'\n');
     
%%%print pattern element details of background only%%%%

fprintf(fid3,'ZoneID');
fprintf(fid3,';HueAngle');
fprintf(fid3,';Chroma');
fprintf(fid3,';Luminance');
fprintf(fid3,';SWccq');
fprintf(fid3,';MWccq');
fprintf(fid3,';LWccq');
fprintf(fid3,';DBLccq');
fprintf(fid3,';RelativeSize');
fprintf(fid3,';TotalSize\n');

for i=1:nz_background
    fprintf(fid3,'%6f',i);
    fprintf(fid3,';%6.4f',hue_background(i));
    fprintf(fid3,';%6.4f',chr_background(i));
    fprintf(fid3,';%6.4f',lum_background(i));
    fprintf(fid3,';%6.4f',zones_ID_background{i,4});
    fprintf(fid3,';%6.4f',zones_ID_background{i,6});
    fprintf(fid3,';%6.4f',zones_ID_background{i,8});
    fprintf(fid3,';%6.4f',zones_ID_background{i,10});
    fprintf(fid3,';%6.4f',zf_background(i));
    fprintf(fid3,';%6.4f',zones_ID_background{i,12});
    fprintf(fid3,'\n');
end

%%%%Safe animal only outputs in animal vs. background analysis%%%%%%%%%%%%

if blr==1
    oname1=[num2str(dist) '_mm_PatternStatistics_animal_only_blurred.csv']; fid1=fopen(oname1,'wt');
    oname2=[num2str(dist) '_mm' '_transitionmatrix_animal_only_blurred.csv']; fid2=fopen(oname2,'wt');
    oname3=[num2str(dist) '_mm' '_colour_pattern_details_animal_only_blurred.csv']; fid3=fopen(oname3,'wt');
    
    saveas(figure(4),[num2str(dist) '_mm' '_animal_blurred_zones.png']);  %saves zone map as .png file
    [num2str(dist) '_mm' 'animal_blurred__zones.png'];
    
    saveas(figure(1),[num2str(dist) 'mm' '_animal_blurred_zones_maxwell.png']);  %saves zone map as .png file
    [num2str(dist) '_mm' 'animal_blurred__zones_maxwell.png'];
    
    saveas(figure(5),[num2str(dist) 'mm' 'animal_blurred_zones_png.png']);  %saves zone map as .png file
    [num2str(dist) '_mm' 'animal_unblurred_zones_maxwell.png'];
    
    save([num2str(dist) 'mm' '_animal_blurred_zones.mat'],'zones_animal'); %saves zone map as matrix
    [num2str(dist) '_mm''animal_blurred_zones.mat'];
    
end;
   
if blr==0
     oname1='PatternStatistics_animal_only_unblurred.csv'; fid1=fopen(oname1,'wt');
     oname2='Transitionmatrix_animal_only_unblurred.csv'; fid2=fopen(oname2,'wt');
     oname3='Colour_pattern_element_details_animal_only_unblurred.csv'; fid3=fopen(oname3,'wt');
    
     
    saveas(figure(4),['animal_unblurred_zones.png']);  %saves zone map as .png file
    ['animal_unblurred_zones.png'];
    
    saveas(figure(5),['animal_unblurred_png.png']);  %saves zone map as .png file
    ['animal_unblurred_zones_maxwell.png'];
    
    saveas(figure(1),['animal_unblurred_zones_maxwell.png']);  %saves zone map as .png file
    ['animal_unblurred_zones_maxwell.png'];
    
    save(['animal_unblurred_zones.mat'],'zones_animal'); %saves zone map as matrix
    ['animal_unblurred_zones.mat'];
    
end;

fprintf(fid1,'%s\n',ohd);


for j=1 fprintf(fid1,'%6.4f',contrasts_animal(j)); end;
for j=2:6 fprintf(fid1,';%6.4f',contrasts_animal(j)); end;
for j=1:3 fprintf(fid1,';%6.4f',adjacency_animal(j)); end;
for j=4:5 fprintf(fid1,';%6.2f',adjacency_animal(j)); end;
for j=6 fprintf(fid1,';%6.4f',adjacency_animal(j)); end;
for j=8 fprintf(fid1,';%6f',adjacency_animal(j));end;
for j=9 fprintf(fid1,';%6.4f',normTDiv_animal); end;
for j=10 fprintf(fid1,';%6.4f',normCDiv_animal); end;
for j=11 fprintf(fid1,';%6.4f',ClrCmplx_animal); end;

fprintf(fid1,'\n\n\n\n\n'); 
fprintf(fid1,'Meaning of columns\n');
fprintf(fid1,'MnHue,mean Hue angle (0 to 1=360degrees) weighted by patch relative area ignoring low chroma\n');
fprintf(fid1,'MnChr,mean Chroma, weighted by patch relative areas\n');
fprintf(fid1,'MnLum,mean Luminance, weighted by patch relative areas\n');
fprintf(fid1,'CVhue,CV (coefficient of variation) of Hue angle, weighted as above\n');
fprintf(fid1,'CVChr,CV of Chroma\n');
fprintf(fid1,'CVLum,CV of luminance\n');
fprintf(fid1,'CDiv,colour diversity\n');
fprintf(fid1,'Tdiv,transition diversity (diversity of off diagonals)\n');
fprintf(fid1,'Cmplx,transitions/sample (animal) a measure of pattern complexity\n');
fprintf(fid1,'PxUD,distance between transitions (pixels) up and down (vertical axis)\n');
fprintf(fid1,'PxLR,distance between transitions left and right (horizontal axis)\n');
fprintf(fid1,'asptr,transition aspect ratio (larger if longer in vertical axis)\n');

fclose(fid1);

fprintf(1,'\n');
fprintf(1,'Pattern stats in %s\n',oname1);
fprintf(1,'\n');
fprintf(1,'Transition Matrix in %s\n',oname2);
fprintf(1,'\n');
fprintf(1,'Pattern element details in %s\n',oname3);

%%%%print transition matrix of animal only%%%%%%%%%%%%%%%%%%

fprintf(fid2, 'Adjacency_Transition_Matrix_animal_only\n\n');
fprintf(fid2,' ');
for j=1:nz_animal fprintf(fid2,'\t%15.0f',zlist_animal(j)); end; fprintf(fid2,'\n');
     for k=1:nz_animal
       fprintf(fid2,'%15.0f',zlist_animal(k));
       for j=1:nz_animal fprintf(fid2,'\t%15.0f',seq_animal(k,j)); end; fprintf(fid2,'\n');
     end;
     fprintf(fid2,'\n');
     
%%%print pattern element details of background only%%%%

fprintf(fid3,'ZoneID');
fprintf(fid3,';HueAngle');
fprintf(fid3,';Chroma');
fprintf(fid3,';Luminance');
fprintf(fid3,';SWccq');
fprintf(fid3,';MWccq');
fprintf(fid3,';LWccq');
fprintf(fid3,';DBLccq');
fprintf(fid3,';RelativeSize');
fprintf(fid3,';TotalSize\n');

for i=1:nz_animal
    fprintf(fid3,'%6f',i);
    fprintf(fid3,';%6.4f',hue_animal(i));
    fprintf(fid3,';%6.4f',chr_animal(i));
    fprintf(fid3,';%6.4f',lum_animal(i));
    fprintf(fid3,';%6.4f',zones_ID_animal{i,4});
    fprintf(fid3,';%6.4f',zones_ID_animal{i,6});
    fprintf(fid3,';%6.4f',zones_ID_animal{i,8});
    fprintf(fid3,';%6.4f',zones_ID_animal{i,10});
    fprintf(fid3,';%6.4f',zf_animal(i));
    fprintf(fid3,';%6.4f',zones_ID_animal{i,12});
    fprintf(fid3,'\n');
end

% safe donut only outputs in animal vs. background analysis

if blr==1
    oname1=[num2str(dist) '_mm_PatternStatistics_Donut_only_blurred.csv']; fid1=fopen(oname1,'wt');
    oname2=[num2str(dist) '_mm' '_transitionmatrix_Donut_only_blurred.csv']; fid2=fopen(oname2,'wt');
    oname3=[num2str(dist) '_mm' '_colour_pattern_details_Donut_only_blurred.csv']; fid3=fopen(oname3,'wt');
    
    saveas(figure(3),[num2str(dist) '_mm' '_donut_blurred_zones.png']);  %saves zone map as .png file
    [num2str(dist) '_mm' 'donut_blurred_zones.png'];
    
    saveas(figure(2),[num2str(dist) 'mm' '_donut_blurred_zones_maxwell.png']);  %saves zone map as .png file
    [num2str(dist) '_mm' 'donut_blurred_zones_maxwell.png'];
    
    saveas(figure(6),[num2str(dist) 'mm' 'donut_blurred_zones_png.png']);  %saves zone map as .png file
    [num2str(dist) '_mm' 'donut_unblurred_zones_maxwell.png'];
    
    save([num2str(dist) 'mm' '_donut_blurred_zones.mat'],'zones_donut'); %saves zone map as matrix
    [num2str(dist) '_mm''donut_blurred_zones.mat'];
    
end;
   
if blr==0
     oname1='PatternStatistics_donut_only_unblurred.csv'; fid1=fopen(oname1,'wt');
     oname2='Transitionmatrix_donut_only_unblurred.csv'; fid2=fopen(oname2,'wt');
     oname3='Colour_pattern_element_details_donut_only_unblurred.csv'; fid3=fopen(oname3,'wt');
    
     
    saveas(figure(7),['donut_unblurred_zones.png']);  %saves zone map as .png file
    ['donut_unblurred_zones.png'];
    
    saveas(figure(8),['donut_unblurred_png.png']);  %saves zone map as .png file
    ['donut_unblurred_zones_maxwell.png'];
    
    saveas(figure(9),['donut_unblurred_zones_maxwell.png']);  %saves zone map as .png file
    ['donut_unblurred_zones_maxwell.png'];
    
    save(['donut_unblurred_zones.mat'],'zones_donut'); %saves zone map as matrix
    ['donut_unblurred_zones.mat'];
    
end;

fprintf(fid1,'%s\n',ohd);


for j=1 fprintf(fid1,'%6.4f',contrasts_donut(j)); end;
for j=2:6 fprintf(fid1,';%6.4f',contrasts_donut(j)); end;
for j=1:3 fprintf(fid1,';%6.4f',adjacency_donut(j)); end;
for j=4:5 fprintf(fid1,';%6.2f',adjacency_donut(j)); end;
for j=6 fprintf(fid1,';%6.4f',adjacency_donut(j)); end;
for j=8 fprintf(fid1,';%6f',adjacency_donut(j));end;
for j=9 fprintf(fid1,';%6.4f',normTDiv_donut); end;
for j=10 fprintf(fid1,';%6.4f',normCDiv_donut); end;
for j=11 fprintf(fid1,';%6.4f',ClrCmplx_donut); end;

fprintf(fid1,'\n\n\n\n\n'); 
fprintf(fid1,'Meaning of columns\n');
fprintf(fid1,'MnHue,mean Hue angle (0 to 1=360degrees) weighted by patch relative area ignoring low chroma\n');
fprintf(fid1,'MnChr,mean Chroma, weighted by patch relative areas\n');
fprintf(fid1,'MnLum,mean Luminance, weighted by patch relative areas\n');
fprintf(fid1,'CVhue,CV (coefficient of variation) of Hue angle, weighted as above\n');
fprintf(fid1,'CVChr,CV of Chroma\n');
fprintf(fid1,'CVLum,CV of luminance\n');
fprintf(fid1,'CDiv,colour diversity\n');
fprintf(fid1,'Tdiv,transition diversity (diversity of off diagonals)\n');
fprintf(fid1,'Cmplx,transitions/sample (animal) a measure of pattern complexity\n');
fprintf(fid1,'PxUD,distance between transitions (pixels) up and down (vertical axis)\n');
fprintf(fid1,'PxLR,distance between transitions left and right (horizontal axis)\n');
fprintf(fid1,'asptr,transition aspect ratio (larger if longer in vertical axis)\n');

fclose(fid1);

fprintf(1,'\n');
fprintf(1,'Pattern stats in %s\n',oname1);
fprintf(1,'\n');
fprintf(1,'Transition Matrix in %s\n',oname2);
fprintf(1,'\n');
fprintf(1,'Pattern element details in %s\n',oname3);

%%%%print transition matrix of background only%%%%%%%%%%%%%%%%%%

fprintf(fid2, 'Adjacency_Transition_Matrix_donut_only\n\n');
fprintf(fid2,' ');
for j=1:nz_donut fprintf(fid2,'\t%15.0f',zlist_donut(j)); end; fprintf(fid2,'\n');
     for k=1:nz_donut
       fprintf(fid2,'%15.0f',zlist_donut(k));
       for j=1:nz_donut fprintf(fid2,'\t%15.0f',seq_donut(k,j)); end; fprintf(fid2,'\n');
     end;
     fprintf(fid2,'\n');
     
%%%print pattern element details of donut only%%%%

fprintf(fid3,'ZoneID');
fprintf(fid3,',HueAngle');
fprintf(fid3,',Chroma');
fprintf(fid3,',Luminance');
fprintf(fid3,',SWccq');
fprintf(fid3,',MWccq');
fprintf(fid3,',LWccq');
fprintf(fid3,',DBLccq');
fprintf(fid3,',RelativeSize');
fprintf(fid3,',TotalSize\n');

for i=1:nz_donut
    fprintf(fid3,'%6f',i);
    fprintf(fid3,';%6.4f',hue_donut(i));
    fprintf(fid3,';%6.4f',chr_donut(i));
    fprintf(fid3,';%6.4f',lum_donut(i));
    fprintf(fid3,';%6.4f',zones_ID_donut{i,4});
    fprintf(fid3,';%6.4f',zones_ID_donut{i,6});
    fprintf(fid3,';%6.4f',zones_ID_donut{i,8});
    fprintf(fid3,';%6.4f',zones_ID_donut{i,10});
    fprintf(fid3,';%6.4f',zf_donut(i));
    fprintf(fid3,';%6.4f',zones_ID_donut{i,12});
    fprintf(fid3,'\n');
end

% Add titanic diagramms as output


    end;
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%ADDONS and FUNCTIONS%%%%%%%%%%%%%%%%%%%%%%%%

function zones_out = modifyzonemap(zones_in)

global H zonesG
zonesG = zones_in;
%Code to manually modify zone map
scr=groot; scr=scr.ScreenSize;
H.fig = figure(1);clf;

H.ax = axes;

H.im = imshow(zonesG,'parent',H.ax);set(gcf,'Position',scr);
hold on
H.pp = rectangle('position',[0 0 1 1],'EdgeColor','r');
hold off
set(H.im,'buttondownfcn',@buttondown)
set(H.fig,'keypressfcn',@keypress,...
	'windowbuttonmotionfcn',[])
H.ttext = {'Modify: Rclick selects zone, Lclick paints over, +- adjusts brush size, z = toggle zoom, d = done'};
H.t = title(H.ttext);
H.zp=0;H.bs = 10; 
H.zoomed=0;
H.endprog=0;
while H.endprog==0
	waitforbuttonpress;
	
end
zones_out = zonesG;
close(H.fig)
end

%%%%%%%%%%%%%

function zones_out = modifyzonemap2(zones_in)

global H zonesG
zonesG = zones_in;
%Code to manually modify zone map
scr=groot; scr=scr.ScreenSize;
H.fig = figure(1);clf;

H.ax = axes;

H.im = imshow(zonesG,'parent',H.ax,'DisplayRange',[]);

set(gcf,'Position',scr);
hold on
H.pp = rectangle('position',[0 0 1 1],'EdgeColor','r');
hold off
set(H.im,'buttondownfcn',@buttondown)
set(H.fig,'keypressfcn',@keypress,...
	'windowbuttonmotionfcn',[])
H.ttext = {'Modify: Rclick selects zone, Lclick paints over, +- adjusts brush size, z = toggle zoom, d = done'};
H.t = title(H.ttext);
H.zp=0;H.bs = 10; 
H.zoomed=0;
H.endprog=0;
while H.endprog==0
	waitforbuttonpress;
	
end
zones_out = zonesG;
close(H.fig)
end

%%%%%%%%%%%%%

function buttondown(A,B)
global H zonesG

switch B.Button
	case 1 %draw over zone
		set(H.im,'buttondownfcn',[])
		set(H.fig,'windowbuttonmotionfcn',@dragmouse2,...
			'windowbuttonupfcn',@buttonup)
	case 3 %select zone for palate
		xy = get(H.ax,'CurrentPoint'); xy = round(xy(1,1:2));
		H.zp = zonesG(xy(2),xy(1));
        
end
end

%%%%%%%%%%%%%

function dragmouse2(~,~)
global H zonesG
xy = get(H.ax,'CurrentPoint');xy = round(xy(1,1:2));
set(H.pp,'position',[xy(1)-H.bs/2 xy(2)-H.bs/2 H.bs H.bs]);
zonesG(xy(2)-H.bs/2:xy(2)+H.bs/2,xy(1)-H.bs/2:xy(1)+H.bs/2,:)=H.zp;
set(H.im,'cdata',zonesG)
end

%%%%%%%%%%%%%

function buttonup(~,~)
global H
set(H.im,'buttondownfcn',@buttondown)
set(H.fig,'windowbuttonmotionfcn',[],...
	'windowbuttonupfcn',[])
set(H.pp,'position',[0 0 1 1])
end

%%%%%%%%%%%%%

function keypress(~,B)
global H zonesG
switch B.Character
	case '+' %Increase brush size
		H.bs = H.bs+2;
		disp(['Brush size: ' num2str(H.bs)]);
	case '-' %Decrease brush size
		H.bs = H.bs-2; 
		disp(['Brush size: ' num2str(H.bs)]);
		if H.bs==0, H.bs=2;end
	case 'z' %Zoom
		if H.zoomed==0
			set(H.im,'buttondownfcn',[])
			set(H.ax,'xlim',[0 size(zonesG,2)],...
				'ylim',[0 size(zonesG,1)])
			set(H.t,'string','Drag rectangle over zoom area')
			k = waitforbuttonpress;
			point1 = get(H.ax,'CurrentPoint');    % button down detected
			finalRect = rbbox;pause(0.01);                 % return figure units
			point2 = get(H.ax,'CurrentPoint');    % button up detected
			point1 = point1(1,1:2);point2 = point2(1,1:2);
			p1 = min(point1,point2);           % calculate locations
			offset = abs(point1-point2);
			set(H.ax,'xlim',[p1(1) p1(1)+offset(1)],'ylim',[p1(2) p1(2)+offset(2)]);
			set(H.im,'buttondownfcn',@buttondown)
			set(H.t,'string',H.ttext)
			H.zoomed=1;
		else
			set(H.ax,'xlim',[0 size(zonesG,2)],...
				'ylim',[0 size(zonesG,1)])
			H.zoomed=0;
		end
	case 'd' %Done
		H.endprog=1;
end
end

%%%%%%%%%%%%
   
function zoomcenter(varargin)
%ZOOMCENTER Zoom in and out of a specifeid point on a 2-D plot.
% ZOOMCENTER(X,Y) zooms the current axis on the point (X,Y) by a
%factor of 2.5.
% ZOOMCENTER(X,Y,FACTOR) zooms the current axis on the point (X,Y) by
%FACTOR.
%
% ZOOMCENTER(AX,...) zooms on the specified axis
%
% Example:
% line
% zoomcenter(.5, .5, 10)
%
% line
% zoomcenter(.7, .3, .5)

nin = nargin;
if nin==0
 error('ZOOMCENTER requires at least 2 inputs');
end
if ishandle(varargin{1})
 ax = varargin{1};
 varargin = varargin(2:end);
 nin = nin-1;
else
 ax = gca;
end
if nin<2
 error('ZOOMCENTER requires specifying both X and Y');
else
 x = varargin{1};
 y = varargin{2};
end
if nin==3
 factor = varargin{3};
else
 factor = 2.5;
end

cax = axis(ax);
daxX = (cax(2)-cax(1))/factor(1)/2;
daxY = (cax(4)-cax(3))/factor(end)/2;
axis(ax,[x+[-1 1]*daxX y+[-1 1]*daxY]); 
   end
    
%%%%%%%%%%%%

function [ex,ey,eh]=DrawFromXY(fignum,xo,yo)
% [ex,ey,eh]=DrawFromXY(fignum,xo,yo);
% given a figure number and a starting point (including outside axes)
%   click and drag endpoint of line to new place on figure
%   will then store x,y endpoint and handle to the line (for erasing it)
global ah xdat ydat;
tx=0; ty=0; th=0; ex=0; ey=0; eh=0;
set(fignum,'WindowButtonDownFcn',@wbdf); ah=gca;  
set(ah,'DrawMode','fast');
ankh=plot(xo,yo,'.k'); %to make this the pivot (current) point
uiwait; %*********************************************************
% title('Click and drag line, release to anchor')
  function wbdf(src,~)
     set(src,'pointer','crosshair')
     th = line('xdata',xo,'ydata',yo,'Marker','+','color','b');
     cp = get(ah,'CurrentPoint');
     xdat = [xo,cp(1,1)];
     ydat = [yo,cp(1,2)];
     set(th,'xdata',xdat,'ydata',ydat); drawnow     
     set(src,'WindowButtonMotionFcn',@wbmf) 
     set(src,'WindowButtonUpFcn',@wbuf)
     %subfunctions
     function wbmf(~,~)
        cp = get(ah,'CurrentPoint');
        xdat = [xo,cp(1,1)];
        ydat = [yo,cp(1,2)];
        set(th,'xdata',xdat,'ydata',ydat); drawnow
     end %wbmf
     function wbuf(src,~)
        tx=xdat(2); ty=ydat(2);
        % fprintf(1,'end-point %6.3f %6.3f\n',tx,ty); 
        set(src,'WindowButtonMotionFcn','');
        set(src,'WindowButtonUpFcn','')
        set(src,'WindowButtonDownFcn','')
        set(src,'Pointer','arrow')
        uiresume;
        ex=tx; ey=ty; eh=th; delete(ankh);
        %delete(findobj(fignum,'type','uicontrol'));   %remove all uicontrols!
        return
end %wbuf
  end %wbdf
end

%%%%%%%%%%%%

function output=rotateAround(image, pointY, pointX, angle, varargin)
% ROTATEAROUND rotates an image.
%   ROTATED=ROTATEAROUND(IMAGE, POINTY, POINTX, ANGLE) rotates IMAGE around
%   the point [POINTY, POINTX] by ANGLE degrees. To rotate the image
%   clockwise, specify a negative value for ANGLE.
%
%   ROTATED=ROTATEAROUND(IMAGE, POINTY, POINTX, ANGLE, METHOD) rotates the
%   image with specified method:
%       'nearest'       Nearest-neighbor interpolation
%       'bilinear'      Bilinear interpolation
%       'bicubic'       Bicubic interpolation
%    The default is fast 'nearest'. Switch to 'bicubic' for nicer results.
%
%   Example
%   -------
%       imshow(rotateAround(imread('eight.tif'), 1, 1, 10));
%
%   See also IMROTATE, PADARRAY.

%   Contributed by Jan Motl (jan@motl.us)
%   $Revision: 1.2 $  $Date: 2014/05/01 12:08:01 $

% Parameter checking.
numvarargs = length(varargin);
if numvarargs > 1
    error('myfuns:somefun2Alt:TooManyInputs', ...
        'requires at most 1 optional input');
end
optargs = {'nearest'};    % Set defaults for optional inputs
optargs(1:numvarargs) = varargin;
[method] = optargs{:};    % Place optional args in memorable variable names

% Initialization.
[imageHeight imageWidth ~] = size(image);
centerX = floor(imageWidth/2+1);
centerY = floor(imageHeight/2+1);

dy = centerY-pointY;
dx = centerX-pointX;

% How much would the "rotate around" point shift if the 
% image was rotated about the image center. 
[theta, rho] = cart2pol(-dx,dy);
[newX, newY] = pol2cart(theta+angle*(pi/180), rho);
shiftX = round(pointX-(centerX+newX));
shiftY = round(pointY-(centerY-newY));

% Pad the image to preserve the whole image during the rotation.
padX = abs(shiftX);
padY = abs(shiftY);

padded = padarray(image, [padY padX]);

% Rotate the image around the center.
rot = imrotate(padded, angle, method, 'loose');

% Crop the image.
output = rot(padY+1-shiftY:end-padY-shiftY, padX+1-shiftX:end-padX-shiftX, :);
end

%%%%%%%%%%%%

function [mn,sd]=WeightedMnSD(x,w)
% [mn,sd]=WeightedMnSD(x,w);
% Calculates the mean and SD for x data with weights w
% INPUT: x values
%        w weights for each value
% Both x and w must be the same length 
% Will automatically remove any rows with NaN in x
t=isnan(x);
if sum(t)>0  %remove NaN rows
 xx=x(t==0); ww=w(t==0); 
 x=xx; w=ww; w=w/sum(w);
end;
n=length(x); n2=length(w);
if n~=n2 
    mn=NaN; sd=NaN; 
    fprintf(1,'X and weights do not have same n\n');
    return; 
end;
sw=sum(w);     %sum of weights
swx=sum(x.*w); %sum of weights times x
mn=swx/sw;     %weighted mean
nnz=sum(w(w>0)>0); %number of nonzero weights
if nnz>1
  s=0;
  for i=1:n
    s=s+w(i)*(x(i)-mn)^2;
  end;
  sd=sqrt(nnz*s/((nnz-1)*sw));    %math.stackexchange.com
else
  sd=0;
end;
end  %weighted MN-sd

%%%%%%%%%%%%%

function u=n_interp(u)

nanLocations = isnan(u);
nanLinearIndexes = find(nanLocations);
nonNanLinearIndexes = setdiff(1:numel(u), nanLinearIndexes);

% Get the x,y,z of all other locations that are non nan.
[xGood, yGood, zGood] = ind2sub(size(u), nonNanLinearIndexes);
for index = 1 : length(nanLinearIndexes);
	thisLinearIndex = nanLinearIndexes(index);
	% Get the x,y,z location
	[x,y,z] = ind2sub(size(u), thisLinearIndex);
	% Get distances of this location to all the other locations
 	distances = sqrt((x-xGood).^2 + (y - yGood) .^ 2 + (z - zGood) .^ 2);
	[sortedDistances, sortedIndexes] = sort(distances, 'ascend');
	% The closest non-nan value will be located at index sortedIndexes(1)
	indexOfClosest = sortedIndexes(1);
	% Get the u value there.
	goodValue = u(xGood(indexOfClosest), yGood(indexOfClosest), zGood(indexOfClosest));
	% Replace the bad nan value in u with the good value.
	u(x,y,z) = goodValue;
    
end;
end

%%%%%%%%%%%%%

function DrawTri(vertnames)
% DrawTri draws triangle for triangular coordinates
% call as drawtri(vertnames) where vertnames holds names of x,y,z
% or      drawtri            which will put in X Y and Z
% coordinates are trix=[0 1 0.5 0]; triy=[0 0 sqrt(3)/2 0];
if nargin==0 vertnames={'X','Z','Y',' '}; end;
ht=sqrt(3)/2; trix=[0 1 0.5 0]; triy=[0 0 ht 0];
plot(trix,triy,'Color','black'); hold on; text(trix,triy,vertnames); axis equal
end

%%%%%%%%%%%%%

function [X,Y] = polyout(x1,y1,delta,join,info)
% POLYOUT  Outset (expand/inflate) a polygon.
% 
%  [X,Y] = POLYOUT(X1,Y1,DELTA,JOIN,JOININFO) outsets the polygon given by
%  points (X1(i),Y1(i)).  The corners are joined based on the (optional) JOIN:
%    JOIN = 'm' or 'miter' ==> exact corners but square at small angles (default)
%               Optional JOININFO is the miter limit, a multiple of DELTA;
%               if the corner point would be moved more than JOININFO*DELTA,
%               then it is squared off instead.  The default value, as well
%               as the minimum allowed, is 2.
%    JOIN = 's' or 'square' ==> square off corners
%               JOININFO is ignored.
%    JOIN = 'r' or 'round'  ==> round corners
%               Required JOININFO sets the precision of points along the arc
%               (smaller JOININFO ==> more points along the arc); for a 180deg
%               arc, the number of points is (pi/acos(1-JOININFO/DELTA) using
%               the same scaling as the polygon points.
%  X1 and Y1 must be vectors of the same size. JOININFO must be a scalar double.
%
%  X and Y are cell arrays (because the result may be multiple polygons)
%  containing the x and y coordinates of the resulting polygon(s).
%  
%  XY = POLYOUT([X1 Y1],DELTA,JOIN,JOININFO) does the same for column vectors
%  X1 and Y1, and XY is {X Y}.
%
% See also CLIPPER, POLYCLIP.

% Copyright (c)2015-17, Prof. Erik A. Johnson <JohnsonE@usc.edu>, 01/28/17

% 09/13/15  EAJ  Initial code
% 01/28/17  EAJ  Update for newer MATLAB versions
% 01/29/17  EAJ  Correct [X1 Y1] code (change x & y to x1 & y1)

narginchk(2,5)
nargoutchk(0,2)
if nargin==2 || (nargin>=3&&ischar(delta))  % POLYOUT([X1 Y1],DELTA,...)
	narginchk(2,4)
	if nargin>=3,
		if nargin==4
			info = join;
		end;
		join = delta;
	else
		join = [];
		info = [];
	end
	delta = y1;
	y1 = x1(:,end/2+1:end);
	x1(:,end/2+1:end) = [];
else
	if nargin<4, join=[]; end;
	if nargin<5, info=[]; end;
end
if isempty(join), join='miter'; end; %default join method
if ~ischar(join), error('JOIN must be a string'); end;
if isempty(info)
	infoArg = {};
else
	infoArg = {info};
end

scale=2^32;
% leave this in array format just in case we later adapt this
pack   = @(p) arrayfun(@(x) struct('x', int64(x.x*scale),'y', int64(x.y*scale)),p);
unpack = @(p) arrayfun(@(x) struct('x',double(x.x)/scale,'y',double(x.y)/scale),p);
packdelta   = @(d) d*scale;
packarc = packdelta;
packmiter = @(m) m;

if lower(join(1))=='r'
	if isempty(info)
		info = (abs(delta)+(delta==0))*(1-cosd(5)); % about every 5 degrees
		% error('A ''round'' JOIN must provice the ArcTolerance as the fourth argument.')
	end
	infoArg = {packarc(info)};
end

if ~isnumeric(x1) || ~isnumeric(y1)
	error('The polygon coordinates must be numeric matrices.')
elseif numel(x1)~=length(x1) || numel(y1)~=length(y1)
	error('Function only handles a single polygon');
elseif numel(x1)~=numel(y1)
	error('X1 must be the same size as Y1, and X2 the same size as Y2.');
end;

x1=x1(:); y1=y1(:); % ensure column

poly1 = struct('x',num2cell(x1,1),'y',num2cell(y1,1));
assert(numel(poly1)==1, 'Only single polygon outsets are allowed in the current code');
poly3 = unpack(clipper(pack(poly1),packdelta(delta),join,infoArg{:}));

if isempty(poly3)
	x = cell(1,0);
	y = cell(1,0);
else
	x = {poly3.x};
	y = {poly3.y};
end

if nargout>=2
	X=x; Y=y;
else
	X={x y};
end
end

%%%%%%%%%%%%

function [tx,ty]=ToTriangle(xyz)
% [tx,ty]=ToTriangle(xyz)
%    converts data matrix xyz=[x y z] of x,y,z coordinates
%    to triangular coordinates tx,ty (both column vectors)
%    with sides=1, height sqrt(3)/2
% Will temporarily make xyz rows sum to 1
% Coordinates of triangle are trix=[0 1 0.5 0]; triy=[0 0 sqrt(3)/2 0];
% Same as plottri but only saves triangular coordinates
[N,C]=size(xyz); ht=sqrt(3)/2; 
su=sum(xyz')'; for i=1:N xyz(i,:)=xyz(i,:)/su(i); end; %row sums=1
sx=(sqrt(3)/2)*xyz; % sx=xyz to make height=1 instead of sides
for i=1:N tx(i)=(sx(i,2)+2*sx(i,3))/sqrt(3); ty(i)=sx(i,2); end;
tx=tx'; ty=ty'; %convert to column vectors
end

%%%%%%%%%%%

function [pclr,cdiv,odfr,tdiv,cmplx,pxud,pxlr,asptr,inout,seq]=DoAdjacency(zones,pixsp,nrds,drw)
% [pclr,cdiv,odfr,tdiv,cmplx,pxud,pxlr,asptr,inout,X,P,seq]=DoAdjacency(zones,pixsp,nrds,drw);
%Calculates adjacency parameters on a label matrix
% This is a newer version which also gives the transition matrix
%INPUT: 
%  zones  matrix of image size, points are class numbers 1...nz.  No 0, -, NaN
%         NOTE: Assumes last zone number (nz) is the background, rest is animal
%         NOTE: zones should be 1,2,3,4...nz; it will insert zero rows 
%                 and columns if some zone numbers are missing
%  pixsp  pixel spacing, number of pixels between transects
%  nrds   number of randomizations for test of deviation from expected trans.
%  drw    not present or 0--do not draw, 1--draw permutated X histogram
%   [...]=DoAdjaceny(z,p,nrds) or (z,p,nrds,0) will not draw, (z,p,nrds,1) will
%   [...]=DoAdjacency(z,p) assumes nrds=20,000 and will not draw histogram
%OUTPUT:
% pclr      color fractions, excluding background
% cdiv      color diversity
% odfr      transition (off diagonal) frequencies
% tdiv      transition diversity
% cmplx     transitions/sample (animal) a measure of pattern diversity
% pxud,pxlr updown, leftright distance between transitions (patch size pixels)
% asptr     transition aspect ratio (larger if longer up-down)
% inout     fraction of transitions within animal
%             if >0.5 more transitons in animal than between animal & background
%             (complex animal on simpler background); <0.5 fancier background
% X,P       test for random adjacency pattern (sum abs difference & permute P)
% seq       transition matrix including background (first row & column)
%              diagonal and upper off-diagonal as sum of both off-diagonals
showhist=0; 
if nargin==4 && drw==1 showhist=1; end;
if nargin==3 drw=0; end;
if nargin==2 nrds=20000; drw=0; end;
z2=zones; [nbr,nbc]=size(z2); 
nz=max(max(z2)); zone=z2;
% nz = number of color classes, assume last is background 
%first get rid of the odd 0 in zones (unclassified)--assign to adjacent value
% John A. Endler, April 2012.  Please let me know of problems and errors
%    John.Endler@deakin.edu.au
%*********************************************************************************************
for g=1:3 %number of attempts to remove 0s
  %find unclassified pixels (0) in zone map
  z3=zone(:); [tsz,i]=size(z3); unc=uint8(zeros(tsz,1)); unc(z3==0)=1;
  unc=reshape(unc,nbr,nbc); %figure; imshow(unc,[]); %unc=1 at unclassed pixels
  cc=bwconncomp(unc); stats=regionprops(cc,'PixelList'); [nzr,i]=size(stats);
  clear z3 tsz i unc cc;
  zone2=zone;
    for zp=1:nzr
     pts=stats(zp).PixelList; px=pts(:,1); py=pts(:,2); [n,i]=size(px);
     for i=1:n
      x=px(i); y=py(i);
      if zone(y,x)<nz
        %7x7 pixel neighborhood centered on x,y [nx(5),ny(5)]
        nx=zeros(25,1); ny=nx;
        nx(1)=x-3;  nx(2)=x-2;  nx(3)=x-1;  nx(4)=x;  nx(5)=x+1;  nx(6)=x+2; 
        nx(7)=x+3;  ny(1:7)=y-3;
        nx(8)=x-3;  nx(9)=x-2;  nx(10)=x-1; nx(11)=x; nx(12)=x+1; nx(13)=x+2; 
        nx(14)=x+3; ny(8:14)=y-2;
        nx(15)=x-3; nx(16)=x-2; nx(17)=x-1; nx(18)=x; nx(19)=x+1; nx(20)=x+2; 
        nx(21)=x+3; ny(15:21)=y-1;
        nx(22)=x-3; nx(23)=x-2; nx(24)=x-1; nx(25)=x; nx(26)=x+1; nx(27)=x+2; 
        nx(28)=x+3; ny(22:28)=y;
        nx(29)=x-3; nx(30)=x-2; nx(31)=x-1; nx(32)=x; nx(33)=x+1; nx(34)=x+2; 
        nx(35)=x+3; ny(29:35)=y+1;
        nx(36)=x-3; nx(37)=x-2; nx(38)=x-1; nx(39)=x; nx(40)=x+1; nx(41)=x+2; 
        nx(42)=x+3; ny(36:42)=y+2;
        nx(43)=x-3; nx(44)=x-2; nx(45)=x-1; nx(46)=x; nx(47)=x+1; nx(48)=x+2; 
        nx(49)=x+3; ny(43:49)=y+3;
        nx(nx<1)=1; ny(ny<1)=1; nx(nx>nbc)=nbc; ny(ny>nbr)=nbr; %positions
        zn=zeros(9,1);
        for j=1:9
          zn(j)=zone(ny(j),nx(j));
        end;
        %new value from most common of neighbor labels
        [cts,edg]=DoHistogram(zn,-0.5,1,3.5,0);
        pcin=sum(cts)/49;
        if pcin>0.05
          idx=(0:3); mxh=max(cts); idl=idx(cts==mxh)'; idl=idl(idl>0);
          [nm,k]=size(idl);
          if nm>0 
            if nm==1 
              zone2(y,x)=idl;
            else %if modes are tied, as in line border
              pr=rand(1); 
              if pr>0.5 
                  zone2(y,x)=idl(1);
              else
                  zone2(y,x)=idl(2);
              end;
            end;
          end;
        end; %if <40% outside
      end; %ignore background
     end; %0 pixel n
    end; %0 pixel group zp;
    zone=zone2; clear zone2 idx mxh idl nm k pr;
end; %g cycles of filling in gaps
z2=zone; 
clear zone g zp nzr pts stats px py n i j x y nx ny zn cts edg pcin;
%*******************************************************************
[mx,my]=meshgrid(2:pixsp:nbc,2:pixsp:nbr);  %pixsp pixel grid
[nr,nc]=size(mx); %rows and columns of samples
%   px=mx(:); py=my(:); 
%   imshow(z2,[]); hold on; plot(px,py,'.'); clear px py; 
%   ttl=[ifn ', ' num2str(pixsp) ' pixel spacing'];
%   title(ttl); 
%take transects and accumulate patch code transition matrix
seq=zeros(nz,nz); nrm=nr-1; ncm=nc-1; 
trr=zeros(nr,1); trc=zeros(1,nc); 
ttrr=trr; ttrc=trc; %mean horizontal & vertical bw/wb transitions

%horizontal transects
for i=1:nr %rows
  t=0; tt=0;
  for j=1:ncm %columns
    p1=z2(my(i,j),mx(i,j)); p2=z2(my(i,j+1),mx(i,j+1));
    if isnan(p1) p1=z2(my(i,j-1),mx(i,j-1)); end;
    if isnan(p2) p2=z2(my(i,j+2),mx(i,j+2)); end;
    if ~isnan(p1) && ~isnan(p2) && p1>0 && p2>0
      seq(p1,p2)=seq(p1,p2)+1; 
      if p1~=p2 
        tt=tt+1;           %all transitions (to get start and end of animal)
        if p1<nz && p2<nz  %just in animal, exclude background (zone nz)
          t=t+1; 
        end;
      end; %transects within animal and to background
    end; %neither NaN
  end; %transect columns
  trr(i)=t; ttrr(i)=tt; %transitions per row (horizontal); animal and all
end; %transect rows

%vertical transects
tseq=seq;
seq=zeros(nz,nz);
for j=1:nc %columns
  t=0; tt=0;
  for i=1:nrm %rows
    p1=z2(my(i,j),mx(i,j)); p2=z2(my(i+1,j),mx(i+1,j));
    if isnan(p1) p1=z2(my(i-1,j),mx(i-1,j)); end;
    if isnan(p2) p2=z2(my(i+2,j),mx(i+2,j)); end;
    if ~isnan(p1) && ~isnan(p2) && p1>0 && p2>0 
      seq(p1,p2)=seq(p1,p2)+1;
      if p1~=p2 
        tt=tt+1;
        if p1<nz && p2<nz %just in animal, exclude background (zone nz)
          t=t+1;
        end;
      end;  %transects within animal and to background
    end; %neither NaN
  end;
  trc(j)=t; ttrc(j)=tt;
end;
seq=seq+tseq; clear tseq i j p1 p2 t tt ncm nrm;  %complete transition matrix
%group off-diagonals (doesn't matter which colour comes first
wseq=zeros(nz,nz); 
for i=1:nz
  wseq(i,i)=seq(i,i);
  for j=(i+1):nz
    wseq(i,j)=seq(i,j)+seq(j,i);
  end;
end;
seq=wseq; clear wseq;
%get lists of diagonal and off-diagonal numbers (including off animal edge)
diag=seq(1,1); for i=2:nz diag=[diag; seq(i,i)]; end; 
offdg=[];
for i=1:nz
  for j=(i+1):nz
    offdg=[offdg; seq(i,j)];
  end;
end; 
%get list of diags and off-diagonals within animal only
nza=nz-1;
adiag=seq(1,1); for i=2:nza adiag=[adiag; seq(i,i)]; end; 
aoffdg=[];
for i=1:nza
  for j=(i+1):nza
    aoffdg=[aoffdg; seq(i,j)];
  end;
end; 
clear i j mx my; 
%VARIABLES *****************************************************************
% z2             zone map
% nbr,nbc        number of rows and columns in image and zone map z2
% nz, nza        number of color classes with and without (nza) background
% seq(1:nz,1:nz) transition matrix for pixsp pixel spacing, upper off-dagonal for both
% diag,offdg     list of diagonals and off diagonals of seq
% adiag,aoffdg   list of diagonals and off-diags of seq excluding background
% pixsp          number of pixels between transects
% nr,nc          numbers of transect rows and columns
% trr(1:nr),trc(1:nc)   number of transitions per row or column animal only
% ttrr(1:nr),ttrc(1:nc) number of transitions per row or column animal+backgnd
%***************************************************************************
%get start and end rows and columns for body using ttrr>0 and ttrc>0
% use this to collect transition lists per row and per column (tr,tc)
tr=trr(ttrr>0); [ntr,j]=size(tr);  %animal transitions/row (horizontal left-right)
tc=trc(ttrc>0); [~,ntc]=size(tc);  %animal transitions/column (vertical up-down)
cmplx=(sum(tr)+sum(tc))/(ntr*ntc); %animal transitions/sample (complexity)
rp=1+pixsp*ntr; cp=1+pixsp*ntc; %size of animal in pixels
% aasp=rp/cp; %animal size aspect ratio (pixels)
% drow=mean(tr)/cp; dcol=mean(tc)/rp; % row and column transition densities
pxud=rp/mean(tc); pxlr=cp/mean(tr); %rougly average patch dimensions
  %pixels between transitions up-down or left-right
asptr=pxud/pxlr;  %aspect ratio of distance between transitions 
  %larger if pattern longer in up-down direction (elongated up-down)
inout=sum(aoffdg)/sum(offdg); %fraction transitions in animal vs to backgrnd
pclr=adiag/sum(adiag);       %animal color fractions
cdiv=SimpsonDiversity(pclr); %animal color diversity
odfr=aoffdg/sum(aoffdg);     %transition frequencies, animal only
tdiv=SimpsonDiversity(odfr); %transition diversity, animal only
%calculate Chi-squared for transition frequences
N=sum(aoffdg); %observed number of offdiags in animal only
etm=[]; c=0; 
for i=1:(nza-1)
  for j=(i+1):nza
    c=c+1; etm(c)=2*pclr(i)*pclr(j);
  end;
end;
etm=etm/sum(etm); %fractions expected off diagonals animal only
eod=(N*etm)'; clear i j k nzm; %expected off diagonal numbers
[nt,i]=size(aoffdg);
obsx=0;
for i=1:nt  
   obsx=obsx+abs(aoffdg(i)-eod(i));
end;
P=1;
% %now do the random subsampling of the off-diagonals from the 
% %  expecteds (etm), number of expected classes (nt) and numbers (N)
% cf=zeros(1,nt+1); sm=0; %make cumulative transition frequency vector
% for i=1:nt
%   sm=sm+etm(i); cf(i+1)=sm;
% end;
% cs=zeros(nrds,1);
% fprintf(1,'Doing randomization test  ');
% for r=1:nrds
%   rd=rand(N,1); ex=histc(rd,cf); ex(nt)=ex(nt)+ex(nt+1); %include =1
%   ex=ex(1:nt)'; %expected numbers  This will handle two cf with same value
%   cst=0;
%   for i=1:nt
%     cst=cst+abs(aoffdg(i)-ex(i));
%   end;
%   cs(r)=cst;
% end;
% if showhist==1
%   figure; hist(cs); my=max(hist(cs));
%   hold on; 
%   plot([obsx obsx],[0,my],'w','LineWidth',3);
%   plot([obsx obsx],[0,my],'k:');
%   fprintf(1,'obs %6.3f\n',obsx);
% end;
% X=obsx;
% P=sum(cs<obsx)/nrds; %randomized P of getting the observed X^2
% fprintf(1,' finished\n');
end

function [counts,edges]=DoHistogram(x,starth,width,endh,drawit,color)
% [counts,edges]=DoHistogram(x,starth,width,endh,drawit,color);
% function to draw histograms with defined bin limits
%INPUT:
%  x          the variable (vector) to be drawn as a histogram
%  starth     start of histogram range
%  width      width of each bar
%  endh       end of histogram range, may be reset by start and width
%  drawit     =1 place on a new figure; =0 do not draw
%  color      color: 'b' 'g' 'r' 'y' 'c' 'm'  as a character
%OUTPUT
%  counts     heights of bars (y axis)
%  edges      start and stop points for each bar (x axis)
ed=(starth:width:endh)'; edges=ed; [n,i]=size(ed); n=n-1;  
counts=zeros(n,1);
for i=1:n
  if i==1 counts(i)=sum(edges(i)<=x & x<=edges(i+1));
  else counts(i)=sum(edges(i)<x & x<=edges(i+1));
  end;
end;
if drawit==1
  hold on;
  for i=1:n
    if counts(i)>0 
        rectangle('Position',[edges(i),0,width,counts(i)],'FaceColor',color); 
    else
        plot([edges(i) edges(i+1)],[0 0],'k');
    end;
  end;
  hold off;
end;
end

function h=SimpsonDiversity(dvec)
% h=SimpsonDiversity(dvec);
%  Calculates Inverse Simpson's diversity on dvec (column vector)
%   which gives the equivalent number of equally common classes.
% dvec contains relative frequencies (sum to 1.0)
%   will convert a row vector to a column vector first
%   will ensure that sum of dvec = 1.
[~,i]=size(dvec); if i>1 dvec=dvec'; [~,i]=size(dvec); end;
tot=sum(dvec);
if tot>0
   dvec=dvec/tot; %must add to 1
   h=1/sum(dvec.^2);
else
   h=0;
end;
end

%%%%%%%%%%%
