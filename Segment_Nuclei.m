%% Nuclear segmentation demo

% This demonstration walks through a simple method for nuclei segmentation
% from a Whole Slide Image (WSI). Additional folders and files in this
% folder may be operating system dependent. This was originally written on
% a Windows computer, however, if openslide is not needed (for WSI formats
% like svs, tif

% Clean up command window, workspace, etc.
clc
close all
clear all


% Defining image name (works if the image is in the current working
% directory)
img_name = 'TY02 17309 TG26 VEH M PAS.svs';

% Since this is a large file, we want to pick a region of the image to read
% and apply our nuclear segmentation methods too.
% The bounding box here defines minimum x, maximum x, minimum y, and
% maximum y coordinates which define a rectangular region in the WSI.
bbox = [67960,68726,28348,29032];

% Since this image has an SVS file extension, we can natively read this
% region using Matlab's imread command.
image_region = imread(img_name,'Index',1,'PixelRegion',{bbox(3:4),bbox(1:2)});

% Showing the extracted image region
figure, imshow(image_region), axis image, title('Extracted Image Region')

% Using color deconvolution to separate stains
[h_stain,pas_stain,residual] = colour_deconvolution(image_region,'H PAS');

% Showing each of the stains separately
figure
subplot(1,3,1), imshow(imcomplement(h_stain)), axis image, title('Hematoxylin Stain Channel')
subplot(1,3,2), imshow(imcomplement(pas_stain)), axis image, title('PAS Stain Channel')
subplot(1,3,3), imshow(imcomplement(residual)), axis image, title('Residual (leftover)')
sgtitle('Output of Color Deconvolution')

% We can use the Hematoxylin channel to extract nuclei. This command
% thresholds the image by a value that is the result of a modified Otsu's
% threshold.
nuclei_mask = imbinarize(adapthisteq(imcomplement(h_stain)),2.3*graythresh(imcomplement(h_stain)));

% Showing results
figure, imshow(nuclei_mask), axis image, title('Initial Nuclei Mask')

% Initial results are pretty good, but we want to increase the specificity
% of our nuclear segmentations using other image analysis methods.

% Morphological closing (dilation followed by erosion), allows us to
% connect some of the broken up objects in the image.
morph_image = bwmorph(nuclei_mask,'close',Inf);

% Filling in holes in the nuclei
morph_image = imfill(morph_image,'holes');

% Removing small objects which aren't nuclei
morph_image = bwareaopen(morph_image,200);

% Showing the results
figure, imshow(morph_image), axis image, title('Output of Morphological Operations')

%% Now let's write the results to an XML file for viewing in Aperio ImageScope.

% Aperio uses XML formatted text files for storing annotations. We can
% write that in Matlab using structs followed by a struct2xml conversion.

% Iterating through each nucleus in the mask by first defining an index for
% each.
labeled_image = bwlabel(morph_image);

% The number of nuclei in the image is the maximum label index.
n_nuclei = max(labeled_image,[],'all');

% Initializing cell array for storing nuclei regions.
nuc_regions = {};

for i = 1:n_nuclei
    
    % Individual nucleus mask
    individual_image = labeled_image==i;

    % Getting the coordinates of the boundary
    [y_coords,x_coords] = find(bwperim(individual_image));
    
    % Modifying these coords so they are in the WSI coordinate system as
    % opposed to the extracted image region. This is done by just adding
    % the minimum x and y coordinate used to extract the region.
    y_coords = y_coords+bbox(3);
    x_coords = x_coords+bbox(1);

    % Adding vertex coordinates to a struct.
    new_region = [];
    for j = 1:length(x_coords)

        new_region.Vertices.Vertex{1,j}.Attributes.X = num2str(round(x_coords(j)));
        new_region.Vertices.Vertex{1,j}.Attributes.Y = num2str(round(y_coords(j)));

    end
    % Adding a region to the XML file.
    nuc_regions.Region{i} = new_region;

end

% Adding all of this to the main struct
annotations_xml.Annotations.Annotation.Regions = nuc_regions;

% Defining XML file path
xml_save_path = strrep(img_name,'svs','xml');

% Converting struct and saving as xml
struct2xml(annotations_xml,xml_save_path)

