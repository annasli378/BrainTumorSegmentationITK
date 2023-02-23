# BrainTumorSegmentationITK

A project aimed at semi-automatic segmentation of brain tumors on magnetic resonance images of the head

## Input Data
1. 3D Images as VTK - scaled to intensity 0-255
2. Coordinates (x, y, z) to 2 pixels in area of the tumor

## Command arguments:
Properties -> Debugging -> Command Arguments:
"pth_to_coordinates_txt_file" "path_to_img_vtk" "path_to_save_img_vtk"

You can also add those paths inside code.

## Project structure:
- src - folder with source code and cmake file
- dane - img data
- wyniki - folder for saving results

