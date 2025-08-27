/*
 * Macro template to process multiple images in a folder
 */

#@ File (label = "Input directory (mask_structures)", style = "directory") input
#@ File (label = "Output directory (mask_structures_smoothed)", style = "directory") output_mask_smoothed
#@ File (label = "Output directory (mesh_structures)", style = "directory") output_mesh
//#@ String (label = "File suffix", value = ".tif") suffix

// See also Process_Folder.py for a version of this code
// in the Python scripting language.

processFolder(input);

// function to scan folders/subfolders/files to find files with correct suffix
function processFolder(input) {
	list = getFileList(input);
	list = Array.sort(list);
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + File.separator + list[i]))
			processFolder(input + File.separator + list[i]);
		//if(endsWith(list[i], suffix))
		processFile(input, output_mask_smoothed, output_mesh, list[i]);
	}
}

function processFile(input, output_mask_smoothed, output_mesh, file) {
	// Do the processing here by adding your own code.
	// Leave the print statements until things work, then remove them.
	print("Processing: " + input + File.separator + file);
		
	open(input + File.separator + file);
	run("Smooth (3D)", "method=Gaussian sigma=1 use");
	
	//run("Options...", "iterations=1 count=1 black do=Nothing");
	//setBackgroundColor(0, 0, 0);
	//setOption("BlackBackground", true);
	//run("Make Binary", "method=Otsu background=Dark black");
	selectImage("Smoothed");
	run("Properties...", "unit=micron pixel_width=50 pixel_height=50 voxel_depth=50");
	
	filename_smoothed = output_mask_smoothed + File.separator + file + "_smoothed.tiff";
	saveAs("Tiff", filename_smoothed);
	
	filename_mesh = output_mesh + File.separator + file + "_smoothed.obj";
	
	commandLine = "stack=" + file + "_smoothed.tiff" + " threshold=50 resampling=1 red green blue save=" + filename_mesh;
	print(commandLine);
	
	run("Wavefront .OBJ ...", commandLine);
	//run("Wavefront .OBJ ...", "stack=sample3_anno_structures_X32Y32Z32_nearest2_smoothed.tiff threshold=50 resampling=4 red green blue save=/run/media/admin/LeGAO/Project_Human_Hypothalamus/task_0708_sample3_render_structures/sample3_anno_structures_X32Y32Z32_nearest2_smoothed.tiff.obj");
	
	close("*");
	
	print("Saving finished" + input);
}
