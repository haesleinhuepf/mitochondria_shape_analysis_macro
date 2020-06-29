// measure_and_count.ijm
// 
// This script thresholds bright objects on dark background and measures shape 
// of objects with a area above a given value.
// Furthermore, the script runs over a folder, its subfolders and analyses all images.
// If it finds an image with multiple z-planes, it applies a maximum projection before
// analysing the image.
// 
// Author: Robert Haase, rhaase@mpi-cbg.de.
// June 2020
// License: BSD3
// 
// Copyright 2020 Robert Haase, Max Planck Institute for Molecular Cell Biology and Genetics Dresden
// 
// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
// 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
// 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//

// Set a threshold for size:
minimum_size_in_squared_microns = 1.0;

// Where are the images located?
foldername = "/path/to/image/files/";
output_csv = "/path/to/result/files/output.csv";

if (File.exists(output_csv)) {
	File.delete(output_csv);
}
File.append("filename,aspect_ratio,circularity,roundness,solidity,count,number_of_nuclei", output_csv);


// go through the given folder
analyse_folder(foldername);


function analyse_folder(foldername) {
	mainfolder = getFileList(foldername);
	
	for (i = 0; i < lengthOf(mainfolder); i++) {
		filename = foldername + mainfolder[i];
		if (File.isDirectory(filename)) {
			analyse_folder(filename + "/");
		} else if (endsWith(filename, ".lsm")) {
			analyse_image(filename);
		}
	}
}




function analyse_image(filename) {
	roiManager("reset");

	open(filename);

	getDimensions(width, height, channels, slices, frames);
	if (slices > 1) {
		run("Z Project...", "projection=[Max Intensity]");
	}

	// go through all annotated nuclei and multiply pixels with 0
	number_of_nuclei = 0;
	multi_roi_filename = replace(filename, ".lsm", "_nucl.lsm.zip");
	if (File.exists(multi_roi_filename)) {
		roi_filename = multi_roi_filename;
	} else {
		roi_filename = replace(filename, ".lsm", "_nucl.lsm.roi");
	}
	if (File.exists(roi_filename)) {
		roiManager("open", roi_filename);
		
		for (r = 0; r < roiManager("count"); r++) {
			roiManager("select", r);
			run("Multiply...", "value=0");
			number_of_nuclei++;
		}
		run("Select None");
		roiManager("reset");
	}
	

	// denoise the image
	run("Median...", "radius=2 slice");
	run("Duplicate...", "duplicate channels=1");

	// threshold it
	setAutoThreshold("Moments dark");
	run("Convert to Mask");
	//run("Watershed");
	run("Erode");

	// measure
	roiManager("reset");
	run("Clear Results");
	run("Set Measurements...", "area shape redirect=None decimal=3");
	run("Analyze Particles...", "size=" + minimum_size_in_squared_microns + "-Infinity clear display add");
	
	roiManager("save", filename + "_roi.zip");

	// get mean shape descriptors over all objects in one image.
	// handle with care. These underlying measurements may not be 
	// normal distributed. For deeper studying shape of objects, it 
	// might make sense to count objects of a given shape instead 
	// of determining mean shape. See also: https://youtu.be/IaY2zNKrtYg?t=122  
	count = nResults();
	if (nResults() > 1) {
		run("Summarize");
		aspect_ratio = getResult("AR", nResults() - 4);
		circularity = getResult("Circ.", nResults() - 4);
		roundness = getResult("Round", nResults() - 4);
		solidity = getResult("Solidity", nResults() - 4);
	} else {
		aspect_ratio = getResult("AR", nResults() - 1);
		circularity = getResult("Circ.", nResults() - 1);
		roundness = getResult("Round", nResults() - 1);
		solidity = getResult("Solidity", nResults() - 1);
	}

	File.append(filename + "," + aspect_ratio + "," + circularity + "," + roundness + "," + solidity + "," + count + "," + number_of_nuclei, output_csv);
	
	// save results as image
	close();
	run("Duplicate...", "duplicate channels=1");
	run("Enhance Contrast", "saturated=0.35");
	for (r = 0; r < roiManager("count"); r++) {
		roiManager("select", r);
		Roi.setStrokeColor("magenta");
	}
	run("Select None");
	if (File.exists(roi_filename)) {
		roiManager("open", roi_filename);
		num_rois_before = roiManager("count");
		for (r = num_rois_before; r < roiManager("count"); r++) {
			roiManager("select", r);
			Roi.setStrokeColor("yellow");
		}
		run("Select None");
	}
	roiManager("show all with labels");
	run("Flatten");
	saveAs("PNG", filename + "_results.png");
	
	// clean up
	run("Close All");
}