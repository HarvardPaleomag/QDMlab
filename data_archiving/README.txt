# Archiving QDM data in MagIC


## Notes on conventions for storing QDM data in MagIC format

First, a few notes about conventions for storing QDM data in MagIC. Although MagIC's organization of data into sites, samples, and specimens is straightforward for conventional paleomagnetic core samples, QDM analyses characterize magnetization at a sub-core scale. 

For analyses where individual source net magnetic moments are reported, each source is a specimen and each inversion at a treatment step is a measurement. Typically, each thick section used for mapping is considered a sample and the bulk rock (or meteorite chip/subsample) from which the thick section is produced is considered a site.

Compared to a typical 2G measurement each QDM measurement should also use the columns “meas_pos_x”, “meas_pos_y”, “inversion_height”, “inversion_residuals”, “QDM_GF_factor”, and “inversion_type” to describe the fitting. The columns “QDM_context_image”, “QDM_B111_map”, “QDM_Bz_map” (optional), and “QDM_fitted_map” (optional) should be used to name the appropriate raw data maps. B111 and Bz maps should be uploaded in condensed MagIC format while context visible light maps and fitted/cropped maps should be uploaded as jpegs. QDM maps without net magnetic moment analyses can be uploaded as a set of “QDM_context_image”, “QDM_B111_map”, and “QDM_Bz_map.”

## Converting B111 and Bz data into condensed MagIC format

For simple experiments, or experiments where no net moment analysis is performed, the script Bfile2MagIC can be used to convert individual magnetic field maps into condensed MagIC format. (Note: by default, this script will not produce a header for the condensed MagIC file--a pre-made header must be explicitly provided to the 'header' keyword argument or added after the fact.)

## Process for converting QDM experiments into MagIC format

There are four steps to preparing QDM data to add to MagIC:
1.	Converting QDM maps into condensed MagIC format
2.	Building the measurements table for any relevant fit data
3.	Collecting any additional context images/fitted maps for upload
4.	Assembling the final contribution file

When fitting is involved, it is typically easiest to perform the first three of these tasks simultaneously, processing one experiment and one specimen at a time.

The script results2MagIC can be used to automate much of this process by using the data stored in the results structures produced by QDMLab. By default, it writes the fit results in a results structure into the MagIC contribution file, converts the Bz and B111 maps into condensed MagIC format, makes a copy of the context image, and saves the context image and the converted Bz and B111 maps to a specified location. If desired, cropped/fitted maps can also be produced from the information in the results structure and saved for MagIC upload.

results2MagIC takes as input:
1.	An ordered list of file paths to data for the experiment
2.	File paths for produced MagIC data files and measurements table
3.	A fit results file (output of fit_source_series)
4.	Output preferences (as keyword arguments) 
5.	Information about the experiment
	•	Experiment name
	•	Specimen name
	•	Method codes
	•	Measurement orientation
	•	Applied AC field, temperature, etc. for each file
	•	Software used
	•	Analysts
	•	… many other measurement fields in the MagIC data model!

(Note: Most measurement fields in the data model that are likely to be relevant to QDM data are included already as keyword arguments. It is straightforward to modify the supplied script to include additional measurement fields, if necessary, or manually add additional columns to the measurements table.)

Running the script once will add lines to the measurements table for each step of the experiment for the supplied source. To assemble the full measurements table, one need only repeat this process while varying the source names and experiment details. For well-organized data, this process is typically simple to automate.
