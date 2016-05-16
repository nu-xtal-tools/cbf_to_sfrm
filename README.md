# cbf_to_sfrm 2016 (alpha version) by N.T.Johnson and M.R.Probert, Newcastle University, UK. :microscope:
Code to convert .cbf diffraction frames into .sfrm format. (Currently only supports data from Diamond Light Source, UK, Beamline i19, EH1 &amp; EH2)

Tested with: python (2.7.5), numpy (1.7.1) (on scientific linux).

#IMPORTANT:
- Bruker is not associated with this software and will not support this. (Please direct any queries to N.T.Johnson (N.Johnson5@ncl.ac.uk) or M.R.Probert (Michael.Probert@ncl.ac.uk).)
- Bruker reserves the right to make modifications to their frame format in the future which may affect this software's operation.
- Code is still in alpha version. Use entirely at your own risk!!!

##Notes:
- Currently frames are padded; masks for i19, eh1 and eh2 are available.
- Further modifications will be made to test a version with no padding for Apex3 (where square detectors are not required).

## Required files:
- `cbf_to_sfrm_2016.py`
- `pilatus_header_c2sedit.py` - *a minor edit of pilatus_header.py, provided with kind permission by _Marcus Mueller, Dectris Ltd._ *
- Default text file - *default files for i19-eh1/eh2 are included in repository*
- *Mask files are also currently advised for processing of data, as frames have been padded*

## To Run:
**Useage of this file:** cbf_to_sfrm_2016.py, cbf_file, name_string of cbf file, default file, output directory, number of runs, number of processes to be run, number of frames to add together into the final dataset

**cbf_file:** full file location of one of the files to be converted

**name_string:** file name up to (but not including) the run number

**default file:** text file containing information required by .sfrm format that is not present in .cbf header (as well as any detector offsets and the angle conversion required)

**output directory:** location of folder to which the frames will be saved (currently: previous directories have to be present, the code will make a folder in the directory above output directory given (unless it already exists)).

**number of runs:** number of runs to be converted (*currently: will start numbering from 1 - given number, can't only convert from 3-5 for example*)

**number of processes to be run:** code is parallelised. give number of processes you wish to set away.

**number of frames to add together in the final dataset:** the code can add up the images from a number of consequtive cbf frames into one sfrm frame. e.g. 1 for no addition, 5 for 5 cbf frames to be added to make 1 sfrm frame. *Note: this feature has not been extensively tested*

e.g. `python cbf_to_sfrm_2016-compress.py /home/Documents/cbfs/cbf_100K_02_00004.cbf cbf_100K_ default_file.txt /home/Documents/converted_cbf/cbf_100K 4 7 2`

for runs 1-4 of cbf_100K... to be converted with 7 processes being set away and two cbf frames being added together into one sfrm frame.

### Bugs to be fixed:
*This will be updated as bugs are reported*

### Updates:
*This will be updated as program is developed*
