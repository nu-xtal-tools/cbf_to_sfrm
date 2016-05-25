### Notes ####
#
# cbf to sfrm - 2016 (alpha version) by N.T.Johnson and M.R.Probert, Newcastle University, UK.
# created for use with data from Diamond Light Source, UK, Beamline i19.
#
#IMPORTANT:
# - Bruker is not associated with this software and will not support this. (Please direct any queries to N.T.Johnson (N.Johnson5@ncl.ac.uk) or M.R.Probert (Michael.Probert@ncl.ac.uk))
# - Bruker reserves the right to make modifications to their frame format in the future which may affect this software's operation.
# - Code is still in alpha version. Use entirely at your own risk!!!
#
# - currently frames are padded, masks for i19 eh1 and eh2 are available.
# - further modifications will be made to test no padding for apex3.
#
# Angles (in default_file):
# - if kappa - angle conversion from EH2 kappa geometry to APEX conventions
# - if euleurian - angle conversion from EH1 angles to APEX conventions
# - no conversion, just takes values from cbf files as is (and adds in any offsets from default files)
#
# -Tested with: python (2.7.5), numpy (1.7.1)
#
#
#
#    'cbf_to_sfrm_2016' converts .cbf diffraction frames to .sfrm format.
#    Copyright (C) 2016, N.T.Johnson & M.R.Probert
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#


import datetime as datetime
import numpy as np
import sys, glob, os, multiprocessing, time
from struct import unpack
from pilatus_header_c2sedit import PilatusHeader #pilatus_header.py from DECTRIS (very slightly edited)

print datetime.datetime.now().time()

def Conversion(cbf_file, name_string, default_file, dir_out, runs, cores, compression):
        """ Entire conversion is run from this function, ensures it has the correct files, then converts them."""
        #if code is not working comment out multiprocessing and should see error messages
        global output_dir
        #global name_string
        from_default(default_file) #gets all values from default file
        file_list = get_file_list(cbf_file, name_string,runs, compression) #gets list of files to be accessed
        output_dir = make_dir(dir_out) # output_dir is the name of the directory the file will be output in
        ###multiprocessing###
        p = multiprocessing.Pool(cores)
        p.map(run_conversion, file_list, chunksize=1)
        ###
        ####regular processing###
        # if having problems with multiprocessing, comment out multiprocessing and uncomment regular processing below to help figure out problem
        #for i in file_list:
        #        print i
        #        run_conversion(i)
        #run_conversion(file_list)
        ####
        print datetime.datetime.now().time()
        #self.make_dm(self.some_marker) 
        #self.write_dm() #writes a mask file for cbf. only works for eh1 and eh2 of I19 (2015)

def cbf_name(starter,r_no, f_no): #list contains (starter,r_no, f_no)
        """creates cbf name from run and file number"""
        l = [starter, r_no, "_", f_no, ".cbf"]
        return "".join(l)

def get_file_list(cbf_file, name_string, runs, compression):
        """creates list of all frames to be converted"""
        global name
        name = name_string #name_string = everything before the run 
        end_index = cbf_file.rindex(name_string) #in
        starter = cbf_file[:len(name_string)+end_index]
        run_string = cbf_file[end_index+len(name_string):] #file location up to run number and file name
        run_strings = run_string.split(".")[-2]
        r = run_strings.split("_")
        run_no_z = len(r[0])
        frame_no_z = len(r[1])
        run_len = []
        #need number of frames in each run
        for i in range(1,runs+1):
            run_search = starter + str(i).rjust(run_no_z,"0") + "*"
            run_len.append(len(glob.glob(run_search)))
        file_list = []
        counter = 0
        for j in run_len:
            #j = number of runs in each run
            counter += 1
            for i in range(1,j+1,compression):
                l = [starter, str(counter).rjust(run_no_z,"0"), "_", str(i).rjust(frame_no_z,"0"), ".cbf"]
                file_list.append("".join(l))
                #print "".join(l)
        print " %d files found to convert" % len(file_list)
        return file_list

def make_dir(output_dir):
        """makes output directory, wheere the files will be saved."""
        #global dir_string
        top_dir = output_dir
        try:
            dir_string = top_dir + "/"
            os.chdir(dir_string)
        except:
            dir_string = top_dir + "/"
            os.mkdir(dir_string)
        os.chdir(top_dir)
        return dir_string

def kappa_to_euleurian(kappa, k_omega, k_phi, k_2theta, k_offset, omega_offset):
        """converts kappa geometry to euleurian geometry needed for .sfrm header - currently only designed for i19, EH2, Diamond Light Source"""
        deg2rad = float(np.pi)/180.0
        rad2deg = float(180.0/np.pi)

        rad_kappa = kappa * deg2rad
        rad_k_offset = k_offset * deg2rad
        rad_k_omega = k_omega * deg2rad
        rad_k_phi = k_phi * deg2rad

        delta = np.arctan(np.tan(rad_kappa/2) * np.cos(rad_k_offset))

        rad_e_omega = rad_k_omega + delta
        rad_e_chi = 2*np.arcsin(np.sin(rad_kappa/2) * np.sin(rad_k_offset))
        rad_e_phi = rad_k_phi + delta

        e_omega =  (omega_offset + (rad_e_omega * rad2deg) + 180 )
        e_chi = - (rad_e_chi * rad2deg)
        e_phi = -(rad_e_phi*rad2deg + 90) #returned to normal
        e_2theta = k_2theta

        return e_omega, e_chi, e_phi, e_2theta

def default_value(info, chars):
        """chars = maximum number of characters in line, if info is greater, only first/last chars will appear"""
        info = info.strip()
        value = info.split("=")[1]
        value = value.lstrip()
        value = value[:chars]
        return value
            
def from_default(default_file):
        """takes values from detector default file to use in the .sfrm header"""
        global default_vals
        open_def = open(default_file, "r")
        items = open_def.readlines()
        default_vals = {}
        default_keys = ["Data_type","Site","User_name","Sample_ID","Setname","User_comment_1","User_comment_2","User_comment_3","User_comment_4","User_comment_5","User_comment_6","Original_frame_filename",
                        "Flood_field_correction_filename","Brass_plate_correction_filename","Xray_target_material","Voltage","Current","Filter_monochromator_setting","Low_temp_flag","exp_temp_in_hundredths","ccd_temp_in_hundredths","X_c","Y_c",
                        "Mag","dX","dY","dDist","pitch","roll","yaw","Display_lookup","Display_limit1","Display_limit2","Det_type","Brass_hole_spacing","Phosphor_distance","Be_window_thickness","Det_4","Det_7",
                        "no_of_exposures","per_eposure_bias","baseline_offset","orientation","over_scan_flag","readnoise","e_ADU","e_photon","bias","full_scale","Chemical_formula","Crystal_morphology",
                        "Crystal_colour","Crystal_dimension_1","Crystal_dimension_2","Crystal_dimension_3","Density_measurement_method","gain","high_speed_time","scale","offset","auto_full_scale","mono_2theta",
                        "mono_roll","mono_tilt","attenuation_factor","kappa_offset","omega_offset","geometry","chi_offset"]

        for k in default_keys:
            for line in items:
                if k == "omega_offset":
                    if line.startswith(k):
                        offset = default_value(line, 72)
                        default_vals["omega_offset"] = float(offset)
                elif k == "kappa_offset" :
                    if line.startswith(k):
                        offset = default_value(line, 72)
                        default_vals["kap_offset"] = float(offset)
                elif k == "chi_offset":
                    if line.startswith(k):
                        offset = default_value(line, 72)
                        default_vals["chi_offset"] = float(offset)
                else:
                    if line.startswith(k):
                        default_vals[k] = default_value(line, 72)
           #return default_vals 

def run_conversion(file_name):
        """the function from which the conversion code is run"""
        global name
        global cbf_list
        try: 
                g = open(file_name,"r")
                g.close()
                opens = True
        except:
                print file_name
                opens = False
                print "not there"
        if opens == True:
                cbf_list = get_names(file_name, name) #gets list of cbf files to be compressed together
                img_string, over_1_string, over_2_string = from_cbf(file_name) #reads values from .cbf header (and assigns self.some_marker), then goes through unpack_cbf which unpacks the cbf file and converts it to bruker format, saves it as a
                sfrm_name = get_sfrm_name(file_name, name)#generates name for output in standard .sfrm convention
                header_string = compile_header(some_marker) #creates header string, some_marker - whether its 2048 or 1024.
                write_sfrm(header_string, img_string,  sfrm_name, over_1_string, over_2_string) #writes .sfrm file
        else:
                #print "file still not there"
                pass



def get_names(file_name, cbf_string):
        """creates list of all frames to be converted"""
        # assumes cbf format is textRun_frameno.cbf
        # name_starter is everything up to the run number, including the _
        end_index = file_name.rindex(name_string)
        starter = file_name[:len(name_string)+end_index]
        run_string = file_name[end_index+len(name_string):] #file location up to run number and file name
        run_strings = run_string.split(".")[-2]
        r = run_strings.split("_")
        run_no_z = len(r[0])
        frame_no_z = len(r[1])
        run = int(r[0])
        frame = int(r[1])
        cbf_list = []
        for i in range(frame,frame+compression):
                l = [starter, str(run).rjust(run_no_z,"0"), "_", str(i).rjust(frame_no_z,"0"), ".cbf"]
                cbf_list.append("".join(l))
        return cbf_list

def get_sfrm_name(file_name, cbf_string):
        """changes epseranto file name to sfrm file name; when mutliple frames are added up .sfrm files are numbered differently, need to change so apex2 can recognise the sequence of files"""
        global output_dir
        global header_vals
        global compression
        end_index = file_name.rindex(cbf_string)
        run_string = file_name[len(cbf_string)+end_index:]
        info = run_string.split(".")[0]
        run = info.split("_")[0]
        frame = info.split("_")[1]
        f_calc = (int(frame[1:])+ (compression-1))/compression #should calculate new frame number based on current frame number and compression
        f_no = str(f_calc)
        header_vals["run"] = run.lstrip("0")
        run = run.rjust(2,"0")
        header_vals["frame_no"] = f_no
        return output_dir + cbf_string.strip("_") + "_" + run + "_" + f_no.rjust(4,"0") + ".sfrm"

def cbf_matrix(filename, nvals, pointer, some_marker):
        """uncompresses cbf image and converts to .sfrm compatible format: byte_offset compression only"""
        n_count = 0
        bpv = 0
        ints = []
        fil = open(filename,"rb")
        c = fil.read()
        while n_count < nvals:
            val = unpack("b",c[pointer:pointer+1])
            if val[0] > -128:
                bpv += val[0]
                ints.append(bpv)
                pointer += 1
                n_count += 1
            else:
                val = unpack("h",c[pointer+1:pointer+3])
                if  val[0] > -32768:
                    bpv += val[0]
                    ints.append(bpv)
                    pointer += 3
                    n_count += 1
                else:
                    val = unpack("i",c[pointer+3:pointer+7])
                    if  val[0] > -2147483648: #changed these as it shouldn't hit + 2147483648, same with lower vals.
                        bpv += val[0]
                        ints.append(bpv)
                        pointer += 7
                        n_count += 1
                    else:
                        print "help, very high numbers"
                        sys.exit()
        c = ""
        img = np.array(ints)
        reset = np.where(img < 0)
        img[reset] = 0
        return img
        

def padd_final_matrix(img, some_marker):
        """creates the image binary (in .sfrm compression) and padds matrix to square. However, Apex3 may not require square detectors, so this may be unneccesary.
            Further information on format of image binary can be found in SAINT manual from BRUKER AXS."""
        global header_vals
        max_tuple = np.unravel_index(img.argmax(),img.shape)
        if img[max_tuple] < 255: #no overflow table needed
                dt = np.uint8
                img = img.astype(dt)
                img_string = img.tostring(dt)
                header_vals["count_1"] = 0
                header_vals["count_2"] = 0
                over_1_string = ""
                over_2_string = ""
        elif 255 <= img[max_tuple] <65535: # only need 1 overflow table
                lst = np.where(img>254)
                over_1 = []
                for i in range(0,len(lst[0])):
                        over_1.append(img[lst[0][i]]) #create list of numbers that are 255 and above
                img[lst] = 255 # set values above 255 to 255
                max_tuple = np.unravel_index(img.argmax(),img.shape)
                ### overflow table now ###
                count_1 = len(over_1)
                header_vals["count_1"] = count_1
                if count_1%8 != 0: #table length needs to be padded up to multiples of 16 bytes
                    padd = 8 - (count_1%8)
                    while padd > 0:
                        over_1.append(0)
                        padd -= 1
                ddt = np.uint16 # for overflow table 1 (16 bit/2 byte binary)           
                over_string = np.array(over_1)
                over_1 = over_string.astype(ddt)
                over_1_string = over_1.tostring(ddt)
                count_2 = 0 # no second overflow table needed
                over_2_string = ""
                header_vals["count_2"] = count_2
                #sys.exit()
        elif 65535 <= img[max_tuple] < 4294967295: #need 2 overflow tables
                lst = np.where(img>254)
                over_1 = []
                for i in range(0,len(lst[0])):
                        over_1.append(img[lst[0][i]]) #create list of numbers that are 255 and above
                img[lst] = 255 # set values above 255 to 255
                count_1 = len(over_1)
                over_1 = np.array(over_1)
                lst_over = np.where(over_1>65534)
                over_2 = []
                for i in range(0,len(lst_over[0])):
                        over_2.append(over_1[lst_over[0][i]]) #create list of numbers that are 255 and above
                over_1[lst_over[0]] = 65535 # set values above 65535 to 65535
                ### overflow table 1 now ###
                header_vals["count_1"] = count_1
                padd_1 = []
                if count_1%8 != 0: #table length needs to be padded up to multiples of 16 bytes, each number in table == 2 bytes
                    padd = 8 - (count_1%8)
                    while padd > 0:
                        padd_1.append(0)
                        padd -= 1
                ddt = np.uint16 # for overflow table 1 (16 bit/2 byte binary)           
                padd_a = np.array(padd_1)
                overs = np.append(over_1, padd_a)
                over_string = np.array(overs)
                over_1 = over_string.astype(ddt)
                over_1_string = over_1.tostring(ddt)
                ### overflow table 2 now ###
                count_2 = len(over_2)
                header_vals["count_2"] = count_2
                if count_2%8 != 0: #table length needs to be padded up to multiples of 16 bytes
                    padd = 4 - (count_2%4)
                    while padd > 0:
                        over_2.append(0)
                        padd -= 1
                dddt = np.uint32 # for overflow table 1 (16 bit/2 byte binary)           
                over_string = np.array(over_2)
                over_2 = over_string.astype(dddt)
                over_2_string = over_2.tostring(dddt)
                len(over_1_string)
                len(over_2_string)
        xsize = header_vals["xsize"]
        ysize = header_vals["ysize"]
        ints = np.resize(img, (xsize,ysize))
        img = ""
        threshold_list = np.where(ints>64000)
        header_vals['nover64'] = len(threshold_list)
        header_vals["min"] =np.amin(ints)
        header_vals["max"] =np.amax(ints)
        header_vals["ncounts"] = str(ints.sum())
        if some_marker == 1024:
            padding_x = 1024-xsize #padding to 1024 x 1024 size
            padding_y = xsize-ysize # to make square matrix (padded in y first, then both x and y padded again)
        elif some_marker == 2048:
            padding_x = 2048-xsize
            padding_y = xsize-ysize # to make square matrix (padded in y first, then both x and y padded again)
        img = np.pad(ints, [(0, 0),(padding_y,0)], mode = "constant", constant_values=0)
        img = np.transpose(img)
        img = np.pad(img,[(padding_x, 0),(padding_x,0)], mode = "constant", constant_values=0)
        ##making image here
        dt = np.uint8
        img = img.astype(dt)
        img_string = img.tostring(dt)
        ####
        max_tuple = np.unravel_index(img.argmax(),img.shape)
        header_vals["max_col"] = str(max_tuple[1]) + ".000000" # for use when writing frame header
        header_vals["max_row"] = str(max_tuple[0]) + ".000000"# for use when writing frame header
        header_vals["nrows"] =len(img[0]) # for use when writing frame header
        header_vals["ncols"] =len(img) # for use when writing frame header
        header_vals["npixelb"] = 1 # for use when writing frame header
        header_vals["padding_x"] = some_marker-xsize #for use when working out frame beam centre
        header_vals["padding_y"] = xsize-ysize #for use when working out frame beam centre
        img = ""
        return img_string, over_1_string, over_2_string

def from_cbf(file_name):
        """takes values from the .cbf header to be used in when writing sfrm header, uses slight modification of dectris_header.py"""
        global cbf_matrix
        global cbf_list
        global some_marker
        global header_vals
        header = PilatusHeader(file_name)
        header_vals = header.header_dict
        header_vals["hist_string"] =  header.get_date_time()
        op = open(file_name,"r")
        d = op.readlines()
        indicies = [i for i, s in enumerate(d) if "--CIF-BINARY-FORMAT-SECTION--" in s] #pilatus_header.py does not take the binary format from the cbf header
        compression_dic_old = {"     conversions" : "?",
                               "Content-Type:" : "?",
                               "X-Binary-Size:" : 0,
                               "X-Binary-ID:" : 0,
                               "X-Binary-Element-Type:" : "?",
                               "X-Binary-Element-Byte-Order:" : "?",
                               "X-Binary-Number-of-Elements:" : 0,
                               "X-Binary-Size-Fastest-Dimension:" :0 ,
                               "X-Binary-Size-Second-Dimension:" : 0,
                               "X-Binary-Size-Padding:" : 0,
                               }
        compression_dic = {}
        for line in d:
            for key in compression_dic_old:
                if line.startswith(key):
                    if key == "     conversions":
                        l = line.strip()
                        l = l.split("=")
                        compression_dic[l[0]] = l[1]
                    else:
                        l = line.strip()
                        l = l.split()
                        compression_dic[l[0]] = l[1]
        op.close()
        if compression_dic["conversions"] != '"x-CBF_BYTE_OFFSET"':
            print "compression unsupported. programme exiting"
            sys.exit()
        else:
            #print "hi"
            op = open(file_name,"rb")
            c = op.read()
            start_bin = "\x0c\x1a\x04\xd5"
            pointer = c.index(start_bin) + 4
            nvals = int(compression_dic["X-Binary-Number-of-Elements:"])
            xsize = int(compression_dic["X-Binary-Size-Second-Dimension:"])
            header_vals["xsize"] = xsize
            ysize = int(compression_dic["X-Binary-Size-Fastest-Dimension:"])
            header_vals["ysize"] = ysize           
            op.close()
            if xsize > 1024 or ysize > 1024:
                some_marker = 2048
                some_marker = 2048
            else:
                some_marker = 1024
                some_marker = 1024
            #this goes into sizing matrix
            img = np.zeros((xsize*ysize),dtype=np.uint32)
            #print "hi"
            for item in cbf_list:
                im = cbf_matrix(item,nvals, pointer, some_marker) # uses self.pointer to iterate through binary, uncompressing the image
                imgs = np.add(img,im)
                img = imgs
        reset = np.where(img < 0)
        img[reset] = 0
        image, over_1_string, over_2_string = padd_final_matrix(img,some_marker) #img string
        return image, over_1_string, over_2_string

def make_dm(self,some_marker):
        """makes dynamic masks for use in integration. currently not used by code, but retained for completeness"""
        if some_marker == 1024:
            dm_mat = np.ones((619,487),dtype=np.uint32)
            for j in range(0,487): #adds in the blank lines in dynamic mask
                for i in range(195,212):
                    dm_mat[i,j] = 0
                for i in range(407,424):
                    dm_mat[i,j] = 0
            reset = np.where(dm_mat>0)
            dm_mat[reset]=64
            to_square = 619 - 487
            padding = 1024 - 619
            dm_mat = np.pad(dm_mat,[(0,0),(to_square,0)],mode="constant",constant_values=0)
            dm_mat = np.transpose(dm_mat)
            dm_mat = np.pad(dm_mat,[(padding,0),(padding,0)],mode="constant",constant_values=0)
            dmt = np.uint8
            dm = dm_mat.astype(dmt)
            self.dm = dm.tostring(dmt)
        elif some_marker == 2048:
            dm_mat = np.ones((1679,1475),dtype=np.uint32)
            for j in range(0,1475): #adds in the blank horizontal lines in dynamic mask
                for i in range(195,212):
                    dm_mat[i,j] = 0
                for i in range(407,424):
                    dm_mat[i,j] = 0
                for i in range(619,636):
                    dm_mat[i,j] = 0
                for i in range(831,848):
                    dm_mat[i,j] = 0
                for i in range(1043,1060):
                    dm_mat[i,j] = 0
                for i in range(1255,1272):
                    dm_mat[i,j] = 0
                for i in range(1467,1484):
                    dm_mat[i,j] = 0
            for i in range(0,1679): #vertical lines, therefore every thing in row will be 0
                for j in range(487,494):
                    dm_mat[i,j] = 0
                for j in range(981,988):
                    dm_mat[i,j] = 0
            reset = np.where(dm_mat>0)
            dm_mat[reset]=64
            to_square = 1679 - 1475
            padding = 2048 - 1679
            dm_mat = np.pad(dm_mat,[(0,0),(to_square,0)],mode="constant",constant_values=0)
            dm_mat = np.transpose(dm_mat)
            dm_mat = np.pad(dm_mat,[(padding,0),(padding,0)],mode="constant",constant_values=0)
            dmt = np.uint8
            dm = dm_mat.astype(dmt)
            self.dm = dm.tostring(dmt)
            


def write_sfrm(header_string, img_string,sfrm_name, over_1_string, over_2_string):
        """writes .sfrm file"""
        global header_vals
        #dm making currently commented out - no longer needed, but kept for completion
        output_name = sfrm_name
        #dm_name = self.dm_name
        output = open(output_name,'wb')
        #dm_out = open(dm_name,"wb")
        sfrm_output = header_string + img_string
        #sfrm_dm_output = self.dm_header + self.dm
        #handling overflows, if they are present
        if header_vals["count_1"] != 0:
            sfrm_output += over_1_string
            #print "adding string one"
        if header_vals["count_2"] != 0: #if no overflows, no table
            sfrm_output += over_2_string
            #print "adding string two"
        output.write(sfrm_output)
        #dm_out.write(sfrm_dm_output)
        output.close()
        #dm_out.close()

def write_dm(self):
        """No longer needed, as can use static mask when processing data. Retained for completion. Writes mask file - currently only done using the sfrm from the very last frame converted"""
        dm_name = self.dm_name
        dm_out = open(dm_name,"wb")
        sfrm_dm_output = self.dm_header + self.dm
        dm_out.write(sfrm_dm_output)
        dm_out.close()
        
def compile_header(some_marker):
        """creates string with all header items in it;
           .sfrm header; 80 byte lines =  8 byte description + 72 bytes for info
           header block - 512 bytes (blocks present in multiples of 5)

           Some information required by .sfrm is not included in .cbf file, therefore a defaults file is used to supply additional information.
           Documentation of header information required in .sfrm format can be found from Bruker corporation"""
        global default_vals
        global header_vals
        global compression
        header_list = []
        header_list.append("FORMAT :" + "100".ljust(72)) #code works for 100 format, not 85
        header_list.append("VERSION:" + "15".ljust(72))
        header_list.append("HDRBLKS:" + "15".ljust(72)) #header size in 512 byte blocks
        header_list.append("TYPE   :" + str(default_vals["Data_type"]).ljust(72)) #all default vals items are imported from the default values text file
        header_list.append("SITE   :" + str(default_vals["Site"]).ljust(72))
        header_list.append("MODEL  :" + str(header_vals["Detector_identifier"]).ljust(72))
        header_list.append("USER   :" + str(default_vals["User_name"]).ljust(72))
        header_list.append("SAMPLE :" + str(default_vals["Sample_ID"]).ljust(72))
        header_list.append("SETNAME:" + str(default_vals["Setname"]).ljust(72))
        header_list.append("RUN    :" + str(header_vals["run"]).ljust(72)) #run number within the data set
        header_list.append("SAMPNUM:" + "0".ljust(72)) #specimen number within data set
        header_list.append("TITLE  :" + "Converted from cbf by cbf_to_sfrm.py (N Johnson M Probert ALPHA VERSION)".ljust(72)) #user comments (8 lines)
        original_string = "original: %s %s" % (header_vals["Detector_identifier"], header_vals["hist_string"])
        header_list.append("TITLE  :" + original_string.ljust(72))
        header_list.append("TITLE  :" + str(default_vals["User_comment_1"]).ljust(72))#user comments
        header_list.append("TITLE  :" + str(default_vals["User_comment_2"]).ljust(72))#user comments
        header_list.append("TITLE  :" + str(default_vals["User_comment_3"]).ljust(72))#user comments
        header_list.append("TITLE  :" + str(default_vals["User_comment_4"]).ljust(72))#user comments
        header_list.append("TITLE  :" + str(default_vals["User_comment_5"]).ljust(72))#user comments
        header_list.append("TITLE  :" + str(default_vals["User_comment_6"]).ljust(72))#user comments
        header_list.append("NCOUNTS:" + str(header_vals["ncounts"]).ljust(36) + "0".ljust(36)) #total number of counts in frame, number of counts in underflow
        header_list.append("NOVERFL:" + "-1".ljust(24) +str(header_vals["count_1"]).ljust(24) + str(header_vals["count_2"]).ljust(24))#number of underflows, number in first overflow, number in second overflow
        header_list.append("MINIMUM:" + str(header_vals["min"]).ljust(72))#minimum counts in a pixel
        header_list.append("MAXIMUM:" + str(header_vals["max"]).ljust(72))#maximum counts in a pixel
        header_list.append("NONTIME:" + "-1".ljust(72))
        header_list.append("NLATE  :" + "0".ljust(72))
        header_list.append("FILENAM:" + str(default_vals["Original_frame_filename"]).ljust(72))
        #creating date-time string
        dt = datetime.datetime.today()
        creates = dt.timetuple() # 0 = year, 1 = month, 2 = day, hour = 3, min = 4, sec = 5
        dstring1 = "%s/%s/%s" %(creates[1],creates[2],creates[0])
        dstring2 = "%s:%s:%s" %(creates[3],creates[4],creates[5])
        date_string = dstring1.ljust(36) + dstring2.ljust(36)
        header_list.append("CREATED:" + date_string.ljust(72)) #date and time file created
        #exposure time calculator
        t = float(header_vals["Exposure_time"])*compression
        t = str(t).ljust(9,"0")
        t = round(float(t),6)

        t2 = float(header_vals["Exposure_period"])*compression
        t2 = str(t2).ljust(9,"0")
        t2 = round(float(t2),6)
        
        header_list.append("CUMULAT:" + str(t).ljust(72))#accumulated exposure time in seconds
        header_list.append("ELAPSDR:" + str(t2).ljust(72))#requested time to last exposure (s)
        header_list.append("ELAPSDA:" + str(t).ljust(72))#actual time for last exposure(s)
        header_list.append("OSCILLA:" + "0".ljust(72)) #nonzero if acquired by oscillation
        header_list.append("NSTEPS :" + "1".ljust(72)) #number of steps or oscillations in this frame
        #calculating angles - need to be converted from kappa to eulerian geometry if required
        axis_assigned = False
        oscillation_axis = header_vals["Oscillation_axis"]
        if oscillation_axis == "X CW" : #currently what i19 EH1 has set as their oscillation axis
            #assuming it is either a fixed phi or fixed omega rotation
            p = float(header_vals["Phi_increment"])
            o = float(header_vals["Omega_increment"])
            if p > 0:
                #print "phi"
                oscillation_axis = "PHI"
            elif o > 0:
                #print "omega"
                oscillation_axis = "OMEGA"
            if p > 0 and o > 0:
                print "oh dear"
                print header_vals["Oscillation_axis"]
                print "oscillation axis currently not supported"
                sys.exit()
            #sys.exit()
        #This angle conversion is currently only set up to handle detectors at Beamline i19 at Diamond Light Source (EH1 and EH2)
        if oscillation_axis == "OMEGA":
            axis_val = "2"
            if default_vals["geometry"] == "kappa":
            
                kappa = float(header_vals["Kappa"])
                k_omega = float(header_vals["Start_angle"])
                k_phi = float(header_vals["Phi"])
                k_2theta = float(header_vals["Detector_2theta"])
                es_omega, es_chi, es_phi, es_2theta = kappa_to_euleurian(kappa, k_omega, k_phi, k_2theta, default_vals["kap_offset"], default_vals["omega_offset"])

                ke_omega = k_omega + float(header_vals["Angle_increment"])*compression
                ke_phi = k_phi + float(header_vals["Phi_increment"])
                kappa_e = kappa + float(header_vals["Chi_increment"])
                #assumes 2theta doesn't change
                ee_omega, ee_chi, ee_phi, ee_2theta = kappa_to_euleurian(kappa_e, ke_omega, ke_phi, k_2theta,default_vals["kap_offset"], default_vals["omega_offset"])
                dif = np.absolute(float(header_vals["Angle_increment"]))*compression
                range_ = ee_omega - es_omega
            elif default_vals["geometry"] == "eulerian":
                es_omega = -(float(header_vals["Omega"])+ default_vals["omega_offset"])
                es_phi = float(header_vals["Phi"])
                es_2theta = -(float(header_vals["Detector_2theta"]))
                es_chi = -(float(header_vals["Chi"]) + default_vals["chi_offset"])
                #
                ee_2theta = -(float(header_vals["Detector_2theta"]))
                ee_omega = es_omega - float(header_vals["Omega_increment"])*compression
                ee_phi = (es_phi +float(header_vals["Phi_increment"]))
                ee_chi = es_chi + float(header_vals["Chi_increment"])
                dif = np.absolute(float(header_vals["Angle_increment"]))*compression
                range_ = ee_omega - es_omega
            elif default_vals["geometry"] == "no_convert":
                es_omega = float(header_vals["Omega"])+ default_vals["omega_offset"]
                es_phi = float(header_vals["Phi"])
                es_2theta = float(header_vals["Detector_2theta"])
                es_chi = float(header_vals["Chi"]) + default_vals["chi_offset"]
                #
                ee_2theta = float(header_vals["Detector_2theta"])
                ee_omega = es_omega + float(header_vals["Omega_increment"])*compression
                ee_phi = (es_phi + float(header_vals["Phi_increment"]))
                ee_chi = es_chi + float(header_vals["Chi_increment"])
                dif = np.absolute(float(header_vals["Angle_increment"]))*compression
                range_ = ee_omega - es_omega
            header_list.append("RANGE  :" + str(round(np.absolute(range_),6)).ljust(72)) #unsigned scan range in decimal degrees
            header_list.append("START  :" + str(round(es_omega,6)).ljust(72)) #start scan angle value, decimal degrees
            header_list.append("INCREME:" + str(round(range_,6)).ljust(72)) #signed scan angle increment between frames

        elif oscillation_axis == "PHI":
            axis_val = "3"
            if default_vals["geometry"] == "kappa":
                kappa = float(header_vals["Kappa"])
                k_phi = float(header_vals["Phi"])
                k_omega = float(header_vals["Omega"])
                k_2theta = float(header_vals["Detector_2theta"])
                es_omega, es_chi, es_phi, es_2theta = kappa_to_euleurian(kappa, k_omega, k_phi, k_2theta, default_vals["kap_offset"], default_vals["omega_offset"])

                ke_omega = k_omega + float(header_vals["Omega_increment"])
                ke_phi = k_phi + float(header_vals["Phi_increment"])*compression
                kappa_e = kappa + float(header_vals["Chi_increment"])
                #assumes 2theta doesn't change
                ee_omega, ee_chi, ee_phi, ee_2theta = kappa_to_euleurian(kappa_e, ke_omega, ke_phi, k_2theta,default_vals["kap_offset"], default_vals["omega_offset"])
                dif = np.absolute(float(header_vals["Angle_increment"]))*compression
                range_ = ee_phi - es_phi
            elif default_vals["geometry"] == "eulerian":
                es_omega = -(float(header_vals["Omega"])+ default_vals["omega_offset"])
                es_phi = (float(header_vals["Phi"]))
                es_2theta = -(float(header_vals["Detector_2theta"]))
                es_chi = -(float(header_vals["Chi"]) + default_vals["chi_offset"])  # 
                #
                ee_2theta = -(float(header_vals["Detector_2theta"]))
                ee_omega = es_omega - float(header_vals["Omega_increment"])
                ee_phi = es_phi + float(header_vals["Phi_increment"]*compression)
                ee_chi = es_chi - float(header_vals["Chi_increment"])
                dif = np.absolute(float(header_vals["Angle_increment"]))*compression
                #sys.exit()
                range_ = ee_phi - es_phi
            elif default_vals["geometry"] == "no_convert":
                es_omega = -(float(header_vals["Omega"])+ default_vals["omega_offset"])
                es_phi = float(header_vals["Phi"])
                es_2theta = -(float(header_vals["Detector_2theta"]))
                es_chi = float(header_vals["Chi"]) + default_vals["chi_offset"]  # 
                #
                ee_2theta = float(header_vals["Detector_2theta"])
                ee_omega = es_omega - float(header_vals["Omega_increment"])
                ee_phi = es_phi + float(header_vals["Phi_increment"]*compression)
                ee_chi = es_chi + float(header_vals["Chi_increment"])
                dif = np.absolute(float(header_vals["Angle_increment"]))*compression
                range_ = ee_phi - es_phi
            header_list.append("RANGE  :" + str(round(np.absolute(range_),6)).ljust(72)) #unsigned scan range in decimal degrees
            header_list.append("START  :" + str(round(es_phi,6)).ljust(72)) #start scan angle value, decimal degrees
            header_list.append("INCREME:" + str(round(range_,6)).ljust(72)) #signed scan angle increment between frames
            
        else:
            print header_vals["Oscillation_axis"]
            print "oscillation axis currently not supported"
            sys.exit()
        #comment
        header_list.append("NUMBER :" + str(header_vals["frame_no"]).ljust(72)) #sequence number of this frame in series
        header_list.append("NFRAMES:" + str(header_vals["N_oscillations"]).ljust(72)) #total number of frames in the series
        header_list.append("ANGLES :" + str(round(es_2theta,6)).ljust(18) + str(round(es_omega,6)).ljust(18) + str(round(es_phi,6)).ljust(18) + str(round(es_chi,6)).ljust(18)) #2Th, Omg, Phi, Chi
        #calculating number of pixels > 64000
        header_list.append("NOVER64:" +str(header_vals["nover64"]).ljust(24)+ "0".ljust(24)*2) #no of pixels over 64k
        header_list.append("NPIXELB:" + str(header_vals["npixelb"]).ljust(36) + "1".ljust(36)) #bytes/pixel in main image, bytes/pixel in underflow table
        header_list.append("NROWS  :" + str(header_vals["nrows"]).ljust(72)) #rasters in frame
        header_list.append("NCOLS  :" + str(header_vals["ncols"]).ljust(72)) #pixels/raster
        header_list.append("WORDORD:" + "0".ljust(72)) #order of bytes in word (0=LSB first)
        header_list.append("LONGORD:" + "0".ljust(72)) #order of words in a longword
        header_list.append("TARGET :" + str(default_vals["Xray_target_material"]).ljust(72)) #x-ray target material
        header_list.append("SOURCEK:" + str(default_vals["Voltage"]).ljust(72)) #x-ray source kV
        header_list.append("SOURCEM:" + str(default_vals["Current"]).ljust(72)) #source milliamps
        header_list.append("FILTER :" + str(default_vals["Filter_monochromator_setting"]).ljust(72)) #filter/monochromator setting
        header_list.append("CELL   :" + "0.000000".ljust(14)*4 + "0.000000".ljust(16)) #unit cell; A,B,C,ALPHA,BETA,GAMMA
        header_list.append("CELL   :" + "0.000000".ljust(72))
        header_list.append("MATRIX :" + "0.000000".ljust(14)*4 + "0.000000".ljust(16)) #orientation matrix (P3 conventions)
        header_list.append("MATRIX :" + "0.000000".ljust(18)*4)
        header_list.append("LOWTEMP:" + str(default_vals["Low_temp_flag"]).ljust(24) + str(default_vals["exp_temp_in_hundredths"]).ljust(24) + str(default_vals["ccd_temp_in_hundredths"]).ljust(24))#low temp flag, experimental temp in hundredths C, ccd temp in hundredths C
        header_list.append("ZOOM   :" + str(default_vals["X_c"]).ljust(24)+ str(default_vals["Y_c"]).ljust(24)+ str(default_vals["Mag"]).ljust(24)) #zoom Xc, Yc, Mag
        #calculating beam center
        cent = header_vals["Beam_xy"]
        cent_x = cent[0]
        cent_y = cent[1] #if cbf is measured from bottom left corner, values 'should' remain the same after padding
        padding = header_vals["padding_x"] + header_vals["padding_y"]
        cent_x = cent_x + padding   #need to add padding in pixels to x, as measured from bottom left corner of detector (where the padding is put) 
        cent_x = round(cent_x,6)
        cent_y = round(cent_y,6)
        cent_x = str(cent_x).ljust(18)
        cent_y = str(cent_y).ljust(18)
        header_list.append("CENTER :" + cent_x + cent_y + cent_x + cent_y) #x,y, (in pixels) of direct beam at 2-theta=0
        #detector distance
        distanc = float(header_vals["Detector_distance"])*100 #converts to cm
        distanc_string = str(np.round(distanc,6))
        header_list.append("DISTANC:" + distanc_string.ljust(36)*2) #sample-detector distance (cm)       
        header_list.append("TRAILER:" + "0".ljust(72)) #byte pointer to trailer info
        header_list.append("COMPRES:" + "NONE".ljust(72)) #compression scheme if any
        header_list.append("LINEAR :" + "1.000000".ljust(36) + "0.000000".ljust(36)) #linear scale, offset for pixel values
        header_list.append("PHD    :" + "0.000000".ljust(36)*2) #pulse height settings
        header_list.append("PREAMP :" + "0".ljust(36)*2) #preamp gain settings
        header_list.append("CORRECT:" + str(default_vals["Flood_field_correction_filename"]).ljust(72)) #flood field correction filename
        header_list.append("WARPFIL:" + str(default_vals["Brass_plate_correction_filename"]).ljust(72)) #brassplate correction filename
        wavelen = str(header_vals["Wavelength"])
        wavelen = wavelen.ljust(9,"0")
        wav = round(float(wavelen),6)
        wav = str(wav).ljust(18)
        header_list.append("WAVELEN:" + wav*4) #wavelengths
        row = str(header_vals["max_row"]).ljust(36) #should remain the same
        col = str(header_vals["max_col"]).ljust(36)
        header_list.append("MAXXY  :" + col + row) #x,y pixel number of maximum counts
        header_list.append("AXIS   :" + axis_val.ljust(72)) #scan axis (1-4 for 2Th, Omg, Phi, Chi)
        header_list.append("ENDING :" + str(round(ee_2theta,6)).ljust(18) + str(round(ee_omega,6)).ljust(18) + str(round(ee_phi,6)).ljust(18) + str(round(ee_chi,6)).ljust(18)) #actual goniometer angles at end of frame
        #theta_ang_in_deg = self.start[0] #2 theta
        #theta_rad = np.radians(theta_ang_in_deg)
        #dist = -(self.det[2]) #d_det in cm, needs to be in mm for this.
        #x_correct = (np.tan(theta_rad)*dist)/self.pix_size
        #x_correct = round(float(x_correct),6)
        header_list.append("DETPAR :"+ str(default_vals["dX"]).ljust(14)+ str(default_vals["dY"]).ljust(14)+ str(default_vals["dDist"]).ljust(14) + str(default_vals["pitch"]).ljust(14) + str(default_vals["roll"]).ljust(16)) #detector position corrections (dX, dY, dDist, Pitch, Roll, Yaw(below))
        header_list.append("DETPAR :" + str(default_vals["yaw"]).ljust(72) )#Yaw
        header_list.append("LUT    :" + str(default_vals["Display_lookup"]).ljust(72)) #recommended display lookup table
        header_list.append("DISPLIM:" + str(default_vals["Display_limit1"]).ljust(36) + str(default_vals["Display_limit2"]).ljust(36)) #recomended display limits
        header_list.append("PROGRAM:" + "ALPHA VERSION: cbf_to_sfrm by N Johnson M Probert (Newcastle University)".ljust(72)) #nature and version of programme writing frame
        header_list.append("ROTATE :" + "0".ljust(72)) #nonzero if acquired by rotation
        header_list.append("BITMASK:" + "$NULL".ljust(72)) #filename of active pixel mask
        header_list.append("OCTMASK:" + "0".ljust(12)*3 + "1023".ljust(12)*2 + "2046".ljust(12)) #min x, min x+y, min y, max x-y, max x, max x+y,
        header_list.append("OCTMASK:" + "1023".ljust(36)*2) #max y, max y-x
        header_list.append("ESDCELL:" + "0.000000".ljust(14)*4 + "0.000000".ljust(16)) #unit cell parameter standard deviations
        header_list.append("ESDCELL:" + "0.000000".ljust(72))
        pix_size_tup = (header_vals["Pixel_size"])
        #calculating pixels per cm, way .sfrm currently requires pixel size to be input
        pix_size = float(pix_size_tup[0])*1000
        pix_mm = (10/pix_size)*(512/float(some_marker))  #assumes detector size in 512x512 format
        pix_mm = round(pix_mm,6)
        header_list.append("DETTYPE:" + str(default_vals["Det_type"]).ljust(21) + str(pix_mm).ljust(12) + str(default_vals["Phosphor_distance"]).ljust(12) + str(default_vals["Det_4"]).ljust(2) +  str(default_vals["Brass_hole_spacing"]).ljust(12) + str(default_vals["Be_window_thickness"]).ljust(12) + str(default_vals["Det_7"]).ljust(1)) #det-type (needs a name), pix per cm (scaled), phosphor, ?, brass, Be, ?)
        #det-type (needs a name), pix per cm (scaled), phosphor, ?, brass, Be, ?)
        #currently APEX2 needs a recognized detector (which could cause problems if values that are built in for the detector affect the final refinement).
        header_list.append("NEXP   :" + str(default_vals["no_of_exposures"]).ljust(14) + str(default_vals["per_eposure_bias"]).ljust(14) + str(default_vals["baseline_offset"]).ljust(14) + str(default_vals["orientation"]).ljust(14) + str(default_vals["over_scan_flag"]).ljust(16)) #no. of exposures, per-exposure bias level in ADU, baseline offset, orientation, overscan flag (1 if overscan used)
        header_list.append("CCDPARM:" + str(default_vals["readnoise"]).ljust(14)+ str(default_vals["e_ADU"]).ljust(14)+ str(default_vals["e_photon"]).ljust(14)+ str(default_vals["bias"]).ljust(14) + str(default_vals["full_scale"]).ljust(16)) #ccd parameters; readnoise, e/ADU, e/photon, bias, full scale
        #these five can all be changed, not required values, but can be entered for reference/completeness
        #can be edited in default value files
        header_list.append("CHEM   :" + str(default_vals["Chemical_formula"] ).ljust(72)) #chemical formula
        header_list.append("MORPH  :" + str(default_vals["Crystal_morphology"]).ljust(72)) #crystal morphology
        header_list.append("CCOLOR :" + str(default_vals["Crystal_colour"]).ljust(72)) #crystal colour
        header_list.append("CSIZE  :" + str(default_vals["Crystal_dimension_1"]).ljust(24)+ str(default_vals["Crystal_dimension_2"]).ljust(24)+ str(default_vals["Crystal_dimension_3"]).ljust(24) ) #crystal dimensions (3 ea)
        header_list.append("DNSMET :" + str(default_vals["Density_measurement_method"]).ljust(72)) #density measurement method
        #
        header_list.append("DARK   :" + "NONE".ljust(72)) #name of dark correction
        header_list.append("AUTORNG:" + str(default_vals["gain"]).ljust(14) + str(default_vals["high_speed_time"]).ljust(14) + str(default_vals["scale"]).ljust(14) + str(default_vals["offset"]).ljust(14) + str(default_vals["auto_full_scale"]).ljust(16)) # auto-ranging: gain, high-speed time, scale, offset, full set
        header_list.append("ZEROADJ:"+ "0.000000".ljust(18)*4) #goniometer zero corections (refined in least squares)
        header_list.append("XTRANS :"+ "0.000000".ljust(24)*3) #sample transitions (refined in least squares)
        header_list.append("HKL&XY :"+ "0.000000".ljust(14)*4 + "0.000000".ljust(16)) #hkl and pixel xy for reciprocal space scan
        header_list.append("AXES2  :"+ "0.000000".ljust(18)*4) # diffractometer setting linear axes
        header_list.append("ENDING2:"+ "0.000000".ljust(18)*3 + "0.000000".ljust(18)) #actual gonometer axes at end of frame (used by GADDS)
        header_list.append("FILTER2:"+ str(default_vals["mono_2theta"]).ljust(18)+ str(default_vals["mono_roll"]).ljust(18)+ str(default_vals["mono_tilt"]).ljust(18)+ str(default_vals["attenuation_factor"]).ljust(18)) #monochromator 2theta angle, monochromator roll angle, beam tilt angle, attenuation factor
        #
        header_list.append("LEPTOS :" + "".ljust(72))
        header_list.append("CFR: HDR:  IMG: " + "".ljust(62,".") + "\x1a" + "\x04" ) #final lines of header
        header_string = "".join(header_list) # joins everything in header list into a string
        #pixelb_index = [i for i, s in enumerate(header_list) if 'NPIXELB:' in s] #changes bytes/pixel value in dynamic mask header value
        #pix_string = "NPIXELB:" + "1".ljust(36) + "1".ljust(36)
        #header_list[pixelb_index[0]] = pix_string #resets pixel count to 1 bperp for dynamic mask
        #self.dm_header = "".join(header_list)
        return header_string 

if len(sys.argv) != 8:
    print "useage of this file: %s, cbf_file, name_string of cbf file, default file, output directory, number of runs, number of processes to be run, number of frames to add" %sys.argv[0]
    print "cbf_file = full file location of one of the files to be converted, name_string = file name up to (but not including) the run number"
    print "e.g. python %s /home/Documents/cbfs/cbf_100K_02_00004.cbf cbf_100K_ default_file.txt /home/Documents/converted_cbf/cbf_100K 4 7 2" %sys.argv[0]
    print "for runs 1-4 of cbf_100K to be converted with 7 processes being set away and two cbf frames being added together into one sfrm frame."
    sys.exit()

cbf_file = sys.argv[1]
name_string = sys.argv[2]
default_file = sys.argv[3]
dir_out = sys.argv[4]
runs = int(sys.argv[5])
cores = int(sys.argv[6])
compression = int(sys.argv[7])

Conversion(cbf_file, name_string, default_file, dir_out, runs, cores, compression)

print datetime.datetime.now().time()
