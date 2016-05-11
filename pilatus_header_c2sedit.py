#==========================================================================
# == pilatus_header.py
# == Parsing detector and experiment parameters from Pilatus file headers
# == 19.09.2011
# == Marcus Mueller (marcus.mueller@dectris.com)
# == Dectris Ltd.
#==========================================================================

#minor edit made for cbf_to_sfrm_2016.py compatibility, Newcastle Crystallography

import re
import sys
SPACE_EQUIVALENT_CHARACTERS = '#:=,()'
NON_OPTIONAL_KEYWORDS = {
#key: (pattern, [value_indeces], type)
'Detector_identifier': ('Detector ', [slice(1, None)], str),
'Pixel_size': ('Pixel_size', [1, 4], float),
'Silicon': ('Silicon', [3], float),
'Exposure_time': ('Exposure_time', [1], float),
'Exposure_period': ('Exposure_period', [1], float),
'Tau': ('Tau', [1], float),
'Count_cutoff': ('Count_cutoff', [1], int),
'Threshold_setting': ('Threshold_setting', [1], float),
'Gain_setting': ('Gain_setting', [1, 2], str),
'N_excluded_pixels': ('N_excluded_pixels', [1], int),
'Excluded_pixels': ('Excluded_pixels', [1], str),
'Flat_field': ('Flat_field', [1], str),
'Trim_file': ('Trim_file', [1], str),
'Image_path': ('Image_path', [1], str),
}
OPTIONAL_KEYWORDS = {
'Wavelength': ('Wavelength', [1], float),
'Energy_range': ('Energy_range', [1, 2], float),
'Detector_distance': ('Detector_distance', [1], float),
'Detector_Voffset': ('Detector_Voffset', [1], float),
'Beam_xy': ('Beam_xy', [1, 2], float),
'Beam_x': ('Beam_xy', [1], float),
'Beam_y': ('Beam_xy', [2], float),
'Flux': ('Flux', [1], str),
'Filter_transmission': ('Filter_transmission', [1], float),
'Start_angle': ('Start_angle', [1], float),
'Angle_increment': ('Angle_increment', [1], float),
'Detector_2theta': ('Detector_2theta', [1], float),
'Polarization': ('Polarization', [1], float),
'Alpha': ('Alpha', [1], float),
'Kappa': ('Kappa', [1], float),
'Phi': ('Phi ', [1], float), #minor edit to this line
'Phi_increment': ('Phi_increment', [1], float),
'Chi': ('Chi ', [1], float), #minor edit to this line
'Chi_increment': ('Chi_increment', [1], float),
'Oscillation_axis': ('Oscillation_axis', [slice(1, None)], str),
'N_oscillations': ('N_oscillations', [1], int),
'Start_position': ('Start_position', [1], float),
'Position_increment': ('Position_increment', [1], float),
'Shutter_time': ('Shutter_time', [1], float),
'Omega': ('Omega ', [1], float), #minor edit to include omega
'Omega_increment' : ('Omega_increment', [1], float) #minor edit to include omega incremenent
}
ALL_KEYWORDS = {}
ALL_KEYWORDS.update(NON_OPTIONAL_KEYWORDS)
ALL_KEYWORDS.update(OPTIONAL_KEYWORDS)
class PilatusHeader(object):
    """
    Class for parsing contents of a Pilatus cbf header from a minimal
    or full cbf file.
    Parsing the Pilatus cbf header populates the header_dict dictionary,
    using Pilatus cbf header keywords as keys.
    """
    def __init__(self, filepath, non_binary_length=4096):
        self.filepath = filepath
        self.non_binary_length = non_binary_length
        self.header_dict = {}
        self.header_lines = []
        self.read_header_lines()
        self.parse_header()

    def _rawheader(self):
        return open(self.filepath, 'rb').read(self.non_binary_length)

    def has_pilatus_cbf_convention(self):
        # Check for the _array_data.header_convention data item
        pilatus_header_pattern = re.compile(
        r'''_array_data.header_convention +["']?(SLS|PILATUS)'''
        r'''_\d+(\.?\d*)*["']?''')
        return bool(pilatus_header_pattern.search(self._rawheader()))
    def read_header_lines(self):
        """
        Populate the self.header_lines list.
        """
        contents_pattern = re.compile(r'''_array_data.header_contents\s+;.*?;''',re.DOTALL)
        contents_match = contents_pattern.search(self._rawheader())
        assert contents_match is not None
        self.header_lines = contents_match.group().splitlines()
        
    def _spaced_header_lines(self):
        """
        Return header_lines with all space equivalent charecters converted
        to space.
        """
        spaced_header_lines = []
        for line in self.header_lines:
            for space_equivalent in SPACE_EQUIVALENT_CHARACTERS:
                line = line.replace(space_equivalent, ' ')
            spaced_header_lines.append(line)
        return spaced_header_lines

    def parse_header(self):
        """
        Populate self.header_dict with contents of Pilatus cbf header
        """
        assert self.has_pilatus_cbf_convention()
        if len(self.header_lines) == 0:
            self.read_header_lines()
        # parse the header lines
        for key, (pattern, valueindices, datatype) in ALL_KEYWORDS.items():
            for line in self._spaced_header_lines():
                if pattern in line:
                    values = []
                    for idx in valueindices:
                        try:
                            # handle multiple or single values
                            if isinstance(idx, slice):
                                values += line.split()[idx]
                            else:
                                values.append(line.split()[idx])
                        except IndexError:
                            print ('No value at index %d on header line:'%s % (idx, line))
                    value = self._datatype_handling(values, key, datatype)
                    if value is not None:
                        self.header_dict[key] = value
                        
    def _datatype_handling(self, values, key, datatype):
        # handle rare cases of value "not set"
        if datatype is float and values[0] == 'not':
            # NON_OPTIONAL_KEYWORDS should always have value, at least NaN
            if key in NON_OPTIONAL_KEYWORDS:
                return float('NaN')
            else:
                return None
        # do the conversion for standard cases
        if len(values) == 1:
            values = datatype(values[0])
        else:
            if datatype is str:
                values = ' '.join(values)
            else:
                values = tuple([datatype(v) for v in values])
        return values
    
    def get_beam_xy_mm(self, factor=1000):
        return tuple([n * size * factor for n, size in zip(
            self.header_dict['Beam_xy'], self.header_dict['Pixel_size'])])

    def get_date_time(self):
        """
        Return date and time of image acquistion.
        Works for format of current camserver versions
        2011-06-04T04:57:02.976
        or format of old camserver versions
        2011/Sep/12 09:21:27.252
        """
        date_time_pattern = re.compile(
            r'(\d{4}-\d{2}-\d{2}T|\d{4}/\D+/\d{2} )\d{2}:\d{2}:\d{2}.\d+')
        return date_time_pattern.search(self._rawheader()).group()

    def get_time(self):
        time_pattern = re.compile(r'\d{2}:\d{2}:\d{2}.\d+')
        return time_pattern.search(self.get_date_time()).group()

    date_time = property(get_date_time)
    time = property(get_time)

if __name__ == '__main__':
    if len(sys.argv) == 3:
        header = PilatusHeader(sys.argv[1], int(sys.argv[2]))
    elif len(sys.argv) == 2:
        header = PilatusHeader(sys.argv[1])
    else:
        print "Usage: %s cbf_file_path [bytes_to_read]" % sys.argv[0]
        sys.exit()
    #header.parse_header()
    header_vals = header.header_dict
    print len(header_vals)
    print type(header_vals)
    print header.get_date_time()
    for k, v in header.header_dict.items():
        print k, v

    print header.header_dict["Flux"]
