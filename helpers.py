
from glob import glob

def get_tle_lines(sat_name, tle_dir='./tles'):
    line0, line1, line2 = '', '', ''
    tles = glob(tle_dir + '/*.txt')
    for tle in tles:
        with open(tle, 'r') as f:
            text = 'abc'
            while text != '':
                text = f.readline()
                if sat_name.lower() in text.lower():
                    line0 = text
                    line1 = f.readline()
                    line2 = f.readline()
                    break
        if line0 != '':
            break
    return line0, line1, line2