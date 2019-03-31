import urllib2

response = urllib2.urlopen('https://www.celestrak.com/NORAD/elements/orbcomm.txt')
html = response.read()

if '21576' in html: # check if NORAD ID of ORBCOMM-X is in text file
	with open('./tles/orbcomm.txt', 'w') as f:
		f.write(html)
	print("Updated Orbcomm TLE file.")
else:
	print("Error updating TLE: {}".format(html))