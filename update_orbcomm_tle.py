try:
	from urllib2 import urlopen 
except:
	from urllib.request import urlopen

response = urlopen('https://celestrak.org/NORAD/elements/gp.php?GROUP=orbcomm&FORMAT=tle')
html = response.read().decode('ascii')

if '21576' in html: # check if NORAD ID of ORBCOMM-X is in text file
	with open('./tles/orbcomm.txt', 'w') as f:
		f.write(html)
	print("Updated Orbcomm TLE file.")
else:
	print("Error updating TLE: {}".format(html))
