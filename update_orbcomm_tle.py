try:
	from urllib2 import urlopen 
except:
	from urllib.request import urlopen

response = urlopen('https://www.celestrak.com/NORAD/elements/orbcomm.txt')
html = str(response.read())

if '21576' in html: # check if NORAD ID of ORBCOMM-X is in text file
	with open('./tles/orbcomm.txt', 'w') as f:
		f.write(html)
	print("Updated Orbcomm TLE file.")
else:
	print("Error updating TLE: {}".format(html))