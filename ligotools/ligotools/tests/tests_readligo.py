import ligotools

def test_read_hdf5(filename, readstrain=True):
	ligotools.read_hdf5(filename, readstrain=True)
	return None

def test_getstrain(): #Checking that the output is of the desired form.
	assert len(ligotools.getstrain(start=1, stop=2, ifo=H1)[1])>=3,'getstrain output is incorrect. Should include 3 keys from output 2.'
