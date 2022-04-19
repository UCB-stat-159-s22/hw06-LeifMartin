import readligo

def test_read_hdf5(filename, readstrain=True):
	read_hdf5(filename, readstrain=True)
	return None

def test_getstrain(m_strain, meta, m_dq): #Checking that the output is of the desired form.
	assert 
