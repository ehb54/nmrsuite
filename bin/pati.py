#!/opt/miniconda2/bin/python
import json, sys, StringIO

if __name__=='__main__':

	argv_io_string = StringIO.StringIO(sys.argv[1])
	json_variables = json.load(argv_io_string)


#	path = base_directory.replace('\/','/') + "/"
#	os.chdir(path)
#	sys.path.append('./')
	
	output = {} 
	output["note"] = "pati.py executable";

	print json.dumps(output)
		
