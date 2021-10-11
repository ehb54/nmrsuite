#!/opt/miniconda2/bin/python
import json, sys, StringIO

if __name__=='__main__':

	argv_io_string = StringIO.StringIO(sys.argv[1])
	json_variables = json.load(argv_io_string)


#	path = base_directory.replace('\/','/') + "/"
#	os.chdir(path)
#	sys.path.append('./')

	output = {} 
	output["note"] = "test.py executable";

        output['_textarea'] = "JSON output from executable:\n" + json.dumps( output, indent=4 ) + "\n\n";
        output['_textarea'] += "JSON input to executable:\n" + json.dumps( json_variables, indent=4 ) + "\n";

	print json.dumps(output)
		
