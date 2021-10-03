#!/usr/bin/python3
import json, sys, os, io, subprocess

if __name__=='__main__':

        argv_io_string = io.StringIO(sys.argv[1])
#	json_variables = json.load(sys.stdin)	
        json_variables = json.load(argv_io_string)

        axes = json_variables['axes']
        axesl = float(json_variables['axesl'])
        bicelle = float(json_variables['bicelle'])
        dockpdb = json_variables['dockpdb']     
        model = int(json_variables['model'])
        pdb = json_variables['pdb'][0]
        pdb_base = os.path.basename(pdb)
        directory = os.path.dirname(pdb)
        os.chdir(directory)
        if 'rdc' in json_variables:
            rdc_base = os.path.basename(json_variables['rdc'][0])

        exec_string = "java -jar /opt/genapp/pati/bin/pati-1.1.jar"
        if axes != "":
            exec_string += " -axes " + str(axes)
        exec_string += " -axesl " + str(axesl)
        exec_string += " -bicelle " + str(bicelle)
        if 'dock' in json_variables:
            exec_string += " -dock"
        if len(dockpdb) > 0:
            exec_string += " -dockpdb " + str(dockpdb)
        if 'draw' in json_variables:
            exec_string += " -draw"
        exec_string += " -model " + str(model)
        if 'nostat' in json_variables:
            exec_string += " -nostat"
        exec_string += " -pdb " + str(pdb_base)
        if 'rdc' in json_variables:
            exec_string += " -rdc " + str(rdc_base)
        if 'robust' in json_variables:
            exec_string += " -robust"

#	path = base_directory.replace('\/','/') + "/"
#	os.chdir(path)
#        sys.path.append('./')
	
        output_string = str(subprocess.check_output(exec_string, shell=True))
        output_string = output_string.replace("\\n", "\n")

        results_file = open("output.txt", "wt")
        print(output_string, file=results_file)
        results_file.close()

        output = {} 
        output['mainoutput'] = str(directory + "/output.txt")
        if len(axes) > 0:
            if 'rdc' in json_variables:
                output['axesoutput'] = str(directory + "/" + axes + ".py")
            output['axesoutputpredicted'] = str(directory + "/" + axes + "_predicted.py")
        
        if 'dock' in json_variables and 'rdc' in json_variabls and len(dockpdb) > 0:
            output['dockoutput'] = str(directory + "/" + dockpdb)

        if 'draw' in json_variables and 'rdc' in json_variables:
            # Make first graph
            graph_file = open("graph.txt", "rt")
            points = json.loads(graph_file.read())
            graph_title = points['title']
            del points['title']
            points['mode'] = 'markers'
            points['name'] = 'Points'
            xmin = min(min(points['x']), min(points['y']))
            xmax = max(max(points['x']), max(points['y']))
            line_data = {'x' : [xmin,xmax], 'y' : [xmin,xmax], 'name' : 'y=x', 'mode' : 'lines'}            

            graph_data = {}
            graph_data['data'] = [points, line_data]
            graph_data['layout'] = {"title" : graph_title, 'xaxis': {'title':'RDC Experimental'}, 'yaxis':{'title':'RDC Computed'}}
            output['graphoutput'] = graph_data

            # Make chain graphs
            graph_file = open("graph2.txt", "rt")
            graphs = json.loads(graph_file.read())
            chains = graphs['graphs']
            num_chains = len(graphs)

            chain_json = {}
            chain_data = []
            for i, chain in enumerate(chains):
                chain_title = chain['title']
                del chain['title']
                chain['type'] = 'bar'
                chain['name'] = chain_title
                chain['xaxis'] = 'x' + str(i)
                chain['yaxis'] = 'y' + str(i)
                chain_data.append(chain)
            chain_layout = {'grid': {'rows' : num_chains, 'columns' : 1, 'pattern' : 'independent', 'roworder' : 'top to bottom'}, 'title' : 'Residuals by Chain', 'xaxis': {'title':'Residue/Nucleotide #'}, 'yaxis':{'title':'Residuals'}}
            chain_json['data'] = chain_data
            chain_json['layout'] = chain_layout
            output['chainoutput'] = chain_json

        print(json.dumps(output))

