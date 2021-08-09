#!/usr/bin/python
import json, sys, time

from io import StringIO

from genapp3 import genapp

from RunMaxEntropy import run_max_entropy


if __name__=='__main__':
    argv_io_string = StringIO(sys.argv[1])
    json_variables = json.load(argv_io_string)

    ### initialize the genapp object
    ga = genapp(json_variables)
    '''
    if none checked:
            stop
        else:
                if max ent checked:
                       run max ent and put results in run folder

                check in run folder to get four files

                if (filter by weight checked):
                        do filter by weight

                if (cluster by rmsd checked):
                        do cluster by rmsd 
    if (max_ent not checked):
            look in run folder to get four files

        if max ent checked:
                run max ent and put results in run folder
        then:
                do the same as above

        so maybe


        


    '''
    output = {}
    output['plotbar'] = {
            "data": [
                    {
                            "x": [
                                    "giraffes",
                                    "orangutans",
                                    "monkeys"
                            ],
                            "y": [
                                    20,
                                    14,
                                    23
                            ],
                            "type": "bar"
                    }
            ]
    }

    output['plotline'] = {
            "data" : [
                    {
                            "x": [1, 2, 3, 4],
                            "y": [10, 15, 13, 17],
                            "mode": "markers",
                            "marker": {
                                    "color": "rgb(219, 64, 82)",
                                    "size": 12
                            }
                    },
                    {
                            "x" : [2, 3, 4, 5],
                            "y" : [16, 5, 11, 9],
                            "mode" : "lines",
                            "line" : {
                                    "color" : "rgb(55, 128, 191)",
                                    "width": 3
                            }
                    },
                    {
                            "x" : [1, 2, 3, 4],
                            "y" : [12, 9, 15, 12],
                            "mode" : "lines+markers",
                            "marker" : {
                                    "color" : "rgb(128, 0, 128)",
                                    "size": 8
                            },
                            "line" : {
                                    "color" : "rgb(128, 0, 128)",
                                    "width" : 1
                            }
                    }
            ],
            "layout" : {
                    "title" : "Line and Scatter Styling"
            }
    }

    # send both plots to the UI via tcpmessage
    ga.tcpmessage( output )

    time.sleep(5)

    # add a plot point to 'plotbar'

    output['plotbar']['data'][0]['x'].append( "llamas" )
    output['plotbar']['data'][0]['y'].append( "17" )

    # update the UI via tcpmessage

    ga.tcpmessage( { "plotbar" : output['plotbar'] } )

    time.sleep(5)

    # add another plot point to 'plotbar'

    output['plotbar']['data'][0]['x'].append( "alpacas" )
    output['plotbar']['data'][0]['y'].append( "22" )

    # change values for the 2nd data trace (index 1) of 'plotline'

    output['plotline']['data'][1]['y'] = [ 9, 7, 3, 6 ];

    # update both plots in a single tcpmessage
    ga.tcpmessage( { "plotbar" : output['plotbar'], "plotline" : output['plotline'] } )

    time.sleep(5)

    # change values for the 2nd data trace of plotline again for the final output

    output['plotline']['data'][1]['y'] = [ 17, 6, 12, 5 ];

    # optional debugging output to the textarea
    output['_textarea'] = "JSON output from executable:\n" + json.dumps( output, indent=4 ) + "\n\n";
    output['_textarea'] += "JSON input to executable:\n" + json.dumps( json_variables, indent=4 ) + "\n";

    # output final object to the UI
    print(json.dumps(output))