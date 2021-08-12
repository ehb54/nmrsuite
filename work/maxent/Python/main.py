#!/usr/bin/python3.6

import json, sys, time
from io import StringIO ## for Python3
from genapp3 import genapp ## python3

if __name__=='__main__':

        argv_io_string = StringIO(sys.argv[1])
        json_variables = json.load(argv_io_string)

        ### initialize the genapp object
        ga = genapp( json_variables )

        ga.udpprogress( 0 );
        ga.udpmessage( { "output3" : 0.0 } );
        ga.udpmessage( { "_textarea" : "udp message to _textarea\n" } );
        ga.tcpmessage( { "_textarea" : "tcp message to _textarea\n" } );
        ga.udpmessage( { "output2" : "udp message to output2" } );
        ga.tcpmessage( { "output2" : "tcp message to output2" } );

        output = {}

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
        ga.udpprogress( 0.3 );
        ga.udpmessage( { "output3" : 0.3 } );

        time.sleep(5)

        # add a plot point to 'plotbar'
 
        output['plotbar']['data'][0]['x'].append( "llamas" )
        output['plotbar']['data'][0]['y'].append( "17" )

        # update the UI via tcpmessage

        ga.tcpmessage( { "plotbar" : output['plotbar'] } )
        ga.udpprogress( 0.5 );
        ga.udpmessage( { "output3" : 0.5 } );

        time.sleep(5)

        # add another plot point to 'plotbar'

        output['plotbar']['data'][0]['x'].append( "alpacas" )
        output['plotbar']['data'][0]['y'].append( "22" )

        # change values for the 2nd data trace (index 1) of 'plotline'

        output['plotline']['data'][1]['y'] = [ 9, 7, 3, 6 ];

        # update both plots in a single tcpmessage
        ga.tcpmessage( { "plotbar" : output['plotbar'], "plotline" : output['plotline'] } )
        ga.udpprogress( 0.7 );
        ga.udpmessage( { "output3" : 0.7 } );

        time.sleep(5)

        # change values for the 2nd data trace of plotline again for the final output

        output['plotline']['data'][1]['y'] = [ 17, 6, 12, 5 ];

        # optional debugging output to the textarea
#        output['_textarea'] = "JSON output from executable:\n" + json.dumps( output, indent=4 ) + "\n\n";
#        output['_textarea'] += "JSON input to executable:\n" + json.dumps( json_variables, indent=4 ) + "\n";

        # output final object to the UI
        output[ '_progress' ] = 1.0
        output[ 'output3' ] = 1.0
        print(json.dumps(output))
