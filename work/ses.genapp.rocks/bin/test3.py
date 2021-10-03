#!/usr/bin/python3.6

import os
import re

server_path = "/var/www/html//results/users/mcasertano/no_project_specified/PDB"
os.chdir( "/opt/genapp/parnmr/results/users/mcasertano/Server_PDB_Path_11_07_20" )
new_server_path = re.sub( r'/[^/]+$', '', os.getcwd() ) + re.sub( r'^.*results/users/[^\/]+', '', server_path )

print( "new_server_path:" + new_server_path )



