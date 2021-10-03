# Old code for the PDB path thing

"""
            OLD CODE FOR THIS METHOD:
            if (os.path.exists("PDB")):
                shutil.rmtree("PDB")

            # The below line gets the full pth of the PDB source directory
            source = os.path.join(os.path.dirname(os.path.realpath(server_path)), os.path.basename(server_path))
            destination = os.path.join(os.getcwd(), "PDB")
            try:
                shutil.move(source, destination)
            except Exception as e:
                printQuit(f"There was an error in obtaining the PDB files on the server. The system failed the following
                    move operation:\n\nSource: {source}\nDestination: {destination}\nError: {e}")
            
             files = os.listdir(server_path)
            files = absoluteFilePaths(server_path)
            #"results/users/mcasertano/no_project_specified/struct00002.pdb", "results/users/mcasertano/no_project_specified/PDB"
            # current: /opt/genapp/parnmr/output/html5/results/users/mcasertano/no_project_specified
            #fullPath = combineOverlappingPaths(os.getcwd(), path)

            # Source: /opt/genapp/parnmr/output/html5/results/users/mcasertano/Local_PDB_Files_11_07_20/PDB
            # Destination: /opt/genapp/parnmr/output/html5/results/users/mcasertano/Local_PDB_Files_11_07_20/PDB
            # results/users/[username]/[__]

            
            path1 = os.path.join(project_prefix, server_path)
            path2 = "../Local_PDB_Directory_11_07_20/PDB"
            printQuit(f"Path 1: {path1} ({os.path.exists(path1)})\n
            Path 2: {path2} ({os.path.exists(path2)})\n
            {os.getcwd()}")

            if (os.path.exists("PDB")):
                shutil.rmtree("PDB")
            os.mkdir("PDB")
            try:
                for filename in files:
                    shutil.move(os.path.join(project_prefix, filename), os.path.join("PDB", filename))
            except Exception as e:
                printQuit(f"Error. \n\n Files: {str(files)}. Error: {str(e)}")

"""