'''
def saveL0Solutions(output):
    key_1 = "===Best"
    key_2 = "["

    if os.path.exists(L0_folder_name):
        shutil.rmtree(L0_folder_name)
    
    time.sleep(.1) # Sometimes a delay is needed
     
    os.mkdir(L0_folder_name)
    
    os.chdir(L0_folder_name)

    line_numbers_1 = []
    line_numbers_2 = []
    L0_numbers = []
    Chi2L_numbers = []
    rel_error_numbers = []

    #output = output.append(5) # Remove after testing

    for index, line in enumerate(output):
        if (line.lstrip().startswith(key_1)):
            line_numbers_1.append(index)
        if (line.startswith(key_2)):
            line_numbers_2.append(index)

    num_L0s = len(line_numbers_1)

    if (num_L0s != len(line_numbers_2)):
        raise Exception ("Output error. Please make sure the entire output is printing.\nLine numbers 1: {}.\nLine numbers 2: {}.\n Output: {}".format(line_numbers_1, line_numbers_2, output))

    columns = ["Index", "Exp. Data Info", "Exp. Value" "Pred. Value", "Relative Err."]

    # Captures all the output -- check if this 
    for i in range (len(line_numbers_1)):
        table = []
        start_line = line_numbers_1[i]
        end_line = line_numbers_2[i]
        relevant_lines = output[start_line+5:end_line]
        L0_number = int(re.sub("[^0-9]", "", output[start_line]))
        L0_numbers.append(L0_number)

        Chi2L_line = output[start_line + 3]
        colon_locations = [x.start() for x in re.finditer(':', Chi2L_line)]
        try:
            Chi2L = float(Chi2L_line[colon_locations[-1] + 1:])
        except Exception as e:
            output = f"{colon_locations}   {Chi2L_line}    {e}"
            printQuit(output)
        Chi2L_numbers.append(Chi2L)


        rel_error_line = output[start_line + 1]
        colon_locations = [x.start() for x in re.finditer(':', rel_error_line)]
        try:
            rel_error = float(rel_error_line[colon_locations[0] + 1:])
        except Exception as e:
            output = f"{colon_locations}   {rel_error_line}    {e}"
            printQuit(output)
        rel_error_numbers.append(rel_error)


        for line in relevant_lines:
            line = list(map(lambda x: x.strip(), line.split("\t")))
            table.append(line)

        table = np.array(table)

        df = pd.DataFrame(data=table[:, 1:], index=table[:, 0], columns=columns)
        if (L0_number == 0):
            df.to_csv("BestPossibleSolution.txt", sep='\t')
        else:
            df.to_csv("L0={}.txt".format(L0_number), sep='\t')
        #np.savetxt("L0={}.txt".format(L0_number), np.array(table), fmt='%s')
    
    os.chdir("..")

    return (L0_numbers, Chi2L_numbers, rel_error_numbers)
'''