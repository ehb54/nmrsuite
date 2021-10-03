# Place after y is defined (try-except statement), before y is redefined


    """ ===Best L0=1 Solution===
    Relative Error: 0.5592553791798197
    Chi2: 5453.266985593737
    L: 70, Chi2/L: 77.90381407991053 """


    

    """ L0-norm_list = []
    relative_error_list = []
    columns_list = []
    weights_list = [] """
    
"""
[[ 1.962055    4.7564434   3.8757345  ... -2.5194603  -2.5194603
  -2.5194603 ]
 [ 2.8499926   2.1388422  -1.3436634  ... -4.2734486  -5.6024531
  -5.6024529 ]
 [-5.4186163  -4.78244    -3.212591   ...  1.9601804   0.42877677
  -1.2780167 ]
 ...
 [ 3.3110707   3.7154739   1.805861   ...  1.2981889   1.0233878
   2.5849734 ]
 [-1.2021795   3.6883057   3.8412075  ...  2.8213401   3.2715511
   2.2740167 ]
 [ 0.46721389 -5.597108   -6.6658696  ... -7.5647794  -7.9535617
  -7.4620285 ]] 
 
 [6 ] 
 
 [0.4031567605433127 ] 
    """

def getFunction (lines, matrxfile, datafile):
    newOutput = ""
    matrix = np.loadtxt(matrixfile)
    data = np.loadtxt(datafile)

    experimental_values = data[:, 0]
    experimental_errors = data[:, 1]
    L0_solution_lines_start = lines.index("Best solutions by l0-norm (column index starts at 1):")
    L0_solution_lines_end = [lines.index(l) for l in lines if l.startswith("L-curve: ")][0]
    L0_solution_lines = list(range(L0_solution_lines_start + 1, L0_solution_lines_end - 1))
    for i, L0_solution_line in enumerate(L0_solution_lines):
        L0_solution_line = lines[L0_solution_line]
        try:
            m = re.match(r'l0-norm=(.*): Relative Error=(.*), Columns=(.*), Weights=(.*)', L0_solution_line)
            l0_norm = m.group(1)
            rel_error_value = m.group(2)
            columns = np.array(stringToList(m.group(3))).astype(int)-1
            weights = np.array(stringToList(m.group(4))).astype(float)
        except:
            printQuit (f"""L0_solution_line: {L0_solution_line}
            \n\n""")

        output = f"{str(matrix)} \n {matrix.dtype} \n \n {str(columns)} \n {columns.dtype} \n \n {str(weights)} \n {weights.dtype} \n \n"
        #printQuit(output)
        #printQuit(f"{matrix} \n\n {columns} {matrix[:, columns]}")
        #predicted_values = np.multiply(matrix[:, columns], weights)

        try:
            if (columns.shape[0] == 1):
                predicted_values = np.multiply(matrix[:, columns], weights)
            else:
                predicted_values = np.matmul(matrix[:, columns], weights)
        except Exception as e:
            printQuit(f"""\n\nMatrix: {matrix}
            \n\nColumns: {columns}
            \n\nWeights: {weights}
            \n\nMatrix[:, columns] shape: {matrix[:, columns].shape}
            \n\nWeights shape: {weights.shape}
            \n\nException: {e}""")


        try:
            Chi2 = np.sum(np.divide(np.subtract(experimental_values, predicted_values), experimental_errors))
        except:
            printQuit(f"""Experimental values:{experimental_values}
            \n\nExperimental errors: {experimental_errors}
            \n\nMatrix: {matrix}
            \n\nColumns: {columns}
            \n\nWeights: {weights}
            \n\nPredicted values: {predicted_values}
            \n\nMatrix[:, columns] shape: {matrix[:, columns].shape}
            \n\nWeights shape: {weights.shape}""")
        L = np.shape(data)[0] # Length


        indices = np.array(list(range(1, L+1)))
        rel_error_values = np.divide(np.subtract(predicted_values, experimental_values), experimental_errors)

        newOutput += f"===Best L0={i} Solution===\nRelative Error: {rel_error_value}\nChi2: {Chi2}\nL: {L}, Chi2/L: {Chi2/L}"

        #output = f"{tableString} \n\n Predicted Values: {str(predicted_values)} \n\n Experimental Values: {str(experimental_values)}"
        
        columns = ["<index>", "<exp. value>", "<pred. value>", "<relative err.>"]
        data = data[:, 0]
        predicted_values = predicted_values.flatten()
        rel_error_values = rel_error_values.flatten()

        printQuit(f"""Data: {data}
        \n\nIndices: {indices}
        \n\nExperimental Values: {experimental_values}
        \n\nPredicted Values: {predicted_values}
        Relative Error Values: {rel_error_values}""")
        
        data = np.column_stack((indices, experimental_values, predicted_values.flatten(), rel_error_values[0]))
        
 

        row_format ="{:>15}" * (len(columns) + 1)
        newOutput += "\n" + row_format.format("", *columns)
        for row in data:
            newOutput += "\n" + ("{:>15}  " * (len(row) + 1)).format("", *row)
        #for i in range (L):
            #newOutput += "\n" + "%.5f %.5f %.5f %.5f" % (indices[i], experimental_values[i], predicted_values[i], rel_error_value[i])

        '''
        Need:
        Rel Error
        Chi2
        L
        Chi2
        '''

        """ L0-norm_list.append(m.group(0))
        relative_error_list.append(m.group(1))
        columns_list.append(m.group(2))
        weights_list.append(m.group(3)) """

    

    printQuit(newOutput)
    #printQuit(f"{newOutput} \n {outputStr}")