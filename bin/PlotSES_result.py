#!/usr/bin/python3.6

#The way I thought about this: the user will run SES and see the L curve
#output image. Then the user will pick which ensemble size he/she wants
#based on the curve. A new flag of results will be post-analysis, where the user
#will insert the ensemble size he/she picked and the final image result will
#appear.

import numpy as np

from scipy.cluster.hierarchy import fcluster

sortrows = lambda matrix: matrix[np.lexsort(matrix.T[::-1])] # Equivalent of Matlab's sortrows
is_checked = lambda json_variables, checkbox_id : checkbox_id in json_variables.keys()

def split(data_filename, size=2147483647):

    # We are dealing with a file which does not have 
    lines = open(data_filename, "r").readlines()[2:] # Removes first two rows which are not part of data
    num_lines = len(lines)
    maximum_line_length = len(lines[-1].split())-2 # Removes first two columns which are not part of data
    file_read = np.zeros((num_lines, maximum_line_length))

    for i, line in enumerate(lines):
        data = np.array(line.split()[2:], dtype="float")
        
        if (len(data) > 2*size+1):
            break
        
        file_read[i,0:np.size(data)] = data

    file_read_cut = file_read[0:i+1,:] # Only takes the first i lines (which have data)
    
    return file_read_cut


def weight(data, good_solrows, j):

    size_good_solrows = np.size(good_solrows)

    conf = np.zeros((j, size_good_solrows))
    conf_weight = np.zeros((j, size_good_solrows))

    
    for i in range (j):
        conf[i, 0:size_good_solrows] = data[good_solrows[0:size_good_solrows],1+i]
        conf_weight[i,0:size_good_solrows] = data[good_solrows[0:size_good_solrows],j+2+i-1]
    
    #normalizing
    
    for i in range (size_good_solrows):
        conf_weight[0:j,i] = conf_weight[0:j,i]/np.sum(conf_weight[:,i])

    return conf, conf_weight 

def PlotSES_result (json_variables, solutions_output_file, num_of_structs):
    #name_file will be the file saved by SES, it'll just load it here.
    name_file = solutions_output_file
    #name_file =  'C:\Users\raquel\Dropbox\Fushman Internship\SES 2020\parnmr\bin\SES_output_1.txt'

    size = int(json_variables["size"])
    number = int(json_variables["number"])

    data = split(name_file, size) #split is a function which reads the file up to the size picked

    good_solrows = []
    solrows = []

    weights_sol = np.zeros(size)

    for i in range (np.size(data[:,0])):
        
        if data[i,2*size] != 0: #I just want the solutions of fixed picked size 
            for j in range(size):
                weights_sol[j] = data[i,size+1+j] #saving the weights for each solution
            
            weights_total = np.sum(weights_sol)
            
            if weights_total > 0.9 and weights_total < 1.1:
                good_solrows.append(i) #rows of the good solutions that add up to 0.9-1.1.

            solrows.append(i) #all solutions (even if they don't add up to 1)

    if (len(good_solrows) == 0):
        good_solrows = solrows

    #weight is a function that gets the structures and the weights of all conformers from the good_solrows.
    [conf,conf_weight] = weight(data,good_solrows,size)

    figure_1 = { # xlim ([0 length(good_solrows)+1])
            "data": [
                {
                    "x": np.arange(1, len(good_solrows) + 1).tolist(),
                    "y": data[good_solrows[0:], 0].tolist(),
                    "mode": "markers",
                    "marker": {
                        "color": "Red"
                    }
                }
            ],
            "layout": {
                    "title": "Chi2 for different solutions",
                    "xaxis": {
                        "title": "Solution #"
                    },
                    "yaxis": {
                        "title": "Chi2"
                    }
            }
        }

    figure_2 = { # xlim ([0 length(good_solrows)+1])
            "data": [
                {
                    "x": 100*conf_weight[0, 0:np.size(good_solrows)].tolist(),
                    "y": data[good_solrows[0:], 0].tolist(),
                    "mode": "markers",
                    "marker": {
                        "color": "Red"
                    }
                }
            ],
            "layout": {
                    "title": "Weights of the major conformer in different solutions",
                    "xaxis": {
                        "title": "Weight of first conformer (#)"
                    },
                    "yaxis": {
                        "title": "Chi2"
                    }
            }
        }


    # #The 3rd image is the one with the actual result:

    #This is working for n,bin now (output of n and bin is exactly like in Matlab)  - RGLC
    conf_0 = list(np.intc(conf[0,:]))
    res = {}

    for i in conf_0:
        res[i] = conf_0.count(i)

    j = 0
    n = np.zeros(np.size(np.unique(conf[0,:])))
    bin = np.zeros(np.size(np.unique(conf[0,:])))

    for i in res:
        n[j] = i
        bin[j] = res[i]
        j = j+1

    n_bin = np.transpose(np.vstack((n,bin)))
    n_bin = sortrows(n_bin)
    n = n_bin[:,1]
    bin = n_bin[:,0]

    #return data

    if is_checked(json_variables, "cluster_by_rmsd_checkbox"):
        tree_filename = json_variables["tree_input"][0]
        rmsd_cut = json_variables["rmsd_cut"]

        tree = np.loadtxt(tree_filename)

        T = np.zeros((num_of_structs, 2))
        new_T = np.zeros((np.size(bin), 2))
        
        T[:,0] = np.arange(num_of_structs)
        T[:,1] = fcluster(Z=tree, t=rmsd_cut, criterion="distance") # the RMSD clustering is chose to be 4A in here
        
        new_T[:,0] = bin
        new_T[:,1] = T[np.intc(bin)-1,1]
        new_T = sortrows(new_T)
        
        size_clusters = np.size(np.unique(new_T[:,1]))
        clusters = np.unique(new_T[:,1])
        aux = []
        represent = np.zeros((size_clusters,1))
        weight_rep = np.zeros((size_clusters,1))
        chi2_rep = np.zeros((size_clusters,1))
        for i in range (size_clusters):
            aux = np.where(new_T[:,1] == clusters[i])
            chi2 = np.zeros((np.size(aux),1))
            weights_clust = np.zeros((np.size(aux),size))
            for j in range (np.size(aux)):
                aux2 = np.where(data[good_solrows,1] == new_T[aux[0][j],0])
                chi2[j] = data[good_solrows[aux2[0][0]],0] 
                weights_clust[j,:] = data[good_solrows[aux2[0][0]],size+1:2*size+1]
            aux3 = np.where(chi2 == min(chi2))
            represent[i] = new_T[aux[0][aux3[0][0]],0]
            weight_rep[i] = weights_clust[aux3[0][0],0]/np.sum(weights_clust[aux3[0][0],:])
            chi2_rep[i] = chi2[aux3[0][0]]
                
        weight_rep = 100*weight_rep
        leg = np.hstack((np.hstack((chi2_rep,represent)),weight_rep))
        leg = sortrows(leg)

        if np.size(leg[:,1]) > number:
            leg = leg[0:number,:]
        
    else:

        bin_size = np.size(bin)

        lines = np.zeros(bin_size)
        save_chi2 = np.zeros(bin_size)

        for i in range(bin_size): #this loop saves the minimum chi2 for a given first conformer structure # among the solutions
            save_find = np.where(bin[i] == data[good_solrows,1])
            lines[i] = good_solrows[np.min(save_find)]
            save_chi2[i] = data[good_solrows[np.min(save_find)],0]

        aux = np.sort(save_chi2)
        y = 0

        lines = lines.astype(int)

        if number > np.size(aux):
            number = np.size(aux)

        repetition = np.zeros(number)
        confor = np.zeros(number)

        for i in range (np.size(lines)): #this loop selects only the lowest chi2 solutions based on how much solutions the user wants
            if data[lines[i],0] <= aux[number-1]:
                repetition[y] = n[i]
                confor[y] = bin[i]
                y = y + 1

        avg_weight = np.zeros(number)
        chi2 = np.zeros(number)

        for i in range (np.size(confor)):
            a = np.where(conf[0,:] == confor[i])[0]
            avg_weight[i] = np.sum(conf_weight[0,a])/np.size(a)
            chi2[i] = data[good_solrows[a[0]],0] 

        avg_weight = 100*avg_weight
        leg = np.vstack((np.vstack((chi2, confor)), avg_weight)).T
        leg = sortrows(leg)

    if (np.shape(leg)[0] < number):
        number = np.shape(leg)[0]

    ticks = leg[:, 1].tolist() # This produces the right values, but we need to make it recognized as a string
    ticks = list(map(lambda x: str(int(x)) + "‎", ticks)) # Added the empty character


    figure_3 = { # xlim ([0 length(good_solrows)+1])
            "data": [
                {
                    "x": ticks, # Convert to string
                    "y": leg[:,0].tolist(),
                    "text": list(map(lambda x: str(round(x, 5)) + "%", leg[:,2].tolist())),
                    "type": "bar",
                    "marker": {
                        "color": "Red"
                    }
                }
            ],
            "layout": {
                    "title": "Solution SES",
                    "xaxis": {
                        "title": "First Conformer"
                    },
                    "yaxis": {
                        "title": "Chi2",
                        "range": [np.min(leg[:,0])-0.05*np.min(leg[:,0]), np.max(leg[:,0])+0.1*np.max(leg[:,0])]
                    },
            }
        }

    return figure_1, figure_2, figure_3


# json_variables = {
#     "cluster_by_rmsd_checkbox": 1,
#     "tree_input": ["tree.txt"],
#     "rmsd_cut": 4
# }

# f1, f2, f3 = PlotSES_result(json_variables, "SES_output_1.txt")
# print (f1)
# print (f2)
# print (f3)

    # figure(3)
    # x = string(leg(:,2))
    # y = leg(:,1)
    # h = bar(y, 'grouped', 'r')
    # hold on
    # hold off
    # set(gca, 'XTick', 1:length(x),'XTickLabel',x)
    # ax = gca
    # ay = gca
    # hAx=gca            # get a variable for the current axes handle
    # for j=1:length(y)
    #     labels(1,j)={num2str(leg(j,3),'#.0f')+"#" + ' ' + num2str(leg(j,1),'#.1f')}
    # end   
    # labels = labels'
    # hT=[]              # placeholder for text object handles
    # for i=1:length(h)  # iterate over number of bar objects
    #   hT=[hT,text(h(i).XData+h(i).XOffset,h(i).YData,labels(:,i), ...
    #           'VerticalAlignment','bottom','horizontalalign','center')]
    # end
    # txt2 = 'weight (#) \chi^2'
    # text(0.8,max(leg(:,1))+0.05*max(leg(:,1)),txt2)
    # ax.TickDir = 'out'
    # ay.TickDir = 'out'
    # xlabel ('First Conformer')
    # set(gca, 'FontName', 'Arial', 'Fontsize', 14)
    # ylabel ('X^2')
    # ylim([min(leg(:,1))-0.05*min(leg(:,1)) max(leg(:,1))+0.1*max(leg(:,1))])
    # box off
    # #legend(strcat(' ', num2str(leg(:,2))), 'location', 'northwest')
    # title ('Solution SES')
