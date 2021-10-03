#!/usr/bin/python
def vec_plot(resid):
    bond_vec = []
    with open('bond_vector') as in_f:
        for line in in_f: 
            bond_vec.extend(line.strip().split('\n'))
#list of list of strings
    bond_vec = list(items[1:-1].split(',') for items in bond_vec)
    bond_x = []
    bond_y = []
    bond_z = []
    break_vec = [0,None]
    for items in bond_vec:
        bond_x.append(float(items[0]))
        bond_x.append(break_vec[0])
        bond_x.append(break_vec[1])
        bond_y.append(float(items[1]))
        bond_y.append(break_vec[0])
        bond_y.append(break_vec[1])
        bond_z.append(float(items[2]))
        bond_z.append(break_vec[0])
        bond_z.append(break_vec[1])

    label = []
    break_1 = ""
    for ii in range(len(resid)):
        mylab = "residue " + str(int(resid[ii]))
        label.append(mylab)
	label.append(break_1)
        label.append(break_1)

    #assert(len(label) == len(bond_x))

    vec_plotly = {
    "data":[
  {    
    "type": 'scatter3d',
    "mode": 'lines+markers',
    "x": bond_x,
    "y": bond_y,
    "z": bond_z,
    "text": label,
    "line": {
        "width": 6,
        "colorscale": "Viridis"},
        "marker": {
            "size": 3.5,
            "colorscale": "Greens"
  }
}],
    "layout":{
      "title": "Bond Orientations Plot",
      "showlegend": False,
      "width": 550,
      "height": 550 }
}
    return vec_plotly
