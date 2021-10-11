#!/opt/miniconda2/bin/python

def load_f_lists(f,deli):
    tmp = [x.strip().split(deli) for x in f.read().splitlines()]
    return tmp

def chi2_plot():
    dot_x =[]
    dot_y =[]
    dot_z =[]
    with open('chi2_data') as in_f:
        chi2_scatter = load_f_lists(in_f, '\t')
    for dot in chi2_scatter:
        dot_x.append(float(dot[0]))
        dot_y.append(float(dot[1]))
        dot_z.append(float(dot[2]))
    surf_x = []
    surf_y = []
    surf_z = []
    with open('chi2_surface') as in_f:
        chi2_surf_data = load_f_lists(in_f, ' ')
    for surfs in chi2_surf_data:
        surf_x.append(float(surfs[0]))
        surf_y.append(float(surfs[1]))
        surf_z.append(float(surfs[2]))   

    chi2_plotly = { 
   "data":[
  {
    "opacity":0.8,
    "colorscale":[
    [ 0 ,'rgb(255,0,0)'],
    [0.5,'rgb(0, 255, 0)'],
    [1, 'rgb(0, 0, 255)']
],
    "type": 'mesh3d',
    "x": surf_x,
    "y": surf_y,
    "z": surf_z,
    "intensity": surf_z,
    "showlegend": True
  },
  {
    "x": dot_x, 
    "y": dot_y, 
    "z": dot_z,
    "mode": 'markers',
    "marker": {
		"size": 2,
		"line": {
		    "color": 'rgba(217, 217, 217, 0.14)',
		    "width": 0.5},
		"opacity": 0.8},
    "type": 'scatter3d',
    "showlegend": False
  }
],
   "layout": { 
     "title": "Chi-square Plot",
     "scene":{
       "xaxis": {
         "title":"alpha"},
       "yaxis": {
         "title":"beta"},
       "zaxis": {
         "title":"Chi^2"
      }
    },
     "width": 550,
     "height": 550 
  }
}
    return chi2_plotly

