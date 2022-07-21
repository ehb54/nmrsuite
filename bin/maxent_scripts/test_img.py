#!/opt/miniconda3/bin/python

import plotly.graph_objects as go
import os 
import shutil

line_plot = {
        "data": [
            {
                "x": [-1, 0, 1, 2, 3],
                "y": [1, 0, 1, 4, 9],
                "mode": "markers",
                "marker": {
                    "color": "Black"
                },
                "name": "Test"
            },   
        ]
}

line_plot_figure = go.Figure(line_plot)
run_directory = "image"
if (os.path.exists(run_directory)):
    shutil.rmtree(run_directory)
os.mkdir(run_directory)
line_plot_figure.write_image(os.path.join(run_directory, "line_plot.png"))
