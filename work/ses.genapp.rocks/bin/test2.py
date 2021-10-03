import plotly
import plotly.graph_objs as go

x = [0, 1, 2, 3, 4, 5, 6]
y = [150, 35, 5, 1, .1, 0]

from plotly.subplots import make_subplots
import plotly.graph_objects as go

fig = make_subplots(rows=1, cols=2)

fig.add_trace(
    go.Scatter(x=x, y=y),
    row=1, col=2
)

fig.update_layout(xaxis_type="log")

fig.add_trace(
    go.Scatter(x=x, y=y),
    row=1, col=1
)

fig.update_layout(height=600, width=1600, title_text="Side By Side Subplots")
plotly.offline.plot(fig, filename='name.html')

""" plotly.offline.plot({
    "data": [go.Scatter(x=x, y=y)],
        "layout": go.Layout(title="line chart")
	}, auto_open=True)

 """