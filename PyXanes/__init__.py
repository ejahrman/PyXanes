import plotly.offline as py
import plotly.tools as tls
py.init_notebook_mode()

def plotly_show():
    #get fig and convert to plotly
    fig = plt.gcf()
    plotlyfig = tls.mpl_to_plotly(fig, resize=True)
    
    #fix dumb automatic formatting choices
    plotlyfig['layout']['xaxis1']['tickfont']['size']=14
    plotlyfig['layout']['xaxis1']['titlefont']['size']=16
    plotlyfig['layout']['yaxis1']['tickfont']['size']=14
    plotlyfig['layout']['yaxis1']['titlefont']['size']=16
    plotlyfig['layout']['showlegend'] = True
    
    # plot
    py.iplot(plotlyfig)