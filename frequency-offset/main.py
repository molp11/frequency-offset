from os import listdir
from os.path import join, dirname

from io import StringIO
import base64

import pandas as pd
import numpy as np
import pwlf

from sympy import Symbol

from bokeh.io import curdoc
from bokeh.layouts import column, row
from bokeh.models import ColumnDataSource, Slider, TextInput, DataTable, TableColumn, Select, CustomJS
from bokeh.models.widgets import Button
from bokeh.plotting import figure

def poly(x, A):
    
    e, out = 0, 0
    for a in A:
        out += a*x**e
        e += 1
    
    return out

def write_poly(A):
    
    out = ''
    i = 0
    for a in A:
        if i == 0:
            if a < 0:
                out = '- ' + str(np.abs(a))
            else:
                out = '+ ' + str(np.abs(a))
        elif i == len(A) - 1:
            out = str(a) + '*x**' + str(i) + ' ' + out
        else:
            if a < 0:
                out = '- ' + str(np.abs(a)) + "*x**" + str(i) + ' ' + out
            else:
                out = '+ ' + str(np.abs(a)) + "*x**" + str(i) + ' ' + out
        i+=1
        
    return out

def group_reduce(itmp):
    
    grouped = itmp.groupby('x')
    otmp = pd.DataFrame(columns=list(itmp))
    for g in grouped:
        otmp = otmp.append((g[1].max() + g[1].min())/2, ignore_index=True)
    otmp = otmp.sort_values('x')

    return otmp

def get_symbolic_eqn(pwlf_, segment_number):
    x = Symbol('x')
    if pwlf_.degree < 1:
        raise ValueError('Degree must be at least 1')
    if segment_number < 1 or segment_number > pwlf_.n_segments:
        raise ValueError('segment_number not possible')
    # assemble degree = 1 first
    for line in range(segment_number):
        if line == 0:
            my_eqn = pwlf_.beta[0] + (pwlf_.beta[1])*(x-pwlf_.fit_breaks[0])
        else:
            my_eqn += (pwlf_.beta[line+1])*(x-pwlf_.fit_breaks[line])
    # assemble all other degrees
    if pwlf_.degree > 1:
        for k in range(2, pwlf_.degree + 1):
            for line in range(segment_number):
                beta_index = pwlf_.n_segments*(k-1) + line + 1 
                my_eqn += (pwlf_.beta[beta_index])*(x-pwlf_.fit_breaks[line])**k
    return my_eqn.simplify()

def update_poly(attrname, old, new):
    
    itmp = pd.DataFrame(data={'x' : temp_source.data['y'], 't' : temp_source.data['x']})
    otmp = group_reduce(itmp)

    tx = np.linspace(min(temp_source.data['x']), max(temp_source.data['x']), num=len(temp_source.data['x']))
    pt = np.polyfit(otmp.t, otmp.x, torder.value)[::-1]
    xnew = [poly(ti, pt) for ti in temp_source.data['x']]
    ty = [poly(ti, pt) for ti in tx]
    eqn_list = [write_poly(pt)]
    segment_list = ['Temperature Sensor']
    nt = torder.value
    
    temp_source.data.ty = ty
    temp_source.data.tx = tx
    plot_source.data.x = xnew
    
    s = segments.value
    n = order.value
    
    itmp = pd.DataFrame(data={'x' : plot_source.data['x'], 'y' : plot_source.data['y']})
    otmp = group_reduce(itmp)
    
    if s == 1:
        p = np.polyfit(otmp['x'], otmp['y'], n)[::-1]
        fy = [poly(xi, p) for xi in plot_source.data['fx']]
        me = np.max([np.abs(poly(xi, p) - yi) for xi, yi in zip(plot_source.data['x'], plot_source.data['y'])])
        poly_source.data['p'] = [pi for pi in p]
        eqn_list.append(write_poly(p))
        segment_list.append(['(' + str(min(plot_source.data['x'])) + ' ' + str(max(plot_source.data['x'])) + ')'])
    else:
        p = pwlf.PiecewiseLinFit(otmp['x'], otmp['y'], degree=n)
        p.fit(s)
        fy = p.predict(plot_source.data['fx'])
        me = np.max([np.abs(p.predict(xi) - yi) for xi, yi in zip(plot_source.data['x'], plot_source.data['y'])])
        x = Symbol('x')
        for i in range(p.n_segments):
            eqn_list.append(str(get_symbolic_eqn(p, i + 1)))
        for ix in range(1, len(p.fit_breaks)):
            segment_list.append('(' + str(p.fit_breaks[ix-1]) + ' ' + str(p.fit_breaks[ix]) + ')')
    
    plot_source.data['fy'] = fy
    table_source.data = {'segment' : [s for s in segment_list], 'eqn_list' : [eq for eq in eqn_list]}
    plot.title.text = 'Max Frequency Error (ppb): ' + str(me)
        
def load_data(attrname, old, new):
    
    df = pd.read_csv(join(dirname(__file__), 'data/'+data_select.value)).dropna()
    
    load_page(df)
    
'''    xname, yname, tname = df.columns[2], df.columns[1], df.columns[3]
    for name in list(df):
        if 'mbient' in name or 'MBIENT' in name:
            xname = name
        elif 'ensor' in name or 'ENSOR' in name or 'TMP' in name:
            tname = name
        elif 'ppb' in name or 'PPB' in name:
            yname = name

    xsample_select.value = xname
    xsample_select.options = list(df)
    tsample_select.value = tname
    tsample_select.options = list(df)
    ysample_select.value = yname
    ysample_select.options = list(df)

    itmp = pd.DataFrame(data={'x' : x, 'y' : y, 't' : t})
    otmp = group_reduce(itmp)

    tx = np.linspace(min(t), max(t), num=len(t))
    pt = np.polyfit(otmp.t, otmp.x, 1)[::-1]
    xnew = [poly(ti, pt) for ti in t]
    ty = [poly(ti, pt) for ti in tx]
    eqn_list = [write_poly(pt)]
    segment_list = ['Temperature Sensor']

    itmp = pd.DataFrame(data={'x' : xnew, 'y' : y})
    otmp = group_reduce(itmp)

    p = np.polyfit(otmp.x, otmp.y, 5)[::-1]
    fx = np.linspace(min(xnew), max(xnew), num=len(xnew))
    fy = [poly(xi, p) for xi in fx]
    me = np.max([np.abs(poly(xi, p) - yi) for xi, yi in zip(xnew, y)])
    eqn_list.append([write_poly(p)])
    segment_list.append(['(' + str(min(xnew)) + ' ' + str(max(xnew)) + ')'])

    temp_source.data = dict(x=t, y=x, tx=tx, ty=ty)
    plot_source.data = dict(x=xnew, y=y, fx=fx, fy=fy)
    table_source.data = dict(segment=segment_list, eqn_list=eqn_list)
'''
def file_callback(attname, old, new):
    
    raw_contents = file_source.data['file_contents'][0]
    prefix, b64_contents = raw_contents.split(',', 1)
    file_contents = base64.b64decode(b64_contents).decode('utf-8-sig', errors='ignore')
    file_io = StringIO(file_contents)
    
    df = pd.read_csv(file_io).dropna()
    
    load_page(df)
'''    xname, yname, tname = df.columns[2], df.columns[1], df.columns[3]
    for name in list(df):
        if 'mbient' in name or 'MBIENT' in name:
            xname = name
        elif 'ensor' in name or 'ENSOR' in name or 'TMP' in name:
            tname = name
        elif 'ppb' in name or 'PPB' in name:
            yname = name

    xsample_select.value = xname
    xsample_select.options = list(df)
    tsample_select.value = tname
    tsample_select.options = list(df)
    ysample_select.value = yname
    ysample_select.options = list(df)
                                  
    x = df[xname].values
    y = df[yname].values
    t = df[tname].values

    itmp = pd.DataFrame(data={'x' : x, 'y' : y, 't' : t})
    otmp = group_reduce(itmp)

    tx = np.linspace(min(t), max(t), num=len(t))
    pt = np.polyfit(otmp.t, otmp.x, 1)[::-1]
    xnew = [poly(ti, pt) for ti in t]
    ty = [poly(ti, pt) for ti in tx]
    eqn_list = [write_poly(pt)]
    segment_list = ['Temperature Sensor']

    itmp = pd.DataFrame(data={'x' : xnew, 'y' : y})
    otmp = group_reduce(itmp)

    p = np.polyfit(otmp.x, otmp.y, 5)[::-1]
    fx = np.linspace(min(xnew), max(xnew), num=len(xnew))
    fy = [poly(xi, p) for xi in fx]
    me = np.max([np.abs(poly(xi, p) - yi) for xi, yi in zip(xnew, y)])
    eqn_list.append([write_poly(p)])
    segment_list.append(['(' + str(min(xnew)) + ' ' + str(max(xnew)) + ')'])

    temp_source.data = dict(x=t, y=x, tx=tx, ty=ty)
    plot_source.data = dict(x=xnew, y=y, fx=fx, fy=fy)
    table_source.data = dict(segment=segment_list, eqn_list=eqn_list)
'''
def choose_column(attrname, old, new):
    
    xname = xsample_select.value
    tname = tsample_select.value
    yname = ysample_select.value
    
    x = df[xname].values
    y = df[yname].values
    t = df[tname].values

    itmp = pd.DataFrame(data={'x' : x, 'y' : y, 't' : t})
    otmp = group_reduce(itmp)

    tx = np.linspace(min(t), max(t), num=len(t))
    pt = np.polyfit(otmp.t, otmp.x, 1)[::-1]
    xnew = [poly(ti, pt) for ti in t]
    ty = [poly(ti, pt) for ti in tx]
    eqn_list = [write_poly(pt)]
    segment_list = ['Temperature Sensor']

    itmp = pd.DataFrame(data={'x' : xnew, 'y' : y})
    otmp = group_reduce(itmp)

    p = np.polyfit(otmp.x, otmp.y, 5)[::-1]
    fx = np.linspace(min(xnew), max(xnew), num=len(xnew))
    fy = [poly(xi, p) for xi in fx]
    me = np.max([np.abs(poly(xi, p) - yi) for xi, yi in zip(xnew, y)])
    eqn_list.append([write_poly(p)])
    segment_list.append(['(' + str(min(xnew)) + ' ' + str(max(xnew)) + ')'])

    temp_source.data = dict(x=t, y=x, tx=tx, ty=ty)
    plot_source.data = dict(x=xnew, y=y, fx=fx, fy=fy)
    table_source.data = dict(segment=segment_list, eqn_list=eqn_list)

def load_page(df):
    
    xname, yname, tname = df.columns[2], df.columns[1], df.columns[3]
    for name in list(df):
        if 'mbient' in name or 'MBIENT' in name:
            xname = name
        elif 'ensor' in name or 'ENSOR' in name or 'TMP' in name:
            tname = name
        elif 'ppb' in name or 'PPB' in name:
            yname = name

    x = df[xname].values
    y = df[yname].values
    t = df[tname].values

    itmp = pd.DataFrame(data={'x' : x, 'y' : y, 't' : t})
    otmp = group_reduce(itmp)

    tx = np.linspace(min(t), max(t), num=len(t))
    pt = np.polyfit(otmp.t, otmp.x, 1)[::-1]
    xnew = [poly(ti, pt) for ti in t]
    ty = [poly(ti, pt) for ti in tx]
    eqn_list = [write_poly(pt)]
    segment_list = ['Temperature Sensor']

    itmp = pd.DataFrame(data={'x' : xnew, 'y' : y})
    otmp = group_reduce(itmp)

    p = np.polyfit(otmp.x, otmp.y, 5)[::-1]
    fx = np.linspace(min(xnew), max(xnew), num=len(xnew))
    fy = [poly(xi, p) for xi in fx]
    me = np.max([np.abs(poly(xi, p) - yi) for xi, yi in zip(xnew, y)])
    eqn_list.append([write_poly(p)])
    segment_list.append(['(' + str(min(xnew)) + ' ' + str(max(xnew)) + ')'])

    global xsample_select
    xsample_select = Select(title='Select Temperature Sensor', options=list(df), value=xname)
    xsample_select.on_change('value', choose_column)
    
    global tsample_select
    tsample_select = Select(title='Select Ambient Temperature', options=list(df), value=tname)
    tsample_select.on_change('value', choose_column)
    
    global ysample_select
    ysample_select = Select(title='Select Frequency Offset', options=list(df), value=yname)
    ysample_select.on_change('value', choose_column)

    global temp_source
    temp_source = ColumnDataSource(data=dict(x=t, y=x, tx=tx, ty=ty))
    global plot_source
    plot_source = ColumnDataSource(data=dict(x=xnew, y=y, fx=fx, fy=fy))
    global poly_source
    poly_source = ColumnDataSource(data=dict(p=p))
    global table_source
    table_source = ColumnDataSource(data=dict(segment=segment_list, eqn_list=eqn_list))

    global tplot
    tplot = figure(plot_height=400, plot_width=500, tools="crosshair,pan,reset,save,wheel_zoom")
    tplot.xaxis.axis_label = 'Temperature Sensor'
    tplot.yaxis.axis_label = 'Ambient Temperature'
    tplot.circle('x', 'y', source=temp_source, fill_alpha=0.6, color='grey')
    tplot.line('tx', 'ty', source=temp_source, line_width=3, line_alpha=0.8, color='black')

    global plot
    plot = figure(plot_height=400, plot_width=500, tools="crosshair,pan,reset,save,wheel_zoom",
                 title='Max Frequency Error (ppb): ' + str(me))
    plot.xaxis.axis_label = 'Temperature Estimate'
    plot.yaxis.axis_label = 'Frequency Offset (ppb)'
    plot.circle('x', 'y', source=plot_source, fill_alpha=0.6, color='grey')
    plot.line('fx', 'fy', source=plot_source, line_width=3, line_alpha=0.8, color='black')

    columns = [
        TableColumn(field="segment", title="Temp. Range"),
        TableColumn(field="eqn_list", title="Equation")
    ]
    
    global data_table
    data_table = DataTable(source=table_source, columns=columns, width=1000)

    #temp = Slider(title="Temperature", value=25., start=np.min(x), end=np.max(x), step=1., width=800)
    #temp.on_change('value', update_temp)

    global torder
    torder = Slider(title="Temperature: Polynomial Order", value=1, start=1, end=5, step=1, width=200)
    torder.on_change('value', update_poly)

    global order
    order = Slider(title="Frequency: Polynomial Order", value=5, start=1, end=5, step=1, width=200)
    order.on_change('value', update_poly)
    
    global segments
    segments = Slider(title="Frequency: # Segments", value=1, start=1, end=5, step=1, width=200)
    segments.on_change('value', update_poly)

    global download_button
    download_button = Button(label="Download Table to CSV", button_type="primary", width=200)
    download_button.callback = CustomJS(args=dict(source=table_source, file_name='output.csv'),
                                   code=open(join(dirname(__file__), "download.js")).read())

    global file_source
    file_source = ColumnDataSource({'file_contents':[], 'file_name':[]})
    file_source.on_change('data', file_callback)

    global upload_button
    upload_button = Button(label="Upload Local File", button_type="success", width=200)
    upload_button.callback = CustomJS(args=dict(file_source=file_source),
                                   code=open(join(dirname(__file__), "upload.js")).read())
    
    curdoc().clear()
    curdoc().add_root(column(row(data_select, upload_button, download_button), row(column(xsample_select, tsample_select, ysample_select, torder, order, segments), column(tplot), column(plot)), data_table))
    curdoc().title = "frequency-offset"
    
#def update_temp(attrname, old, new):
#
#    offset = poly(temp.value, poly_source.data['p'])
#    plot_source.data['y'] -= offset    
#    p = poly_source.data['p']
#    p[0] = p[0] - offset
#    poly_source.data['p'] = [pi for pi in p]
#    plot_source.data['fy'] = [poly(xi, poly_source.data['p']) for xi in plot_source.data['fx']]

sample_data = listdir(join(dirname(__file__), 'data'))
parts = []
for sd in sample_data:
    parts.append(sd.split('.')[0])
    
data_select = Select(title='Select Sample Data File', options=sample_data, value=sample_data[0])
data_select.on_change('value', load_data)

df = pd.read_csv(join(dirname(__file__), 'data/'+data_select.value))

load_page(df)