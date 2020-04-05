from os import listdir
from os.path import join, dirname

from io import StringIO
import base64

import pandas as pd
import numpy as np
import pwlf

from sympy import Symbol, Number

from bokeh.io import curdoc
from bokeh.layouts import column, row
from bokeh.models import ColumnDataSource, Slider, TextInput, DataTable, TableColumn, Select, CustomJS, CheckboxButtonGroup
from bokeh.models.widgets import Button
from bokeh.plotting import figure


# general methods
def round_expr(expr, num_digits):
    
    '''
    reformat numbers in sympy expression
    '''
    
    return expr.xreplace({n : round(n, num_digits) for n in expr.atoms(Number)})


def hz2ppb(x, a):
    
    '''
    convert frequency in Hz to frequency offset in ppb
    '''

    return (x - a) * 1e9 / a


def poly(x, A):
    
    '''
    return N-degree polynomial given number 'x' and coefficient array 'A'
    '''
    
    e, out = 0, 0
    for a in A:
        out += a*x**e
        e += 1

    return out


def write_poly(A):
    
    '''
    return polynomial equation as text given coefficient array 'A'
    '''
    
    out = ''
    i = 0
    for a in A:
        if i == 0:
            if a < 0:
                out = '- ' + str("{:.8e}".format(np.abs(a)))
            else:
                out = '+ ' + str("{:.8e}".format(np.abs(a)))
        elif i == len(A) - 1:
            out = str("{:.8e}".format(a)) + '*x**' + str(i) + ' ' + out
        else:
            if a < 0:
                out = '- ' + str("{:.8e}".format(np.abs(a))) + "*x**" + str(i) + ' ' + out
            else:
                out = '+ ' + str("{:.8e}".format(np.abs(a))) + "*x**" + str(i) + ' ' + out
        i+=1
        
    return out


def group_reduce(df):
    
    '''
    average multiple frequency (offset) values measured at each temperature 
    '''
    
    x = df[xsample_select.value].values
    y = df[ysample_select.value].values
    t = df[tsample_select.value].values    
    inp = pd.DataFrame(data={'x' : x, 'y' : y, 't' : t})
    grouped = inp.groupby('x')
    out = pd.DataFrame(columns=list(inp))
    for g in grouped:
        out = out.append((g[1].max() + g[1].min())/2, ignore_index=True)
    out = out.sort_values('x')

    return out


def get_symbolic_eqn(pwlf_, segment_number):
    
    '''
    convert pwlf output to equation in text format
    implementation from https://github.com/cjekel/piecewise_linear_fit_py/issues/44
    '''
    
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
    
    out = round_expr(my_eqn.simplify(), 8)
    
    return out


def fit_temp(df):
    
    '''
    fit sensor measurements to estimate ambient temperature
    '''
    
    t = df[xsample_select.value].values # sensor temperature readings    
    n = int(torder.value)
    reduced = group_reduce(df)
    tx = np.linspace(min(t), max(t), num=len(t))
    pt = np.polyfit(reduced.x, reduced.t, n)[::-1]
    xnew = [poly(ti, pt) for ti in t]
    ty = [poly(ti, pt) for ti in tx]
    eqn = write_poly(pt)

    return  tx, ty, xnew, eqn


def fit_poly(df, xnew):
    
    '''
    fit relationship between frequency and estimated temperature
    '''
    
    eqn_list, segment_list = [], []
    s = int(segments.value)
    n = int(order.value)
    y = df[ysample_select.value].values
    reduced = group_reduce(df)
    fx = np.linspace(min(xnew), max(xnew), num=len(xnew))
    
    if s == 1:
        p = np.polyfit(reduced['t'], reduced['y'], n)[::-1]
        fy = [poly(xi, p) for xi in fx]
        me = np.max([np.abs(poly(xi, p) - yi) for xi, yi in zip(xnew, y)])
        eqn_list.append(write_poly(p))
        segment_list.append('(' + str(np.round(min(xnew), 2)) + ' ' + str(np.round(max(xnew), 2)) + ')')
    
    else:
        p = pwlf.PiecewiseLinFit(reduced['t'], reduced['y'], degree=n)
        p.fit(s)
        fy = p.predict(plot_source.data['fx'])
        me = np.max([np.abs(p.predict(xi) - yi) for xi, yi in zip(xnew, y)])
        x = Symbol('x')
        for i in range(p.n_segments):
            eqn_list.append(str(get_symbolic_eqn(p, i + 1)))
        for ix in range(1, len(p.fit_breaks)):
            segment_list.append('(' + str(np.round(p.fit_breaks[ix-1], 2)) + ' ' + str(np.round(p.fit_breaks[ix], 2)) + ')')
    
    x = fx
    if 'x' in pcustom.value:
        cy = eval(pcustom.value)
    
    else:
        cy = fy
        
    return fx, fy, cy, eqn_list, segment_list, me
    
    
# callback methods
def update(attrname, old, new):
    
    '''
    update plots based on widget events
    '''
    
    # temporary title
    plot.title.text = 'calculating...'
    
    # fit ambient temperature based on sensor readings
    df_copy = df.copy()
    if len(transform.active) > 0:
        df_copy[ysample_select.value] = hz2ppb(df[ysample_select.value], float(basef.value))
        plot.yaxis.axis_label = 'Frequency Variation from Nominal (ppb)'
    else:
        plot.yaxis.axis_label = 'Frequency (Hz)'
    tx, ty, xnew, eqn0 = fit_temp(df_copy)
    
    # fit relationship between frequency and estimated temperature
    fx, fy, cy, eqn1, seg1, me = fit_poly(df_copy, xnew)
    
    # update sources
    x = df_copy[xsample_select.value].values
    y = df_copy[ysample_select.value].values
    t = df_copy[tsample_select.value].values
    temp_source.data = dict(x=x, y=t, tx=tx, ty=ty)
    plot_source.data = dict(x=xnew, y=y, fx=fx, fy=fy, cy=cy)
    segment_list = ['Temperature Sensor'] + seg1
    eqn_list = [eqn0] + eqn1
    table_source.data = dict(eqn_list=eqn_list, segment=segment_list)
    
    # update title
    if len(transform.active) > 0:
        plot.title.text = 'Max error of fit (ppb): ' + str(me)
    else:
        plot.title.text = 'Max error of fit (Hz): ' + str(me)
        
        
def load_data(attrname, old, new):
    
    '''
    load data when new file source is selected from sample data dropdown list
    '''
    
    global df
    df = pd.read_csv(join(dirname(__file__), 'data/'+data_select.value)).dropna()

    load_page(df)

    
def file_callback(attname, old, new):
    
    '''
    read uploaded data and reload page
    '''

    raw_contents = file_source.data['file_contents'][0]
    prefix, b6_contents = raw_contents.split(',', 1)
    file_contents = base64.b64decode(b64_contents).decode('utf-8-sig', errors='ignore')
    file_io = StringIO(file_contents)

    global df
    df = pd.read_csv(file_io).dropna()

    load_page(df)
    

# method to (re-)generate page
def load_page(df):
    
    '''
    load page with new data selected from sample files or local file
    '''
    
    # take initial guess at column contents
    xname, yname, tname = df.columns[1], df.columns[0], df.columns[2]
    for name in list(df):
        if 'mbient' in name or 'MBIENT' in name:
            tname = name
        elif 'ensor' in name or 'ENSOR' in name or 'TMP' in name:
            xname = name
        elif 'requ' in name or 'REQU' in name:
            yname = name

    # select temperature sensor data
    global xsample_select
    xsample_select = Select(title='Temperature Sensor Measurements', options=list(df), value=xname)
    xsample_select.on_change('value', update)

    # select ambient temperature data
    global tsample_select
    tsample_select = Select(title='Chamber Ambient Temperature', options=list(df), value=tname)
    tsample_select.on_change('value', update)

    # select frequency data
    global ysample_select
    ysample_select = Select(title='Frequency Measurements', options=list(df), value=yname)
    ysample_select.on_change('value', update)

    # select polynomial order for temperature estimation
    global torder
    torder = TextInput(title="Temperature Estimate: Polynomial Order", value='1', width=250)
    torder.on_change('value', update)

    # select order for frequency vs temperature polynomial
    global order
    order = TextInput(title="Frequency Fit: Polynomial Order", value='5', width=250)
    order.on_change('value', update)

    # select number of segments to optimize
    global segments
    segments = TextInput(title="Frequency Fit: # Segments", value='1', width=250)
    segments.on_change('value', update)
    
    # set base frequency
    global basef
    basef = TextInput(title="Set Nominal Frequency", value='1e7', width=250)
    basef.on_change('value', update)
    
    # convert frequency to ppb offset
    global transform
    transform = CheckboxButtonGroup(labels=["Convert Hz to ppb"], active=[], width=250)
    transform.on_change('active', update)
    
    # calculate user-defined polynomial fit
    global pcustom
    pcustom = TextInput(title="Fit User-Defined Equation", value='', width=250)
    pcustom.on_change('value', update)
    
    # fit ambient temperature based on sensor readings
    tx, ty, xnew, eqn0 = fit_temp(df)
    
    # fit relationship between frequency and estimated temperature
    fx, fy, cy, eqn1, seg1, me = fit_poly(df, xnew)
    
    # column data for left plot
    global temp_source
    temp_source = ColumnDataSource(data=dict(x=df[xname].values, y=df[tname].values, tx=tx, ty=ty))
    
    # column data for right plot
    global plot_source
    plot_source = ColumnDataSource(data=dict(x=xnew, y=df[ysample_select.value].values, fx=fx, fy=fy, cy=cy))
    
    # column data for equation table
    segment_list = ['Temperature Sensor'] + seg1
    eqn_list = [eqn0] + eqn1
    global table_source
    table_source = ColumnDataSource(data=dict(segment=segment_list, eqn_list=eqn_list))

    # set up left plot
    global tplot
    tplot = figure(plot_height=400, plot_width=500, tools="crosshair,pan,reset,save,wheel_zoom")
    tplot.xaxis.axis_label = 'Temperature Sensor'
    tplot.yaxis.axis_label = 'Ambient Temperature'
    tplot.circle('x', 'y', source=temp_source, fill_alpha=0.6, color='grey')
    tplot.line('tx', 'ty', source=temp_source, line_width=3, line_alpha=0.8, color='black')
    
    # set up right plot
    global plot
    plot = figure(plot_height=400, plot_width=500, tools="crosshair,pan,reset,save,wheel_zoom",
                 title='Max error of fit (Hz): ' + str(me))
    plot.xaxis.axis_label = 'Temperature Estimate'
    plot.yaxis.axis_label = 'Frequency (Hz)'
    plot.circle('x', 'y', source=plot_source, fill_alpha=0.6, color='grey')
    plot.line('fx', 'fy', source=plot_source, line_width=3, line_alpha=1., color='black')
    plot.line('fx', 'cy', source=plot_source, line_width=3, line_alpha=1., color='black', line_dash="4 4")
    
    # set up data table for equations
    columns = [
        TableColumn(field="segment", title="Temp. Range", width=200),
        TableColumn(field="eqn_list", title="Equation", width=800)
    ]
    global data_table
    data_table = DataTable(source=table_source, columns=columns, width=1000, fit_columns=False)

    # download fit table to csv
    global download_button
    download_button = Button(label="Download Table to CSV", button_type="primary", width=200)
    download_button.callback = CustomJS(args=dict(source=table_source, file_name='output.csv'),
                                   code=open(join(dirname(__file__), "download.js")).read())

    # reload page when new sample data is selected
    global file_source
    file_source = ColumnDataSource({'file_contents':[], 'file_name':[]})
    file_source.on_change('data', file_callback)

    # reload page wehen new data file is uploaded
    global upload_button
    upload_button = Button(label="Upload Local File", button_type="success", width=200)
    upload_button.callback = CustomJS(args=dict(file_source=file_source),
                                   code=open(join(dirname(__file__), "upload.js")).read())
    
    # clear page and display new
    curdoc().clear()
    curdoc().add_root(column(row(data_select, upload_button, download_button), 
                             row(column(xsample_select, tsample_select, ysample_select, torder, order, segments, basef, transform, pcustom), 
                                 column(row(tplot, plot), data_table))))
    curdoc().title = "frequency-offset"


# assemble list of files in sample data directory
sample_data = listdir(join(dirname(__file__), 'data'))
parts = [sd.split('.')[0] for sd in sample_data]

# create dropdown menu to select different sample data file
data_select = Select(title='Select Sample Data File', 
                     options=sample_data, value='data-180601001-2123.csv')
data_select.on_change('value', load_data)

# read sample file, dropping rows containing any empty column values
df = pd.read_csv(join(dirname(__file__), 'data/'+data_select.value)).dropna()

load_page(df)